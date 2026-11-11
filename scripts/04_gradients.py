# -*- coding: utf-8 -*-
"""04 表型 -> 坏死-感染进程：三库结局梯度 + 调整OR + 超越早期严重度 + 连续风险敏感性
   输出 results/:
     rates_by_pheno.csv      各库x结局x表型 n/events/rate/Wilson CI + chi2 p + trend p
     adj_or.csv              逐结局 logistic (P1参照) 多套调整: crude / +age+sex / +etiology / +earlySAP / +CT(本院)
     beyond_severity_lr.csv  表型在 (age+sex+etiology+earlySAP) 之上的 LR 检验
     rates_nonsevere.csv     非重症亚组内梯度
     risk_continuous.csv     六变量连续风险 risk6 / 简模型 risk3: 每+1SD OR 与 四分位率
     rates_sens.csv          敏感性: 剔除70例全插补 / MIMIC坏死未指明按无 / 完整病例
"""
import pandas as pd
import numpy as np
import os, sys, warnings
import statsmodels.api as sm
from scipy.stats import chi2_contingency, chi2 as chi2dist, norm
warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = lambda var, sub: os.environ.get(var, os.path.join(ROOT, "inputs", sub))
P = ROOT
RES, DAT = os.path.join(P, "results"), os.path.join(P, "data")

loc = pd.read_csv(os.path.join(DAT, "local.csv"))
mim = pd.read_csv(os.path.join(DAT, "mimic.csv"))
eic = pd.read_csv(os.path.join(DAT, "eicu.csv"))
mim_icd = mim[mim.necro_specified == 1].copy()          # 坏死状态已编码 1649
mim_icd10 = mim[mim.icd10_ap == 1].copy()               # ICD-10 子集 1963 (未指明按无, 敏感性)
mim_icd10["necro_coded0"] = mim_icd10.necro_coded.fillna(0)
mim_icd10["ipn_coded0"] = mim_icd10.ipn_coded.fillna(0)

def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    den = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
    return (c - h, c + h)

def rates(df, y, cohort, outcome, extra=None):
    rows = []
    tab = pd.crosstab(df.pheno, df[y])
    chi_p = chi2_contingency(tab).pvalue if tab.shape == (4, 2) else np.nan
    try:
        tr = sm.Logit(df[y].values, sm.add_constant(df.pheno.values.astype(float))).fit(disp=0)
        trend_p = tr.pvalues[1]
    except Exception:
        trend_p = np.nan
    for ph in [1, 2, 3, 4]:
        s = df[df.pheno == ph]
        k, n = int(s[y].sum()), len(s)
        lo, hi = wilson(k, n)
        rows.append(dict(cohort=cohort, outcome=outcome, phenotype=f"P{ph}", n=n, events=k,
                         rate=k / n if n else np.nan, lo=lo, hi=hi, chi2_p=chi_p, trend_p=trend_p,
                         **(extra or {})))
    return rows

# ---------- 结局清单 ----------
LOCAL_OUT = [("necroCT", "CT necrosis"),
             ("intervention", "Drainage/necrosectomy"), ("pcd", "Percutaneous drainage"),
             ("LOS14", "LOS > 14 d"), ("persOF", "Persistent OF"), ("icu", "ICU admission"),
             ("death", "In-hospital death"), ("mctsi_high", "MCTSI >= 6 (admission CT)"),
             ("apfc", "Peripancreatic collection (CT)"), ("SAP", "SAP (recorded)")]
MIMIC_OUT = [("necro_coded", "Coded necrosis"), ("ipn_coded", "Coded infected necrosis")]
MIMIC_FULL = [("sepsis", "Sepsis code"), ("epi_any", "Extra-pancreatic infection"),
              ("pneumonia", "Pneumonia"), ("bacteremia", "Bacteremia"), ("uti", "UTI"),
              ("cholangitis", "Cholangitis"), ("cdiff", "C. difficile"),
              ("LOS14", "LOS > 14 d"), ("death", "In-hospital death"),
              ("support_persistent", "Organ support >= 48 h"), ("icu", "ICU admission"),
              ("SAP_proxy", "SAP proxy")]
EICU_OUT = [("necro_dx", "Necrotizing (diagnosis string)"), ("abscess_dx", "Abscess/infected collection"),
            ("sepsis_dx", "Sepsis"), ("drain_dx", "Percutaneous drainage"), ("epi_any", "Extra-pancreatic infection"),
            ("LOS14", "LOS > 14 d"), ("death", "In-hospital death")]

rows = []
for y, lab in LOCAL_OUT:
    df = loc
    rows += rates(df, y, "Local", lab)
for y, lab in MIMIC_OUT:
    d_ = mim_icd.copy(); d_[y] = d_[y].astype(int)
    rows += rates(d_, y, "MIMIC-IV (ICD-10, necrosis specified)", lab)
for y, lab in MIMIC_FULL:
    rows += rates(mim, y, "MIMIC-IV", lab)
for y, lab in EICU_OUT:
    rows += rates(eic, y, "eICU", lab)
rt = pd.DataFrame(rows)
rt.to_csv(os.path.join(RES, "rates_by_pheno.csv"), index=False)
print("=== 结局梯度 ===")
piv = rt.assign(txt=lambda d: (d.rate * 100).round(1).astype(str) + "% (" + d.events.astype(str) + "/" + d.n.astype(str) + ")")
print(piv.pivot_table(index=["cohort", "outcome"], columns="phenotype", values="txt", aggfunc="first").to_string())
print(rt.groupby(["cohort", "outcome"]).trend_p.first().map(lambda p: f"{p:.2e}").to_string())

# ---------- 调整 OR ----------
def fit_or(df, y, covars, cohort, outcome, model):
    X = pd.get_dummies(df.pheno.astype(int), prefix="P").astype(float).iloc[:, 1:]
    for c in covars:
        if c == "etiology":
            X = pd.concat([X, pd.get_dummies(df.etiology, prefix="E", drop_first=True).astype(float)], axis=1)
        else:
            X[c] = df[c].astype(float).values
    X = sm.add_constant(X)
    yv = df[y].astype(int).values
    out = []
    try:
        m = sm.Logit(yv, X).fit(disp=0, maxiter=200, method="bfgs")
        for cn in ["P_2", "P_3", "P_4"]:
            b, se = m.params[cn], m.bse[cn]
            sep = abs(b) > 8 or se > 8
            out.append(dict(cohort=cohort, outcome=outcome, model=model, pheno=cn.replace("P_", "P"),
                            OR=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se),
                            p=m.pvalues[cn], n=len(df), events=int(yv.sum()), separation=int(sep)))
    except Exception as ex:
        for cn in ["P2", "P3", "P4"]:
            out.append(dict(cohort=cohort, outcome=outcome, model=model, pheno=cn, OR=np.nan, lo=np.nan,
                            hi=np.nan, p=np.nan, n=len(df), events=int(yv.sum()), separation=1))
    return out

def lr_beyond(df, y, covars, cohort, outcome, label):
    """表型在 covars 之上的 LR 检验"""
    def design(with_pheno):
        X = pd.DataFrame(index=df.index)
        for c in covars:
            if c == "etiology":
                X = pd.concat([X, pd.get_dummies(df.etiology, prefix="E", drop_first=True).astype(float)], axis=1)
            else:
                X[c] = df[c].astype(float).values
        if with_pheno:
            X = pd.concat([X, pd.get_dummies(df.pheno.astype(int), prefix="P").astype(float).iloc[:, 1:]], axis=1)
        return sm.add_constant(X)
    yv = df[y].astype(int).values
    m0 = sm.Logit(yv, design(False)).fit(disp=0, maxiter=200, method="bfgs")
    m1 = sm.Logit(yv, design(True)).fit(disp=0, maxiter=200, method="bfgs")
    lr = 2 * (m1.llf - m0.llf)
    return dict(cohort=cohort, outcome=outcome, base_model=label, LR=lr, df=3, p=chi2dist.sf(lr, 3),
                n=len(df), events=int(yv.sum()), auc_base=np.nan, auc_pheno=np.nan)

from sklearn.metrics import roc_auc_score
orr, lrr = [], []
# 本院
loc["Sex"] = loc.Sex_male
L_COV = {"crude": [], "age_sex": ["Age", "Sex"], "age_sex_etio": ["Age", "Sex", "etiology"],
         "age_sex_etio_SAP": ["Age", "Sex", "etiology", "SAP"],
         "age_sex_etio_SAP_MCTSI": ["Age", "Sex", "etiology", "SAP", "MCTSI"]}
for y, lab in LOCAL_OUT:
    if y in ("SAP", "mctsi_high"):
        continue
    df = loc
    df = df.assign(Sex=df.Sex_male)
    for mname, cov in L_COV.items():
        if y in ("necroCT", "apfc") and "MCTSI" in cov:
            continue
        if y == "persOF" and "SAP" in cov:
            continue
        orr += fit_or(df, y, cov, "Local", lab, mname)
    if y not in ("persOF",):
        lrr.append(lr_beyond(df, y, ["Age", "Sex", "etiology", "SAP"], "Local", lab, "age+sex+etiology+SAP"))
        if y not in ("necroCT", "apfc"):
            lrr.append(lr_beyond(df, y, ["Age", "Sex", "etiology", "SAP", "MCTSI"], "Local", lab, "age+sex+etiology+SAP+MCTSI"))
# MIMIC
mim["Sex"] = mim.Sex_male; mim["Age"] = mim.Age_imp
mim_icd = mim[mim.necro_specified == 1].copy()
M_COV = {"crude": [], "age_sex": ["Age", "Sex"], "age_sex_etio": ["Age", "Sex", "etiology"],
         "age_sex_etio_SAP": ["Age", "Sex", "etiology", "SAP_proxy"]}
for y, lab in MIMIC_OUT:
    d_ = mim_icd.copy(); d_[y] = d_[y].astype(int)
    for mname, cov in M_COV.items():
        orr += fit_or(d_, y, cov, "MIMIC-IV (ICD-10, necrosis specified)", lab, mname)
    lrr.append(lr_beyond(d_, y, ["Age", "Sex", "etiology", "SAP_proxy"], "MIMIC-IV (ICD-10, necrosis specified)", lab, "age+sex+etiology+SAPproxy"))
for y, lab in MIMIC_FULL:
    if y == "SAP_proxy":
        continue
    for mname, cov in M_COV.items():
        if y == "support_persistent" and "SAP_proxy" in cov:
            continue
        orr += fit_or(mim, y, cov, "MIMIC-IV", lab, mname)
    if y != "support_persistent":
        lrr.append(lr_beyond(mim, y, ["Age", "Sex", "etiology", "SAP_proxy"], "MIMIC-IV", lab, "age+sex+etiology+SAPproxy"))
# eICU (APACHE 调整, 缺失指示)
eic["Sex"] = eic.Sex_male; eic["Age"] = eic.Age_imp
eic["apache_miss"] = eic.apachescore.isna().astype(int)
eic["apache_f"] = eic.apachescore.fillna(eic.apachescore.median())
E_COV = {"crude": [], "age_sex": ["Age", "Sex"], "age_sex_apache": ["Age", "Sex", "apache_f", "apache_miss"]}
for y, lab in EICU_OUT:
    for mname, cov in E_COV.items():
        orr += fit_or(eic, y, cov, "eICU", lab, mname)
    lrr.append(lr_beyond(eic, y, ["Age", "Sex", "apache_f", "apache_miss"], "eICU", lab, "age+sex+APACHE"))
orf = pd.DataFrame(orr); orf.to_csv(os.path.join(RES, "adj_or.csv"), index=False)
lrf = pd.DataFrame(lrr); lrf.to_csv(os.path.join(RES, "beyond_severity_lr.csv"), index=False)
print("\n=== P4 vs P1 调整OR (age_sex_etio_SAP 或对应) ===")
sel = orf[(orf.pheno == "P4") & orf.model.isin(["age_sex_etio_SAP", "age_sex_apache"])]
print(sel[["cohort", "outcome", "model", "OR", "lo", "hi", "p", "events", "separation"]].round(3).to_string(index=False))
print("\n=== 表型超越早期严重度 LR ===")
print(lrf[["cohort", "outcome", "base_model", "LR", "p", "events"]].round(4).to_string(index=False))

# ---------- 非重症亚组梯度 ----------
rows = []
for y, lab in [("necroCT", "CT necrosis"),
               ("intervention", "Drainage/necrosectomy"), ("LOS14", "LOS > 14 d")]:
    df = loc
    rows += rates(df[df.SAP == 0], y, "Local non-SAP", lab)
for y, lab in MIMIC_OUT:
    d_ = mim_icd[mim_icd.SAP_proxy == 0].copy(); d_[y] = d_[y].astype(int)
    rows += rates(d_, y, "MIMIC-IV non-SAPproxy (ICD-10, necrosis specified)", lab)
for y, lab in [("epi_any", "Extra-pancreatic infection"), ("sepsis", "Sepsis code"), ("LOS14", "LOS > 14 d")]:
    rows += rates(mim[mim.SAP_proxy == 0], y, "MIMIC-IV non-SAPproxy", lab)
ns = pd.DataFrame(rows); ns.to_csv(os.path.join(RES, "rates_nonsevere.csv"), index=False)
print("\n=== 非重症亚组 ===")
nsp = ns.assign(txt=lambda d: (d.rate * 100).round(1).astype(str) + "% (" + d.events.astype(str) + "/" + d.n.astype(str) + ")")
print(nsp.pivot_table(index=["cohort", "outcome"], columns="phenotype", values="txt", aggfunc="first").to_string())
print(ns.groupby(["cohort", "outcome"]).trend_p.first().map(lambda p: f"{p:.2e}").to_string())

# ---------- 连续风险敏感性 ----------
def cont(df, y, score, cohort, outcome):
    z = (df[score] - df[score].mean()) / df[score].std()
    yv = df[y].astype(int).values
    m = sm.Logit(yv, sm.add_constant(z.values)).fit(disp=0)
    b, se = m.params[1], m.bse[1]
    q = pd.qcut(df[score], 4, labels=["Q1", "Q2", "Q3", "Q4"])
    out = dict(cohort=cohort, outcome=outcome, score=score, OR_perSD=np.exp(b), lo=np.exp(b - 1.96 * se),
               hi=np.exp(b + 1.96 * se), p=m.pvalues[1], auc=roc_auc_score(yv, df[score]), n=len(df), events=int(yv.sum()))
    for qq in ["Q1", "Q2", "Q3", "Q4"]:
        s = yv[q == qq]
        out[f"rate_{qq}"] = s.mean(); out[f"n_{qq}"] = len(s); out[f"ev_{qq}"] = int(s.sum())
    return out
rows = []
for y, lab in LOCAL_OUT:
    if y in ("SAP", "mctsi_high"):
        continue
    df = loc
    for sc in ["risk6", "risk3"]:
        rows.append(cont(df, y, sc, "Local", lab))
for y, lab in MIMIC_OUT:
    d_ = mim_icd.copy(); d_[y] = d_[y].astype(int)
    for sc in ["risk6", "risk3"]:
        rows.append(cont(d_, y, sc, "MIMIC-IV (ICD-10, necrosis specified)", lab))
for y, lab in MIMIC_FULL:
    if y == "SAP_proxy":
        continue
    for sc in ["risk6", "risk3"]:
        rows.append(cont(mim, y, sc, "MIMIC-IV", lab))
for y, lab in EICU_OUT:
    for sc in ["risk6", "risk3"]:
        rows.append(cont(eic, y, sc, "eICU", lab))
rc = pd.DataFrame(rows); rc.to_csv(os.path.join(RES, "risk_continuous.csv"), index=False)
print("\n=== 连续风险 每+1SD OR ===")
print(rc[["cohort", "outcome", "score", "OR_perSD", "lo", "hi", "p", "auc", "rate_Q1", "rate_Q4"]].round(3).to_string(index=False))

# ---------- 敏感性 ----------
rows = []
l2 = loc[~loc.allimp]
for y, lab in [("necroCT", "CT necrosis"), ("intervention", "Drainage/necrosectomy"), ("LOS14", "LOS > 14 d")]:
    rows += rates(l2, y, "Local excl. fully-imputed (n=1273)", lab)
rows += rates(loc[loc.cohort == "dev"], "intervention", "Local derivation 2017-2020", "Drainage/necrosectomy")
rows += rates(loc[loc.cohort == "val"], "intervention", "Local temporal 2021", "Drainage/necrosectomy")
rows += rates(loc[loc.cohort == "dev"], "necroCT", "Local derivation 2017-2020", "CT necrosis")
rows += rates(loc[loc.cohort == "val"], "necroCT", "Local temporal 2021", "CT necrosis")
d_ = mim_icd10.copy(); d_["necro_coded0"] = d_.necro_coded0.astype(int); d_["ipn_coded0"] = d_.ipn_coded0.astype(int)
rows += rates(d_, "necro_coded0", "MIMIC-IV ICD-10 all (unspecified = no, n=1963)", "Coded necrosis")
rows += rates(d_, "ipn_coded0", "MIMIC-IV ICD-10 all (unspecified = no, n=1963)", "Coded infected necrosis")
cc = mim_icd[mim_icd.complete6 == 1].copy()
for y, lab in MIMIC_OUT:
    cc[y] = cc[y].astype(int); rows += rates(cc, y, "MIMIC-IV complete cases (ICD-10, necrosis specified)", lab)
rows += rates(mim[mim.complete6 == 1], "epi_any", "MIMIC-IV complete cases", "Extra-pancreatic infection")
rows += rates(mim[mim.transfer_in == 0], "epi_any", "MIMIC-IV non-transfer", "Extra-pancreatic infection")
rows += rates(mim_icd[mim_icd.transfer_in == 0].assign(necro_coded=lambda d: d.necro_coded.astype(int)), "necro_coded", "MIMIC-IV non-transfer (ICD-10, necrosis specified)", "Coded necrosis")
rows += rates(eic[eic.complete6 == 1], "epi_any", "eICU complete cases", "Extra-pancreatic infection")
rows += rates(eic[eic.complete6 == 1], "necro_dx", "eICU complete cases", "Necrotizing (diagnosis string)")
ss = pd.DataFrame(rows); ss.to_csv(os.path.join(RES, "rates_sens.csv"), index=False)
print("\n=== 敏感性 ===")
ssp = ss.assign(txt=lambda d: (d.rate * 100).round(1).astype(str) + "% (" + d.events.astype(str) + "/" + d.n.astype(str) + ")")
print(ssp.pivot_table(index=["cohort", "outcome"], columns="phenotype", values="txt", aggfunc="first").to_string())
print(ss.groupby(["cohort", "outcome"]).trend_p.first().map(lambda p: f"{p:.2e}").to_string())
print("\n已落盘: rates_by_pheno / adj_or / beyond_severity_lr / rates_nonsevere / risk_continuous / rates_sens")
