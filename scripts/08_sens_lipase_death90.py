# -*- coding: utf-8 -*-
"""08 MIMIC 补充敏感性: (a) K85 + 脂肪酶 >=3xULN (>=180 U/L, 入院48h内最大值) 子集重复坏死/感染梯度
   (b) 30/90/365 天死亡 (patients.dod) 按表型梯度 + 年龄性别调整 OR
   输出: results/rates_sens.csv 追加行 (标记 source=08), results/death_longterm.csv"""
import pandas as pd
import numpy as np
import os, sys, warnings
import statsmodels.api as sm
from scipy.stats import chi2_contingency
warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = lambda var, sub: os.environ.get(var, os.path.join(ROOT, "inputs", sub))
P = ROOT
RES, DAT = os.path.join(P, "results"), os.path.join(P, "data")
mim = pd.read_csv(os.path.join(DAT, "mimic.csv"))
lp = pd.read_csv(os.path.join(DAT, "mimic_lipase_dod.csv"))
_lp_cols = ["Lipase", "Lipase_max48", "death30", "death90", "death365", "dod_days"]
mim = mim.drop(columns=[c for c in _lp_cols if c in mim.columns])  # idempotent re-runs
mim = mim.merge(lp[["hadm_id"] + _lp_cols], on="hadm_id", how="left")
mim.to_csv(os.path.join(DAT, "mimic.csv"), index=False)

def wilson(k, n, z=1.96):
    p = k / n; den = 1 + z**2 / n; c = (p + z**2 / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
    return c - h, c + h
def rates(df, y, cohort, outcome):
    rows = []
    tab = pd.crosstab(df.pheno, df[y]); chi_p = chi2_contingency(tab).pvalue if tab.shape == (4, 2) else np.nan
    tr = sm.Logit(df[y].astype(int).values, sm.add_constant(df.pheno.values.astype(float))).fit(disp=0)
    for ph in [1, 2, 3, 4]:
        s = df[df.pheno == ph]; k, n = int(s[y].sum()), len(s); lo, hi = wilson(k, n)
        rows.append(dict(cohort=cohort, outcome=outcome, phenotype=f"P{ph}", n=n, events=k, rate=k / n, lo=lo, hi=hi,
                         chi2_p=chi_p, trend_p=tr.pvalues[1]))
    return rows
lip = mim[(mim.icd10_ap == 1) & (mim.necro_specified == 1) & (mim.Lipase_max48 >= 180)].copy()
rows = []
for y, lab in [("necro_coded", "Coded necrosis"), ("ipn_coded", "Coded infected necrosis")]:
    lip[y] = lip[y].astype(int); rows += rates(lip, y, "MIMIC-IV ICD-10 + lipase >=3xULN (necrosis specified)", lab)
lip2 = mim[mim.Lipase_max48 >= 180].copy()
for y, lab in [("epi_any", "Extra-pancreatic infection"), ("sepsis", "Sepsis code"), ("LOS14", "LOS > 14 d")]:
    rows += rates(lip2, y, "MIMIC-IV lipase >=3xULN", lab)
ss = pd.read_csv(os.path.join(RES, "rates_sens.csv"))
ss = ss[~ss.cohort.str.contains("lipase")]
ss = pd.concat([ss, pd.DataFrame(rows)], ignore_index=True); ss.to_csv(os.path.join(RES, "rates_sens.csv"), index=False)
print("脂肪酶子集 n =", len(lip), "/", len(lip2))
print(pd.DataFrame(rows).assign(txt=lambda d: (d.rate * 100).round(1).astype(str) + "% (" + d.events.astype(str) + "/" + d.n.astype(str) + ")")
      .pivot_table(index=["cohort", "outcome"], columns="phenotype", values="txt", aggfunc="first").to_string())

# 长期死亡
rows = []
for y, lab in [("death", "In-hospital death"), ("death30", "30-day death"), ("death90", "90-day death"), ("death365", "1-year death")]:
    rows += rates(mim, y, "MIMIC-IV", lab)
dl = pd.DataFrame(rows)
ors = []
for y, lab in [("death30", "30-day death"), ("death90", "90-day death"), ("death365", "1-year death")]:
    X = pd.get_dummies(mim.pheno.astype(int), prefix="P").astype(float).iloc[:, 1:]
    X["Age"] = mim.Age_imp; X["Sex"] = mim.Sex_male; X["SAP_proxy"] = mim.SAP_proxy
    m = sm.Logit(mim[y].values, sm.add_constant(X)).fit(disp=0)
    for cn in ["P_2", "P_3", "P_4"]:
        b, se = m.params[cn], m.bse[cn]
        ors.append(dict(outcome=lab, pheno=cn.replace("P_", "P"), OR=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=m.pvalues[cn], model="age_sex_SAPproxy"))
dl.to_csv(os.path.join(RES, "death_longterm.csv"), index=False)
pd.DataFrame(ors).to_csv(os.path.join(RES, "death_longterm_or.csv"), index=False)
print(dl.assign(txt=lambda d: (d.rate * 100).round(1).astype(str) + "% (" + d.events.astype(str) + "/" + d.n.astype(str) + ")")
      .pivot_table(index="outcome", columns="phenotype", values="txt", aggfunc="first").to_string())
print(pd.DataFrame(ors).round(3).to_string(index=False))
print("已落盘: rates_sens.csv(追加) / death_longterm.csv / death_longterm_or.csv; data/mimic.csv 已并入脂肪酶与 dod")
