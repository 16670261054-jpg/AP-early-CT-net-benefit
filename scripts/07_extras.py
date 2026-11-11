# -*- coding: utf-8 -*-
"""07 文献建议并入的增补分析
   a E-value: P4 vs P1 调整 OR (age_sex_etio_SAP) -> 常见结局用 sqrt(OR) 近似 RR 后算 E-value
   b P4 (vs 非P4) 对各结局的 敏感度/特异度/LR+/LR- (Fagan 素材) + 验前/验后概率
   c 亚组一致性: P4 vs P1 (或表型序数每+1) 在 病因/年龄/性别/早期SAP 分层内的 OR (本院干预+CT坏死; MIMIC 编码坏死+胰外感染)
   d 院内死亡累积发生 (Aalen-Johansen, 出院存活为竞争事件) 按表型 (MIMIC/eICU) + 28d 时点 CIF
   e MIMIC 诊断序号: 胰外感染码 vs 感染坏死码 的 seq_num 比较 (顺序代理)
   f 低风险层内 CT 规则 (MCTSI>=6) 的 DCA
   g 表型分类器一致性 (models/pheno_clf_xgb.pkl) 在本院数据
   输出 results/: evalue.csv, lr_fagan.csv, subgroup_forest.csv, cif_death.csv, cif_death_28d.csv,
                  dx_seq_order.csv, dca_lowrisk_ct.csv
"""
import pandas as pd
import numpy as np
import os, sys, warnings, joblib
import statsmodels.api as sm
warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = lambda var, sub: os.environ.get(var, os.path.join(ROOT, "inputs", sub))
P = ROOT
RES, DAT = os.path.join(P, "results"), os.path.join(P, "data")
loc = pd.read_csv(os.path.join(DAT, "local.csv")); loc["Sex"] = loc.Sex_male
mim = pd.read_csv(os.path.join(DAT, "mimic.csv")); mim["Sex"] = mim.Sex_male; mim["Age"] = mim.Age_imp
eic = pd.read_csv(os.path.join(DAT, "eicu.csv")); eic["Sex"] = eic.Sex_male; eic["Age"] = eic.Age_imp
mim_icd = mim[mim.necro_specified == 1].copy()

# ---------- a E-value ----------
def evalue(rr):
    if pd.isna(rr):
        return np.nan
    if rr < 1:
        rr = 1 / rr
    return rr + np.sqrt(rr * (rr - 1)) if rr > 1 else 1.0
orf = pd.read_csv(os.path.join(RES, "adj_or.csv"))
rt = pd.read_csv(os.path.join(RES, "rates_by_pheno.csv"))
prev = rt.groupby(["cohort", "outcome"]).apply(lambda d: d.events.sum() / d.n.sum()).rename("prev").reset_index()
sel = orf[(orf.pheno == "P4") & orf.model.isin(["age_sex_etio_SAP", "age_sex_apache"]) & (orf.separation == 0)].merge(prev, on=["cohort", "outcome"], how="left")
rows = []
for _, r in sel.iterrows():
    common = r.prev > 0.15
    rr = np.sqrt(r.OR) if common else r.OR
    # CI 的 E 值用靠近 1 的界: OR>1 用 lo, OR<1 用 hi; CI 跨 1 则 E=1
    ci_bound = r.lo if r.lo > 1 else (r.hi if r.hi < 1 else np.nan)
    ci_rr = (np.sqrt(ci_bound) if common else ci_bound) if not np.isnan(ci_bound) else np.nan
    rows.append(dict(cohort=r.cohort, outcome=r.outcome, OR=r.OR, lo=r.lo, hi=r.hi, prevalence=r.prev,
                     rr_approx=rr, rr_method="sqrt(OR)" if common else "OR~RR",
                     evalue_point=evalue(rr), evalue_ci=evalue(ci_rr) if not np.isnan(ci_rr) else 1.0))
ev = pd.DataFrame(rows); ev.to_csv(os.path.join(RES, "evalue.csv"), index=False)
print("=== E-value ===")
print(ev.round(3).to_string(index=False))

# ---------- b LR+/LR- (P4 vs others) ----------
def lr_stats(df, y, cohort, outcome):
    d = df[df[y].notna()]; yv = d[y].astype(int).values; t = (d.pheno == 4).astype(int).values
    tp = ((t == 1) & (yv == 1)).sum(); fp = ((t == 1) & (yv == 0)).sum(); fn = ((t == 0) & (yv == 1)).sum(); tn = ((t == 0) & (yv == 0)).sum()
    se = tp / (tp + fn); sp = tn / (tn + fp)
    lrp = se / (1 - sp) if sp < 1 else np.inf; lrn = (1 - se) / sp if sp > 0 else np.inf
    pre = yv.mean(); odds = pre / (1 - pre)
    post_pos = odds * lrp / (1 + odds * lrp); post_neg = odds * lrn / (1 + odds * lrn)
    return dict(cohort=cohort, outcome=outcome, n=len(yv), events=int(yv.sum()), P4_n=int(t.sum()), sens=se, spec=sp,
                LR_pos=lrp, LR_neg=lrn, pretest=pre, posttest_P4=post_pos, posttest_nonP4=post_neg, PPV=tp / (tp + fp), NPV=tn / (tn + fn))
rows = []
for y, lab in [("necroCT", "CT necrosis"), ("intervention", "Drainage/necrosectomy"),
               ("persOF", "Persistent OF"), ("LOS14", "LOS > 14 d")]:
    rows.append(lr_stats(loc, y, "Local", lab))
for y, lab in [("necro_coded", "Coded necrosis"), ("ipn_coded", "Coded infected necrosis")]:
    rows.append(lr_stats(mim_icd, y, "MIMIC-IV (ICD-10, necrosis specified)", lab))
for y, lab in [("sepsis", "Sepsis code"), ("epi_any", "Extra-pancreatic infection"), ("LOS14", "LOS > 14 d"), ("death", "In-hospital death")]:
    rows.append(lr_stats(mim, y, "MIMIC-IV", lab))
for y, lab in [("necro_dx", "Necrotizing (diagnosis string)"), ("LOS14", "LOS > 14 d"), ("death", "In-hospital death")]:
    rows.append(lr_stats(eic, y, "eICU", lab))
lrf = pd.DataFrame(rows); lrf.to_csv(os.path.join(RES, "lr_fagan.csv"), index=False)
print("\n=== P4 LR ===")
print(lrf.round(3).to_string(index=False))

# ---------- c 亚组一致性 (表型序数每 +1 级 的 OR, 分层内; 及 P4 vs P1) ----------
def sub_or(df, y, strata, cohort, outcome):
    out = []
    d0 = df[df[y].notna()].copy()
    for svar, slab, groups in strata:
        for gname, mask in groups(d0):
            d = d0[mask]
            yv = d[y].astype(int).values
            if len(d) < 40 or yv.sum() < 5:
                out.append(dict(cohort=cohort, outcome=outcome, subgroup_var=slab, subgroup=gname, n=len(d), events=int(yv.sum()),
                                OR_ordinal=np.nan, lo=np.nan, hi=np.nan, p=np.nan)); continue
            X = sm.add_constant(pd.DataFrame({"ph": d.pheno.astype(float).values, "Age": d.Age.astype(float).values, "Sex": d.Sex.astype(float).values}))
            if svar == "Age":
                X = X.drop(columns="Age")
            if svar == "Sex":
                X = X.drop(columns="Sex")
            try:
                m = sm.Logit(yv, X).fit(disp=0, maxiter=200, method="bfgs")
                b, se = m.params["ph"], m.bse["ph"]
                out.append(dict(cohort=cohort, outcome=outcome, subgroup_var=slab, subgroup=gname, n=len(d), events=int(yv.sum()),
                                OR_ordinal=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=m.pvalues["ph"]))
            except Exception:
                out.append(dict(cohort=cohort, outcome=outcome, subgroup_var=slab, subgroup=gname, n=len(d), events=int(yv.sum()),
                                OR_ordinal=np.nan, lo=np.nan, hi=np.nan, p=np.nan))
    # 交互 p (表型序数 x 分层) 仅对二分层
    return out

def strata_local(d):
    return [("Age", "Age", lambda d: [("< 60 y", d.Age < 60), (">= 60 y", d.Age >= 60)]),
            ("Sex", "Sex", lambda d: [("Male", d.Sex == 1), ("Female", d.Sex == 0)]),
            ("etiology", "Etiology", lambda d: [(e, d.etiology == e) for e in ["Biliary", "Alcoholic", "Hypertriglyceridemic", "Other"]]),
            ("SAP", "Early severity", lambda d: [("Non-SAP", d.SAP == 0), ("SAP", d.SAP == 1)]),
            ("cohort", "Period", lambda d: [("2017-2020", d.cohort == "dev"), ("2021", d.cohort == "val")]),
            ("BMI", "BMI", lambda d: [("< 28", d.BMI < 28), (">= 28", d.BMI >= 28)])]
def strata_mimic(d):
    return [("Age", "Age", lambda d: [("< 60 y", d.Age < 60), (">= 60 y", d.Age >= 60)]),
            ("Sex", "Sex", lambda d: [("Male", d.Sex == 1), ("Female", d.Sex == 0)]),
            ("etiology", "Etiology", lambda d: [(e, d.etiology == e) for e in ["Biliary", "Alcoholic", "Unspecified", "Other"]]),
            ("SAP_proxy", "Early severity", lambda d: [("Non-SAP proxy", d.SAP_proxy == 0), ("SAP proxy", d.SAP_proxy == 1)]),
            ("transfer_in", "Admission source", lambda d: [("Direct", d.transfer_in == 0), ("Transfer", d.transfer_in == 1)]),
            ("complete6", "Six variables", lambda d: [("Complete", d.complete6 == 1), ("Imputed", d.complete6 == 0)])]
rows = []
for y, lab in [("intervention", "Drainage/necrosectomy"), ("necroCT", "CT necrosis"), ("LOS14", "LOS > 14 d")]:
    rows += sub_or(loc, y, strata_local(loc), "Local", lab)
for y, lab in [("necro_coded", "Coded necrosis")]:
    rows += sub_or(mim_icd, y, strata_mimic(mim_icd), "MIMIC-IV (ICD-10, necrosis specified)", lab)
for y, lab in [("epi_any", "Extra-pancreatic infection"), ("sepsis", "Sepsis code"), ("LOS14", "LOS > 14 d")]:
    rows += sub_or(mim, y, strata_mimic(mim), "MIMIC-IV", lab)
sf = pd.DataFrame(rows); sf.to_csv(os.path.join(RES, "subgroup_forest.csv"), index=False)
print("\n=== 亚组内 表型序数每+1 OR ===")
print(sf.round(3).to_string(index=False))

# ---------- d 累积发生 (Aalen-Johansen) ----------
def aj_cif(t, ev):
    """ev: 1=死亡(目标), 2=出院存活(竞争), 0=删失. 返回 DataFrame(time, cif_death, cif_disch)"""
    df = pd.DataFrame({"t": t, "e": ev}).sort_values("t")
    times = np.unique(df.t[df.e > 0])
    n_at = len(df); S = 1.0; cif1 = 0.0; cif2 = 0.0; out = [(0.0, 0.0, 0.0)]
    for tt in times:
        d1 = ((df.t == tt) & (df.e == 1)).sum(); d2 = ((df.t == tt) & (df.e == 2)).sum()
        n_at = (df.t >= tt).sum()
        if n_at == 0:
            break
        cif1 += S * d1 / n_at; cif2 += S * d2 / n_at
        S *= (1 - (d1 + d2) / n_at)
        out.append((tt, cif1, cif2))
    return pd.DataFrame(out, columns=["time", "cif_death", "cif_discharge"])
rows, pts = [], []
for name, d_ in [("MIMIC-IV", mim), ("eICU", eic)]:
    if name == "MIMIC-IV":
        t = np.where(d_.death == 1, d_.death_day, d_.los); e = np.where(d_.death == 1, 1, 2)
    else:
        t = d_.los.values; e = np.where(d_.death == 1, 1, 2)
    ok = np.isfinite(t) & (t > 0)
    for ph in [1, 2, 3, 4]:
        m = ok & (d_.pheno.values == ph)
        c = aj_cif(np.minimum(t[m], 60), np.where(t[m] > 60, 0, e[m]))
        c["cohort"] = name; c["phenotype"] = f"P{ph}"; c["n"] = int(m.sum())
        rows.append(c)
        for day in [7, 14, 28]:
            cc = c[c.time <= day]
            pts.append(dict(cohort=name, phenotype=f"P{ph}", day=day, cif_death=cc.cif_death.iloc[-1], cif_discharge=cc.cif_discharge.iloc[-1], n=int(m.sum())))
cif = pd.concat(rows); cif.to_csv(os.path.join(RES, "cif_death.csv"), index=False)
cp = pd.DataFrame(pts); cp.to_csv(os.path.join(RES, "cif_death_28d.csv"), index=False)
print("\n=== 累积死亡发生 (竞争: 出院) ===")
print(cp.pivot_table(index=["cohort", "phenotype"], columns="day", values="cif_death").round(4).to_string())

# ---------- e 诊断序号 ----------
M = IN("MIMIC_HOSP", "mimic-iv-3.1/hosp")
dx = pd.read_csv(os.path.join(M, "diagnoses_icd.csv.gz"), usecols=["hadm_id", "seq_num", "icd_code", "icd_version"])
dx = dx[dx.hadm_id.isin(set(mim.hadm_id))].copy(); dx["icd_code"] = dx.icd_code.str.strip()
ipn_codes = dx[(dx.icd_version == 10) & dx.icd_code.str.match(r"^K85\d2$")]
epi_codes = dx[dx.icd_code.str.startswith(("J13", "J14", "J15", "J16", "J17", "J18", "J690", "R7881", "N390", "K830", "A047")) & (dx.icd_version == 10)]
both = set(ipn_codes.hadm_id) & set(epi_codes.hadm_id)
rows = []
for h in both:
    rows.append(dict(hadm_id=h, seq_ipn=ipn_codes[ipn_codes.hadm_id == h].seq_num.min(), seq_epi=epi_codes[epi_codes.hadm_id == h].seq_num.min()))
so = pd.DataFrame(rows)
summ = pd.DataFrame([dict(n_both=len(so), epi_ranked_before_ipn=int((so.seq_epi < so.seq_ipn).sum()) if len(so) else 0,
                          ipn_ranked_before_epi=int((so.seq_ipn < so.seq_epi).sum()) if len(so) else 0,
                          median_seq_ipn=so.seq_ipn.median() if len(so) else np.nan, median_seq_epi=so.seq_epi.median() if len(so) else np.nan,
                          median_seq_epi_all=epi_codes.groupby("hadm_id").seq_num.min().median(),
                          median_seq_ipn_all=ipn_codes.groupby("hadm_id").seq_num.min().median(),
                          note="seq_num = billing rank, not chronology; descriptive only")])
summ.to_csv(os.path.join(RES, "dx_seq_order.csv"), index=False)
print("\n=== 诊断序号 ===\n", summ.to_string(index=False))

# ---------- f 低风险层 CT 规则 DCA ----------
rows = []
for stratum, mask in [("Low-risk (P1/P2, non-SAP)", loc.pheno.isin([1, 2]) & (loc.SAP == 0)), ("P4 non-SAP", (loc.pheno == 4) & (loc.SAP == 0)),
                      ("SAP", loc.SAP == 1), ("All", loc.pheno > 0)]:
    s = loc[mask]
    for y, lab in [("intervention", "Drainage/necrosectomy"), ("persOF", "Persistent OF")]:
        yv = s[y].astype(int).values; rule = (s.MCTSI >= 6).astype(int).values; prev = yv.mean()
        for t in np.arange(0.02, 0.41, 0.01):
            tp = np.sum((rule == 1) & (yv == 1)) / len(yv); fp = np.sum((rule == 1) & (yv == 0)) / len(yv)
            rows.append(dict(stratum=stratum, outcome=lab, threshold=t, n=len(s), events=int(yv.sum()), prevalence=prev,
                             nb_ct_rule=tp - fp * t / (1 - t), nb_treat_all=prev - (1 - prev) * t / (1 - t), nb_none=0.0))
pd.DataFrame(rows).to_csv(os.path.join(RES, "dca_lowrisk_ct.csv"), index=False)

# ---------- g 分类器一致性 ----------
clf = joblib.load(os.path.join(ROOT, "models", "pheno_clf_xgb.pkl"))
R = ["Ca", "LYMPH_pct", "Age", "RBC", "SCr", "ALT"]
pr = clf.predict(loc[R].values) + 1
agree = pd.DataFrame([dict(cohort="Local (all)", n=len(loc), agreement=(pr == loc.pheno.values).mean())])
agree.to_csv(os.path.join(RES, "clf_agreement_local.csv"), index=False)
print("\n分类器一致性:", agree.to_string(index=False))
print("\n已落盘: evalue / lr_fagan / subgroup_forest / cif_death(+28d) / dx_seq_order / dca_lowrisk_ct / clf_agreement_local")
