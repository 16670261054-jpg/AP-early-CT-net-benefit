# -*- coding: utf-8 -*-
"""06 次要结局与描述：
   a 死亡时相 (MIMIC: <=7 / 8-14 / >14 d) 按表型 + 感染坏死/胰外感染死亡者的时相
   b 器官支持时长分布 (MIMIC 313 例) 按表型: 一过性(<48h) vs 持续(>=48h)
   c LOS 分布 (三库) 按表型: 中位数/IQR + 每患者长表 (雨云图)
   d 胰外感染 vs 感染坏死: MIMIC 各表型内二者频数与共现; 各感染部位死亡率
   e 表型 x MCTSI x 结局 三维表 (本院)
   f Table1: 三库基线按表型 (中位数 IQR / n %)
   g 病因分布按表型 (本院 + MIMIC ICD-10)
   输出 results/: death_phase.csv, support_duration.csv, los_long.csv, los_summary.csv,
                  epi_vs_ipn.csv, infection_site_mortality.csv, pheno_mctsi_outcome.csv,
                  table1_local.csv, table1_mimic.csv, table1_eicu.csv, etiology_by_pheno.csv
"""
import pandas as pd
import numpy as np
import os, sys, warnings
from scipy.stats import kruskal, chi2_contingency
warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = lambda var, sub: os.environ.get(var, os.path.join(ROOT, "inputs", sub))
P = ROOT
RES, DAT = os.path.join(P, "results"), os.path.join(P, "data")
loc = pd.read_csv(os.path.join(DAT, "local.csv"))
mim = pd.read_csv(os.path.join(DAT, "mimic.csv"))
eic = pd.read_csv(os.path.join(DAT, "eicu.csv"))

# ---------- a 死亡时相 ----------
dd = mim[mim.death == 1].copy()
rows = []
for ph in [1, 2, 3, 4]:
    s = dd[dd.pheno == ph]
    for phase in ["<=7d", "8-14d", ">14d"]:
        rows.append(dict(cohort="MIMIC-IV", phenotype=f"P{ph}", phase=phase, deaths=int((s.death_phase == phase).sum()),
                         total_deaths=len(s), n_pheno=int((mim.pheno == ph).sum())))
for grp, mask in [("Coded infected necrosis", dd.ipn_coded == 1), ("Coded necrosis (any)", dd.necro_coded == 1),
                  ("Extra-pancreatic infection", dd.epi_any == 1), ("Sepsis code", dd.sepsis == 1),
                  ("No infection code", (dd.epi_any == 0) & (dd.sepsis == 0) & (dd.ipn_coded != 1))]:
    s = dd[mask]
    for phase in ["<=7d", "8-14d", ">14d"]:
        rows.append(dict(cohort="MIMIC-IV", phenotype=grp, phase=phase, deaths=int((s.death_phase == phase).sum()),
                         total_deaths=len(s), n_pheno=np.nan))
dp = pd.DataFrame(rows); dp.to_csv(os.path.join(RES, "death_phase.csv"), index=False)
print("=== 死亡时相 ===")
print(dp.pivot_table(index="phenotype", columns="phase", values="deaths", aggfunc="first").to_string())
print("死亡日中位数:", dd.death_day.median().round(1), " IQR:", dd.death_day.quantile([.25, .75]).round(1).tolist())

# ---------- b 器官支持时长 ----------
sp = mim[mim.support_any == 1].copy()
sp["max_h"] = sp[["vent_h", "vaso_h"]].max(axis=1)
sp["category"] = np.where(sp.max_h >= 48, "Persistent (>=48 h)", "Transient (<48 h)")
rows = []
for ph in [1, 2, 3, 4]:
    s = sp[sp.pheno == ph]; n_ph = int((mim.pheno == ph).sum())
    for cat in ["Transient (<48 h)", "Persistent (>=48 h)"]:
        ss = s[s.category == cat]
        rows.append(dict(phenotype=f"P{ph}", category=cat, n=len(ss), n_pheno=n_ph, pct_of_pheno=len(ss) / n_ph,
                         deaths=int(ss.death.sum()), median_h=ss.max_h.median() if len(ss) else np.nan))
su = pd.DataFrame(rows); su.to_csv(os.path.join(RES, "support_duration.csv"), index=False)
sp[["hadm_id", "pheno", "vent_h", "vaso_h", "max_h", "category", "death"]].to_csv(os.path.join(RES, "support_long.csv"), index=False)
print("\n=== 器官支持 ===")
print(su.round(3).to_string(index=False))

# ---------- c LOS ----------
lo = pd.concat([
    pd.DataFrame({"cohort": "Local", "phenotype": "P" + loc.pheno.astype(str), "los": loc.LOS}),
    pd.DataFrame({"cohort": "MIMIC-IV", "phenotype": "P" + mim.pheno.astype(str), "los": mim.los}),
    pd.DataFrame({"cohort": "eICU", "phenotype": "P" + eic.pheno.astype(str), "los": eic.los})])
lo = lo[lo.los.notna() & (lo.los > 0)]
lo.to_csv(os.path.join(RES, "los_long.csv"), index=False)
rows = []
for co, s in lo.groupby("cohort"):
    kw = kruskal(*[g.los.values for _, g in s.groupby("phenotype")]).pvalue
    for ph, g in s.groupby("phenotype"):
        rows.append(dict(cohort=co, phenotype=ph, n=len(g), median=g.los.median(), q1=g.los.quantile(.25),
                         q3=g.los.quantile(.75), pct_gt14=(g.los > 14).mean(), kruskal_p=kw))
ls = pd.DataFrame(rows); ls.to_csv(os.path.join(RES, "los_summary.csv"), index=False)
print("\n=== LOS ===")
print(ls.round(3).to_string(index=False))

# ---------- d 胰外感染 vs 感染坏死 ----------
mi = mim[mim.necro_specified == 1]
rows = []
for ph in [1, 2, 3, 4]:
    s = mi[mi.pheno == ph]
    rows.append(dict(phenotype=f"P{ph}", n=len(s), ipn=int(s.ipn_coded.sum()), epi=int(s.epi_any.sum()),
                     both=int(((s.ipn_coded == 1) & (s.epi_any == 1)).sum()),
                     epi_only=int(((s.ipn_coded == 0) & (s.epi_any == 1)).sum()),
                     ipn_only=int(((s.ipn_coded == 1) & (s.epi_any == 0)).sum()),
                     sepsis=int(s.sepsis.sum()), ratio_epi_to_ipn=s.epi_any.sum() / max(s.ipn_coded.sum(), 1)))
rows.append(dict(phenotype="All", n=len(mi), ipn=int(mi.ipn_coded.sum()), epi=int(mi.epi_any.sum()),
                 both=int(((mi.ipn_coded == 1) & (mi.epi_any == 1)).sum()),
                 epi_only=int(((mi.ipn_coded == 0) & (mi.epi_any == 1)).sum()),
                 ipn_only=int(((mi.ipn_coded == 1) & (mi.epi_any == 0)).sum()),
                 sepsis=int(mi.sepsis.sum()), ratio_epi_to_ipn=mi.epi_any.sum() / max(mi.ipn_coded.sum(), 1)))
ev = pd.DataFrame(rows); ev.to_csv(os.path.join(RES, "epi_vs_ipn.csv"), index=False)
print("\n=== 胰外感染 vs 感染坏死 (MIMIC ICD-10 坏死已指明) ===")
print(ev.round(2).to_string(index=False))
rows = []
for co, df in [("MIMIC-IV", mim), ("eICU", eic)]:
    for site, col in [("Any extra-pancreatic infection", "epi_any"), ("Pneumonia", "pneumonia"), ("Bacteremia", "bacteremia"),
                      ("UTI", "uti"), ("Cholangitis", "cholangitis"), ("C. difficile", "cdiff"),
                      ("Coded infected necrosis", "ipn_coded" if co == "MIMIC-IV" else "abscess_dx"),
                      ("Sepsis", "sepsis" if co == "MIMIC-IV" else "sepsis_dx")]:
        s1 = df[df[col] == 1]; s0 = df[df[col] == 0]
        rows.append(dict(cohort=co, site=site, n=len(s1), deaths=int(s1.death.sum()), mort=s1.death.mean() if len(s1) else np.nan,
                         mort_without=s0.death.mean(), los_median=s1.los.median() if len(s1) else np.nan,
                         los_median_without=s0.los.median(),
                         p4_share=(s1.pheno == 4).mean() if len(s1) else np.nan))
ism = pd.DataFrame(rows); ism.to_csv(os.path.join(RES, "infection_site_mortality.csv"), index=False)
print("\n=== 感染部位死亡率 ===")
print(ism.round(3).to_string(index=False))

# ---------- e 表型 x MCTSI x 结局 ----------
loc["mctsi_cat"] = pd.cut(loc.MCTSI, [-1, 2, 4, 6, 10], labels=["0-2", "4", "6", "8-10"]).astype(str)
rows = []
for (ph, mc), s in loc.groupby(["pheno", "mctsi_cat"]):
    rows.append(dict(phenotype=f"P{ph}", mctsi=mc, n=len(s), necroCT=int(s.necroCT.sum()), intervention=int(s.intervention.sum()),
                     intervention_rate=s.intervention.mean(), persOF=int(s.persOF.sum()), persOF_rate=s.persOF.mean(),
                     LOS14_rate=s.LOS14.mean(),
                     death=int(s.death.sum())))
pm = pd.DataFrame(rows); pm.to_csv(os.path.join(RES, "pheno_mctsi_outcome.csv"), index=False)
print("\n=== 表型 x MCTSI: 干预率 ===")
print(pm.pivot(index="phenotype", columns="mctsi", values="intervention_rate").round(3).to_string())
print(pm.pivot(index="phenotype", columns="mctsi", values="n").to_string())

# ---------- f Table 1 ----------
def t1(df, cont, cat, cohort):
    rows = []
    for v, lab in cont:
        s = df[v]
        kw = kruskal(*[g.dropna().values for _, g in s.groupby(df.pheno)]).pvalue
        r = dict(cohort=cohort, variable=lab, type="cont", p=kw, missing=int(s.isna().sum()))
        for ph in [1, 2, 3, 4]:
            g = s[df.pheno == ph].dropna()
            r[f"P{ph}"] = f"{g.median():.1f} [{g.quantile(.25):.1f}-{g.quantile(.75):.1f}]" if len(g) else ""
        r["All"] = f"{s.median():.1f} [{s.quantile(.25):.1f}-{s.quantile(.75):.1f}]"
        rows.append(r)
    for v, lab in cat:
        s = df[v]
        ok = s.notna()
        try:
            p = chi2_contingency(pd.crosstab(df.pheno[ok], s[ok])).pvalue
        except Exception:
            p = np.nan
        r = dict(cohort=cohort, variable=lab, type="cat", p=p, missing=int(s.isna().sum()))
        for ph in [1, 2, 3, 4]:
            g = s[(df.pheno == ph) & ok]
            r[f"P{ph}"] = f"{int(g.sum())} ({g.mean()*100:.1f})" if len(g) else ""
        r["All"] = f"{int(s[ok].sum())} ({s[ok].mean()*100:.1f})"
        rows.append(r)
    r = dict(cohort=cohort, variable="n", type="n", p=np.nan, missing=0, All=str(len(df)))
    for ph in [1, 2, 3, 4]:
        r[f"P{ph}"] = str(int((df.pheno == ph).sum()))
    return pd.DataFrame([r] + rows)

loc["Biliary"] = (loc.etiology == "Biliary").astype(int); loc["Alcoholic"] = (loc.etiology == "Alcoholic").astype(int)
loc["HTG"] = (loc.etiology == "Hypertriglyceridemic").astype(int)
t1l = t1(loc, [("Age", "Age, years"), ("Ca", "Calcium, mmol/L"), ("LYMPH_pct", "Lymphocyte, %"), ("RBC", "RBC, 10^12/L"),
               ("SCr", "Creatinine, umol/L"), ("ALT", "ALT, U/L"), ("Albumin", "Albumin, g/L"), ("BUN", "BUN, mmol/L"),
               ("WBC", "WBC, 10^9/L"), ("BMI", "BMI, kg/m2"), ("BISAP", "BISAP"), ("MCTSI", "MCTSI (admission CT)"), ("LOS", "Length of stay, d")],
          [("Sex_male", "Male"), ("Biliary", "Biliary etiology"), ("Alcoholic", "Alcoholic etiology"), ("HTG", "Hypertriglyceridemic etiology"),
           ("SAP", "SAP (recorded)"), ("MSAP_or_SAP", "MSAP or SAP"), ("mctsi_high", "MCTSI >= 6"), ("necroCT", "CT necrosis"),
           ("apfc", "Peripancreatic collection"), ("pcd", "Percutaneous drainage"),
           ("necrosectomy", "Necrosectomy/debridement"), ("intervention", "Any local intervention"), ("persOF", "Persistent OF"),
           ("icu", "ICU admission"), ("LOS14", "LOS > 14 d"), ("death", "In-hospital death")], "Local")
t1l.to_csv(os.path.join(RES, "table1_local.csv"), index=False)
mim["Biliary"] = (mim.etiology == "Biliary").astype(int); mim["Alcoholic"] = (mim.etiology == "Alcoholic").astype(int)
mim["ICD10"] = mim.icd10_ap
t1m = t1(mim, [("Age", "Age, years"), ("Ca", "Calcium, mmol/L"), ("LYMPH_pct", "Lymphocyte, %"), ("RBC", "RBC, 10^12/L"),
               ("SCr", "Creatinine, umol/L"), ("ALT", "ALT, U/L"), ("Albumin", "Albumin, g/L"), ("los", "Length of stay, d")],
          [("Sex_male", "Male"), ("ICD10", "ICD-10 coded"), ("Biliary", "Biliary etiology (ICD-10)"), ("Alcoholic", "Alcoholic etiology (ICD-10)"),
           ("transfer_in", "Transfer from other hospital"), ("complete6", "All six variables measured"), ("icu", "ICU admission"),
           ("SAP_proxy", "SAP proxy"), ("support_persistent", "Organ support >= 48 h"), ("necro_coded", "Coded necrosis (ICD-10 specified)"),
           ("ipn_coded", "Coded infected necrosis (ICD-10 specified)"), ("sepsis", "Sepsis code"), ("epi_any", "Extra-pancreatic infection code"),
           ("pneumonia", "Pneumonia"), ("bacteremia", "Bacteremia"), ("uti", "UTI"), ("cholangitis", "Cholangitis"), ("cdiff", "C. difficile"),
           ("LOS14", "LOS > 14 d"), ("death", "In-hospital death")], "MIMIC-IV")
t1m.to_csv(os.path.join(RES, "table1_mimic.csv"), index=False)
t1e = t1(eic, [("Age", "Age, years"), ("Ca", "Calcium, mmol/L"), ("LYMPH_pct", "Lymphocyte, %"), ("RBC", "RBC, 10^12/L"),
               ("SCr", "Creatinine, umol/L"), ("ALT", "ALT, U/L"), ("ALB_gL", "Albumin, g/L"), ("apachescore", "APACHE IVa score"), ("los", "Length of stay, d")],
          [("Sex_male", "Male"), ("complete6", "All six variables measured"), ("necro_dx", "Necrotizing pancreatitis (diagnosis)"),
           ("abscess_dx", "Abscess/infected collection"), ("sepsis_dx", "Sepsis"), ("drain_dx", "Percutaneous drainage"),
           ("epi_any", "Extra-pancreatic infection"), ("LOS14", "LOS > 14 d"), ("death", "In-hospital death")], "eICU")
t1e.to_csv(os.path.join(RES, "table1_eicu.csv"), index=False)
print("\n=== Table1 本院 ===")
print(t1l.to_string(index=False))

# ---------- g 病因 ----------
rows = []
for co, df in [("Local", loc), ("MIMIC-IV (ICD-10)", mim[mim.icd10_ap == 1])]:
    ct = pd.crosstab(df.pheno, df.etiology)
    for ph in ct.index:
        for et in ct.columns:
            rows.append(dict(cohort=co, phenotype=f"P{ph}", etiology=et, n=int(ct.loc[ph, et]), frac=ct.loc[ph, et] / ct.loc[ph].sum()))
pd.DataFrame(rows).to_csv(os.path.join(RES, "etiology_by_pheno.csv"), index=False)
print("\n已落盘: death_phase / support_duration(+long) / los_long+summary / epi_vs_ipn / infection_site_mortality / pheno_mctsi_outcome / table1_x3 / etiology_by_pheno")
