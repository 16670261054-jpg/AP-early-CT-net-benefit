# -*- coding: utf-8 -*-
"""02 MIMIC-IV 队列构建（AP_NECRO）
   基础: inputs/pheno_inputs/data/mimic.csv (4440, 六变量+插补+SAP_proxy+death) + labels_mimic.csv (表型)
         + mimic_severity.csv (死亡/出院时点, ICU, 乳酸, 早期支持)
   新提取: diagnoses_icd -> ICD-10 K85 第4位病因 / 第5位坏死状态; 脓毒症; 胰外感染五类
           admissions -> admission_location (转院)
           mimic_ap_labevents.csv.gz -> 白蛋白首值 [-6,+24h]
           inputs/probes/mimic_support_durations_probe.csv -> 通气/升压小时数
   输出: data/mimic.csv + results/audit_mimic.csv
"""
import pandas as pd
import numpy as np
import os, sys, joblib
import statsmodels.api as sm
sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = lambda var, sub: os.environ.get(var, os.path.join(ROOT, "inputs", sub))

P = ROOT
PH = IN("PHENO_INPUTS", "pheno_inputs")
AP = IN("AP_SAP_INPUTS", "sap_inputs")
XT = IN("AP_PROBES", "probes")
M = IN("MIMIC_HOSP", "mimic-iv-3.1/hosp")
R = ["Ca", "LYMPH_pct", "Age", "RBC", "SCr", "ALT"]

mc = pd.read_csv(os.path.join(PH, "data", "mimic.csv"))
lab = pd.read_csv(os.path.join(PH, "results", "labels_mimic.csv"))[["hadm_id", "pheno"]]
sev = pd.read_csv(os.path.join(PH, "data", "mimic_severity.csv"))
mc = mc.merge(lab, on="hadm_id").merge(sev, on="hadm_id", how="left")
assert len(mc) == 4440
hset = set(mc.hadm_id)

# ---------- 诊断编码 ----------
dx = pd.read_csv(os.path.join(M, "diagnoses_icd.csv.gz"),
                 usecols=["hadm_id", "seq_num", "icd_code", "icd_version"])
dx = dx[dx.hadm_id.isin(hset)].copy()
dx["icd_code"] = dx.icd_code.astype(str).str.strip()
k85 = dx[(dx.icd_version == 10) & dx.icd_code.str.startswith("K85")].copy()
k85["etio_digit"] = k85.icd_code.str[3:4]
k85["necro_digit"] = k85.icd_code.str[4:5]
# 每住院取最严重坏死状态 (2 感染坏死 > 1 无菌坏死 > 0 无)
k85["necro_rank"] = k85.necro_digit.map({"0": 0, "1": 1, "2": 2}).fillna(-1)
g = k85.groupby("hadm_id")
necro_rank = g.necro_rank.max()
etio_first = k85.sort_values("seq_num").groupby("hadm_id").etio_digit.first()
mc["icd10_ap"] = mc.hadm_id.isin(k85.hadm_id).astype(int)
mc["necro_rank"] = mc.hadm_id.map(necro_rank)
mc["necro_coded"] = np.where(mc.icd10_ap == 1, (mc.necro_rank >= 1).astype(float), np.nan)
mc["ipn_coded"] = np.where(mc.icd10_ap == 1, (mc.necro_rank == 2).astype(float), np.nan)
mc["necro_specified"] = np.where(mc.icd10_ap == 1, (mc.necro_rank >= 0).astype(float), np.nan)
etio_map = {"0": "Idiopathic", "1": "Biliary", "2": "Alcoholic", "3": "Drug", "8": "Other", "9": "Unspecified"}
mc["etiology"] = mc.hadm_id.map(etio_first).map(etio_map).fillna("ICD9_unknown")

def has(prefixes, version=None):
    m = dx.icd_code.str.startswith(tuple(prefixes))
    if version is not None:
        m &= dx.icd_version == version
    return mc.hadm_id.isin(dx.loc[m, "hadm_id"]).astype(int)

# 脓毒症: ICD-10 A40/A41/R652x ; ICD-9 038/99591/99592/78552
mc["sepsis"] = (has(["A40", "A41", "R6520", "R6521"], 10) | has(["038", "99591", "99592", "78552"], 9)).astype(int)
# 胰外感染 (与预实验一致的五类)
mc["pneumonia"] = (has(["J13", "J14", "J15", "J16", "J17", "J18", "J690"], 10) |
                   has(["481", "482", "483", "485", "486", "5070"], 9)).astype(int)
mc["bacteremia"] = (has(["R7881"], 10) | has(["7907"], 9)).astype(int)
mc["uti"] = (has(["N390"], 10) | has(["5990"], 9)).astype(int)
mc["cholangitis"] = (has(["K830"], 10) | has(["5761"], 9)).astype(int)
mc["cdiff"] = (has(["A047"], 10) | has(["00845"], 9)).astype(int)
mc["epi_any"] = mc[["pneumonia", "bacteremia", "uti", "cholangitis", "cdiff"]].max(axis=1)

# ---------- 入院来源 ----------
adm = pd.read_csv(os.path.join(M, "admissions.csv.gz"), usecols=["hadm_id", "admission_location"])
mc = mc.merge(adm, on="hadm_id", how="left")
mc["transfer_in"] = mc.admission_location.astype(str).str.startswith("TRANSFER").astype(int)

# ---------- 白蛋白 ----------
le = pd.read_csv(os.path.join(AP, "data", "mimic_ap_labevents.csv.gz"), parse_dates=["charttime"])
alb = le[(le.item == "Albumin") & le.hadm_id.isin(hset)].copy()
admt = pd.read_csv(os.path.join(M, "admissions.csv.gz"), usecols=["hadm_id", "admittime"], parse_dates=["admittime"])
alb = alb.merge(admt, on="hadm_id")
alb["hrs"] = (alb.charttime - alb.admittime).dt.total_seconds() / 3600
alb = alb[(alb.hrs >= -6) & (alb.hrs <= 24)].sort_values("charttime").groupby("hadm_id").valuenum.first()
mc["Albumin"] = mc.hadm_id.map(alb) * 10          # g/dL -> g/L
# 缺失以 dev 中位数填补 (同六变量策略)
med = pd.read_csv(os.path.join(ROOT, "models", "impute_medians.csv"), index_col=0)["0"]
mc["Albumin_imp"] = mc.Albumin.fillna(med["Albumin"])

# ---------- 支持时长 / 结局时点 ----------
sd = pd.read_csv(os.path.join(XT, "mimic_support_durations_probe.csv"))
mc = mc.merge(sd, on="hadm_id", how="left")
mc["vent_h"] = mc.vent_h.fillna(0); mc["vaso_h"] = mc.vaso_h.fillna(0)
mc["support_any"] = ((mc.vent_h > 0) | (mc.vaso_h > 0)).astype(int)
mc["support_persistent"] = ((mc.vent_h >= 48) | (mc.vaso_h >= 48)).astype(int)
mc["los"] = mc.disch_hrs / 24.0
mc["LOS14"] = (mc.los > 14).astype(int)
mc["death_day"] = mc.death_hrs / 24.0
mc["death_phase"] = pd.cut(mc.death_day, [-np.inf, 7, 14, np.inf], labels=["<=7d", "8-14d", ">14d"]).astype(str)
mc.loc[mc.death != 1, "death_phase"] = "alive"
mc["icu"] = mc.icu_in_hrs.notna().astype(int)
mc["nolab6"] = (mc[R].isna().all(axis=1)).astype(int)

# ---------- 连续风险 ----------
gb = joblib.load(os.path.join(ROOT, "models", "GBDT_reduced.pkl"))
mc["risk6"] = gb.predict_proba(mc[[f + "_imp" for f in R]].values)[:, 1]
coef = pd.read_csv(os.path.join(P, "results", "simple_model_coef.csv")).set_index("term").beta
lin = coef["const"] + coef["Age"] * mc.Age_imp + coef["Albumin"] * mc.Albumin_imp + coef["LYMPH_pct"] * mc.LYMPH_pct_imp
mc["risk3"] = 1 / (1 + np.exp(-lin))

# ---------- 剔除负/零 LOS ----------
neg = mc.los <= 0
audit = pd.DataFrame([
    dict(item="n_total", value=len(mc)),
    dict(item="n_los_nonpositive_excluded", value=int(neg.sum())),
    dict(item="n_icd10_subset", value=int(mc.icd10_ap.sum())),
    dict(item="n_icd10_necro_specified", value=int((mc.necro_specified == 1).sum())),
    dict(item="necro_coded", value=int((mc.necro_coded == 1).sum())),
    dict(item="ipn_coded", value=int((mc.ipn_coded == 1).sum())),
    dict(item="sepsis", value=int(mc.sepsis.sum())),
    dict(item="epi_any", value=int(mc.epi_any.sum())),
    dict(item="pneumonia", value=int(mc.pneumonia.sum())),
    dict(item="bacteremia", value=int(mc.bacteremia.sum())),
    dict(item="uti", value=int(mc.uti.sum())),
    dict(item="cholangitis", value=int(mc.cholangitis.sum())),
    dict(item="cdiff", value=int(mc.cdiff.sum())),
    dict(item="transfer_in", value=int(mc.transfer_in.sum())),
    dict(item="nolab6", value=int(mc.nolab6.sum())),
    dict(item="complete6", value=int(mc.complete6.sum())),
    dict(item="albumin_available", value=int(mc.Albumin.notna().sum())),
    dict(item="death", value=int(mc.death.sum())),
    dict(item="death_cohort", value=int(mc.loc[~neg, "death"].sum())),
    dict(item="epi_any_cohort", value=int(mc.loc[~neg, "epi_any"].sum())),
    dict(item="sepsis_cohort", value=int(mc.loc[~neg, "sepsis"].sum())),
    dict(item="support_any", value=int(mc.support_any.sum())),
    dict(item="support_persistent", value=int(mc.support_persistent.sum())),
    dict(item="icu", value=int(mc.icu.sum())),
])
mc = mc[~neg].copy()
mc.to_csv(os.path.join(P, "data", "mimic.csv"), index=False)
audit.to_csv(os.path.join(P, "results", "audit_mimic.csv"), index=False)
print(audit.to_string(index=False))
print("\n病因分布:\n", mc.etiology.value_counts().to_string())
print("\n死亡时相:\n", mc.death_phase.value_counts().to_string())
print("\n已落盘 data/mimic.csv (n=%d), results/audit_mimic.csv" % len(mc))
