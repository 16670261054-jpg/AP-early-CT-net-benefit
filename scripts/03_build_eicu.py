# -*- coding: utf-8 -*-
"""03 eICU 队列构建（AP_NECRO, 方向性验证）
   基础: inputs/pheno_inputs/data/eicu.csv (941) + labels_eicu.csv
   合并: inputs/probes/eicu_late_outcomes_probe.csv (坏死/脓肿/脓毒症/引流/LOS)
         inputs/probes/eicu_epi_infection_probe.csv (胰外感染五类)
         apachePatientResult -> APACHE IV 分数 (调整用)
   输出: data/eicu.csv + results/audit_eicu.csv
"""
import pandas as pd
import numpy as np
import os, sys, joblib
sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = lambda var, sub: os.environ.get(var, os.path.join(ROOT, "inputs", sub))

P = ROOT
PH = IN("PHENO_INPUTS", "pheno_inputs")
AP = IN("AP_SAP_INPUTS", "sap_inputs")
XT = IN("AP_PROBES", "probes")
E = IN("EICU_DIR", "eicu-crd-2.0")
R = ["Ca", "LYMPH_pct", "Age", "RBC", "SCr", "ALT"]

e = pd.read_csv(os.path.join(PH, "data", "eicu.csv"))
lab = pd.read_csv(os.path.join(PH, "results", "labels_eicu.csv"))[["patientunitstayid", "pheno"]]
lo = pd.read_csv(os.path.join(XT, "eicu_late_outcomes_probe.csv")).drop(columns=["pheno"])
ep = pd.read_csv(os.path.join(XT, "eicu_epi_infection_probe.csv")).drop(columns=["pheno"])
e = e.merge(lab, on="patientunitstayid").merge(lo, on="patientunitstayid", how="left") \
     .merge(ep, on="patientunitstayid", how="left")
assert len(e) == 941
e = e.rename(columns={"necro": "necro_dx", "abscess": "abscess_dx", "sepsis": "sepsis_dx",
                      "drain": "drain_dx", "any_epi": "epi_any"})
ap = pd.read_csv(os.path.join(E, "apachePatientResult.csv.gz"),
                 usecols=["patientunitstayid", "apachescore", "apacheversion", "predictedhospitalmortality"])
ap = ap[ap.apacheversion == "IVa"].drop_duplicates("patientunitstayid")
e = e.merge(ap[["patientunitstayid", "apachescore"]], on="patientunitstayid", how="left")
e.loc[e.apachescore < 0, "apachescore"] = np.nan
e["LOS14"] = (e.los > 14).astype(int)
e["death"] = e.hosp_mort.astype(int)
e["transfer_in"] = e.hospitaladmitsource.astype(str).str.contains("Other Hospital", case=False).astype(int)

gb = joblib.load(os.path.join(ROOT, "models", "GBDT_reduced.pkl"))
e["risk6"] = gb.predict_proba(e[[f + "_imp" for f in R]].values)[:, 1]
med = pd.read_csv(os.path.join(ROOT, "models", "impute_medians.csv"), index_col=0)["0"]
e["Albumin_imp"] = e.ALB_gL.fillna(med["Albumin"])
coef = pd.read_csv(os.path.join(P, "results", "simple_model_coef.csv")).set_index("term").beta
lin = coef["const"] + coef["Age"] * e.Age_imp + coef["Albumin"] * e.Albumin_imp + coef["LYMPH_pct"] * e.LYMPH_pct_imp
e["risk3"] = 1 / (1 + np.exp(-lin))

audit = pd.DataFrame([
    dict(item="n_total", value=len(e)),
    dict(item="hospitals", value=int(e.hospitalid.nunique())),
    dict(item="necro_dx", value=int(e.necro_dx.sum())),
    dict(item="abscess_dx", value=int(e.abscess_dx.sum())),
    dict(item="sepsis_dx", value=int(e.sepsis_dx.sum())),
    dict(item="drain_dx", value=int(e.drain_dx.sum())),
    dict(item="epi_any", value=int(e.epi_any.sum())),
    dict(item="LOS_gt14", value=int(e.LOS14.sum())),
    dict(item="death", value=int(e.death.sum())),
    dict(item="apache_available", value=int(e.apachescore.notna().sum())),
    dict(item="complete6", value=int(e.complete6.sum())),
    dict(item="transfer_in", value=int(e.transfer_in.sum())),
])
e.to_csv(os.path.join(P, "data", "eicu.csv"), index=False)
audit.to_csv(os.path.join(P, "results", "audit_eicu.csv"), index=False)
print(audit.to_string(index=False))
print("已落盘 data/eicu.csv, results/audit_eicu.csv")
