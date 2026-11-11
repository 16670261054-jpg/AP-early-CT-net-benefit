# -*- coding: utf-8 -*-
"""01 本院队列构建（AP_NECRO：入院表型 -> 坏死-感染进程）
   输入: inputs/probes/local_late_outcomes_probe.csv (1343例, dev1096+val247, 含表型/晚期结局/MCTSI)
         models/GBDT_reduced.pkl (six-variable continuous risk; derivation patients receive cross-fitted predictions)
   输出: data/local.csv + results/audit_local.csv + results/simple_model_coef.csv
   审计规则(定题记录):
     - allimp: 六项化验全插补者 70 例 -> 标记, 主分析保留, 敏感性剔除
     - 干预 = 经皮引流(pcd) 或 手术字符串含 坏死/脓肿/清创/感染 (排除 LC/ERCP/胆道手术)
"""
import pandas as pd
import numpy as np
import os, sys, re, joblib
import statsmodels.api as sm
sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = lambda var, sub: os.environ.get(var, os.path.join(ROOT, "inputs", sub))

P = ROOT
XT = IN("AP_PROBES", "probes")
AP = IN("AP_SAP_INPUTS", "sap_inputs")
R = ["Ca", "LYMPH_pct", "Age", "RBC", "SCr", "ALT"]

d = pd.read_csv(os.path.join(XT, "local_late_outcomes_probe.csv"))
assert len(d) == 1343 and d.pheno.between(1, 4).all()

# ---------- 结局 ----------
d["necroCT"] = d.necroCT.astype(int)
surg = d.surg.fillna("").astype(str)
d["necrosectomy"] = surg.str.contains(r"坏死|脓肿|清创|感染|假性囊肿", regex=True).astype(int)
d["intervention"] = ((d.pcd == 1) | (d.necrosectomy == 1)).astype(int)
d["LOS14"] = (d.LOS > 14).astype(int)
d["persOF"] = d.persOF.astype(int)
d["death"] = d.death.astype(int)
d["icu"] = d.icu.astype(int)
d["MCTSI"] = d.MCTSI.astype(int)
d["mctsi_high"] = (d.MCTSI >= 6).astype(int)    # 重度 CT 表现 (Mortele 2004: 8-10 重; 我们用>=6 定义中重)
d["sev3"] = d.severity3
d["MSAP_or_SAP"] = d.severity3.isin(["MSAP", "SAP"]).astype(int)


# ---------- 六变量连续风险 (JTM 文 GBDT) 与三变量简模型 ----------
gb = joblib.load(os.path.join(ROOT, "models", "GBDT_reduced.pkl"))
d["risk6"] = gb.predict_proba(d[R].values)[:, 1]
# derivation-period patients: cross-fitted predictions (5 x 5-fold, same hyperparameters, averaged),
# so that no patient is scored by a model trained on that patient; 2021 keeps the full-derivation model
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
_dev = (d.cohort == "dev").values
_X, _y = d[R].values.astype(float)[_dev], d.SAP.values[_dev]
for _rep in range(5):
    _o = np.zeros(_dev.sum())
    for _tr, _te in StratifiedKFold(5, shuffle=True, random_state=_rep).split(_X, _y):
        _m = GradientBoostingClassifier(**gb.get_params()).fit(_X[_tr], _y[_tr])
        _o[_te] = _m.predict_proba(_X[_te])[:, 1]
    _sum = _o if _rep == 0 else _sum + _o
d.loc[_dev, "risk6"] = _sum / 5
dev = d[d.cohort == "dev"]
Xs = sm.add_constant(dev[["Age", "Albumin", "LYMPH_pct"]].astype(float))
m3 = sm.Logit(dev.SAP.values, Xs).fit(disp=0)
d["risk3"] = m3.predict(sm.add_constant(d[["Age", "Albumin", "LYMPH_pct"]].astype(float)))
coef = pd.DataFrame({"term": m3.params.index, "beta": m3.params.values, "se": m3.bse.values,
                     "OR": np.exp(m3.params.values), "p": m3.pvalues.values})
coef.to_csv(os.path.join(P, "results", "simple_model_coef.csv"), index=False)

# ---------- 病因 ----------
d["etiology"] = np.select([d.Etio_Alcoholic == 1, d.Etio_Hypertriglyceridemic == 1, d.Etio_Other == 1],
                          ["Alcoholic", "Hypertriglyceridemic", "Other"], default="Biliary")

out_cols = R + ["Sex_male", "SAP", "MSAP_or_SAP", "sev3", "etiology", "BUN", "Albumin", "WBC", "HCT", "PLT",
                "BMI", "BISAP", "row_idx", "hid", "cohort", "pheno", "risk6", "risk3", "allimp",
                "MCTSI", "mctsi_high", "necroCT", "apfc", "pleural", "pcd",
                "necrosectomy", "intervention", "LOS", "LOS14", "persOF", "icu", "ntube", "death"]
d[out_cols].to_csv(os.path.join(P, "data", "local.csv"), index=False)

# ---------- 审计汇总 ----------
summ = pd.DataFrame([
    dict(item="n_total", value=len(d)),
    dict(item="n_dev", value=int((d.cohort == "dev").sum())),
    dict(item="n_val", value=int((d.cohort == "val").sum())),
    dict(item="n_allimp6", value=int(d.allimp.sum())),
    dict(item="necroCT", value=int(d.necroCT.sum())),
    dict(item="pcd", value=int(d.pcd.sum())),
    dict(item="necrosectomy", value=int(d.necrosectomy.sum())),
    dict(item="intervention_any", value=int(d.intervention.sum())),
    dict(item="LOS_gt14", value=int(d.LOS14.sum())),
    dict(item="persOF", value=int(d.persOF.sum())),
    dict(item="death", value=int(d.death.sum())),
    dict(item="MCTSI_ge6", value=int(d.mctsi_high.sum())),
    dict(item="LOS_min", value=round(d.LOS.min(), 3)),
])
summ.to_csv(os.path.join(P, "results", "audit_local.csv"), index=False)
print(summ.to_string(index=False))
print("\n简模型 (dev, SAP ~ Age+Albumin+LYMPH_pct):")
print(coef.round(4).to_string(index=False))
print("\n已落盘: data/local.csv, results/audit_local.csv, simple_model_coef.csv")
