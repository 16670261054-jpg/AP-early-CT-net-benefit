# -*- coding: utf-8 -*-
"""05c 低风险层终版定义

背景：原低风险层 = (pheno in {1,2}) & (SAP == 0)。SAP 依修订版 Atlanta 需持续器官
衰竭 >48h 才能判定，入院当天不可知，故原层不是一条可在入院执行的规则。

终版定义（完全入院可得、阈值不看验证数据）：
  低风险 = risk6 < c，其中 c = derivation set (2017-2020, n=1096) risk6 的第 33.3 百分位。
  阈值只由 derivation set 决定，不在全队列上重新 qcut，避免阈值借用了验证期信息。

输出 results/:
  lowrisk_final_threshold.csv   阈值本身与来源
  lowrisk_final_table.csv       各层事件数/率（全队列、derivation、2021 时间验证）
  lowrisk_final_dca.csv         "MCTSI>=6 才干预" vs treat-all vs treat-none 的净获益
  lowrisk_final_sens.csv        敏感性：P1/P2 纯表型层、原 P1/P2 非SAP 层
"""
import pandas as pd
import numpy as np
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = lambda var, sub: os.environ.get(var, os.path.join(ROOT, "inputs", sub))

P = ROOT
RES, DAT = os.path.join(P, "results"), os.path.join(P, "data")
loc = pd.read_csv(os.path.join(DAT, "local.csv"))

# ---------- 1. 阈值：只从 derivation set 取 ----------
dev = loc[loc.cohort == "dev"]
CUT = float(np.quantile(dev.risk6, 1.0 / 3.0))
pd.DataFrame([dict(threshold=CUT, source="derivation set 2017-2020 lowest tertile of risk6",
                   n_derivation=len(dev), quantile=1 / 3,
                   note="threshold fixed in derivation data only; applied unchanged to 2021 validation period")]
             ).to_csv(os.path.join(RES, "lowrisk_final_threshold.csv"), index=False)
print(f"低风险阈值 risk6 < {CUT:.6f}  (derivation n={len(dev)} 的最低三分位)")

loc["lowrisk_final"] = (loc.risk6 < CUT).astype(int)

OUTS = [("necroCT", "CT necrosis"), ("intervention", "Drainage/necrosectomy"),
        ("persOF", "Persistent OF"), ("LOS14", "LOS>14d"), ("death", "In-hospital death")]


def describe(s, label, scope):
    r = dict(stratum=label, scope=scope, n=len(s), pct=len(s) / len(loc[loc.cohort.isin(s.cohort.unique())]) if len(s) else np.nan,
             n_SAP=int(s.SAP.sum()), pct_SAP=s.SAP.mean(),
             mctsi_ge6_n=int((s.MCTSI >= 6).sum()), mctsi_ge6_pct=(s.MCTSI >= 6).mean())
    for y, lab in OUTS:
        v = s[y].dropna()
        r[f"{lab}_n_eval"] = len(v)
        r[f"{lab}_ev"] = int(v.sum())
        r[f"{lab}_rate"] = v.mean() if len(v) else np.nan
        hi = s.loc[(s.MCTSI >= 6) & s[y].notna(), y]
        lo = s.loc[(s.MCTSI < 6) & s[y].notna(), y]
        r[f"{lab}_ev_ct_ge6"] = int(hi.sum())
        r[f"{lab}_ev_ct_lt6"] = int(lo.sum())
    return r


# ---------- 2. 主表：全队列 + 两个时期 ----------
rows = []
for scope, sub in [("Full ward cohort", loc), ("Derivation 2017-2020", loc[loc.cohort == "dev"]),
                   ("Temporal validation 2021", loc[loc.cohort == "val"])]:
    rows.append(describe(sub[sub.lowrisk_final == 1], "Low risk (risk6 below derivation tertile-1 cutoff)", scope))
    rows.append(describe(sub[sub.lowrisk_final == 0], "Higher risk (risk6 at or above cutoff)", scope))
tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(RES, "lowrisk_final_table.csv"), index=False)
print("\n=== 终版低风险层 ===")
cols = ["stratum", "scope", "n", "pct", "n_SAP", "mctsi_ge6_n",
        "Drainage/necrosectomy_ev", "Drainage/necrosectomy_rate",
        "Drainage/necrosectomy_ev_ct_ge6", "Drainage/necrosectomy_ev_ct_lt6",
        "Persistent OF_ev", "In-hospital death_ev", "CT necrosis_ev"]
print(tab[cols].round(4).to_string(index=False))

# ---------- 3. 敏感性层 ----------
sens = [describe(loc[loc.pheno.isin([1, 2])], "Sensitivity A: phenotype P1/P2 (admission-only)", "Full ward cohort"),
        describe(loc[loc.pheno.isin([1, 2]) & (loc.SAP == 0)], "Sensitivity B: P1/P2 non-SAP (original, needs 48h)", "Full ward cohort"),
        describe(loc[(loc.lowrisk_final == 1) & (loc.SAP == 0)], "Sensitivity C: low risk & non-SAP", "Full ward cohort")]
pd.DataFrame(sens).to_csv(os.path.join(RES, "lowrisk_final_sens.csv"), index=False)
print("\n=== 敏感性层 ===")
print(pd.DataFrame(sens)[cols].round(4).to_string(index=False))

# ---------- 4. 决策曲线 ----------
dca = []
STRATA = [("Low risk (risk6 < cutoff)", loc.lowrisk_final == 1),
          ("Higher risk (risk6 >= cutoff)", loc.lowrisk_final == 0),
          ("Sensitivity A: P1/P2", loc.pheno.isin([1, 2])),
          ("Sensitivity B: P1/P2 non-SAP", loc.pheno.isin([1, 2]) & (loc.SAP == 0)),
          ("Whole ward cohort", loc.pheno > 0)]
for name, mask in STRATA:
    s = loc[mask]
    for y, lab in [("intervention", "Drainage/necrosectomy"), ("persOF", "Persistent OF")]:
        yv = s[y].astype(int).values
        rule = (s.MCTSI >= 6).astype(int).values
        prev = yv.mean()
        for t in np.arange(0.02, 0.41, 0.01):
            tp = np.sum((rule == 1) & (yv == 1)) / len(yv)
            fp = np.sum((rule == 1) & (yv == 0)) / len(yv)
            w = t / (1 - t)
            dca.append(dict(stratum=name, outcome=lab, threshold=round(float(t), 3), n=len(s),
                            events=int(yv.sum()), prevalence=prev,
                            nb_ct_rule=tp - fp * w, nb_treat_all=prev - (1 - prev) * w, nb_none=0.0))
d = pd.DataFrame(dca)
d.to_csv(os.path.join(RES, "lowrisk_final_dca.csv"), index=False)

print("\n=== CT 规则 (MCTSI>=6) 同时优于 treat-none 与 treat-all 的阈值窗 · Drainage/necrosectomy ===")
for name, _ in STRATA:
    sub = d[(d.stratum == name) & (d.outcome == "Drainage/necrosectomy")]
    win = sub[(sub.nb_ct_rule > 0) & (sub.nb_ct_rule > sub.nb_treat_all)]
    if len(win):
        print(f"{name:34s} 窗口 {win.threshold.min():.2f}-{win.threshold.max():.2f}"
              f"  max NB {win.nb_ct_rule.max():.5f}  (= 每千人净多识别 {win.nb_ct_rule.max()*1000:.1f} 例真阳)")
    else:
        print(f"{name:34s} 无窗口")

print("\n已落盘: lowrisk_final_threshold / lowrisk_final_table / lowrisk_final_sens / lowrisk_final_dca")
