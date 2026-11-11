# -*- coding: utf-8 -*-
"""Table S3 / S3b：低危层定义、事件数、CT 规则净获益窗口。
全部数字从 results/*.csv 读回后组装，不在本脚本内重算。"""
import pandas as pd
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = lambda var, sub: os.environ.get(var, os.path.join(ROOT, "inputs", sub))

P = ROOT
RES = os.path.join(P, "results")
tab = pd.read_csv(os.path.join(RES, "lowrisk_final_table.csv"))
sens = pd.read_csv(os.path.join(RES, "lowrisk_final_sens.csv"))
dca = pd.read_csv(os.path.join(RES, "lowrisk_final_dca.csv"))
CUT = float(pd.read_csv(os.path.join(RES, "lowrisk_final_threshold.csv")).threshold.iloc[0])
NALL = 1343


def win(st, o="Drainage/necrosectomy"):
    s = dca[(dca.stratum == st) & (dca.outcome == o) & (dca.threshold <= 0.30)]
    w = s[(s.nb_ct_rule > 0) & (s.nb_ct_rule > s.nb_treat_all)]
    if not len(w):
        return "none", "-"
    return f"{w.threshold.min()*100:.0f}-{w.threshold.max()*100:.0f}%", f"{w.nb_ct_rule.max():.4f}"


rows = []


def add(label, defn, adm, r, key):
    ww, mx = win(key)
    d = {}
    d["Stratum"] = label
    d["Definition"] = defn
    d["Executable at admission"] = adm
    d["n"] = int(r["n"])
    d["% of ward cohort"] = f"{int(r['n'])/NALL*100:.1f}%"
    d["CT necrosis, n"] = int(r["CT necrosis_ev"])
    d["Drainage/necrosectomy, n (%)"] = (f"{int(r['Drainage/necrosectomy_ev'])} "
                                         f"({r['Drainage/necrosectomy_rate']*100:.1f}%)")
    d["Persistent OF, n"] = int(r["Persistent OF_ev"])
    d["In-hospital deaths, n"] = int(r["In-hospital death_ev"])
    d["MCTSI >= 6 rule: threshold window"] = ww
    d["MCTSI >= 6 rule: max net benefit"] = mx
    rows.append(d)


full = tab[tab.scope == "Full ward cohort"]
lo = full[full.stratum.str.startswith("Low risk")].iloc[0]
hi = full[full.stratum.str.startswith("Higher risk")].iloc[0]
add("Primary: low risk",
    f"Continuous six-variable risk < {CUT:.4f} (lowest tertile of the 2017-2020 derivation period)",
    "Yes", lo, "Low risk (risk6 < cutoff)")
add("Primary: higher risk", f"Continuous six-variable risk >= {CUT:.4f}",
    "Yes", hi, "Higher risk (risk6 >= cutoff)")
add("Sensitivity A", "Host phenotype P1 or P2", "Yes",
    sens[sens.stratum.str.startswith("Sensitivity A")].iloc[0], "Sensitivity A: P1/P2")
add("Sensitivity B", "Host phenotype P1 or P2 and not severe AP",
    "No - severe AP requires organ failure > 48 h",
    sens[sens.stratum.str.startswith("Sensitivity B")].iloc[0], "Sensitivity B: P1/P2 non-SAP")

allr = {"n": NALL}
for k in ["CT necrosis_ev", "Drainage/necrosectomy_ev", "Persistent OF_ev", "In-hospital death_ev"]:
    allr[k] = int(lo[k] + hi[k])
allr["Drainage/necrosectomy_rate"] = allr["Drainage/necrosectomy_ev"] / NALL
add("Reference: whole ward cohort", "All admissions", "-", allr, "Whole ward cohort")

t = pd.DataFrame(rows)
t.to_csv(os.path.join(RES, "TableS3_lowrisk_strata.csv"), index=False, encoding="utf-8-sig")

v = tab[tab.scope != "Full ward cohort"].copy()
v["pctf"] = (v["pct"] * 100).round(1).astype(str) + "%"
vv = v[["scope", "stratum", "n", "pctf", "CT necrosis_ev", "Drainage/necrosectomy_ev",
        "Persistent OF_ev", "In-hospital death_ev"]]
vv.columns = ["Period", "Stratum", "n", "% of period", "CT necrosis, n",
              "Drainage/necrosectomy, n", "Persistent OF, n", "In-hospital deaths, n"]
vv.to_csv(os.path.join(RES, "TableS3b_lowrisk_by_period.csv"), index=False, encoding="utf-8-sig")

print(t.to_string(index=False))
print()
print(vv.to_string(index=False))
print("\n已落盘 TableS3_lowrisk_strata.csv / TableS3b_lowrisk_by_period.csv")
