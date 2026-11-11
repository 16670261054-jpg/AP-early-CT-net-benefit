# 09 落盘 Figure2 (表型画像) 所需数据: 六变量 z-score 画像 + risk6 连续风险长表
# 规则: 出图只读盘 -> 本脚本先把数字落盘 results/
import os
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
RES = os.path.join(BASE, "results")

VARS = ["Age_imp", "Ca_imp", "LYMPH_pct_imp", "RBC_imp", "SCr_imp", "ALT_imp"]
LAB = {"Age_imp": "Age", "Ca_imp": "Calcium", "LYMPH_pct_imp": "Lymphocyte %",
       "RBC_imp": "RBC", "SCr_imp": "Creatinine", "ALT_imp": "ALT"}

def load(name, cohort):
    df = pd.read_csv(os.path.join(DATA, name))
    if cohort == "Local":
        # local.csv 无 *_imp 列名 -> 原列即分析用列 (构建时已插补)
        ren = {"Age": "Age_imp", "Ca": "Ca_imp", "LYMPH_pct": "LYMPH_pct_imp",
               "RBC": "RBC_imp", "SCr": "SCr_imp", "ALT": "ALT_imp"}
        df = df.rename(columns=ren)
    df = df[["pheno", "risk6"] + VARS].copy()
    df["cohort"] = cohort
    return df

frames = [load("local.csv", "Local"), load("mimic.csv", "MIMIC-IV"), load("eicu.csv", "eICU")]

prof_rows, risk_rows = [], []
for df in frames:
    cohort = df["cohort"].iloc[0]
    d = df.dropna(subset=["pheno"]).copy()
    d["pheno"] = d["pheno"].astype(int)
    # 队列内 z 标准化 (对全队列), 再按表型取均值
    for v in VARS:
        mu, sd = d[v].mean(), d[v].std()
        d[v + "_z"] = (d[v] - mu) / sd
    for p, g in d.groupby("pheno"):
        for v in VARS:
            prof_rows.append({"cohort": cohort, "phenotype": f"P{p}", "variable": LAB[v],
                              "z_mean": g[v + "_z"].mean(),
                              "z_lo": g[v + "_z"].mean() - 1.96 * g[v + "_z"].std() / np.sqrt(len(g)),
                              "z_hi": g[v + "_z"].mean() + 1.96 * g[v + "_z"].std() / np.sqrt(len(g)),
                              "raw_median": g[v].median(), "n": len(g)})
    r = d.dropna(subset=["risk6"])[["pheno", "risk6"]].copy()
    r["cohort"] = cohort
    r["phenotype"] = "P" + r["pheno"].astype(str)
    risk_rows.append(r[["cohort", "phenotype", "risk6"]])

prof = pd.DataFrame(prof_rows)
prof.to_csv(os.path.join(RES, "pheno_profile_z.csv"), index=False)
risk = pd.concat(risk_rows, ignore_index=True)
risk.to_csv(os.path.join(RES, "risk6_by_pheno_long.csv"), index=False)

# risk6 组间 Kruskal-Wallis (每队列) 落盘
from scipy.stats import kruskal
kw = []
for cohort, g in risk.groupby("cohort"):
    groups = [x["risk6"].values for _, x in g.groupby("phenotype")]
    kw.append({"cohort": cohort, "kruskal_p": kruskal(*groups).pvalue,
               **{f"median_{p}": g.loc[g.phenotype == p, "risk6"].median() for p in ["P1", "P2", "P3", "P4"]}})
pd.DataFrame(kw).to_csv(os.path.join(RES, "risk6_by_pheno_kw.csv"), index=False)

print("rows:", len(prof), len(risk))
print(prof.head(8).to_string())
print(pd.DataFrame(kw).to_string())
