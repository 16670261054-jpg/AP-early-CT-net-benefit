# -*- coding: utf-8 -*-
"""m01 MIMIC 补提取: 血脂肪酶/淀粉酶首值 [-6,+48h] (itemid 50956 / 50867,53087) + patients.dod (90 天死亡)
   输出 data/mimic_lipase_dod.csv"""
import pandas as pd
import numpy as np
import os, sys
sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = lambda var, sub: os.environ.get(var, os.path.join(ROOT, "inputs", sub))
M = IN("MIMIC_HOSP", "mimic-iv-3.1/hosp")
DAT = os.path.join(ROOT, "data")
mc = pd.read_csv(os.path.join(DAT, "mimic.csv"), usecols=["hadm_id"])
hset = set(mc.hadm_id)
adm = pd.read_csv(os.path.join(M, "admissions.csv.gz"), usecols=["subject_id", "hadm_id", "admittime", "dischtime", "deathtime"],
                  parse_dates=["admittime", "dischtime", "deathtime"])
adm = adm[adm.hadm_id.isin(hset)].copy()
h2t = dict(zip(adm.hadm_id, adm.admittime))
ITEMS = {50956: "Lipase", 50867: "Amylase", 53087: "Amylase"}
keep = []
for ch in pd.read_csv(os.path.join(M, "labevents.csv.gz"), usecols=["hadm_id", "itemid", "charttime", "valuenum"],
                      parse_dates=["charttime"], chunksize=5_000_000):
    ch = ch[ch.itemid.isin(ITEMS) & ch.hadm_id.notna() & ch.valuenum.notna()]
    ch = ch[ch.hadm_id.astype(int).isin(hset)]
    if len(ch):
        keep.append(ch)
lab = pd.concat(keep, ignore_index=True)
lab["hadm_id"] = lab.hadm_id.astype(int)
lab["hrs"] = (lab.charttime - lab.hadm_id.map(h2t)).dt.total_seconds() / 3600
lab = lab[(lab.hrs >= -6) & (lab.hrs <= 48)].sort_values("charttime")
lab["item"] = lab.itemid.map(ITEMS)
first = lab.groupby(["hadm_id", "item"]).valuenum.first().unstack()
mx = lab.groupby(["hadm_id", "item"]).valuenum.max().unstack().add_suffix("_max48")
out = adm[["subject_id", "hadm_id", "admittime", "dischtime", "deathtime"]].merge(first, on="hadm_id", how="left").merge(mx, on="hadm_id", how="left")
pat = pd.read_csv(os.path.join(M, "patients.csv.gz"), usecols=["subject_id", "dod"], parse_dates=["dod"])
out = out.merge(pat, on="subject_id", how="left")
out["dod_days"] = (out.dod - out.admittime.dt.normalize()).dt.days
out["death90"] = ((out.dod_days.notna()) & (out.dod_days <= 90)).astype(int)
out["death30"] = ((out.dod_days.notna()) & (out.dod_days <= 30)).astype(int)
out["death365"] = ((out.dod_days.notna()) & (out.dod_days <= 365)).astype(int)
out = out.drop(columns=["admittime", "dischtime", "deathtime", "dod"])
out.to_csv(os.path.join(DAT, "mimic_lipase_dod.csv"), index=False)
print(f"n={len(out)} lipase_avail={out.Lipase.notna().sum()} lipase>=3xULN(>=180)={int((out.Lipase_max48>=180).sum())} "
      f"death30={out.death30.sum()} death90={out.death90.sum()} death365={out.death365.sum()}")
print("已落盘 data/mimic_lipase_dod.csv")
