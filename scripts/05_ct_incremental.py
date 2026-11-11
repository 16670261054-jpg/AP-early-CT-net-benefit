# -*- coding: utf-8 -*-
"""05 本院：入院 CT (MCTSI) 在化验模型之上的增量价值
   结局: intervention (引流或清创 81) / persOF (175) / LOS14 (564)
   基线模型 A: 表型哑变量 + age + sex + etiology  ; A': risk6 + age + sex + etiology (连续风险版)
   增量模型 B: A + MCTSI (连续 0-10)
   指标: 10x 重复 5 折交叉验证 out-of-fold 预测 -> AUC(bootstrap CI), DeLong p, 连续NRI, IDI, Brier, 校准斜率;
         决策曲线 (阈值 2-40%); 亚组增量 (按表型 / risk6 三分位 / SAP 与否): 各亚组内 MCTSI 的 AUC 与事件率;
         "CT 无增量人群" = 低风险层 (P1+P2 且 SAP=0) 中 MCTSI>=6 与 <6 的事件率
   输出 results/: ct_incr_metrics.csv, ct_incr_roc.csv, ct_incr_dca.csv, ct_incr_subgroup.csv,
                  ct_lowrisk_table.csv, ct_incr_calib.csv, mctsi_dist.csv
"""
import pandas as pd
import numpy as np
import os, sys, warnings
import statsmodels.api as sm
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss
from scipy.stats import norm
warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = lambda var, sub: os.environ.get(var, os.path.join(ROOT, "inputs", sub))
SEED = 123
P = ROOT
RES, DAT = os.path.join(P, "results"), os.path.join(P, "data")

loc = pd.read_csv(os.path.join(DAT, "local.csv"))
loc["Sex"] = loc.Sex_male

def design(df, base, with_ct):
    X = pd.DataFrame(index=df.index)
    if base == "pheno":
        X = pd.concat([X, pd.get_dummies(df.pheno.astype(int), prefix="P").astype(float).iloc[:, 1:]], axis=1)
    elif base == "risk6":
        X["risk6_logit"] = np.log(df.risk6.clip(1e-4, 1 - 1e-4) / (1 - df.risk6.clip(1e-4, 1 - 1e-4)))
    elif base == "none":
        pass
    X["Age"] = df.Age.astype(float); X["Sex"] = df.Sex.astype(float)
    X = pd.concat([X, pd.get_dummies(df.etiology, prefix="E", drop_first=True).astype(float)], axis=1)
    if with_ct:
        X["MCTSI"] = df.MCTSI.astype(float)
    return X

def oof_pred(X, y, n_rep=10):
    """重复5折 out-of-fold 概率 (每重复一份, 取平均)"""
    pr = np.zeros((n_rep, len(y)))
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=n_rep, random_state=SEED)
    for i, (tr, te) in enumerate(rskf.split(X, y)):
        r = i // 5
        m = LogisticRegression(C=1e6, max_iter=2000).fit(X.iloc[tr], y[tr])
        pr[r, te] = m.predict_proba(X.iloc[te])[:, 1]
    return pr.mean(0)

# ---- DeLong (fastDeLong 标准实现) ----
def compute_midrank(x):
    J = np.argsort(x); Z = x[J]; N = len(x); T = np.zeros(N)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N); T2[J] = T
    return T2

def fast_delong(preds, y):
    order = np.argsort(-y); preds = preds[:, order]; m = int(y.sum()); n = len(y) - m
    pos, neg = preds[:, :m], preds[:, m:]
    k = preds.shape[0]
    tx = np.array([compute_midrank(pos[r]) for r in range(k)])
    ty = np.array([compute_midrank(neg[r]) for r in range(k)])
    tz = np.array([compute_midrank(preds[r]) for r in range(k)])
    aucs = tz[:, :m].sum(1) / m / n - (m + 1.0) / (2.0 * n)
    v01 = (tz[:, :m] - tx) / n; v10 = 1.0 - (tz[:, m:] - ty) / m
    sx = np.cov(v01); sy = np.cov(v10)
    return aucs, sx / m + sy / n

def delong_p(p1, p2, y):
    aucs, cov = fast_delong(np.vstack([p1, p2]), y)
    l = np.array([[1, -1]])
    z = np.abs(np.diff(aucs)) / np.sqrt(l @ cov @ l.T)
    return float(2 * (1 - norm.cdf(z))[0, 0])

def boot_ci(fn, y, *preds, B=1000):
    rng = np.random.RandomState(SEED); vals = []
    for _ in range(B):
        b = rng.choice(len(y), len(y), replace=True)
        if y[b].sum() < 5:
            continue
        vals.append(fn(y[b], *[p[b] for p in preds]))
    return np.nanpercentile(vals, [2.5, 97.5])

def nri_idi(y, p0, p1):
    ev, ne = y == 1, y == 0
    nri = (np.mean(p1[ev] > p0[ev]) - np.mean(p1[ev] < p0[ev])) + (np.mean(p1[ne] < p0[ne]) - np.mean(p1[ne] > p0[ne]))
    idi = (p1[ev].mean() - p0[ev].mean()) - (p1[ne].mean() - p0[ne].mean())
    return nri, idi

def calib_slope(y, p):
    lg = np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))
    m = sm.Logit(y, sm.add_constant(lg)).fit(disp=0)
    return m.params[1], m.params[0]

OUTS = [("intervention", "Drainage/necrosectomy"),
        ("persOF", "Persistent OF"), ("LOS14", "LOS > 14 d")]
met, roc_rows, dca_rows, cal_rows = [], [], [], []
preds_store = {}
for y, lab in OUTS:
    df = loc[loc[y].notna()].reset_index(drop=True)
    yv = df[y].astype(int).values
    P_ = {}
    for base in ["pheno", "risk6", "none"]:
        for ct in [False, True]:
            if base == "none" and not ct:
                continue
            key = (base, ct)
            P_[key] = oof_pred(design(df, base, ct), yv)
    P_[("mctsi_only", True)] = df.MCTSI.values.astype(float) / 10  # 单独 MCTSI (无需拟合, 序数)
    preds_store[y] = (yv, P_)
    for key, p in P_.items():
        base, ct = key
        auc = roc_auc_score(yv, p)
        ci = boot_ci(lambda yy, pp: roc_auc_score(yy, pp), yv, p)
        row = dict(outcome=lab, y=y, base=base, with_ct=int(ct), n=len(yv), events=int(yv.sum()),
                   AUC=auc, AUC_lo=ci[0], AUC_hi=ci[1])
        if base in ("pheno", "risk6", "none"):
            row["Brier"] = brier_score_loss(yv, p)
            sl, ic = calib_slope(yv, p); row["calib_slope"] = sl; row["calib_intercept"] = ic
        if ct and base in ("pheno", "risk6"):
            p0 = P_[(base, False)]
            row["delong_p_vs_base"] = delong_p(p0, p, yv)
            nri, idi = nri_idi(yv, p0, p)
            nci = boot_ci(lambda yy, a, b: nri_idi(yy, a, b)[0], yv, p0, p)
            ici = boot_ci(lambda yy, a, b: nri_idi(yy, a, b)[1], yv, p0, p)
            row.update(NRI=nri, NRI_lo=nci[0], NRI_hi=nci[1], IDI=idi, IDI_lo=ici[0], IDI_hi=ici[1],
                       dAUC=auc - roc_auc_score(yv, p0))
        met.append(row)
        # ROC 曲线点
        if key in [("pheno", False), ("pheno", True), ("risk6", False), ("risk6", True), ("mctsi_only", True)]:
            fpr, tpr, _ = roc_curve(yv, p)
            for f_, t_ in zip(fpr, tpr):
                roc_rows.append(dict(outcome=lab, model=f"{base}{'+CT' if ct and base != 'mctsi_only' else ''}", fpr=f_, tpr=t_))
        # DCA
        if key in [("pheno", False), ("pheno", True), ("risk6", False), ("risk6", True), ("none", True)]:
            prev = yv.mean()
            for t in np.arange(0.02, 0.41, 0.01):
                pred_pos = p >= t
                tp = np.sum(pred_pos & (yv == 1)) / len(yv); fp = np.sum(pred_pos & (yv == 0)) / len(yv)
                nb = tp - fp * t / (1 - t)
                dca_rows.append(dict(outcome=lab, model=f"{base}{'+CT' if ct else ''}", threshold=t, net_benefit=nb,
                                     treat_all=prev - (1 - prev) * t / (1 - t)))
        # 校准 (十分位)
        if key in [("pheno", True), ("risk6", True), ("pheno", False)]:
            q = pd.qcut(pd.Series(p).rank(method="first"), 10, labels=False)
            for g in range(10):
                mask = q == g
                cal_rows.append(dict(outcome=lab, model=f"{base}{'+CT' if ct else ''}", decile=g + 1,
                                     pred_mean=p[mask].mean(), obs_rate=yv[mask].mean(), n=int(mask.sum())))
metf = pd.DataFrame(met); metf.to_csv(os.path.join(RES, "ct_incr_metrics.csv"), index=False)
pd.DataFrame(roc_rows).to_csv(os.path.join(RES, "ct_incr_roc.csv"), index=False)
pd.DataFrame(dca_rows).to_csv(os.path.join(RES, "ct_incr_dca.csv"), index=False)
pd.DataFrame(cal_rows).to_csv(os.path.join(RES, "ct_incr_calib.csv"), index=False)
print("=== CT 增量指标 ===")
print(metf[["outcome", "base", "with_ct", "events", "AUC", "AUC_lo", "AUC_hi", "dAUC", "delong_p_vs_base", "NRI", "NRI_lo", "NRI_hi", "IDI", "Brier", "calib_slope"]].round(3).to_string(index=False))

# ---------- 亚组增量: 各表型 / 风险三分位 / SAP 内 MCTSI 的判别力与事件率 ----------
sub_rows = []
loc["risk6_tert"] = pd.qcut(loc.risk6, 3, labels=["T1 low", "T2 mid", "T3 high"]).astype(str)
loc["lowrisk"] = ((loc.pheno.isin([1, 2])) & (loc.SAP == 0)).astype(int)
for y, lab in OUTS:
    df = loc[loc[y].notna()]
    for gvar, gname in [("pheno", "Phenotype"), ("risk6_tert", "Risk6 tertile"), ("SAP", "SAP")]:
        for g, s in df.groupby(gvar):
            yv = s[y].astype(int).values
            if len(yv) < 20 or yv.sum() < 3 or yv.sum() == len(yv):
                auc, lo, hi = np.nan, np.nan, np.nan
            else:
                auc = roc_auc_score(yv, s.MCTSI)
                lo, hi = boot_ci(lambda yy, pp: roc_auc_score(yy, pp), yv, s.MCTSI.values.astype(float), B=500)
            hi_ct = s[s.MCTSI >= 6]; lo_ct = s[s.MCTSI < 6]
            sub_rows.append(dict(outcome=lab, group_var=gname, group=str(g) if gvar != "pheno" else f"P{g}",
                                 n=len(s), events=int(yv.sum()), rate=yv.mean(),
                                 auc_mctsi=auc, auc_lo=lo, auc_hi=hi,
                                 n_ct_high=len(hi_ct), ev_ct_high=int(hi_ct[y].sum()), rate_ct_high=hi_ct[y].mean() if len(hi_ct) else np.nan,
                                 n_ct_low=len(lo_ct), ev_ct_low=int(lo_ct[y].sum()), rate_ct_low=lo_ct[y].mean() if len(lo_ct) else np.nan))
subf = pd.DataFrame(sub_rows); subf.to_csv(os.path.join(RES, "ct_incr_subgroup.csv"), index=False)
print("\n=== 亚组内 MCTSI 判别 ===")
print(subf[["outcome", "group_var", "group", "n", "events", "rate", "auc_mctsi", "rate_ct_high", "rate_ct_low"]].round(3).to_string(index=False))

# ---------- 低风险层 (P1/P2 且非 SAP) 中的 CT 表现 ----------
rows = []
for stratum, mask in [("Low-risk (P1/P2, non-SAP)", loc.lowrisk == 1), ("P3/P4 or SAP", loc.lowrisk == 0),
                      ("P1 non-SAP", (loc.pheno == 1) & (loc.SAP == 0)), ("P2 non-SAP", (loc.pheno == 2) & (loc.SAP == 0)),
                      ("P3 non-SAP", (loc.pheno == 3) & (loc.SAP == 0)), ("P4 non-SAP", (loc.pheno == 4) & (loc.SAP == 0)),
                      ("SAP", loc.SAP == 1)]:
    s = loc[mask]
    r = dict(stratum=stratum, n=len(s), pct_of_cohort=len(s) / len(loc),
             mctsi_ge6_n=int((s.MCTSI >= 6).sum()), mctsi_ge6_pct=(s.MCTSI >= 6).mean(),
             mctsi_ge8_n=int((s.MCTSI >= 8).sum()))
    for y, lab in [("necroCT", "necroCT"), ("intervention", "intervention"), ("persOF", "persOF"), ("LOS14", "LOS14"), ("death", "death")]:
        r[f"{lab}_events"] = int(s[y].sum()); r[f"{lab}_rate"] = s[y].mean()
        r[f"{lab}_ev_ctlow"] = int(s.loc[s.MCTSI < 6, y].sum()); r[f"{lab}_ev_cthigh"] = int(s.loc[s.MCTSI >= 6, y].sum())
    rows.append(r)
lr = pd.DataFrame(rows); lr.to_csv(os.path.join(RES, "ct_lowrisk_table.csv"), index=False)
print("\n=== 分层 CT 表现 ===")
print(lr.round(3).to_string(index=False))

# ---------- MCTSI 分布按表型 ----------
md = loc.groupby(["pheno", "MCTSI"]).size().reset_index(name="n")
md["phenotype"] = "P" + md.pheno.astype(str)
md["frac"] = md.n / md.groupby("pheno").n.transform("sum")
md.to_csv(os.path.join(RES, "mctsi_dist.csv"), index=False)
print("\n已落盘: ct_incr_metrics / ct_incr_roc / ct_incr_dca / ct_incr_calib / ct_incr_subgroup / ct_lowrisk_table / mctsi_dist")
