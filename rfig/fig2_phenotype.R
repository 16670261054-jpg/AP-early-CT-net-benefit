source("00_config.R")
# Figure 2 表型画像与跨库可迁移性 (ML内核显性化):
# a 六变量 z-score 画像热图 (三库, k-means质心分配 + XGBoost分类器一致率)
# b 六变量连续SAP风险(GBDT)分布 -- 表型捕获但不等同于监督严重度风险
# c 病因构成; d 入院MCTSI分布 (本院)
pz <- fread(file.path(RES, "pheno_profile_z.csv"))
pz[, cohort := factor(cohort, c("Local", "MIMIC-IV", "eICU"), labels = c(L("China ward cohort", "中国病房队列"), "MIMIC-IV", "eICU"))]
vord <- c("Age", "Calcium", "Lymphocyte %", "RBC", "Creatinine", "ALT")
vlab <- c(Age = L("Age", "年龄"), Calcium = L("Calcium", "血钙"), `Lymphocyte %` = L("Lymphocyte %", "淋巴%"),
          RBC = L("RBC", "红细胞"), Creatinine = L("Creatinine", "肌酐"), ALT = "ALT")
pz[, var := factor(vlab[variable], vlab[vord])]
ag <- fread(file.path(RES, "clf_agreement_local.csv"))
pa <- ggplot(pz, aes(var, phenotype, fill = pmax(pmin(z_mean, 1.5), -1.5))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", z_mean)), size = 1.85, family = FONT) +
  facet_wrap(~cohort, ncol = 1) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1.5, 1.5),
                       name = L("Mean z-score\n(within cohort)", "组内z均值")) +
  labs(x = NULL, y = NULL,
       subtitle = sprintf(L("Nearest-centroid assignment;\nXGBoost classifier agreement %.1f%% (n = %s)",
                            "最近质心分配; XGBoost分类器一致率 %.1f%% (n=%s)"),
                          ag$agreement[1] * 100, format(ag$n[1], big.mark = ","))) +
  theme(legend.position = "right", legend.key.height = unit(4, "mm"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.subtitle = element_text(size = 5.3, colour = "grey30", lineheight = 1))

# b risk6 分布 (雨云)
rk <- fread(file.path(RES, "risk6_by_pheno_long.csv"))
rk[, cohort := factor(cohort, c("Local", "MIMIC-IV", "eICU"), labels = levels(pz$cohort))]
rk[, phenotype := factor(phenotype)]
kw <- fread(file.path(RES, "risk6_by_pheno_kw.csv"))
kw[, cohort := factor(cohort, c("Local", "MIMIC-IV", "eICU"), labels = levels(pz$cohort))]
kw[, ptxt := paste0("Kruskal-Wallis ", p_star(kruskal_p))]
set.seed(SEED)
rk[, xj := as.numeric(phenotype) - 0.18 + runif(.N, -0.07, 0.07)]
pb <- ggplot(rk, aes(phenotype, risk6, fill = phenotype)) +
  stat_slab(aes(thickness = after_stat(pdf)), scale = 0.6, alpha = 0.85, position = position_nudge(x = 0.15), colour = NA, density = "bounded") +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.28, position = position_nudge(x = 0.01), alpha = 0.9) +
  geom_point(aes(x = xj), shape = 21, size = 0.25, alpha = 0.25, stroke = 0.1) +
  geom_text(data = kw, aes(x = 0.55, y = Inf, label = ptxt), inherit.aes = FALSE, hjust = 0, vjust = 1.4, size = 1.7, family = FONT, colour = "grey25") +
  facet_wrap(~cohort, nrow = 1) +
  scale_fill_manual(values = COL_PHENO, labels = PHENO_LAB(), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  labs(x = NULL, y = L("Six-variable continuous SAP risk (GBDT)", "六变量连续SAP风险 (GBDT)")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.5))

# c 病因构成
et <- fread(file.path(RES, "etiology_by_pheno.csv"))
et[, coh := factor(cohort, c("Local", "MIMIC-IV (ICD-10)"), labels = c(L("China ward", "中国病房"), "MIMIC-IV"))]
eord <- c("Biliary", "Hypertriglyceridemic", "Alcoholic", "Drug", "Idiopathic", "Other", "Unspecified")
elab <- c(Biliary = L("Biliary", "胆源性"), Hypertriglyceridemic = L("Hypertriglyceridemic", "高甘油三酯"),
          Alcoholic = L("Alcoholic", "酒精性"), Drug = L("Drug", "药物性"), Idiopathic = L("Idiopathic", "特发性"),
          Other = L("Other", "其他"), Unspecified = L("Unspecified", "未指明"))
et[, etio := factor(elab[etiology], elab[eord])]
COL_ETIO <- c("#4393C3", "#F4A582", "#B2182B", "#92C5DE", "grey70", "grey50", "grey85")
names(COL_ETIO) <- elab[eord]
pc <- ggplot(et, aes(phenotype, frac * 100, fill = etio)) +
  geom_col(width = 0.7, colour = "grey25", linewidth = 0.2) +
  facet_wrap(~coh, nrow = 1) +
  scale_fill_manual(values = COL_ETIO, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  guides(fill = guide_legend(nrow = 3)) +
  labs(x = NULL, y = L("Etiology, %", "病因构成, %")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.5))

# d 入院MCTSI分布 (本院)
mc <- fread(file.path(RES, "mctsi_dist.csv"))
mc[, MCTSI := factor(MCTSI, c(0, 2, 4, 6, 8, 10))]
COL_MCTSI <- c(`0` = "#DEEBF7", `2` = "#9ECAE1", `4` = "#4393C3", `6` = "#F4A582", `8` = "#D6604D", `10` = "#B2182B")
pd <- ggplot(mc, aes(phenotype, frac * 100, fill = MCTSI)) +
  geom_col(width = 0.7, colour = "grey25", linewidth = 0.2) +
  geom_text(data = mc[frac > 0.05], aes(label = n), position = position_stack(vjust = 0.5), size = 1.7, family = FONT) +
  scale_fill_manual(values = COL_MCTSI, name = L("Admission\nMCTSI", "入院\nMCTSI")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = L("Admission MCTSI, % (China ward cohort)", "入院MCTSI构成, % (中国病房队列)")) +
  theme(legend.position = "right", legend.key.height = unit(3, "mm"), legend.text = element_text(size = 5.5))
p <- pa | (pb / (pc | pd))
p <- p + plot_layout(widths = c(0.52, 1)) + plot_annotation(tag_levels = "a")
save_pub(p, "Figure2_phenotype_profile", 190, 150)
