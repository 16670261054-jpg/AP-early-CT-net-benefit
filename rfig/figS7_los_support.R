source("00_config.R")
# Figure 7 住院时长雨云 (三库) + MIMIC 器官支持一过性/持续 + 长期死亡
lo <- fread(file.path(RES, "los_long.csv"))
lo[, cohort := factor(cohort, c("Local", "MIMIC-IV", "eICU"), labels = c(L("China ward cohort", "中国病房队列"), "MIMIC-IV", "eICU"))]
lo[, phenotype := factor(phenotype)]
lo[, los_c := pmin(los, 60)]
set.seed(SEED)
lo[, xj := as.numeric(phenotype) - 0.18 + runif(.N, -0.07, 0.07)]
ls <- fread(file.path(RES, "los_summary.csv"))
ls[, cohort := factor(cohort, c("Local", "MIMIC-IV", "eICU"), labels = levels(lo$cohort))]
ls[, txt := sprintf("%.1f", median)]
ann <- unique(ls[, .(cohort, kruskal_p)])
ann[, ptxt := paste0("Kruskal-Wallis ", p_star(kruskal_p))]
pa <- ggplot(lo, aes(phenotype, los_c, fill = phenotype)) +
  stat_slab(aes(thickness = after_stat(pdf)), scale = 0.6, alpha = 0.85, position = position_nudge(x = 0.15), colour = NA, density = "bounded") +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.28, position = position_nudge(x = 0.01), alpha = 0.9) +
  geom_point(aes(x = xj), shape = 21, size = 0.25, alpha = 0.3, stroke = 0.1) +
  geom_text(data = ls, aes(phenotype, 58, label = txt), inherit.aes = FALSE, size = 1.8, family = FONT) +
  geom_text(data = ann, aes(x = 0.55, y = Inf, label = ptxt), inherit.aes = FALSE, hjust = 0, vjust = 1.4, size = 1.75, family = FONT, colour = "grey25") +
  geom_hline(yintercept = 14, linetype = 3, linewidth = 0.3, colour = "grey45") +
  facet_wrap(~cohort, nrow = 1) +
  scale_fill_manual(values = COL_PHENO, labels = PHENO_LAB(), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = L("Hospital length of stay, days (capped at 60; median shown)", "住院天数 (截断60; 标注中位数)")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6))

# b 器官支持
su <- fread(file.path(RES, "support_duration.csv"))
su[, category := factor(category, c("Transient (<48 h)", "Persistent (>=48 h)"), labels = c(L("Transient (< 48 h)", "一过性(<48h)"), L("Persistent (≥ 48 h)", "持续(≥48h)")))]
pb <- ggplot(su, aes(phenotype, pct_of_pheno * 100, fill = category)) +
  geom_col(width = 0.7, colour = "grey25", linewidth = 0.2) +
  geom_text(aes(label = ifelse(n > 0, n, "")), position = position_stack(vjust = 0.5), size = 1.8, family = FONT) +
  scale_fill_manual(values = c("#92C5DE", "#B2182B"), name = L("Ventilation or vasopressor", "通气或升压药")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = L("Patients with organ support, % (MIMIC-IV)", "器官支持患者, % (MIMIC-IV)")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6))

# c 长期死亡 (MIMIC, 粗率)
dl <- fread(file.path(RES, "death_longterm.csv"))
dl[, outcome := factor(outcome, c("In-hospital death", "30-day death", "90-day death", "1-year death"),
                       labels = c(L("In-hospital", "院内"), L("30-day", "30天"), L("90-day", "90天"), L("1-year", "1年")))]
pc <- ggplot(dl, aes(outcome, rate * 100, colour = phenotype, group = phenotype)) +
  geom_line(linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo * 100, ymax = hi * 100), width = 0.15, linewidth = 0.3) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = COL_PHENO, labels = PHENO_LAB(), name = NULL) +
  labs(x = NULL, y = L("Mortality, % (MIMIC-IV, unadjusted)", "死亡率, % (MIMIC-IV, 未调整)")) +
  theme(legend.position = "none")
p <- pa / (pb | pc) + plot_layout(heights = c(1, 0.8)) + plot_annotation(tag_levels = "a")
save_pub(p, "FigureS7_LOS_support", 190, 150)
