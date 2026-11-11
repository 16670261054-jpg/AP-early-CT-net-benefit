source("00_config.R")
# Figure 5 入院 CT (MCTSI) 在化验模型之上的增量 (本院): a ROC; b dAUC/NRI/IDI; c DCA; d 表型 x MCTSI 干预率 (de-implementation)
mod_lab <- c(pheno = "Phenotype", `pheno+CT` = "Phenotype + CT", risk6 = "Risk score", `risk6+CT` = "Risk score + CT", mctsi_only = "CT (MCTSI) alone")
roc <- fread(file.path(RES, "ct_incr_roc.csv"))
roc <- roc[outcome %in% c("Drainage/necrosectomy", "Persistent OF")]
roc[, m := factor(mod_lab[model], mod_lab)]
met <- fread(file.path(RES, "ct_incr_metrics.csv"))
met[, m := factor(mod_lab[fifelse(with_ct == 1 & base %in% c("pheno", "risk6"), paste0(base, "+CT"), base)], mod_lab)]
auc_txt <- met[!is.na(m) & outcome %in% c("Drainage/necrosectomy", "Persistent OF"), .(outcome, m, AUC, AUC_lo, AUC_hi)]
auc_txt[, txt := sprintf("%s %.2f (%.2f-%.2f)", m, AUC, AUC_lo, AUC_hi)]
auc_txt <- auc_txt[order(outcome, m)]
auc_txt[, y := 0.30 - 0.055 * (as.integer(m) - 1)]
out_lab <- c(`Drainage/necrosectomy` = L("Drainage or necrosectomy\n(n = 1,343)", "引流或清创\n(n=1,343)"),
             `Persistent OF` = L("Persistent organ failure\n(n = 1,343)", "持续器官衰竭\n(n=1,343)"))
pa <- ggplot(roc, aes(fpr, tpr, colour = m)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3, colour = "grey60") +
  geom_step(linewidth = 0.5) +
  geom_text(data = auc_txt, aes(x = 0.98, y = y, label = txt, colour = m), hjust = 1, size = 1.65, family = FONT, show.legend = FALSE) +
  facet_wrap(~outcome, labeller = labeller(outcome = out_lab)) +
  scale_colour_manual(values = COL_MODEL, name = NULL) +
  coord_equal() +
  scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = L("1 - specificity", "1-特异度"), y = L("Sensitivity", "敏感度")) +
  theme(legend.position = "none", strip.text = element_text(size = 6, lineheight = 0.95), panel.spacing.x = unit(4, "mm"))

# b 增量指标
mb <- met[with_ct == 1 & base %in% c("pheno", "risk6")]
mb[, base_lab := fifelse(base == "pheno", L("+ CT over phenotype", "CT加于表型"), L("+ CT over risk score", "CT加于风险分"))]
mb[, out := factor(outcome, c("Drainage/necrosectomy", "Persistent OF", "LOS > 14 d"),
                   labels = c(L("Drainage/necrosectomy", "引流/清创"), L("Persistent OF", "持续OF"), L("LOS > 14 d", "住院>14天")))]
ml <- rbind(mb[, .(out, base_lab, metric = "dAUC", est = dAUC, lo = NA_real_, hi = NA_real_, p = delong_p_vs_base)],
            mb[, .(out, base_lab, metric = "NRI", est = NRI, lo = NRI_lo, hi = NRI_hi, p = NA_real_)],
            mb[, .(out, base_lab, metric = "IDI", est = IDI, lo = IDI_lo, hi = IDI_hi, p = NA_real_)])
ml[, metric := factor(metric, c("dAUC", "NRI", "IDI"), labels = c("ΔAUC", L("NRI", "NRI"), "IDI"))]
ml[, txt := fifelse(!is.na(p), sprintf("%.3f %s", est, p_star(p)), sprintf("%.3f", est))]
pb <- ggplot(ml, aes(est, out, colour = base_lab)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4, position = position_dodge(0.6), na.rm = TRUE) +
  geom_point(size = 1.5, position = position_dodge(0.6)) +
  geom_text(aes(label = txt), position = position_dodge(0.6), vjust = -0.55, size = 1.45, family = FONT, lineheight = 0.8, show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = COL_TWO, name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.18), breaks = scales::breaks_extended(n = 3)) +
  scale_y_discrete(limits = rev, expand = expansion(add = c(0.5, 0.9))) +
  labs(x = L("Increment from adding MCTSI (10x5-fold out-of-fold; 95% bootstrap CI)", "加入MCTSI的增量 (10x5折袋外; 95% bootstrap CI)"), y = NULL) +
  theme(legend.position = "bottom", panel.spacing.x = unit(4, "mm"))

# c DCA
dca <- fread(file.path(RES, "ct_incr_dca.csv"))
dca <- dca[outcome %in% c("Drainage/necrosectomy", "Persistent OF") & model %in% c("pheno", "pheno+CT", "risk6", "risk6+CT")]
dca[, m := factor(mod_lab[model], mod_lab)]
ta <- unique(dca[, .(outcome, threshold, treat_all)])
pc <- ggplot(dca, aes(threshold * 100, net_benefit, colour = m)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_line(data = ta, mapping = aes(threshold * 100, treat_all), colour = "grey60", linetype = 3, linewidth = 0.4, inherit.aes = FALSE) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~outcome, scales = "free_y", labeller = labeller(outcome = out_lab)) +
  scale_colour_manual(values = COL_MODEL, name = NULL) +
  coord_cartesian(ylim = c(-0.02, NA)) +
  labs(x = L("Threshold probability (%)", "阈值概率 (%)"), y = L("Net benefit", "净获益")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6), strip.text = element_text(size = 6, lineheight = 0.95))

# d 表型 x MCTSI 干预率 tile (de-implementation)
pm <- fread(file.path(RES, "pheno_mctsi_outcome.csv"))
pm[, mctsi := factor(mctsi, c("0-2", "4", "6", "8-10"))]
pm[, txt := sprintf("%.1f%%\n%d/%d", intervention_rate * 100, intervention, n)]
pd <- ggplot(pm, aes(mctsi, phenotype, fill = intervention_rate * 100)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = txt), size = 1.7, family = FONT, lineheight = 0.85) +
  scale_fill_gradient(low = "#F7FBFF", high = "#B2182B", name = L("Drainage or\nnecrosectomy (%)", "引流/清创 (%)")) +
  scale_y_discrete(limits = rev) +
  labs(x = L("Admission MCTSI", "入院MCTSI"), y = NULL) +
  theme(legend.position = "right", legend.key.height = unit(4, "mm"))
p <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1, 0.95), widths = c(1, 1)) + plot_annotation(tag_levels = "a")
save_pub(p, "Figure5_CT_increment", 190, 165)
