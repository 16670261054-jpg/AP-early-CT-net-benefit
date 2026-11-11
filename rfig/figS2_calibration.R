source("00_config.R")
# Figure S2 校准 (10x5折袋外预测, 按预测概率十分位): 观察率 vs 预测均值
cal <- fread(file.path(RES, "ct_incr_calib.csv"))
mod_lab <- c(pheno = "Phenotype", `pheno+CT` = "Phenotype + CT", `risk6+CT` = "Risk score + CT")
cal[, m := factor(mod_lab[model], mod_lab)]
out_ord <- c("Drainage/necrosectomy", "Persistent OF", "LOS > 14 d")
out_lab <- c(L("Drainage or necrosectomy", "引流或清创"),
             L("Persistent OF", "持续器官衰竭"), L("LOS > 14 d", "住院>14天"))
cal[, out := factor(outcome, out_ord, labels = out_lab)]
p <- ggplot(cal, aes(pred_mean, obs_rate, colour = m)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3, colour = "grey60") +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_point(aes(size = n), shape = 21, stroke = 0.5, fill = NA) +
  facet_wrap(~out, scales = "free", nrow = 1) +
  scale_colour_manual(values = COL_MODEL[names(COL_MODEL) %in% mod_lab], name = NULL) +
  scale_size_area(max_size = 2.2, guide = "none") +
  labs(x = L("Mean predicted probability (decile of out-of-fold prediction)", "预测概率均值 (袋外预测十分位)"),
       y = L("Observed event rate", "观察事件率")) +
  theme(legend.position = "bottom", strip.text = element_text(size = 6))
save_pub(p, "FigureS2_calibration", 170, 75)
