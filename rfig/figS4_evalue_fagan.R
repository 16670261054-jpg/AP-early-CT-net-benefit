source("00_config.R")
# Figure S4 未测量混杂稳健性与床旁验后概率: a E-value 棒棒糖; b P4 vs 非P4 验前->验后概率
ev <- fread(file.path(RES, "evalue.csv"))[!(outcome %in% c("Percutaneous drainage")) & !(cohort == "Local" & outcome == "ICU admission")]
ev[, coh := fifelse(grepl("MIMIC", cohort), "MIMIC-IV", fifelse(cohort == "Local", L("China ward", "中国病房"), "eICU"))]
ev[, coh := factor(coh, c(L("China ward", "中国病房"), "MIMIC-IV", "eICU"))]
ev[, out2 := fifelse(outcome == "Sepsis code", "Sepsis", outcome)]
setorder(ev, coh, evalue_point)
ev[, out2 := factor(out2, unique(out2))]
pa <- ggplot(ev, aes(evalue_point, out2)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey45") +
  geom_segment(aes(x = 1, xend = evalue_point, yend = out2), linewidth = 0.35, colour = "grey60") +
  geom_point(aes(x = evalue_ci), shape = 21, size = 1.3, stroke = 0.5, colour = "#4393C3") +
  geom_point(size = 1.6, colour = "#B2182B") +
  facet_wrap(~coh, scales = "free_y", ncol = 1) +
  labs(x = L("E-value for P4 vs P1 association (filled = point estimate, open = CI bound)",
             "P4 vs P1 关联的E值 (实心=点估计, 空心=CI界)"), y = NULL) +
  theme(axis.text.y = element_text(size = 5.5), strip.text = element_text(size = 6))

# b Fagan 风格: 验前 -> P4验后 / 非P4验后
fg <- fread(file.path(RES, "lr_fagan.csv"))
fg[, coh := fifelse(grepl("MIMIC", cohort), "MIMIC-IV", fifelse(cohort == "Local", L("China ward", "中国病房"), "eICU"))]
fg[, out2 := fifelse(outcome == "Sepsis code", "Sepsis", outcome)]
fg[, row := factor(paste0(coh, ": ", out2), rev(unique(paste0(coh, ": ", out2))))]
fl <- melt(fg, id.vars = c("row", "LR_pos", "LR_neg"), measure.vars = c("pretest", "posttest_P4", "posttest_nonP4"),
           variable.name = "what", value.name = "prob")
fl[, what := factor(what, c("posttest_nonP4", "pretest", "posttest_P4"),
                    labels = c(L("Post-test, non-P4", "验后(非P4)"), L("Pre-test", "验前"), L("Post-test, P4", "验后(P4)")))]
pb <- ggplot(fl, aes(prob * 100, row)) +
  geom_line(aes(group = row), colour = "grey60", linewidth = 0.35) +
  geom_point(aes(colour = what), size = 1.7) +
  geom_text(data = fg, aes(x = 100, y = row, label = sprintf("LR+ %.1f / LR- %.2f", LR_pos, LR_neg)),
            hjust = 1, size = 1.5, family = FONT, colour = "grey30") +
  scale_colour_manual(values = c("#4393C3", "grey40", "#B2182B"), name = NULL) +
  scale_x_continuous(limits = c(0, 102), expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = L("Probability of outcome, % (P4 as bedside test)", "结局概率, % (P4作为床旁检验)"), y = NULL) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 5.5), legend.text = element_text(size = 5.5))
p <- pa | pb
p <- p + plot_layout(widths = c(1, 1.15)) + plot_annotation(tag_levels = "a")
save_pub(p, "FigureS4_evalue_fagan", 190, 130)
