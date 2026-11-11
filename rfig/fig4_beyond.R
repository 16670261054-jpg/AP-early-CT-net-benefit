source("00_config.R")
# Figure 4 表型不是早期严重度换皮: a 非重症亚组内梯度; b 超越 age+sex+etiology+早期严重度 的 LR 检验; c 六变量连续风险四分位
ns <- fread(file.path(RES, "rates_nonsevere.csv"))
ns <- ns[!(grepl("MIMIC", cohort) & outcome == "Coded infected necrosis")]
ns[, coh := fifelse(grepl("Local", cohort), L("China non-SAP", "中国 非SAP"), L("MIMIC non-SAP", "MIMIC 非SAP"))]
lab_map <- c(`CT necrosis` = L("CT necrosis", "CT坏死"), `Drainage/necrosectomy` = L("Drainage or necrosectomy", "引流/清创"), `LOS > 14 d` = L("LOS > 14 d", "住院>14天"),
             `Coded necrosis` = L("Coded necrosis", "编码坏死"), `Extra-pancreatic infection` = L("Extra-pancreatic infection", "胰外感染"),
             `Sepsis code` = L("Sepsis", "脓毒症"))
ns[, flab := str_wrap(paste0(coh, ": ", lab_map[outcome]), 12)]
ns[, flab := factor(flab, unique(flab))]
ns[, lab := sprintf("%.1f%%\n%d/%d", rate * 100, events, n)]
ns[, vj := fifelse(phenotype %in% c("P2", "P4"), -0.25, -1.75)]
ann <- ns[, .(ptxt = paste0("trend ", p_star(trend_p[1]))), by = flab]
pa <- ggplot(ns, aes(phenotype, rate * 100, fill = phenotype)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.25) +
  geom_errorbar(aes(ymin = lo * 100, ymax = hi * 100), width = 0.2, linewidth = 0.3) +
  geom_text(aes(y = hi * 100, label = lab, vjust = vj), size = 1.25, family = FONT, lineheight = 0.8) +
  geom_text(data = ann, aes(x = 0.55, y = Inf, label = ptxt), inherit.aes = FALSE, hjust = 0, vjust = 1.35, size = 1.7, family = FONT, colour = "grey25") +
  facet_wrap(~flab, scales = "free_y", nrow = 2) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = COL_PHENO, labels = PHENO_LAB(), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.68))) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = NULL, y = L("Event rate (%)", "事件率 (%)")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.5),
        strip.text = element_text(size = 5.1, lineheight = 0.95), strip.clip = "off")

# b LR 检验
lr <- fread(file.path(RES, "beyond_severity_lr.csv"))
lb <- lr[(cohort == "Local" & base_model == "age+sex+etiology+SAP" & outcome %in% c("CT necrosis", "Peripancreatic collection (CT)", "Drainage/necrosectomy", "LOS > 14 d")) |
         (cohort == "MIMIC-IV (ICD-10, necrosis specified)" & outcome %in% c("Coded necrosis", "Coded infected necrosis")) |
         (cohort == "MIMIC-IV" & outcome %in% c("Sepsis code", "Extra-pancreatic infection", "Pneumonia", "LOS > 14 d", "ICU admission")) |
         (cohort == "eICU" & outcome %in% c("Necrotizing (diagnosis string)", "Abscess/infected collection", "LOS > 14 d"))]
lb[, coh := fifelse(grepl("MIMIC", cohort), "MIMIC-IV", fifelse(cohort == "Local", L("China ward", "中国病房"), "eICU"))]
lb[, coh := factor(coh, c(L("China ward", "中国病房"), "MIMIC-IV", "eICU"))]
lab3 <- c(lab_map, `Peripancreatic collection (CT)` = L("Peripancreatic collection", "胰周积液"), `Coded infected necrosis` = L("Coded infected necrosis", "编码感染坏死"),
          Pneumonia = L("Pneumonia", "肺炎"), `ICU admission` = L("ICU admission", "入ICU"), `Necrotizing (diagnosis string)` = L("Necrotizing pancreatitis", "坏死性胰腺炎"),
          `Abscess/infected collection` = L("Abscess/infected collection", "脓肿/感染积液"))
lb[, out := lab3[outcome]]
lb[, out := factor(out, rev(unique(out)))]
lb[, sig := fifelse(p < 0.05, L("p < 0.05", "p<0.05"), L("p >= 0.05", "p>=0.05"))]
lb[, ptxt := p_star(p)]
pb <- ggplot(lb, aes(LR, out, fill = sig)) +
  geom_col(width = 0.65, colour = "grey25", linewidth = 0.25) +
  geom_vline(xintercept = qchisq(0.95, 3), linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_text(aes(label = ptxt), hjust = -0.1, size = 1.7, family = FONT) +
  facet_wrap(~coh, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("#B2182B", "grey80"), name = NULL, breaks = c(L("p < 0.05", "p<0.05"), L("p >= 0.05", "p>=0.05"))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(x = L("Likelihood-ratio statistic (3 df) for phenotype added to\nage + sex + etiology + early severity (eICU: APACHE IVa)",
             "表型加入 年龄+性别+病因+早期严重度 后的LR统计量(3自由度)"), y = NULL) +
  theme(legend.position = "bottom")

# c 连续六变量风险四分位
rc <- fread(file.path(RES, "risk_continuous.csv"))
rc <- rc[score == "risk6" & ((cohort == "Local" & outcome %in% c("CT necrosis", "Drainage/necrosectomy", "LOS > 14 d")) |
                              (cohort == "MIMIC-IV (ICD-10, necrosis specified)" & outcome == "Coded necrosis") |
                              (cohort == "MIMIC-IV" & outcome %in% c("Extra-pancreatic infection", "Sepsis code")) |
                              (cohort == "eICU" & outcome %in% c("Necrotizing (diagnosis string)", "LOS > 14 d")))]
rl <- melt(rc, id.vars = c("cohort", "outcome", "OR_perSD", "lo", "hi", "p"), measure.vars = patterns(rate = "^rate_Q", n = "^n_Q", ev = "^ev_Q"), variable.name = "Q")
rl[, Q := paste0("Q", Q)]
rl[, coh := fifelse(grepl("MIMIC", cohort), "MIMIC-IV", fifelse(cohort == "Local", L("China ward", "中国病房"), "eICU"))]
rl[, out := lab3[outcome]]
rl[, panel := paste0(coh, ": ", out)]
rl[, panel := factor(panel, unique(panel))]
ann2 <- unique(rl[, .(panel, txt = sprintf("OR/SD %.2f (%.2f-%.2f)", OR_perSD, lo, hi))])
pc <- ggplot(rl, aes(Q, rate * 100, group = 1)) +
  geom_line(colour = "grey50", linewidth = 0.35) +
  geom_point(aes(colour = Q), size = 1.6) +
  geom_text(data = ann2, aes(x = 0.6, y = Inf, label = txt), inherit.aes = FALSE, hjust = 0, vjust = 1.4, size = 1.65, family = FONT, colour = "grey25") +
  facet_wrap(~panel, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = c(Q1 = "#4393C3", Q2 = "#92C5DE", Q3 = "#F4A582", Q4 = "#B2182B"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(x = L("Quartile of six-variable continuous SAP risk", "六变量连续SAP风险四分位"), y = L("Event rate (%)", "事件率 (%)")) +
  theme(strip.text = element_text(size = 5.8))
p <- (pa | pb) / pc + plot_layout(heights = c(1, 0.85), widths = c(1.75, 1)) + plot_annotation(tag_levels = "a")
save_pub(p, "Figure4_beyond_severity", 190, 175)
