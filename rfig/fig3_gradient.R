source("00_config.R")
# Figure 3 坏死-感染进程的三库表型梯度 (a 点-CI 网格; b 调整OR森林)
# 核心结论: 从P1到P4, 坏死/感染坏死/局部干预/长住院在三个独立队列中单调上升
r <- fread(file.path(RES, "rates_by_pheno.csv"))
sel <- rbind(
  r[cohort == "Local" & outcome %in% c("CT necrosis", "Drainage/necrosectomy", "Persistent OF", "LOS > 14 d")],
  r[cohort == "MIMIC-IV (ICD-10, necrosis specified)" & outcome %in% c("Coded necrosis", "Coded infected necrosis")],
  r[cohort == "MIMIC-IV" & outcome %in% c("Sepsis code", "Extra-pancreatic infection", "LOS > 14 d")],
  r[cohort == "eICU" & outcome %in% c("Necrotizing (diagnosis string)", "Abscess/infected collection", "Sepsis", "Extra-pancreatic infection", "LOS > 14 d")])
sel[, coh2 := fifelse(grepl("MIMIC", cohort), "MIMIC-IV", fifelse(cohort == "Local", L("China ward", "中国病房"), "eICU"))]
lab_map <- c(`CT necrosis` = L("CT necrosis", "CT坏死"),
             `Drainage/necrosectomy` = L("Drainage or necrosectomy", "引流或清创"),
             `Persistent OF` = L("Persistent OF", "持续器官衰竭"),
             `LOS > 14 d` = L("LOS > 14 d", "住院>14天"),
             `Coded necrosis` = L("Coded necrosis (ICD-10)", "编码坏死(ICD-10)"),
             `Coded infected necrosis` = L("Coded infected necrosis", "编码感染坏死"),
             `Sepsis code` = L("Sepsis", "脓毒症"), `Sepsis` = L("Sepsis", "脓毒症"),
             `Extra-pancreatic infection` = L("Extra-pancreatic infection", "胰外感染"),
             `Necrotizing (diagnosis string)` = L("Necrotizing pancreatitis", "坏死性胰腺炎"),
             `Abscess/infected collection` = L("Abscess / infected collection", "脓肿/感染积液"))
sel[, flab := str_wrap(paste0(coh2, ": ", lab_map[outcome]), 18)]
sel[, flab := factor(flab, unique(flab))]
sel[, ptxt := p_star(trend_p)]
ann <- sel[, .(ptxt = ptxt[1]), by = flab]
pa <- ggplot(sel, aes(phenotype, rate * 100, colour = phenotype)) +
  geom_line(aes(group = 1), colour = "grey55", linewidth = 0.35) +
  geom_errorbar(aes(ymin = lo * 100, ymax = hi * 100), width = 0.18, linewidth = 0.35) +
  geom_point(size = 1.7) +
  geom_text(data = ann, aes(x = 0.6, y = Inf, label = paste0(L("trend ", "趋势 "), ptxt)), inherit.aes = FALSE,
            hjust = 0, vjust = 1.4, size = 1.75, family = FONT, colour = "grey25") +
  facet_wrap(~flab, scales = "free_y", nrow = 3) +
  scale_colour_manual(values = COL_PHENO, labels = PHENO_LAB(), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
  labs(x = NULL, y = L("Event rate, % (95% Wilson CI)", "事件率, % (95% Wilson CI)")) +
  theme(legend.position = "bottom", strip.text = element_text(size = 5.8, face = "bold", lineheight = 0.95),
        axis.text.x = element_text(size = 6))

# b 调整OR森林 (age+sex+etiology+early severity; eICU: age+sex+APACHE)
o <- fread(file.path(RES, "adj_or.csv"))
ob <- rbind(o[cohort == "Local" & model == "age_sex_etio_SAP" & outcome %in% c("CT necrosis", "Drainage/necrosectomy", "LOS > 14 d", "Peripancreatic collection (CT)")],
            o[cohort == "MIMIC-IV (ICD-10, necrosis specified)" & model == "age_sex_etio_SAP" & outcome == "Coded necrosis"],
            o[cohort == "MIMIC-IV" & model == "age_sex_etio_SAP" & outcome %in% c("Sepsis code", "Extra-pancreatic infection", "LOS > 14 d")],
            o[cohort == "eICU" & model == "age_sex_apache" & outcome %in% c("Necrotizing (diagnosis string)", "LOS > 14 d")])
ob[, coh := fifelse(grepl("MIMIC", cohort), "MIMIC-IV", cohort)]
ob[, coh := factor(coh, c("Local", "MIMIC-IV", "eICU"), labels = c(L("China ward cohort", "中国病房队列"), "MIMIC-IV", "eICU"))]
lab2 <- c(lab_map, `Peripancreatic collection (CT)` = L("Peripancreatic collection (CT)", "胰周积液(CT)"))
ob[, out := factor(lab2[outcome], rev(unique(lab2[outcome])))]
ob[, pheno := factor(pheno, c("P2", "P3", "P4"))]
pb <- ggplot(ob, aes(OR, out, colour = pheno)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4, position = position_dodge(0.6)) +
  geom_point(size = 1.5, position = position_dodge(0.6)) +
  facet_wrap(~coh, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = COL_PHENO[2:4], labels = PHENO_LAB()[2:4], name = L("vs P1", "对照P1")) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16), labels = c("0.25", "0.5", "1", "2", "4", "8", "16")) +
  coord_cartesian(xlim = c(0.2, 20)) +
  labs(x = L("Adjusted OR vs P1 (age, sex, etiology, early severity; eICU: age, sex, APACHE IVa), log scale",
             "调整OR vs P1 (年龄/性别/病因/早期严重度; eICU: 年龄/性别/APACHE), 对数轴"), y = NULL) +
  theme(legend.position = "right", legend.text = element_text(size = 6))
p <- pa / pb + plot_layout(heights = c(1.55, 0.7)) + plot_annotation(tag_levels = "a")
save_pub(p, "Figure3_gradient", 190, 185)
