source("00_config.R")
# Figure 1 三队列流程图（数字读自 results/audit_*.csv）
al <- fread(file.path(RES, "audit_local.csv")); am <- fread(file.path(RES, "audit_mimic.csv")); ae <- fread(file.path(RES, "audit_eicu.csv"))
g <- function(d, k) d[item == k, value]
box <- function(x, y, w, h, lab, fill = "grey97", col = "grey30", size = 2.0) {
  list(annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2,
                fill = fill, colour = col, linewidth = 0.3),
       annotate("text", x = x, y = y, label = lab, size = size, family = FONT, lineheight = 0.95))
}
arrow_ <- function(x1, y1, x2, y2)
  annotate("segment", x = x1, y = y1, xend = x2, yend = y2, linewidth = 0.3,
           colour = "grey40", arrow = arrow(length = unit(1.4, "mm"), type = "closed"))

p <- ggplot() +
  # ---- Local ----
  box(0.17, 0.93, 0.30, 0.10, L(sprintf("Chinese surgical ward, 2017-2021\n1,386 adult AP admissions; 43 excluded\n(no admission date); n = %d, CT in all", g(al, "n_total")),
                                sprintf("中国外科病房队列\n%d例连续成人AP, 2017-2021\n全员入院CT", g(al, "n_total"))), fill = "#EAF2F8") +
  box(0.17, 0.76, 0.30, 0.10, L(sprintf("Phenotype P1-P4 assigned from six admission variables\n(derivation 2017-2020 n = %d; temporal 2021 n = %d)\n%d fully imputed flagged for sensitivity",
                                        g(al, "n_dev"), g(al, "n_val"), g(al, "n_allimp6")),
                                sprintf("六变量入院表型P1-P4\n(推导2017-2020 n=%d; 时间验证2021 n=%d)\n%d例全插补做敏感性", g(al, "n_dev"), g(al, "n_val"), g(al, "n_allimp6"))), fill = "#D6EAF8") +
  box(0.17, 0.565, 0.30, 0.13, L(sprintf("Outcomes\nCT necrosis %d | peripancreatic collection\nDrainage/necrosectomy %d | persistent OF %d\nLOS > 14 d %d | death %d",
                                         g(al, "necroCT"),
                                         g(al, "intervention_any"), g(al, "persOF"), g(al, "LOS_gt14"), g(al, "death")),
                                 sprintf("结局\nCT坏死 %d | 胰周积液\n引流/清创 %d | 持续OF %d\nLOS>14d %d | 死亡 %d",
                                         g(al, "necroCT"),
                                         g(al, "intervention_any"), g(al, "persOF"), g(al, "LOS_gt14"), g(al, "death"))), fill = "white", size = 1.85) +
  box(0.17, 0.40, 0.30, 0.07, L("Incremental value of admission CT (MCTSI)\nover phenotype / risk score", "入院CT(MCTSI)在表型/风险分之上的增量"), fill = "#FDEBD0") +
  # ---- MIMIC ----
  box(0.52, 0.93, 0.30, 0.10, L(sprintf("MIMIC-IV v3.1 (US, 2008-2019)\nAP in first 3 diagnosis positions\nn = %d adults, first admission", g(am, "n_total")),
                                sprintf("MIMIC-IV v3.1 (美国)\n前3诊断位AP\nn=%d 成人首次住院", g(am, "n_total"))), fill = "#EAF2F8") +
  box(0.52, 0.76, 0.30, 0.10, L(sprintf("Full cohort n = %d (1 non-positive LOS excluded)\nExtra-pancreatic infection codes %d | sepsis %d\nLOS > 14 d | in-hospital death %d | 90-day death",
                                        g(am, "n_total") - g(am, "n_los_nonpositive_excluded"), g(am, "epi_any_cohort"), g(am, "sepsis_cohort"), g(am, "death_cohort")),
                                sprintf("全队列 n=%d (剔1例LOS非正)\n胰外感染编码 %d | 脓毒症 %d\nLOS>14d | 院内死亡 %d | 90天死亡",
                                        g(am, "n_total") - g(am, "n_los_nonpositive_excluded"), g(am, "epi_any_cohort"), g(am, "sepsis_cohort"), g(am, "death_cohort"))), fill = "#D6EAF8") +
  box(0.52, 0.565, 0.30, 0.13, L(sprintf("ICD-10 subset n = %d\nNecrosis status specified (5th character) n = %d\nCoded necrosis %d | coded infected necrosis %d\nOrgan support >= 48 h %d",
                                         g(am, "n_icd10_subset"), g(am, "n_icd10_necro_specified"), g(am, "necro_coded"), g(am, "ipn_coded"), g(am, "support_persistent")),
                                 sprintf("ICD-10子集 n=%d\n坏死状态已编码(第5位) n=%d\n编码坏死 %d | 编码感染坏死 %d\n器官支持>=48h %d",
                                         g(am, "n_icd10_subset"), g(am, "n_icd10_necro_specified"), g(am, "necro_coded"), g(am, "ipn_coded"), g(am, "support_persistent"))), fill = "white", size = 1.85) +
  # ---- eICU ----
  box(0.85, 0.93, 0.28, 0.10, L(sprintf("eICU-CRD v2.0 (US, 2014-2015)\nAP by diagnosis table or APACHE dx\nn = %d, %d hospitals", g(ae, "n_total"), g(ae, "hospitals")),
                                sprintf("eICU-CRD v2.0 (美国)\n诊断表/APACHE入科诊断AP\nn=%d, %d家医院", g(ae, "n_total"), g(ae, "hospitals"))), fill = "#EAF2F8") +
  box(0.85, 0.76, 0.28, 0.10, L(sprintf("Directional outcomes (diagnosis strings)\nNecrotizing %d | abscess/infected collection %d\nExtra-pancreatic infection %d | drainage %d\nLOS > 14 d %d | death %d",
                                        g(ae, "necro_dx"), g(ae, "abscess_dx"), g(ae, "epi_any"), g(ae, "drain_dx"), g(ae, "LOS_gt14"), g(ae, "death")),
                                sprintf("方向性结局(诊断串)\n坏死 %d | 脓肿/感染积液 %d\n胰外感染 %d | 引流 %d\nLOS>14d %d | 死亡 %d",
                                        g(ae, "necro_dx"), g(ae, "abscess_dx"), g(ae, "epi_any"), g(ae, "drain_dx"), g(ae, "LOS_gt14"), g(ae, "death"))), fill = "#D6EAF8", size = 1.85) +
  # ---- bottom ----
  box(0.52, 0.27, 0.66, 0.10, L(paste0("Phenotype assigned by nearest derivation centroid in all cohorts\n",
                                       "Outcome gradients | adjustment for age, sex, etiology and early severity | likelihood-ratio tests beyond severity\n",
                                       "Continuous six-variable risk as parallel exposure | competing-risk cumulative incidence of death | E-values"),
                                paste0("外部队列按最近推导质心赋表型\n结局梯度 | 年龄/性别/病因/早期严重度调整 | 超越严重度的LR检验\n",
                                       "六变量连续风险并行 | 死亡累积发生(竞争) | E值")), fill = "#FDEBD0", size = 2.0) +
  arrow_(0.17, 0.88, 0.17, 0.81) + arrow_(0.17, 0.71, 0.17, 0.63) + arrow_(0.17, 0.50, 0.17, 0.435) +
  arrow_(0.52, 0.88, 0.52, 0.81) + arrow_(0.52, 0.71, 0.52, 0.63) +
  arrow_(0.85, 0.88, 0.85, 0.81) +
  arrow_(0.17, 0.365, 0.30, 0.32) + arrow_(0.52, 0.50, 0.52, 0.32) + arrow_(0.85, 0.71, 0.72, 0.32) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.20, 0.99)) + theme_void()
save_pub(p, "Figure1_flow", 190, 125)
