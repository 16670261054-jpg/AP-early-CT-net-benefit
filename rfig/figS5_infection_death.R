source("00_config.R")
# Figure 6 感染部位与死亡时相: a 胰外感染部位 x 表型 (MIMIC); b 胰外感染 vs 编码感染坏死 (MIMIC ICD-10 坏死已指明);
#          c 院内死亡累积发生 (竞争: 出院) 按表型; d 死亡时相 (MIMIC)
r <- fread(file.path(RES, "rates_by_pheno.csv"))
sites <- c("Pneumonia", "Bacteremia", "UTI", "Cholangitis", "C. difficile")
ra <- r[cohort == "MIMIC-IV" & outcome %in% sites]
ra[, out := factor(outcome, sites, labels = c(L("Pneumonia", "肺炎"), L("Bacteremia", "菌血症"), L("Urinary tract", "泌尿道"), L("Cholangitis", "胆管炎"), L("C. difficile", "艰难梭菌")))]
pa <- ggplot(ra, aes(out, rate * 100, fill = phenotype)) +
  geom_col(position = position_dodge(0.75), width = 0.7, colour = "grey25", linewidth = 0.2) +
  geom_errorbar(aes(ymin = lo * 100, ymax = hi * 100), position = position_dodge(0.75), width = 0.25, linewidth = 0.3) +
  scale_fill_manual(values = COL_PHENO, labels = PHENO_LAB(), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = NULL, y = L("Coded infection, % (MIMIC-IV, n = 4,439)", "编码感染, % (MIMIC-IV, n=4,439)"), title = L("Extra-pancreatic infection sites", "胰外感染部位")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.5))

# b
ev <- fread(file.path(RES, "epi_vs_ipn.csv"))[phenotype != "All"]
el <- melt(ev, id.vars = c("phenotype", "n"), measure.vars = c("epi_only", "both", "ipn_only"), variable.name = "cat", value.name = "k")
el[, cat := factor(cat, c("ipn_only", "both", "epi_only"), labels = c(L("Infected necrosis only", "仅感染坏死"), L("Both", "两者"), L("Extra-pancreatic infection only", "仅胰外感染")))]
el[, pct := k / n * 100]
pb <- ggplot(el, aes(phenotype, pct, fill = cat)) +
  geom_col(width = 0.7, colour = "grey25", linewidth = 0.2) +
  geom_text(data = ev, aes(phenotype, (epi + ipn - both) / n * 100, label = sprintf("%d/%d", epi + ipn - both, n)), inherit.aes = FALSE, vjust = -0.3, size = 1.8, family = FONT) +
  scale_fill_manual(values = c("#B2182B", "#F4A582", "#92C5DE"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  guides(fill = guide_legend(nrow = 3)) +
  labs(x = NULL, y = L("Patients with any coded infection, %", "任何编码感染, %"), title = L("Infected necrosis vs extra-pancreatic infection\n(MIMIC-IV ICD-10, necrosis status specified, n = 1,649)", "感染坏死 vs 胰外感染 (MIMIC-IV ICD-10, n=1,649)")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.5))

# c CIF
cif <- fread(file.path(RES, "cif_death.csv"))
cif <- cif[time <= 60]
cp <- fread(file.path(RES, "cif_death_28d.csv"))[day == 28]
cp[, txt := sprintf("%s 28 d: %.1f%%", phenotype, cif_death * 100)]
cp[, y := 0.115 - 0.013 * as.integer(factor(phenotype))]
pc <- ggplot(cif, aes(time, cif_death * 100, colour = phenotype)) +
  geom_step(linewidth = 0.5) +
  geom_text(data = cp, aes(x = 60, y = y * 100, label = txt, colour = phenotype), hjust = 1, size = 1.65, family = FONT, show.legend = FALSE) +
  facet_wrap(~cohort, scales = "free_y") +
  scale_colour_manual(values = COL_PHENO, labels = PHENO_LAB(), name = NULL) +
  labs(x = L("Days from hospital admission", "入院后天数"), y = L("Cumulative incidence of in-hospital death, %\n(discharge alive as competing event)", "院内死亡累积发生, % (出院为竞争事件)")) +
  theme(legend.position = "none")

# d 死亡时相
dp <- fread(file.path(RES, "death_phase.csv"))[phenotype %in% c("P1", "P2", "P3", "P4")]
dp[, phase := factor(phase, c("<=7d", "8-14d", ">14d"), labels = c(L("≤ 7 d", "≤7天"), L("8-14 d", "8-14天"), L("> 14 d", ">14天")))]
dp[, pct := deaths / total_deaths * 100]
pd <- ggplot(dp, aes(phenotype, pct, fill = phase)) +
  geom_col(width = 0.7, colour = "grey25", linewidth = 0.2) +
  geom_text(aes(label = ifelse(deaths > 0, deaths, "")), position = position_stack(vjust = 0.5), size = 1.8, family = FONT) +
  geom_text(data = unique(dp[, .(phenotype, total_deaths)]), aes(phenotype, 102, label = paste0("n=", total_deaths)), inherit.aes = FALSE, size = 1.8, family = FONT, vjust = 0) +
  scale_fill_manual(values = c("#B2182B", "#F4A582", "#92C5DE"), name = L("Death timing", "死亡时相")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = L("Share of in-hospital deaths, % (MIMIC-IV)", "院内死亡构成, % (MIMIC-IV)"), title = L("Timing of in-hospital death", "院内死亡时相")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.5))
p <- (pa | pb) / (pc | pd) + plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "a")
save_pub(p, "FigureS5_infection_death", 190, 160)
