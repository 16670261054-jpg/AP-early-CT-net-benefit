source("00_config.R")
# Figure S6 a 各感染部位患者 vs 无该感染者的院内死亡率 (哑铃); b 长期死亡调整OR (MIMIC)
im <- fread(file.path(RES, "infection_site_mortality.csv"))
im[, site2 := fifelse(site == "Any extra-pancreatic infection", L("Any extra-pancreatic", "任何胰外感染"),
              fifelse(site == "Coded infected necrosis", L("Coded infected necrosis", "编码感染坏死"),
              fifelse(site == "UTI", L("Urinary tract", "泌尿道"), site)))]
ord <- im[cohort == "MIMIC-IV"][order(mort), site2]
im[, site2 := factor(site2, ord)]
iml <- melt(im, id.vars = c("cohort", "site2", "n"), measure.vars = c("mort", "mort_without"), variable.name = "grp", value.name = "m")
iml[, grp := factor(grp, c("mort_without", "mort"), labels = c(L("Without this infection", "无该感染"), L("With this infection", "有该感染")))]
pa <- ggplot(iml, aes(m * 100, site2)) +
  geom_line(aes(group = site2), colour = "grey60", linewidth = 0.35) +
  geom_point(aes(colour = grp), size = 1.7) +
  geom_text(data = im, aes(x = mort * 100, y = site2, label = paste0("n=", n)), vjust = -1, hjust = 0.3, size = 1.5, family = FONT, colour = "grey30") +
  facet_wrap(~cohort, scales = "free_x") +
  scale_colour_manual(values = c("grey45", "#B2182B"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.12))) +
  labs(x = L("In-hospital mortality, %", "院内死亡率, %"), y = NULL) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 5.8))

# b 长期死亡调整OR
dlo <- fread(file.path(RES, "death_longterm_or.csv"))
dlo[, outcome := factor(outcome, c("30-day death", "90-day death", "1-year death"),
                        labels = c(L("30-day death", "30天死亡"), L("90-day death", "90天死亡"), L("1-year death", "1年死亡")))]
dlo[, pheno := factor(pheno, c("P2", "P3", "P4"))]
pb <- ggplot(dlo, aes(OR, outcome, colour = pheno)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4, position = position_dodge(0.6)) +
  geom_point(size = 1.5, position = position_dodge(0.6)) +
  scale_colour_manual(values = COL_PHENO[2:4], labels = PHENO_LAB()[2:4], name = L("vs P1", "对照P1")) +
  guides(colour = guide_legend(nrow = 3)) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = c("0.25", "0.5", "1", "2", "4", "8")) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(xlim = c(0.15, 12)) +
  labs(x = L("Adjusted OR vs P1 (age, sex, SAP proxy; MIMIC-IV), log scale", "调整OR vs P1 (年龄/性别/SAP代理; MIMIC-IV), 对数轴"), y = NULL) +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.5))
p <- pa | pb
p <- p + plot_layout(widths = c(1.25, 1)) + plot_annotation(tag_levels = "a")
save_pub(p, "FigureS6_site_mortality_longterm", 190, 110)
