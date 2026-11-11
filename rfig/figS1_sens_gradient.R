source("00_config.R")
# Figure S1 敏感性分析: 各队列定义下的表型梯度 (剔除全插补/时间段/空白按无/脂肪酶子集/完整病例/非转院)
r <- fread(file.path(RES, "rates_sens.csv"))[!grepl("^Infected necrosis", outcome)]
r[, flab := str_wrap(paste0(cohort, ": ", outcome), 26)]
r[, flab := factor(flab, unique(flab))]
ann <- r[, .(ptxt = paste0(L("trend ", "趋势 "), p_star(trend_p[1]))), by = flab]
p <- ggplot(r, aes(phenotype, rate * 100, colour = phenotype)) +
  geom_line(aes(group = 1), colour = "grey55", linewidth = 0.35) +
  geom_errorbar(aes(ymin = lo * 100, ymax = hi * 100), width = 0.18, linewidth = 0.35) +
  geom_point(size = 1.5) +
  geom_text(data = ann, aes(x = 0.6, y = Inf, label = ptxt), inherit.aes = FALSE,
            hjust = 0, vjust = 1.4, size = 1.7, family = FONT, colour = "grey25") +
  facet_wrap(~flab, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = COL_PHENO, labels = PHENO_LAB(), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(x = NULL, y = L("Event rate, % (95% Wilson CI)", "事件率, % (95% Wilson CI)")) +
  theme(legend.position = "bottom", strip.text = element_text(size = 5.4, lineheight = 0.95),
        axis.text.x = element_text(size = 5.5))
save_pub(p, "FigureS1_sensitivity_gradients", 190, 245)
