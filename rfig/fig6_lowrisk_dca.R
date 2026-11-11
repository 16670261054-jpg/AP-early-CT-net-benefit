source("00_config.R")
CUT_DISK <- fread(file.path(RES, "lowrisk_final_threshold.csv"))$threshold[1]
# Figure S5 De-implementation
#   a 入院可执行低危层 vs 高危层 vs 全队列：MCTSI>=6 才干预 的净获益
#   b 两层内事件被 CT 高/低分捕获情况
#   c 敏感性层（纯表型 P1/P2、原 P1/P2非SAP）的净获益，仅引流/清创
# 数据源：results/lowrisk_final_dca.csv, results/lowrisk_final_table.csv
dca <- fread(file.path(RES, "lowrisk_final_dca.csv"))
tab <- fread(file.path(RES, "lowrisk_final_table.csv"))
tab <- tab[scope == "Full ward cohort"]
N_LO <- tab[grepl("^Low risk", stratum), n]
N_HI <- tab[grepl("^Higher risk", stratum), n]
N_ALL <- N_LO + N_HI

olab <- c(`Drainage/necrosectomy` = L("Drainage or necrosectomy", "引流或清创"),
          `Persistent OF` = L("Persistent OF", "持续器官衰竭"))

XMAX <- 30  # 只画 2-30%，DCA 常规做法；>30% 的阈值在本临床问题里不现实
mk <- function(d) {
  dl <- melt(d[threshold <= XMAX / 100], id.vars = c("strat", "out", "threshold"),
             measure.vars = c("nb_ct_rule", "nb_treat_all"),
             variable.name = "rule", value.name = "nb")
  dl[, rule := factor(rule, c("nb_ct_rule", "nb_treat_all"),
                      labels = c(L("Intervene if MCTSI >= 6", "MCTSI>=6才干预"),
                                 L("Treat all", "全干预")))]
  dl
}
# 每个分面里 CT 规则同时优于 treat-none(0) 与 treat-all 的阈值上界
win_up <- function(d) d[nb_ct_rule > 0 & nb_ct_rule > nb_treat_all,
                        .(xw = max(threshold) * 100), by = .(strat, out)]

# ---- a 主分层 ----
sord <- c("Whole ward cohort", "Low risk (risk6 < cutoff)", "Higher risk (risk6 >= cutoff)")
slab <- c(L(sprintf("All admissions (n = %s)", format(N_ALL, big.mark = ",")), sprintf("全队列(n=%d)", N_ALL)),
          L(sprintf("Low risk, admission-defined (n = %d)", N_LO), sprintf("入院可判低危(n=%d)", N_LO)),
          L(sprintf("Higher risk (n = %d)", N_HI), sprintf("高危(n=%d)", N_HI)))
da <- dca[stratum %in% sord]
da[, strat := factor(stratum, sord, labels = slab)]
da[, out := olab[outcome]]
wa <- win_up(da)
pa <- ggplot(mk(da), aes(threshold * 100, nb, colour = rule)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_vline(data = wa, aes(xintercept = xw), linetype = "22", linewidth = 0.3, colour = "grey30", inherit.aes = FALSE) +
  geom_text(data = wa, aes(x = xw, y = 0.095, label = sprintf("%.0f%%", xw)),
            inherit.aes = FALSE, hjust = -0.15, size = 1.8, family = FONT, colour = "grey20") +
  geom_line(linewidth = 0.5) +
  facet_grid(out ~ strat) +
  coord_cartesian(ylim = c(-0.06, 0.105)) +
  scale_colour_manual(values = c("#B2182B", "grey55"), name = NULL) +
  labs(x = L("Threshold probability (%)", "阈值概率 (%)"), y = L("Net benefit", "净获益")) +
  theme(legend.position = "bottom", strip.text = element_text(size = 5.6))

# ---- b 两层内 CT 高/低分捕获的事件 ----
outs <- c("CT necrosis", "Drainage/necrosectomy", "Persistent OF", "In-hospital death")
sel <- rbindlist(lapply(outs, function(o)
  tab[, .(stratum, outcome = o,
          ctlow = get(paste0(o, "_ev_ct_lt6")), cthigh = get(paste0(o, "_ev_ct_ge6")))]))
sl <- melt(sel, id.vars = c("stratum", "outcome"), variable.name = "ct", value.name = "k")
sl[, ct := factor(ct, c("ctlow", "cthigh"), labels = c(L("MCTSI < 6", "MCTSI<6"), L("MCTSI >= 6", "MCTSI>=6")))]
sl[, strat := factor(ifelse(grepl("^Low risk", stratum), "lo", "hi"), c("lo", "hi"),
                     labels = c(L(sprintf("Low risk (n = %d)", N_LO), sprintf("低危(n=%d)", N_LO)),
                                L(sprintf("Higher risk (n = %d)", N_HI), sprintf("高危(n=%d)", N_HI))))]
olab2 <- c(`CT necrosis` = L("CT necrosis", "CT坏死"),
           `Drainage/necrosectomy` = L("Drainage/\nnecrosectomy", "引流/清创"),
           `Persistent OF` = L("Persistent\nOF", "持续OF"),
           `In-hospital death` = L("In-hospital\ndeath", "院内死亡"))
sl[, out := factor(olab2[outcome], unique(olab2))]
pb <- ggplot(sl, aes(out, k, fill = ct)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.2) +
  geom_text(data = sl[k > 0], aes(label = k), position = position_stack(vjust = 0.5), size = 1.7, family = FONT) +
  facet_wrap(~strat, scales = "free_y") +
  scale_fill_manual(values = c("#92C5DE", "#B2182B"), name = L("Admission CT", "入院CT")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(x = NULL, y = L("Events, n", "事件数")) +
  theme(legend.position = "bottom", strip.text = element_text(size = 5.8),
        axis.text.x = element_text(size = 5.3, lineheight = 0.9))

# ---- c 敏感性层（仅引流/清创）----
sens <- c("Low risk (risk6 < cutoff)", "Sensitivity A: P1/P2", "Sensitivity B: P1/P2 non-SAP")
dc <- dca[stratum %in% sens & outcome == "Drainage/necrosectomy"]
nlab <- dc[, .(n = n[1]), by = stratum]
clab <- setNames(c(L(sprintf("Primary: risk < %.3f (n = %d)", CUT_DISK, nlab[stratum == sens[1], n]), sprintf("主层: risk<%.3f (n=%d)", CUT_DISK, nlab[stratum == sens[1], n])),
                   L(sprintf("Phenotype P1/P2 (n = %d)", nlab[stratum == sens[2], n]), sprintf("表型P1/P2 (n=%d)", nlab[stratum == sens[2], n])),
                   L(sprintf("P1/P2 non-SAP (n = %d)", nlab[stratum == sens[3], n]), sprintf("P1/P2非SAP (n=%d)", nlab[stratum == sens[3], n]))), sens)
dc[, strat := factor(clab[stratum], clab)]
dc[, out := olab["Drainage/necrosectomy"]]
wc <- win_up(dc)
pc <- ggplot(mk(dc), aes(threshold * 100, nb, colour = rule)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_vline(data = wc, aes(xintercept = xw), linetype = "22", linewidth = 0.3, colour = "grey30", inherit.aes = FALSE) +
  geom_text(data = wc, aes(x = xw, y = 0.019, label = sprintf("%.0f%%", xw)),
            inherit.aes = FALSE, hjust = -0.2, size = 1.8, family = FONT, colour = "grey20") +
  geom_line(linewidth = 0.5) +
  facet_wrap(~strat, nrow = 1) +
  coord_cartesian(ylim = c(-0.03, 0.022)) +
  scale_colour_manual(values = c("#B2182B", "grey55"), name = NULL) +
  labs(x = L("Threshold probability (%)", "阈值概率 (%)"),
       y = L("Net benefit\n(drainage/necrosectomy)", "净获益\n(引流/清创)")) +
  theme(legend.position = "none", strip.text = element_text(size = 5.6))

p <- pa / pb / pc + plot_layout(heights = c(1, 0.75, 0.6)) + plot_annotation(tag_levels = "a")
save_pub(p, "Figure6_lowrisk_dca", 190, 205)
cat(sprintf("FigureS5 done | low=%d higher=%d all=%d\n", N_LO, N_HI, N_ALL))
