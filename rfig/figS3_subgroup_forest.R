source("00_config.R")
# Figure S3 亚组一致性森林: 有序表型OR (每级P1->P4) 在预设亚组内 (a 本院; b MIMIC-IV)
sf <- fread(file.path(RES, "subgroup_forest.csv"))
sf[, txt := sprintf("%d/%d", events, n)]
mk <- function(d, xlab) {
  d[, out := factor(out, unique(out))]
  d[, sg := factor(paste0(subgroup_var, ": ", subgroup), rev(unique(paste0(subgroup_var, ": ", subgroup))))]
  ggplot(d, aes(OR_ordinal, sg)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey45") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4, colour = "#4393C3") +
    geom_point(size = 1.3, colour = "#B2182B") +
    geom_text(aes(x = 14, label = txt), hjust = 1, size = 1.5, family = FONT, colour = "grey30") +
    facet_wrap(~out, nrow = 1) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = c("0.25", "0.5", "1", "2", "4", "8")) +
    coord_cartesian(xlim = c(0.2, 16)) +
    labs(x = xlab, y = NULL) +
    theme(strip.text = element_text(size = 6), axis.text.y = element_text(size = 5.3))
}
loc <- sf[cohort == "Local"]
loc[, out := paste0(L("China ward: ", "中国病房: "), outcome)]
pa <- mk(loc, NULL)
mim <- sf[grepl("MIMIC", cohort)]
mim[, out := paste0("MIMIC-IV: ", fifelse(outcome == "Sepsis code", "Sepsis", outcome))]
pb <- mk(mim, L("Ordinal phenotype OR per step P1 -> P4 (adjusted for age and sex), log scale; labels = events/n",
                "有序表型OR (P1->P4每级, 调整年龄性别), 对数轴; 标注=事件/例数"))
p <- pa / pb + plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "a")
save_pub(p, "FigureS3_subgroup_forest", 190, 190)
