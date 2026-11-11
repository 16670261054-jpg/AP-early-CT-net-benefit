# 00_config.R —— AP_NECRO 出版图公共配置（pub-plotting skill 规范: theme_pub / save_pub 四格式 / 全文锁色）
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(patchwork)
  library(scales); library(stringr); library(ggdist)
})

BASE <- dirname(getwd())
RES  <- Sys.getenv("RES_OVERRIDE", file.path(BASE, "results"))
FIGLANG <- Sys.getenv("FIGLANG", "en")
PUBDIR  <- Sys.getenv("PUBDIR_OVERRIDE", file.path(BASE, "figures"))
dir.create(PUBDIR, recursive = TRUE, showWarnings = FALSE)

FONT <- if (FIGLANG == "zh") "Microsoft YaHei" else "Arial"
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"),
               `Microsoft YaHei` = windowsFont("Microsoft YaHei"))
}
L <- function(en, zh) if (FIGLANG == "zh") zh else en
SEED <- 123

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(axis.line = element_line(linewidth = 0.35),
          axis.ticks = element_line(linewidth = 0.35),
          plot.title = element_text(size = base_size + 1, face = "bold"),
          legend.key.size = unit(3.2, "mm"),
          legend.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = base_size, face = "bold"),
          plot.tag = element_text(size = base_size + 2, face = "bold"))
}
theme_set(theme_pub())

# 全文锁色：表型 P1-P4 = 严重度梯度（蓝->红），与表型文一致
COL_PHENO <- c(P1 = "#4393C3", P2 = "#92C5DE", P3 = "#F4A582", P4 = "#B2182B")
COL_TWO   <- c("#4393C3", "#D6604D")
COL_MODEL <- c(`Phenotype` = "#4393C3", `Phenotype + CT` = "#B2182B",
               `Risk score` = "#92C5DE", `Risk score + CT` = "#D6604D", `CT (MCTSI) alone` = "grey40")
COL_COHORT <- c(Local = "#B2182B", `MIMIC-IV` = "#4393C3", eICU = "#F4A582")
PHENO_LAB <- function() c(
  P1 = L("P1 Lymphocyte-preserved", "P1 淋巴保留型"),
  P2 = L("P2 Young hemoconcentrated", "P2 年轻血浓缩型"),
  P3 = L("P3 Hepatobiliary", "P3 肝胆型"),
  P4 = L("P4 Hypocalcemic-lymphopenic", "P4 低钙淋巴减少型"))

# Ann Med 版（全局规则10）：图上只打星号；阈值表写进图注
p_star <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-4, "****", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))))
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))

save_pub <- function(p, name, w, h) {
  ggsave(file.path(PUBDIR, paste0(name, ".png")), p, width = w, height = h,
         units = "mm", dpi = 300, bg = "white")
  ggsave(file.path(PUBDIR, paste0(name, ".pdf")), p, width = w / 25.4, height = h / 25.4,
         device = cairo_pdf, bg = "white")
  ggsave(file.path(PUBDIR, paste0(name, ".svg")), p, width = w / 25.4, height = h / 25.4,
         bg = "white")
  ggsave(file.path(PUBDIR, paste0(name, ".tiff")), p, width = w, height = h,
         units = "mm", dpi = 600, bg = "white", compression = "lzw")
  cat("saved:", name, "->", PUBDIR, "\n")
}
