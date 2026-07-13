# fig_common.R  — shared palette, theme, helpers for all figures
suppressPackageStartupMessages({library(tidyverse); library(patchwork); library(ggrepel)})

PAL <- c(
  Halobacteriota="#eda100", Nanohaloarchaeota="#f0c869", Bacteroidota="#1baf7a",
  Pseudomonadota="#2a78d6", Chloroflexota="#008300", Actinomycetota="#eb6834",
  Myxococcota="#4a3aa7", Deinococcota="#e34948", Bacillota_A="#c98500",
  Gemmatimonadota="#e87ba4", unclassified="#9aa0a6")

SEASON <- c(summer="#eda100", winter="#2a78d6")

theme_pub <- function(base=12) {
  theme_classic(base_size=base) +
    theme(axis.text=element_text(color="black"),
          axis.title=element_text(color="black"),
          legend.key.size=unit(11,"pt"),
          legend.title=element_text(size=base-1),
          legend.text=element_text(size=base-2),
          plot.tag=element_text(face="bold", size=base+3),
          plot.subtitle=element_text(size=base-1))
}
# strip GTDB alphabetic suffixes for display (Halomonas_B -> Halomonas), keep in tables
clean_genus <- function(x) str_replace(x, "_[A-Z]$", "")

# ---- add this helper to fig_common.R ----
# flexible row picker: match a seasonal_function_summary row by keyword, case-insensitive
pick_fn <- function(df, pattern) df[grepl(pattern, df[["function"]], ignore.case=TRUE), ]
