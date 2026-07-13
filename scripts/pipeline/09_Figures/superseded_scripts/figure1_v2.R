suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggrepel)
  library(RColorBrewer)
  library(stringr)
  library(grid)
})

PDF_DEVICE <- cairo_pdf

PROJECT  <- "/data/habib/metagenome/biosurfactant"
BIGSCAPE <- file.path(PROJECT,
  "bigscape_output_all4/network_files/2026-03-21_20-43-03_hybrids_glocal")
DESEQ    <- file.path(PROJECT, "14_deseq2")
OUTDIR   <- file.path(PROJECT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

pub_theme <- theme_classic(base_size = 9, base_family = "sans") +
  theme(
    axis.text        = element_text(size = 7.5, colour = "black"),
    axis.title       = element_text(size = 8.5, colour = "black", face = "bold"),
    axis.line        = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.35, colour = "black"),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 7.5, face = "bold"),
    legend.margin    = margin(0, 2, 0, 2),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.background = element_blank(),
    legend.key       = element_blank(),
    plot.title       = element_text(size = 9, face = "bold", hjust = 0,
                                    margin = margin(b = 3)),
    plot.title.position = "plot",
    plot.subtitle    = element_text(size = 7, colour = "#555555",
                                    margin = margin(b = 4)),
    plot.caption.position = "plot",
    plot.margin      = margin(4, 6, 4, 6),
    panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.25),
    legend.key.size  = unit(0.35, "cm"),
    strip.background = element_rect(fill = "#2C3E50", colour = NA),
    strip.text       = element_text(colour = "white", face = "bold", size = 7.5,
                                    margin = margin(2, 2, 2, 2))
  )

bgc_cols <- c(
  "NRPS"            = "#C0392B",
  "PKS I"           = "#E67E22",
  "PKS other"       = "#F1C40F",
  "PKS-NRP Hybrids" = "#8E44AD",
  "RiPPs"           = "#27AE60",
  "Terpene"         = "#2980B9",
  "Others"          = "#7F8C8D"
)

wrap_class <- function(x) {
  ifelse(x == "PKS-NRP Hybrids", "PKS-NRP\nHybrids", as.character(x))
}

cd /data/habib/metagenome/biosurfactant/scripts/09_Figures

# Add path definitions at the top of the script
python3 << 'EOF'
with open('figure1_v2.R', 'r') as f:
    content = f.read()

header = '''suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggrepel)
  library(RColorBrewer)
  library(stringr)
  library(grid)
})

PDF_DEVICE <- cairo_pdf

PROJECT  <- "/data/habib/metagenome/biosurfactant"
BIGSCAPE <- file.path(PROJECT,
  "bigscape_output_all4/network_files/2026-03-21_20-43-03_hybrids_glocal")
DESEQ    <- file.path(PROJECT, "14_deseq2")
INTEG    <- file.path(PROJECT, "integration")
OUTDIR   <- file.path(PROJECT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

pub_theme <- theme_classic(base_size = 9, base_family = "sans") +
  theme(
    axis.text        = element_text(size = 7.5, colour = "black"),
    axis.title       = element_text(size = 8.5, colour = "black", face = "bold"),
    axis.line        = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.35, colour = "black"),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 7.5, face = "bold"),
    legend.margin    = margin(0, 2, 0, 2),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.background = element_blank(),
    legend.key       = element_blank(),
    plot.title       = element_text(size = 9, face = "bold", hjust = 0,
                                    margin = margin(b = 3)),
    plot.title.position = "plot",
    plot.subtitle    = element_text(size = 7, colour = "#555555",
                                    margin = margin(b = 4)),
    plot.caption.position = "plot",
    plot.margin      = margin(4, 6, 4, 6),
    panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.25),
    legend.key.size  = unit(0.35, "cm"),
    strip.background = element_rect(fill = "#2C3E50", colour = NA),
    strip.text       = element_text(colour = "white", face = "bold", size = 7.5,
                                    margin = margin(2, 2, 2, 2))
  )

bgc_cols <- c(
  "NRPS"            = "#C0392B",
  "PKS I"           = "#E67E22",
  "PKS other"       = "#F1C40F",
  "PKS-NRP Hybrids" = "#8E44AD",
  "RiPPs"           = "#27AE60",
  "Terpene"         = "#2980B9",
  "Others"          = "#7F8C8D"
)

wrap_class <- function(x) {
  ifelse(x == "PKS-NRP Hybrids", "PKS-NRP\\nHybrids", as.character(x))
}

'''

with open('figure1_v2.R', 'w') as f:
    f.write(header + content)
print("Done")
EOF

# Run
Rscript /data/habib/metagenome/biosurfactant/scripts/09_Figures/figure1_v2.R
