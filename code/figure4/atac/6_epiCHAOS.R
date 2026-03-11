# ---- SETUP ----
library(qs2)
library(Signac)
library(Seurat)
library(epiCHAOS)
library(dplyr)
library(ggplot2)

ANALYSIS_PATH = "analysis/atac"
PLOTS_PATH = "plots/atac"
PLOTS_PATH_PDF = "plots/atac/pdf"

CELLS = 300
ITERATIONS = 20

source("code/figure4/helper_functions.R")


# ---- READ SEURAT OBJECT ----
atac <- qs_read(file.path(ANALYSIS_PATH, "3_atac_preprocessed.qs2"))


# ---- EPICHAOS ----
epiCHAOS <- epiCHAOS(
  GetAssayData(atac[["ATAC"]])[VariableFeatures(atac), ],
  meta = atac[[]],
  colname = "sample.name",
  n = CELLS,
  index = NULL,
  plot = F,
  cancer = F,
  subsample = ITERATIONS
)


# ---- SAVE EPICHAOS SCORES ----
qs_save(
  epiCHAOS,
  file.path(
    ANALYSIS_PATH,
    paste0("6_epiCHAOS_", CELLS, "_", ITERATIONS, ".qs2")
  )
)


# ---- PLOT EPICHAOS ----
pl_epi_chaos_score <- epiCHAOS %>%
  mutate(state = sub("group-", "", state)) %>%
  mutate(
    state = factor(
      state,
      levels = c("LOY1", "LOY2", "LOY3", "LOY4", "ROY1", "ROY2", "ROY3", "ROY4")
    )
  ) %>%
  ggplot(aes(state, het.adj, fill = state)) +
  geom_boxplot(outliers = FALSE) +
  # geom_jitter(width = 0.1) +
  scale_fill_manual(values = color.sample.name) +
  scale_y_continuous(lim = c(0, 1)) +
  xlab("Sample") +
  ylab("Heterogeneity Score") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  paste0("6_pl_epi_chaos_score_", CELLS, "_", ITERATIONS, ".png"),
  pl_epi_chaos_score,
  path = PLOTS_PATH
)
ggsave(
  paste0("6_pl_epi_chaos_score_", CELLS, "_", ITERATIONS, ".pdf"),
  pl_epi_chaos_score,
  path = PLOTS_PATH_PDF
)
