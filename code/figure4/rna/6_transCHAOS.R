# ---- SETUP ----
library(dplyr)
library(ggplot2)
library(epiCHAOS)
library(Seurat)
library(qs2)

ANALYSIS_PATH = "analysis/rna"
PLOTS_PATH = "plots/rna"
PLOTS_PATH_PDF = "plots/rna/pdf"

CELLS = 1000
ITERATIONS = 100

source("code/figure4/helper_functions.R")


# ---- READ SEURAT RNA OBJECT ----
rna <- qs_read(file.path(ANALYSIS_PATH, "3_rna_preprocessed.qs2"))


# ---- TRANSCHAOS ----
transCHAOS <- transCHAOS(
  GetAssayData(rna)[VariableFeatures(rna), ],
  meta = rna@meta.data,
  colname = 'sample.name',
  n = CELLS,
  subsample = ITERATIONS
)


# ---- SAVE TRANSCHAOS SCORES ----
qs_save(
  transCHAOS,
  file.path(
    ANALYSIS_PATH,
    paste0("6_trans_chaos_", CELLS, "_", ITERATIONS, ".qs2")
  )
)


# ---- PLOT TRANSCHAOS ----
pl_trans_chaos_score <- transCHAOS %>%
  mutate(state = sub("group-", "", state)) %>%
  mutate(
    state = factor(
      state,
      levels = c("LOY1", "LOY2", "LOY3", "LOY4", "ROY1", "ROY2", "ROY3", "ROY4")
    )
  ) %>%
  ggplot(aes(state, het, fill = state)) +
  geom_boxplot(outliers = FALSE) +
  # geom_jitter(width = 0.1) +
  scale_y_continuous(lim = c(0, 1)) +
  scale_fill_manual(values = color.sample.name) +
  xlab("Sample") +
  ylab("Heterogeneity Score") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  paste0("6_pl_trans_chaos_score_", CELLS, "_", ITERATIONS, ".png"),
  pl_trans_chaos_score,
  path = PLOTS_PATH
)
ggsave(
  paste0("6_pl_trans_chaos_score_", CELLS, "_", ITERATIONS, ".pdf"),
  pl_trans_chaos_score,
  path = PLOTS_PATH_PDF
)
