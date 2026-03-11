# ---- SETUP ----
library(dplyr)
library(ggplot2)


# ---- GGPLOT2 SETUP ----
theme_set(
  theme_classic() +
    theme(
      title = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.title.position = "plot",
      legend.position = "none"
    )
)

color.condition <- c("LOY" = "#0000C0", "ROY" = "#ED7D31")

color.orig.ident <- c(
  'LOY_D11' = '#0000C0',
  'LOY_E5b' = '#0000C0',
  'LOY_E6b' = '#0000C0',
  'LOY_F10' = '#0000C0',
  'WT_A3' = '#ED7D31',
  'WT_C1' = '#ED7D31',
  'WT_C8' = '#ED7D31',
  'WT_D5' = '#ED7D31'
)

color.sample.name <- c(
  'LOY1' = '#0000C0',
  'LOY2' = '#0000C0',
  'LOY3' = '#0000C0',
  'LOY4' = '#0000C0',
  'ROY1' = '#ED7D31',
  'ROY2' = '#ED7D31',
  'ROY3' = '#ED7D31',
  'ROY4' = '#ED7D31'
)

# ---- HELPER FUNCTIONS ----
med_mad <- function(x, mad.cutoff, lower.limit = NULL, upper.limit = NULL) {
  ymin = median(x) - mad.cutoff * mad(x)
  ymax = median(x) + mad.cutoff * mad(x)
  if (!is.null(lower.limit)) {
    ymin = lower.limit
  }
  if (!is.null(upper.limit)) {
    ymax = upper.limit
  }
  return(c(ymin = ymin, ymax = ymax))
}

get_low_qc_cells <- function(meta.data, qc.metric, mad.cutoff = 5, lower.limit = NULL, upper.limit = NULL) {
  meta.data %>%
    mutate(cell_id = rownames(meta.data)) %>%
    group_by(orig.ident) %>%
    mutate(
      ymin = med_mad(!!sym(qc.metric), mad.cutoff, lower.limit, upper.limit)[1],
      ymax = med_mad(!!sym(qc.metric), mad.cutoff, lower.limit, upper.limit)[2]
    ) %>%
    filter(!!sym(qc.metric) < ymin | !!sym(qc.metric) > ymax) %>%
    pull(cell_id)
}

plot_qc_metric <- function(meta.data, qc.metric, mad.cutoff = 5, lower.limit = NULL, upper.limit = NULL) {
  low.qc.cells <- get_low_qc_cells(meta.data, qc.metric, mad.cutoff, lower.limit, upper.limit)
  meta.data %>%
    mutate(cell_id = rownames(meta.data)) %>%
    mutate(low_qc = cell_id %in% low.qc.cells) %>%
    ggplot(aes(sample.name, !!sym(qc.metric), fill = sample.name)) +
    geom_violin() +
    geom_jitter(aes(color = low_qc), width = 0.33, size = 0.1) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "lightgrey")) +
    scale_fill_manual(values = color.sample.name) +
    stat_summary(
      fun = med_mad,
      fun.args = list(mad.cutoff = mad.cutoff, lower.limit = lower.limit, upper.limit = upper.limit),
      geom = "crossbar",
      color = "red"
    ) +
    ggtitle(
      paste0(
        qc.metric,
        "\n",
        if (is.null(lower.limit) | is.null(upper.limit)) paste0("n_mads=", mad.cutoff, " "),
        if (!is.null(lower.limit)) paste0("lower_limit=", round(lower.limit, 1), " ") else "",
        if (!is.null(upper.limit)) paste0("upper_limit=", round(upper.limit, 1)) else ""
      )
    )
}

plot_qc_scatter <- function(meta.data, qc.metric.1, qc.metric.2, mad.cutoff.1 = 5, mad.cutoff.2 = 5) {
  low.qc.cells.1 <- get_low_qc_cells(meta.data, qc.metric.2, mad.cutoff.1)
  low.qc.cells.2 <- get_low_qc_cells(meta.data, qc.metric.1, mad.cutoff.2)
  low.qc.cells <- union(low.qc.cells.1, low.qc.cells.2)
  meta.data %>%
    mutate(cell_id = rownames(meta.data)) %>%
    mutate(low_qc = cell_id %in% low.qc.cells) %>%
    ggplot(aes(!!sym(qc.metric.1), !!sym(qc.metric.2), color = sample.name)) +
    geom_point(aes(alpha = low_qc)) +
    scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.1)) +
    guides(alpha = "none") +
    theme(legend.position = "right")
}
