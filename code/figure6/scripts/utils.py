"""Shared helpers for LOY/ROY classification and statistics."""

from __future__ import annotations

from typing import Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, pearsonr, spearmanr


def add_loy_roy_labels(
    adata,
    y_key: str = "Y_score",
    threshold: float = -0.2,
    label_col: str = "Y_group",
    low_label: str = "loy",
    high_label: str = "roy",
    inplace: bool = True,
) -> pd.Series:
    """Add LOY/ROY labels to adata.obs based on a Y-score threshold.

    Returns the created Series; by default also stores it in adata.obs.
    """
    if y_key not in adata.obs.columns:
        raise KeyError(f"{y_key} not found in adata.obs")

    labels = np.where(adata.obs[y_key] < threshold, low_label, high_label)
    labels = pd.Series(labels, index=adata.obs.index, name=label_col)

    if inplace:
        adata.obs[label_col] = labels.astype("category")

    return labels


def build_loy_roy_table(
    adata,
    group_key: str,
    label_col: str = "Y_group",
    low_label: str = "loy",
    high_label: str = "roy",
) -> pd.DataFrame:
    """Return a counts table of LOY/ROY by group.

    Output columns are guaranteed to include low_label and high_label.
    """
    missing = [c for c in (group_key, label_col) if c not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing columns in adata.obs: {missing}")

    counts = (
        adata.obs.groupby([group_key, label_col], observed=False)
        .size()
        .unstack(fill_value=0)
    )

    for col in (low_label, high_label):
        if col not in counts.columns:
            counts[col] = 0

    return counts[[low_label, high_label]]


def add_loy_roy_proportions(
    counts: pd.DataFrame,
    low_label: str = "loy",
    high_label: str = "roy",
) -> pd.DataFrame:
    """Add total counts and LOY/ROY proportions to a counts table."""
    out = counts.copy()
    out["n_total"] = out[low_label] + out[high_label]
    out["loy_pct"] = out[low_label] / out["n_total"].replace(0, np.nan)
    out["roy_pct"] = out[high_label] / out["n_total"].replace(0, np.nan)
    return out


def chi_square_loy_roy(counts: pd.DataFrame) -> Tuple[float, float, int, np.ndarray]:
    """Run a chi-square test on a LOY/ROY contingency table."""
    if counts.shape[1] < 2:
        raise ValueError("counts must have at least two columns (LOY and ROY)")
    chi2, p_value, dof, expected = chi2_contingency(counts)
    return chi2, p_value, dof, expected


def fisher_pairwise_vs_reference(
    counts: pd.DataFrame,
    reference: str,
    low_label: str = "loy",
    high_label: str = "roy",
) -> pd.DataFrame:
    """Run pairwise Fisher's exact tests vs a reference group.

    Returns a DataFrame with odds_ratio, p_value, and table for each group.
    """
    if reference not in counts.index:
        raise KeyError(f"Reference group '{reference}' not found in counts.index")

    ref = counts.loc[reference]
    results = {}
    for group in counts.index:
        if group == reference:
            continue
        test_table = [
            [ref[low_label], ref[high_label]],
            [counts.loc[group, low_label], counts.loc[group, high_label]],
        ]
        odds_ratio, p_value = fisher_exact(test_table)
        results[group] = {
            "odds_ratio": odds_ratio,
            "p_value": p_value,
            "table": test_table,
        }

    return pd.DataFrame(results).T


def build_s6h_clone_enrichment_table(
    obs: pd.DataFrame,
    *,
    lineage_key: str = "LineageGroup",
    site_key: str = "Sample2",
    label_key: str = "Y_group",
    tree_key: str = "TreeMetRate",
    dataset_key: str = "dataset",
    dataset_value: str = "M5k",
    exclude_lineage_groups: Optional[Sequence[object]] = None,
    exclude_site_prefix: str = "LM0",
    min_cells_per_site: int = 20,
    enrich_pct: float = 0.75,
    weak_cutoff: float = 0.001,
    high_cutoff: float = 0.008,
    tree_summary: str = "median",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Build the canonical Figure S6H clone-level enrichment table.

    Returns `(site_counts_eligible, clone_enriched)` where `clone_enriched`
    contains one row per clone with LOY/ROY enrichment, clone-level TreeMetRate,
    and metastatic phenotype assignments.
    """
    required = [lineage_key, site_key, label_key, tree_key, dataset_key]
    missing = [c for c in required if c not in obs.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    df = obs[required].copy()
    df = df[df[dataset_key].astype(str) == str(dataset_value)]
    df[lineage_key] = df[lineage_key].astype(str)
    df[site_key] = df[site_key].astype(str)

    if exclude_lineage_groups:
        excluded = {str(x) for x in exclude_lineage_groups}
        df = df[~df[lineage_key].isin(excluded)]
    if exclude_site_prefix:
        df = df[~df[site_key].str.startswith(str(exclude_site_prefix))]

    site_counts = (
        df.groupby([lineage_key, site_key, label_key], observed=False)
        .size()
        .unstack(fill_value=0)
    )
    for col in ("LOY", "ROY"):
        if col not in site_counts.columns:
            site_counts[col] = 0
    site_counts = site_counts[["LOY", "ROY"]]
    site_counts["n_total"] = site_counts["LOY"] + site_counts["ROY"]
    site_counts["loy_frac"] = site_counts["LOY"] / site_counts["n_total"].replace(0, np.nan)
    site_counts["roy_frac"] = site_counts["ROY"] / site_counts["n_total"].replace(0, np.nan)
    site_counts_eligible = site_counts[site_counts["n_total"] >= min_cells_per_site].copy()

    def _classify_clone_enrichment(sub: pd.DataFrame) -> Union[float, str]:
        if sub.empty:
            return np.nan
        lf = sub["loy_frac"].to_numpy(dtype=float)
        rf = sub["roy_frac"].to_numpy(dtype=float)
        if np.all(lf >= enrich_pct):
            return "LOY"
        if np.all(rf >= enrich_pct):
            return "ROY"
        return "Mixed"

    clone_enriched = (
        site_counts_eligible.groupby(level=0, observed=False)
        .apply(_classify_clone_enrichment)
        .rename("enriched")
        .to_frame()
    )
    clone_enriched["n_sites_eligible"] = site_counts_eligible.groupby(level=0, observed=False).size()
    clone_enriched = clone_enriched.dropna(subset=["enriched"])

    if tree_summary == "median":
        clone_tree = df.groupby(lineage_key, observed=False)[tree_key].median()
    elif tree_summary == "mean":
        clone_tree = df.groupby(lineage_key, observed=False)[tree_key].mean()
    else:
        raise ValueError("tree_summary must be 'median' or 'mean'")

    clone_enriched = clone_enriched.join(clone_tree.rename(tree_key), how="left")

    def _metapheno(value: float) -> str:
        if value >= high_cutoff:
            return "Highly metastatic"
        if value >= weak_cutoff:
            return "Weakly metastatic"
        return "Non metastatic"

    clone_enriched["metastatic_phenotype"] = (
        clone_enriched[tree_key].fillna(0).map(_metapheno)
    )
    clone_enriched = clone_enriched[clone_enriched["enriched"].isin(["LOY", "ROY"])].copy()

    return site_counts_eligible, clone_enriched


def summarize_s6h_clone_phenotypes(
    clone_enriched: pd.DataFrame,
    *,
    enriched_key: str = "enriched",
    phenotype_key: str = "metastatic_phenotype",
    enriched_order: Sequence[str] = ("LOY", "ROY"),
    phenotype_order: Sequence[str] = (
        "Non metastatic",
        "Weakly metastatic",
        "Highly metastatic",
    ),
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Summarize Figure S6H clone phenotypes and run convenience Fisher tests."""
    summary = (
        clone_enriched.groupby([enriched_key, phenotype_key], observed=False)
        .size()
        .unstack(fill_value=0)
    )
    for phenotype in phenotype_order:
        if phenotype not in summary.columns:
            summary[phenotype] = 0
    summary = summary.reindex(index=list(enriched_order), fill_value=0)
    summary = summary[list(phenotype_order)]

    pct = summary.div(summary.sum(axis=1).replace(0, np.nan), axis=0) * 100

    hi_loy = int(summary.loc["LOY", "Highly metastatic"]) if "LOY" in summary.index else 0
    hi_roy = int(summary.loc["ROY", "Highly metastatic"]) if "ROY" in summary.index else 0
    tot_loy = int(summary.loc["LOY"].sum()) if "LOY" in summary.index else 0
    tot_roy = int(summary.loc["ROY"].sum()) if "ROY" in summary.index else 0
    table_hi = [[hi_loy, tot_loy - hi_loy], [hi_roy, tot_roy - hi_roy]]
    odds_hi, p_hi = fisher_exact(table_hi)

    met_loy = (
        int(summary.loc["LOY", "Weakly metastatic"] + summary.loc["LOY", "Highly metastatic"])
        if "LOY" in summary.index
        else 0
    )
    met_roy = (
        int(summary.loc["ROY", "Weakly metastatic"] + summary.loc["ROY", "Highly metastatic"])
        if "ROY" in summary.index
        else 0
    )
    non_loy = int(summary.loc["LOY", "Non metastatic"]) if "LOY" in summary.index else 0
    non_roy = int(summary.loc["ROY", "Non metastatic"]) if "ROY" in summary.index else 0
    table_any = [[met_loy, non_loy], [met_roy, non_roy]]
    odds_any, p_any = fisher_exact(table_any)

    fisher_rows = pd.DataFrame(
        [
            {
                "test": "Highly_vs_Others",
                "odds_ratio": float(odds_hi),
                "p_value": float(p_hi),
                "table": table_hi,
            },
            {
                "test": "AnyMet_vs_Non",
                "odds_ratio": float(odds_any),
                "p_value": float(p_any),
                "table": table_any,
            },
        ]
    )

    return summary, pct, fisher_rows


def plot_score_correlation(
    adata,
    x_key: str = "Y_score",
    y_key: str = "scTreeMetRate",
    group_key: str = "LineageGroup",
    title: str = "Y_score vs scTreeMetRate",
    save_path: Optional[str] = None,
    dpi: int = 300,
    corr_method: str = "spearman",
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    group_is_categorical: bool = True,
    min_group_size: int = 0,
    groups_to_exclude: Optional[Sequence[str]] = None,
) -> Tuple[float, float]:
    """Scatter plot with correlation (Spearman or Pearson).

    Returns (r, p_value) on the filtered data.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    required = [x_key, y_key, group_key]
    missing = [c for c in required if c not in adata.obs.columns]
    if missing:
        raise KeyError(f"Required columns not found in adata.obs: {missing}")

    plot_df = adata.obs[[x_key, y_key, group_key]].copy()
    plot_df = plot_df.dropna(subset=[x_key, y_key])
    if plot_df.empty:
        raise ValueError("No data available after dropping NaN values")

    if groups_to_exclude:
        plot_df = plot_df[~plot_df[group_key].astype(str).isin([str(g) for g in groups_to_exclude])]

    if min_group_size > 0:
        group_counts = plot_df[group_key].value_counts()
        keep_groups = group_counts[group_counts >= min_group_size].index
        plot_df = plot_df[plot_df[group_key].isin(keep_groups)]

    if plot_df.empty:
        raise ValueError("No data available after filtering by group size/exclusions")

    if group_is_categorical:
        plot_df[group_key] = plot_df[group_key].astype(str)
        legend_title = group_key
    else:
        legend_title = group_key

    if corr_method.lower() == "spearman":
        r, p_value = spearmanr(plot_df[x_key], plot_df[y_key])
        corr_label = "Spearman r"
    elif corr_method.lower() == "pearson":
        r, p_value = pearsonr(plot_df[x_key], plot_df[y_key])
        corr_label = "Pearson r"
    else:
        raise ValueError("corr_method must be 'pearson' or 'spearman'")

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(8, 6), dpi=dpi)

    sns.scatterplot(
        data=plot_df,
        x=x_key,
        y=y_key,
        hue=group_key,
        ax=ax,
        s=10,
        alpha=0.6,
        edgecolor="w",
        linewidth=0.1,
    )

    try:
        sns.regplot(
            data=plot_df,
            x=x_key,
            y=y_key,
            scatter=False,
            color="black",
            line_kws={"linestyle": "--", "alpha": 0.8, "lw": 1.5},
            ax=ax,
        )
    except Exception:
        pass

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    ax.set_title(title)
    ax.set_xlabel(x_key.replace("_", " ").title())
    ax.set_ylabel(y_key.replace("_", " ").title())

    corr_text = f"{corr_label}: {r:.2f}\nP: {p_value:.2e}"
    ax.text(
        0.05,
        0.95,
        corr_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox={"boxstyle": "round,pad=0.5", "fc": "white", "alpha": 0.8},
    )

    ax.legend(title=legend_title, bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")

    plt.show()
    plt.close(fig)

    return r, p_value
