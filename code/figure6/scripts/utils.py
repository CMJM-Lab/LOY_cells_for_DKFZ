"""Shared helpers for LOY/ROY classification and statistics."""

from __future__ import annotations

from typing import Iterable, Optional, Sequence, Tuple

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
