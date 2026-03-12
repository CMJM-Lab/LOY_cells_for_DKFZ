from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Figure6Paths:
    repo_root: Path
    code_dir: Path
    data_dir: Path
    docs_dir: Path
    figures_root: Path
    runtime_dir: Path
    downloads_dir: Path
    raw_dir: Path
    processed_dir: Path


@dataclass(frozen=True)
class NotebookPaths:
    common: Figure6Paths
    save_dir: Path
    figure_dir: Path


def find_repo_root(start: Path | None = None) -> Path:
    cwd = (start or Path.cwd()).resolve()
    for p in [cwd, *cwd.parents]:
        if (p / 'code' / 'figure6' / 'notebooks').exists() and (p / 'data' / 'figure6').exists():
            return p
    raise RuntimeError(
        'Could not locate the repository root. Launch Jupyter from within the repo (or a subdirectory) so code/figure6/notebooks is visible.'
    )


def get_figure6_paths(repo_root: Path | None = None) -> Figure6Paths:
    root = (repo_root or find_repo_root()).resolve()
    data_dir = root / 'data' / 'figure6'
    runtime_dir = data_dir / 'runtime'
    return Figure6Paths(
        repo_root=root,
        code_dir=root / 'code' / 'figure6',
        data_dir=data_dir,
        docs_dir=root / 'docs' / 'figure6',
        figures_root=root / 'figures' / 'figure6',
        runtime_dir=runtime_dir,
        downloads_dir=runtime_dir / 'downloads',
        raw_dir=runtime_dir / 'raw',
        processed_dir=runtime_dir / 'processed',
    )


def get_notebook_paths(name: str, repo_root: Path | None = None) -> NotebookPaths:
    common = get_figure6_paths(repo_root)
    save_dir = common.data_dir / 'derived' / name
    figure_dir = common.figures_root / name
    for d in [common.data_dir, common.runtime_dir, common.downloads_dir, common.raw_dir, common.processed_dir, save_dir, figure_dir]:
        d.mkdir(parents=True, exist_ok=True)
    return NotebookPaths(common=common, save_dir=save_dir, figure_dir=figure_dir)
