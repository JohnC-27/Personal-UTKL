# Overview

Human- and agent-facing documentation for **Personal-UTKL**.

## Docs map

| File | Purpose |
|------|---------|
| [overview.md](overview.md) | What this repo is, who it is for, stack, and this map |
| [structure.md](structure.md) | Directory layout, data/script conventions, key artifacts |
| [goals.md](goals.md) | Research and engineering goals, success criteria, open questions |
| [root_file_catalog.md](root_file_catalog.md) | Catalog of `.root` files and TH2 names/sizes |

For session state, debugging history, and agent rules, see [`../ai/`](../ai/).

## What this repo is

**Personal-UTKL** is a personal working repository for UTKL lab analysis: Python scripts that inspect ROOT histograms, build and tune KDE models of beam/profile distributions, compare model vs data, and produce plots. A second track (`swic/`) covers SWIC scan readout and channel analysis from CSV.

This is **not** a packaged library. Scripts are runnable analysis tools with paths and parameters set near the top of each file.

## Primary workstreams

1. **Nominal / MM xy-position KDE analysis** (`scripts/`, `root_files/`, `plots/`)
   - Fit 1D/2D KDEs (ROOT `TKDE`, RooFit `RooNDKeysPdf`) to nominal and related TH2Ds (e.g. `nominalxyposMM1`–`MM3`, Jan 2026 studies).
   - Optimize bandwidth, mirror mix, amplitude; compare to mz-nominal and other variants.
   - Emit overlays, projections, ratios, and chi² diagnostics as PDF/PNG under `plots/`.

2. **SWIC data readout** (`swic/`)
   - Parse multi-channel scan CSVs, split channel sets, extract peak/minima features, and support figure generation.

## Stack

- **Language:** Python 3 (local `.venv` present)
- **HEP / histograms:** ROOT (PyROOT), RooFit (`RooNDKeysPdf`, `RooDataSet`)
- **Scientific:** NumPy, SciPy, Matplotlib, Pandas (especially SWIC path)
- **Environment:** Activate `.venv` when running non-ROOT-only scripts; ROOT must be available for most `scripts/` tools

## Audience

| Reader | Use this for |
|--------|----------------|
| Human (lab member / self) | Find which script owns a plot, where ROOT inputs live, what goal a change serves |
| Agent | Orient before editing: read `docs/` + `ai/active_task.md` + `ai/instructions.md` first |

## Related paths outside docs

| Path | Role |
|------|------|
| `README.md` | One-line repo identity |
| `ai/` | Agent/human working memory (tasks, context, debug, rules) |
| `.cursor/` | Cursor debug session logs (transient; not the durable `ai/` log) |
