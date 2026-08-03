# Repository structure

```
Personal-UTKL/
├── README.md           # Short identity
├── AGENTS.md           # Top level agent description and instructions
├── docs/               # Stable documentation (this folder)
├── ai/                 # Agent ↔ human working logs and rules
├── scripts/            # Main analysis & plotting scripts (ROOT-heavy)
├── root_files/         # Input/output ROOT (and some CSV/TXT) artifacts
├── plots/              # Generated figures (PDF/PNG)
├── swic/               # SWIC scan data, date folders, readout helpers
├── .venv/              # Local Python virtualenv (not source of truth for physics)
├── .gitignore          # gitignore
└── .cursor/            # IDE debug dumps (ephemeral)
```

## `scripts/`

Runnable analysis tools. Prefer editing here over copying logic into notebooks.

| Pattern | Examples | Role |
|---------|----------|------|
| Fit / build KDE | `2d_kde.py`, `weighted_fixed_kde.py`, `unweighted_fixed_kde.py`, `optimize_kde.py` | Build templates; write ROOT outputs |
| Plot / diagnose | `plot_2d_kde.py`, `plot_kde_chi2_fit.py`, `plot_nominalxyposMM1.py`, `*_projection*` | Overlays, ratios, projections |
| Compare / test | `test_2d_kde_vs_mz_nominal.py`, `test_kde_vs_mz_nominal.py` | Model vs reference comparisons |
| Utilities | `file_inspection.py`, `rebin_TH2.py`, `txt_to_csv.py`, `remove_y_anomaly.py` | Inspect, rebin, convert, clean |
| Exploratory | `both.py`, `nominal_profile.py`, `nominal_profile_2d.py` | Interactive or early sketches; may hardcode absolute paths |

**Conventions agents should expect:**

- Paths often use `os.path.join(..., "..", "root_files"|"plots", ...)` relative to the script.
- Config constants live near the top of the file (`INPUT_ROOT_FILE`, `TARGET_HIST_NAME`, bin/rho settings).
- Indentation in this repo: **two spaces**.
- Batch mode: many scripts set `ROOT.gROOT.SetBatch(True)`.

## `root_files/`

Working store for histograms and fit products. Names encode intent (nominal, run, binning, correction, KDE variant).

Examples:

- Nominal / studies: `nominal.root`, `jan2026studies_nominal.root`, `nominal_corrected.root`
- Rebinned mz-nominal: `mz_nominal_*bin_run*.root`
- KDE products: `jan2026_2d_kde.root`, `weighted_fixed_kde.root`, `kde_chi2_fit.root`
- Systematics-style: `horn_current_p_1sigma.root`, `proton_beam_radius_p_1sigma.root`

Do not treat filenames as a strict schema; when in doubt, open with `file_inspection.py` or `TFile::ls()`.

## `plots/`

Generated figures only. Prefer regenerating via scripts rather than hand-editing plot files. Naming usually mirrors the script or comparison (e.g. `jan2026_2d_kde_*.pdf`, `2d_kde_vs_mz_nominal_*`).

## `swic/`

| Path | Role |
|------|------|
| `swic/2026-07-09/` | Dated scan CSVs (`*-scan_*.csv`) and readout scripts (`readout_functions.py`, testers) |
| `swic/scripts/` | Shared or copied readout testers |
| `swic/Figure_1.pdf` | Figure artifact |

Channel CSVs are wide (≈96 useful columns after skips); `readout_functions.split_channels_numpy` splits into four 24-channel sets.

## `docs/` vs `ai/`

| Folder | Stability | Content |
|--------|-----------|---------|
| `docs/` | Stable, curated | Goals, structure, overview — change when the project direction or layout changes |
| `ai/` | Living | Active task, work/debug/context logs — update every meaningful session |

## What not to treat as source of truth

- `.venv/` — environment only
- `.cursor/debug-*.log` — raw NDJSON debug probes; summarize durable findings into `ai/debugging_log.md`
- `__pycache__/` — ignore
