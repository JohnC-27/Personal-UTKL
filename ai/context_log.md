# Context log

Durable facts and decisions that should **survive chat resets**. Prefer this over re-deriving project knowledge from scratch.

## How to append

Newest entries on top. Tag the kind of fact.

```
### YYYY-MM-DD — <title>
- **Type:** fact | decision | convention | gotcha | pointer
- **Summary:** one or two sentences
- **Details:** optional bullets
- **Source:** script path, human statement, experiment note
```

---

### 2026-07-29 — Repo identity and layout

- **Type:** fact
- **Summary:** Personal-UTKL is a personal UTKL lab working repo: ROOT/PyROOT KDE and profile plotting under `scripts/`, artifacts in `root_files/` and `plots/`, SWIC CSV readout under `swic/`.
- **Details:**
  - Docs live in `docs/`; agent working memory in `ai/`
  - Local `.venv` has NumPy/SciPy/Matplotlib/Pandas; ROOT is external to that venv for most fit/plot scripts
- **Source:** repository inspection + root `README.md`

### 2026-07-29 — Canonical script families (KDE)

- **Type:** pointer
- **Summary:** 2D KDE fit and plot entry points center on `scripts/2d_kde.py` (build/write ROOT) and `scripts/plot_2d_kde.py` (overlays, projections, ratios).
- **Details:**
  - Fit inputs often reference `root_files/jan2026studies_nominal.root` and hist names like `nominal_xypos_1` / `nominalxyposMM1` (verify before use — naming varies)
  - Fit products include `root_files/jan2026_2d_kde.root` (case variants appear in code: `Jan2026_2d_kde.root` vs `jan2026_2d_kde.root` — confirm on disk)
  - RooNDKeysPdf options use strings such as `"a"` / `"am"` (not legacy Mirror enums)
- **Source:** script headers in `scripts/2d_kde.py`, `scripts/plot_2d_kde.py`

### 2026-07-29 — SWIC channel layout

- **Type:** fact
- **Summary:** SWIC CSVs are read as ~96 columns (one index skipped), split into four sets of 24 channels.
- **Details:**
  - Helper: `swic/2026-07-09/readout_functions.py` → `split_channels_numpy`, `get_mins`, `get_max`
  - Dated data: `swic/2026-07-09/*-scan_*.csv`
- **Source:** `readout_functions.py`

### 2026-07-29 — Projection scale is a sensitive area

- **Type:** gotcha
- **Summary:** 1D KDE projection vs TH1 projection comparisons have required explicit scale factors (`x_scale` / `y_scale`); mismatches show up as huge integral ratios even when shapes look plausible.
- **Details:**
  - Debug probes were written to `.cursor/debug-8940bd.log` from `plot_2d_kde.py:plot_projections`
  - Prefer summarizing confirmed root causes in `debugging_log.md` rather than re-reading full NDJSON
- **Source:** `.cursor/debug-8940bd.log`, `scripts/plot_2d_kde.py`
