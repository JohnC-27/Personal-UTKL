# Goals

## Mission

Support UTKL analysis by producing **reproducible, inspectable models and plots** of beam/profile distributions (especially MM xy-position histograms), and by developing **reliable SWIC scan readout** tooling.

## Near-term goals

### KDE / nominal profiles

1. **Faithful 2D KDE templates** of nominal (and study) xy-position TH2Ds using RooFit `RooNDKeysPdf`, including mirror/unmirror mixing and profiled amplitude where needed.
2. **Quantitative agreement** with data via overlays, 1D projections, ratios, and chi² (or chi²-contribution) diagnostics — not just visual similarity.
3. **Correct projection scaling** so KDE marginals match histogram projection normalization (bin width / grid factors must stay consistent).
4. **Comparisons to mz-nominal and other references** across runs and binnings (`75x75`, `100x90`, etc.) so fit quality is not tied to one histogram version.
5. **Keep scripts runnable and parameterizable** from top-of-file constants (or light CLI) so refits and plot regenerations are cheap.

### SWIC

1. **Stable CSV → channel-set pipeline** (96-wide rows → four × 24 channels).
2. **Feature extraction** (per-row extrema / peak indices) that supports scan characterization and figures.
3. **Date-organized raw data** remaining easy to find and reprocess.

## Longer-term / aspirational

- Clearer separation of “fit once → write ROOT” vs “plot from fit ROOT” pipelines.
- Reduce duplicated exploratory scripts; promote a small set of canonical fit + plot entry points.
- Document canonical input ROOT files and hist names per analysis campaign (e.g. Jan 2026).
- Optional: requirements / environment notes so ROOT + venv setup is reproducible across machines.

## Success criteria (practical)

A change is successful when:

- [ ] The intended ROOT output or plot regenerates without manual GUI steps (batch mode).
- [ ] Stats shown on plots (integrals, means, fit params, chi² where applicable) match what the code computes.
- [ ] Projection/ratio diagnostics do not show unexplained overall scale mismatches.
- [ ] `ai/work_log.md` and, if relevant, `ai/debugging_log.md` record what changed and why.

## Open questions

Record lasting unknowns here (not ephemeral bugs — those go in `ai/debugging_log.md`).

| ID | Question | Status |
|----|----------|--------|
| Q1 | Which ROOT file + hist name is canonical for the current Jan 2026 2D KDE campaign? | Open — scripts currently point at `jan2026studies_nominal` / `Jan2026_2d_kde` variants; case and naming differ across files |
| Q2 | Preferred mirror option string / mix strategy for production templates? | Open — see `NDKEYS_*` and mix params in `2d_kde.py` / `plot_2d_kde.py` |
| Q3 | SWIC analysis end product (plots only vs online/offline monitoring hook)? | Open |

## Non-goals (for now)

- Packaging as a public PyPI / shared lab framework
- Replacing official experiment software stacks
- Committing large regenerated plot/ROOT churn without an explicit request
