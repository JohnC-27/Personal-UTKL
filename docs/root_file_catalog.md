# ROOT file catalog (TH2s)

Inventory of `.root` files under `beam_profile_modeling/root_files/`, with TH2 names and sizes from `beam_profile_modeling/scripts/file_inspection.py` logic (refreshed 2026-07-31).

**Naming conventions**

| Suffix / pattern | Meaning |
|------------------|---------|
| `_corrected` | Y-axis trimmed (top anomaly cut; typically y ≲ 80–84 cm) via `remove_y_anomaly.py` — TH2s are not square |
| `_75x75` | Axis range restricted to ≈±75 cm (bin count scales with original density) |
| `_2d_kde` | Fit product from `2d_kde.py` (`target_hist`, `kde_shape`, `kde_template` + meta/graphs) |
| `mz_nominal_*bin_run{1,2}` | Mo’s high-res `NominalxyposMM*` alcove hists; bin count in name; run1 vs run2 |
| `*_p_1sigma` | Systematics-style +1σ variation of the named parameter |

Axis ranges are in the same units as stored in the file (typically cm for alcove xy positions). Where a TH2 has multiple write cycles (`;1`, `;2`), the latest cycle is what `TFile.Get(name)` returns; sizes below match that object.

---

## Index

| File | # TH2s | Description |
|------|--------|-------------|
| [2d_kde.root](#2d_kderoot) | 3 | Early 2D KDE fit on `nominal_corrected` MM1 (25×23) |
| [horn_current_p_1sigma.root](#horn_current_p_1sigmaroot) | 3 | Horn-current +1σ alcove xy (25×25) |
| [jan2026_2d_kde.root](#jan2026_2d_kderoot) | 3 | 2D KDE fit on Jan 2026 nominal MM1 (200×200, uncorrected) |
| [jan2026studies_nominal.root](#jan2026studies_nominalroot) | 3 | Jan 2026 study nominal alcove xy (200×200, ±125 cm) |
| [jan2026studies_nominal_corrected.root](#jan2026studies_nominal_correctedroot) | 3 | Jan 2026 nominal with y anomaly trimmed (200×188) |
| [jan2026studies_nominal_corrected_2d_kde.root](#jan2026studies_nominal_corrected_2d_kderoot) | 3 | Current 2D KDE fit on Jan 2026 corrected MM1 |
| [ml_tracking_hornb_tilt_down_5sigma_beamshift_nX_100um_trackingon.root](#ml_tracking_hornb_tilt_down_5sigma_beamshift_nx_100um_trackingonroot) | 8 | ML tracking: Horn B tilt −5σ + beam −X 100 µm |
| [ml_tracking_hornb_tilt_down_5sigma_beamshift_nY_250um_trackingon.root](#ml_tracking_hornb_tilt_down_5sigma_beamshift_ny_250um_trackingonroot) | 8 | ML tracking: Horn B tilt −5σ + beam −Y 250 µm |
| [ml_tracking_hornb_tilt_down_5sigma_beamshift_pX_100um_trackingon.root](#ml_tracking_hornb_tilt_down_5sigma_beamshift_px_100um_trackingonroot) | 8 | ML tracking: Horn B tilt −5σ + beam +X 100 µm |
| [ml_tracking_hornb_tilt_down_5sigma_beamshift_pY_250um_trackingon.root](#ml_tracking_hornb_tilt_down_5sigma_beamshift_py_250um_trackingonroot) | 8 | ML tracking: Horn B tilt −5σ + beam +Y 250 µm |
| [ml_tracking_hornb_tilt_down_5sigma_nominal_trackingon.root](#ml_tracking_hornb_tilt_down_5sigma_nominal_trackingonroot) | 8 | ML tracking: Horn B tilt −5σ, nominal beam (no shift) |
| [mz_nominal_100bin_run1.root](#mz_nominal_100bin_run1root) | 3 | Mo nominal run1 rebinned to 100×100 |
| [mz_nominal_100bin_run2.root](#mz_nominal_100bin_run2root) | 3 | Mo nominal run2 rebinned to 100×100 |
| [mz_nominal_2000bin_run1.root](#mz_nominal_2000bin_run1root) | 3 | Mo nominal run1 full-res 2000×2000 (±100 cm) |
| [mz_nominal_2000bin_run1_75x75.root](#mz_nominal_2000bin_run1_75x75root) | 3 | Mo run1 restricted to ±75 cm (1500×1500) |
| [mz_nominal_2000bin_run1_corrected.root](#mz_nominal_2000bin_run1_correctedroot) | 3 | Mo run1 y-trimmed (±75 cm, 1500×1500) |
| [mz_nominal_2000bin_run2.root](#mz_nominal_2000bin_run2root) | 3 | Mo nominal run2 full-res 2000×2000 (±100 cm) |
| [mz_nominal_2000bin_run2_75x75.root](#mz_nominal_2000bin_run2_75x75root) | 3 | Mo run2 restricted to ±75 cm (1500×1500) |
| [mz_nominal_2000bin_run2_corrected.root](#mz_nominal_2000bin_run2_correctedroot) | 3 | Mo run2 y-trimmed (±75 cm, 1500×1500) |
| [nominal.root](#nominalroot) | 3 | Coarse lab nominal alcove xy (25×25, ±100 cm) |
| [nominal_75x75.root](#nominal_75x75root) | 3 | Coarse nominal restricted to ≈±76 cm (19×19) |
| [nominal_75x75_2d_kde.root](#nominal_75x75_2d_kderoot) | 3 | 2D KDE fit on `nominal_75x75` MM1 |
| [nominal_corrected.root](#nominal_correctedroot) | 3 | Coarse nominal with y anomaly trimmed (25×23) |
| [nominal_corrected_2d_kde.root](#nominal_corrected_2d_kderoot) | 3 | 2D KDE fit on `nominal_corrected` MM1 |
| [proton_beam_radius_p_1sigma.root](#proton_beam_radius_p_1sigmaroot) | 3 | Proton-beam-radius +1σ alcove xy (25×25) |

Files with **no TH2s** (listed for completeness): `kde_chi2_fit.root`, `unweighted_fixed_kde.root`, `weighted_fixed_kde.root`.

---

## 2d_kde.root

- **Path:** `beam_profile_modeling/root_files/2d_kde.root`
- **Description:** Early `2d_kde.py` fit product trained on `nominalxyposMM1` from `nominal_corrected.root` (25×23 grid). Prefer newer `*_2d_kde.root` files for current campaigns.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `target_hist` | TH2D | 25 × 23 | [-100, 100] | [-100, 84] | 1, 2 | Muon x-y position at Alcove 1 |
| `kde_shape` | TH2D | 25 × 23 | [-100, 100] | [-100, 84] | 1, 2 | mirrored/unmirrored KDE shape |
| `kde_template` | TH2D | 25 × 23 | [-100, 100] | [-100, 84] | 1, 2 | α×KDE(x,y) |

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`).

---

## horn_current_p_1sigma.root

- **Path:** `beam_profile_modeling/root_files/horn_current_p_1sigma.root`
- **Description:** Systematics: horn current at +1σ. Coarse alcove MM1–MM3 xy positions (same binning as `nominal.root`).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `horn_current_p_1sigmaxyposMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 1 |
| `horn_current_p_1sigmaxyposMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 2 |
| `horn_current_p_1sigmaxyposMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 3 |

---

## jan2026_2d_kde.root

- **Path:** `beam_profile_modeling/root_files/jan2026_2d_kde.root`
- **Description:** `2d_kde.py` fit on Jan 2026 study nominal MM1 (`nominal_xypos_1`) before y-correction — full 200×200 (±125 cm). Superseded for corrected campaigns by `jan2026studies_nominal_corrected_2d_kde.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `target_hist` | TH2D | 200 × 200 | [-125, 125] | [-125, 125] | 1, 2 | Muon x-y position at Alcove 1 |
| `kde_shape` | TH2D | 200 × 200 | [-125, 125] | [-125, 125] | 1, 2 | mirrored/unmirrored KDE shape |
| `kde_template` | TH2D | 200 × 200 | [-125, 125] | [-125, 125] | 1, 2 | α×KDE(x,y) |

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`).

---

## jan2026studies_nominal.root

- **Path:** `beam_profile_modeling/root_files/jan2026studies_nominal.root`
- **Description:** Jan 2026 beam-studies nominal muon xy at alcoves 1–3. Hist names use underscore form (`nominal_xypos_*`). Input to `remove_y_anomaly.py`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `nominal_xypos_1` | TH2D | 200 × 200 | [-125, 125] | [-125, 125] | 1 | Muon x-y position at Alcove 1 |
| `nominal_xypos_2` | TH2D | 200 × 200 | [-125, 125] | [-125, 125] | 1 | Muon x-y position at Alcove 2 |
| `nominal_xypos_3` | TH2D | 200 × 200 | [-125, 125] | [-125, 125] | 1 | Muon x-y position at Alcove 3 |

---

## jan2026studies_nominal_corrected.root

- **Path:** `beam_profile_modeling/root_files/jan2026studies_nominal_corrected.root`
- **Description:** Output of `remove_y_anomaly.py` on `jan2026studies_nominal.root` — top-of-y tracking-plane anomaly removed (y ≲ 110 cm). Current training target for Jan 2026 2D KDE (`2d_kde.py`).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `nominal_xypos_1` | TH2D | 200 × 188 | [-125, 125] | [-125, 110] | 1 | Muon x-y position at Alcove 1 |
| `nominal_xypos_2` | TH2D | 200 × 188 | [-125, 125] | [-125, 110] | 1 | Muon x-y position at Alcove 2 |
| `nominal_xypos_3` | TH2D | 200 × 188 | [-125, 125] | [-125, 110] | 1 | Muon x-y position at Alcove 3 |

---

## jan2026studies_nominal_corrected_2d_kde.root

- **Path:** `beam_profile_modeling/root_files/jan2026studies_nominal_corrected_2d_kde.root`
- **Description:** Current `2d_kde.py` fit product on `nominal_xypos_1` from `jan2026studies_nominal_corrected.root`. Default KDE input for `plot_2d_kde.py` and `compare_2d_kde_to_hist.py`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `target_hist` | TH2D | 200 × 188 | [-125, 125] | [-125, 110] | 1, 2 | Muon x-y position at Alcove 1 |
| `kde_shape` | TH2D | 200 × 188 | [-125, 125] | [-125, 110] | 1, 2 | mirrored/unmirrored KDE shape |
| `kde_template` | TH2D | 200 × 188 | [-125, 125] | [-125, 110] | 1, 2 | α×KDE(x,y) |

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`).

---

## ml_tracking_hornb_tilt_down_5sigma_beamshift_nX_100um_trackingon.root

- **Path:** `beam_profile_modeling/root_files/ml_tracking_hornb_tilt_down_5sigma_beamshift_nX_100um_trackingon.root`
- **Description:** ML-tracking simulation: Horn B tilted down 5σ, beam shifted −X by 100 µm, tracking on. Alcove `ShiftxyposMM*` plus five horn/decay/HaDeS tracking-plane `ShiftnewtracksMM*` TH2s (2000×2000).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at HaDeS tracking plane |

---

## ml_tracking_hornb_tilt_down_5sigma_beamshift_nY_250um_trackingon.root

- **Path:** `beam_profile_modeling/root_files/ml_tracking_hornb_tilt_down_5sigma_beamshift_nY_250um_trackingon.root`
- **Description:** Same ML-tracking setup as nX variant, but beam shifted −Y by 250 µm. Same hist names/grid as sibling shift files.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at HaDeS tracking plane |

---

## ml_tracking_hornb_tilt_down_5sigma_beamshift_pX_100um_trackingon.root

- **Path:** `beam_profile_modeling/root_files/ml_tracking_hornb_tilt_down_5sigma_beamshift_pX_100um_trackingon.root`
- **Description:** Same ML-tracking setup, beam shifted +X by 100 µm.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at HaDeS tracking plane |

---

## ml_tracking_hornb_tilt_down_5sigma_beamshift_pY_250um_trackingon.root

- **Path:** `beam_profile_modeling/root_files/ml_tracking_hornb_tilt_down_5sigma_beamshift_pY_250um_trackingon.root`
- **Description:** Same ML-tracking setup, beam shifted +Y by 250 µm.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at HaDeS tracking plane |

---

## ml_tracking_hornb_tilt_down_5sigma_nominal_trackingon.root

- **Path:** `beam_profile_modeling/root_files/ml_tracking_hornb_tilt_down_5sigma_nominal_trackingon.root`
- **Description:** ML-tracking simulation: Horn B tilted down 5σ, nominal beam position (no shift), tracking on. Uses `NominalxyposMM*` / `NominalnewtracksMM*` names. Default comparison hist for `compare_2d_kde_to_hist.py`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `NominalxyposMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 1 |
| `NominalxyposMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 2 |
| `NominalxyposMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 3 |
| `NominalnewtracksMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 1 tracking plane |
| `NominalnewtracksMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 2 tracking plane |
| `NominalnewtracksMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Horn 3 tracking plane |
| `NominalnewtracksMM4` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Decay Pipe tracking plane |
| `NominalnewtracksMM5` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at HaDeS tracking plane |

---

## mz_nominal_100bin_run1.root

- **Path:** `beam_profile_modeling/root_files/mz_nominal_100bin_run1.root`
- **Description:** Mo’s nominal run1 alcove xy, rebinned from 2000-bin source to 100×100 via `rebin_TH2.py`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `NominalxyposMM1` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Muon x-y position at Alcove 1 |
| `NominalxyposMM2` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Muon x-y position at Alcove 2 |
| `NominalxyposMM3` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Muon x-y position at Alcove 3 |

---

## mz_nominal_100bin_run2.root

- **Path:** `beam_profile_modeling/root_files/mz_nominal_100bin_run2.root`
- **Description:** Mo’s nominal run2 alcove xy, rebinned to 100×100 (same pipeline as run1).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `NominalxyposMM1` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Muon x-y position at Alcove 1 |
| `NominalxyposMM2` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Muon x-y position at Alcove 2 |
| `NominalxyposMM3` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Muon x-y position at Alcove 3 |

---

## mz_nominal_2000bin_run1.root

- **Path:** `beam_profile_modeling/root_files/mz_nominal_2000bin_run1.root`
- **Description:** Mo’s full-resolution nominal run1 alcove xy (2000×2000, ±100 cm). Source for rebin and y-trim variants.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `NominalxyposMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 1 |
| `NominalxyposMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 2 |
| `NominalxyposMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 3 |

---

## mz_nominal_2000bin_run1_75x75.root

- **Path:** `beam_profile_modeling/root_files/mz_nominal_2000bin_run1_75x75.root`
- **Description:** Mo run1 restricted to ±75 cm (1500×1500). Used by `test_2d_kde_vs_mz_nominal.py` as high-res comparison target.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `NominalxyposMM1` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 1 |
| `NominalxyposMM2` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 2 |
| `NominalxyposMM3` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 3 |

---

## mz_nominal_2000bin_run1_corrected.root

- **Path:** `beam_profile_modeling/root_files/mz_nominal_2000bin_run1_corrected.root`
- **Description:** Mo run1 with y anomaly trimmed; stored on ±75 cm / 1500×1500 grid (same footprint as `_75x75` sibling).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `NominalxyposMM1` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 1 |
| `NominalxyposMM2` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 2 |
| `NominalxyposMM3` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 3 |

---

## mz_nominal_2000bin_run2.root

- **Path:** `beam_profile_modeling/root_files/mz_nominal_2000bin_run2.root`
- **Description:** Mo’s full-resolution nominal run2 alcove xy (2000×2000, ±100 cm).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `NominalxyposMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 1 |
| `NominalxyposMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 2 |
| `NominalxyposMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 3 |

---

## mz_nominal_2000bin_run2_75x75.root

- **Path:** `beam_profile_modeling/root_files/mz_nominal_2000bin_run2_75x75.root`
- **Description:** Mo run2 restricted to ±75 cm (1500×1500). Sibling comparison target for run2 in `test_2d_kde_vs_mz_nominal.py`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `NominalxyposMM1` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 1 |
| `NominalxyposMM2` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 2 |
| `NominalxyposMM3` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 3 |

---

## mz_nominal_2000bin_run2_corrected.root

- **Path:** `beam_profile_modeling/root_files/mz_nominal_2000bin_run2_corrected.root`
- **Description:** Mo run2 with y anomaly trimmed; ±75 cm / 1500×1500 grid.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `NominalxyposMM1` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 1 |
| `NominalxyposMM2` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 2 |
| `NominalxyposMM3` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 3 |

---

## nominal.root

- **Path:** `beam_profile_modeling/root_files/nominal.root`
- **Description:** Coarse lab nominal alcove xy (25×25, ±100 cm). Hist names lowercase `nominalxyposMM*`. Common input for early 1D KDE scripts (`weighted_fixed_kde.py`, `unweighted_fixed_kde.py`, `optimize_kde.py`).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `nominalxyposMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 1 |
| `nominalxyposMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 2 |
| `nominalxyposMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 3 |

---

## nominal_75x75.root

- **Path:** `beam_profile_modeling/root_files/nominal_75x75.root`
- **Description:** Coarse nominal restricted to ≈±76 cm (19×19). Used with `plot_resized_2d_kde.py`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `nominalxyposMM1` | TH2D | 19 × 19 | [-76, 76] | [-76, 76] | 1 | Muon x-y position at Alcove 1 |
| `nominalxyposMM2` | TH2D | 19 × 19 | [-76, 76] | [-76, 76] | 1 | Muon x-y position at Alcove 2 |
| `nominalxyposMM3` | TH2D | 19 × 19 | [-76, 76] | [-76, 76] | 1 | Muon x-y position at Alcove 3 |

---

## nominal_75x75_2d_kde.root

- **Path:** `beam_profile_modeling/root_files/nominal_75x75_2d_kde.root`
- **Description:** `2d_kde.py` fit product trained on `nominalxyposMM1` from `nominal_75x75.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `target_hist` | TH2D | 19 × 19 | [-76, 76] | [-76, 76] | 1, 2 | Muon x-y position at Alcove 1 |
| `kde_shape` | TH2D | 19 × 19 | [-76, 76] | [-76, 76] | 1, 2 | mirrored/unmirrored KDE shape |
| `kde_template` | TH2D | 19 × 19 | [-76, 76] | [-76, 76] | 1, 2 | α×KDE(x,y) |

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`).

---

## nominal_corrected.root

- **Path:** `beam_profile_modeling/root_files/nominal_corrected.root`
- **Description:** Coarse nominal with top-of-y anomaly trimmed (25×23, y to 84 cm). Training hist for early `2d_kde.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `nominalxyposMM1` | TH2D | 25 × 23 | [-100, 100] | [-100, 84] | 1 | Muon x-y position at Alcove 1 |
| `nominalxyposMM2` | TH2D | 25 × 23 | [-100, 100] | [-100, 84] | 1 | Muon x-y position at Alcove 2 |
| `nominalxyposMM3` | TH2D | 25 × 23 | [-100, 100] | [-100, 84] | 1 | Muon x-y position at Alcove 3 |

---

## nominal_corrected_2d_kde.root

- **Path:** `beam_profile_modeling/root_files/nominal_corrected_2d_kde.root`
- **Description:** `2d_kde.py` fit product on `nominalxyposMM1` from `nominal_corrected.root` (same grid as `2d_kde.root`; newer write of that campaign).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `target_hist` | TH2D | 25 × 23 | [-100, 100] | [-100, 84] | 1, 2 | Muon x-y position at Alcove 1 |
| `kde_shape` | TH2D | 25 × 23 | [-100, 100] | [-100, 84] | 1, 2 | mirrored/unmirrored KDE shape |
| `kde_template` | TH2D | 25 × 23 | [-100, 100] | [-100, 84] | 1, 2 | α×KDE(x,y) |

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`).

---

## proton_beam_radius_p_1sigma.root

- **Path:** `beam_profile_modeling/root_files/proton_beam_radius_p_1sigma.root`
- **Description:** Systematics: proton beam radius at +1σ. Coarse alcove MM1–MM3 xy (25×25).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `proton_beam_radius_p_1sigmaxyposMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 1 |
| `proton_beam_radius_p_1sigmaxyposMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 2 |
| `proton_beam_radius_p_1sigmaxyposMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 3 |

---

## Files with no TH2s

| File | Notes |
|------|-------|
| `kde_chi2_fit.root` | Contains TH1D (`target_hist`, `kde_template`), TF1 (`kde_shape`), TNamed (`fit_meta`) — 1D KDE chi² fit product |
| `unweighted_fixed_kde.root` | 1D unweighted fixed-bandwidth KDE fit product (`unweighted_fixed_kde.py`); no TH2 keys |
| `weighted_fixed_kde.root` | 1D weighted fixed-bandwidth KDE fit product (`weighted_fixed_kde.py`); no TH2 keys |

---

## How to refresh

Point `INPUT_ROOT_FILE` in `beam_profile_modeling/scripts/file_inspection.py` at a file and run:

```bash
python beam_profile_modeling/scripts/file_inspection.py
```

Or re-run inspection over all files in `beam_profile_modeling/root_files/` and update this catalog.
