# ROOT file catalog (TH2s)

Inventory of `.root` files under `beam_profile_modeling/root_files/`, with TH2 names and sizes from `beam_profile_modeling/scripts/file_inspection.py` logic (refreshed 2026-08-07).

**Naming conventions**

| Suffix / pattern | Meaning |
|------------------|---------|
| `_corrected` / `_nom_corrected` | Y-axis trimmed via `remove_y_anomaly.py` (symmetric ±`Y_MAX_CM`). Current ML 100-bin default: ±96 cm → 100×96. ML 25-bin: ±92 cm → 25×23. Jan 2026 100-bin: ±120 cm → 100×96. Older full-res Jan corrected (gone): ±120 cm → 200×192. Older Muon 100-bin corrected (gone): ±90 cm → 100×90. Older coarse corrected may be top-only (asymmetric y) |
| `_75x75` | Axis range restricted to ≈±75 cm (bin count scales with original density) |
| `_100bin` / `_25bin` | Rebinned from a denser square source via `rebin_TH2.py` (or equivalent). ML tracking: 2000→100 or 2000→25 on ±100 cm. Jan 2026: 200→100 on ±125 cm |
| `ml_beamshift_*_{100,25}bin` | Rebinned from the distinct ~279 MB Particle `ml_tracking_…beamshift_*_trackingon` samples |
| `ml_nominal_{100,25}bin` | Rebinned from current Particle `…nominal_trackingon` (keys remain `Shift*`) |
| `_2d_kde` / `*_mm1_2d_kde` | Fit product from `2d_kde.py` (`target_hist`, `kde_shape`, `kde_template` + meta/graphs) |
| `_shifted_(dx,dy)` | Full fit-product ROOT from `shift_kde.py` MODE=`shift` (new file; same key names; `kde_shape`/`kde_template` translated; `fit_meta` gains `shift_x`/`shift_y`) |
| `*_plus_*_sA(...)_sB(...)_c(...)` | `shift_kde.py` MODE=`combine`: shift each leg, then coeff-sum; `fit_meta` gains combine_* / shift_a_* / shift_b_* tags |
| `mz_nominal_*bin_run{1,2}` | Mo’s high-res `NominalxyposMM*` alcove hists; bin count in name; run1 vs run2 |
| `*_p_1sigma` | Systematics-style +1σ variation of the named parameter |

Axis ranges are in the same units as stored in the file (typically cm for alcove xy positions). Where a TH2 has multiple write cycles (`;1`, `;2`), the latest cycle is what `TFile.Get(name)` returns; sizes below match that object.

**Gotchas**

- All five `ml_tracking_hornb_tilt_down_5sigma_*_trackingon.root` files are now distinct ~279 MB **Particle** samples with `Shift*` keys (MM1 integrals ≈2.32–2.33e6). Prefer rebinned `ml_beamshift_*` / `ml_nominal_*` siblings for analysis.
- `ml_tracking_…_nominal_trackingon.root` used to be a smaller (~55 MB) “Muon …” / `Nominal*` file; that content is **gone**. Current nominal trackingon matches the Particle footprint of the beamshift siblings.
- Gone Muon-titled rebinned siblings: `ml_tracking_nominal_100bin(.corrected).root`, `negX_100um_100bin(.corrected).root`, `posX_100um_100bin(.corrected).root` — replaced by `ml_nominal_*` and `ml_beamshift_*`.
- `ml_beamshift_nY_250um_100bin_corrected.root` and `ml_beamshift_pY_250um_100bin_corrected.root` are **bit-identical** (max bin diff 0). Uncorrected nY/pY and 25-bin corrected siblings differ — regenerate the 100-bin corrected pair before trusting ±Y comparisons.
- `jan2026_mm1_2d_kde.root` is trained on `jan2026_100bin_corrected.root` (100×96). Current ML KDE defaults use `ml_nominal_mm1_2d_kde.root` on `ml_nominal_100bin_corrected.root`.

---

## Index

| File | # TH2s | Description |
|------|--------|-------------|
| [2d_kde.root](#2dkderoot) | 3 | Early 2D KDE fit on `nominal_corrected` MM1 (25×23) |
| [horn_current_p_1sigma.root](#horncurrentp1sigmaroot) | 3 | Horn-current +1σ alcove xy (25×25) |
| [jan2026_100bin.root](#jan2026100binroot) | 3 | Jan 2026 nominal rebinned to 100×100 (±125 cm) |
| [jan2026_100bin_corrected.root](#jan2026100bincorrectedroot) | 3 | Jan 2026 100-bin with symmetric y trim ±120 cm (100×96) |
| [jan2026_mm1_2d_kde.root](#jan2026mm12dkderoot) | 3 | 2D KDE fit on Jan 2026 corrected MM1 (100×96) |
| [jan2026studies_nominal.root](#jan2026studiesnominalroot) | 3 | Jan 2026 study nominal alcove xy (200×200, ±125 cm) |
| [ml_beamshift_nX_100um_100bin.root](#mlbeamshiftnx100um100binroot) | 8 | Particle −X 100 µm, rebinned to 100×100 |
| [ml_beamshift_nX_100um_100bin_corrected.root](#mlbeamshiftnx100um100bincorrectedroot) | 3 | −X 100 µm 100-bin; y trim ±96 cm (100×96) |
| [ml_beamshift_nX_100um_25bin.root](#mlbeamshiftnx100um25binroot) | 8 | Particle −X 100 µm, rebinned to 25×25 |
| [ml_beamshift_nX_100um_25bin_corrected.root](#mlbeamshiftnx100um25bincorrectedroot) | 3 | −X 100 µm 25-bin; y trim ±92 cm (25×23) |
| [ml_beamshift_nY_250um_100bin.root](#mlbeamshiftny250um100binroot) | 8 | Particle −Y 250 µm, rebinned to 100×100 |
| [ml_beamshift_nY_250um_100bin_corrected.root](#mlbeamshiftny250um100bincorrectedroot) | 3 | −Y 250 µm 100-bin; y trim ±96 cm (100×96; identical to pY corrected) |
| [ml_beamshift_nY_250um_25bin.root](#mlbeamshiftny250um25binroot) | 8 | Particle −Y 250 µm, rebinned to 25×25 |
| [ml_beamshift_nY_250um_25bin_corrected.root](#mlbeamshiftny250um25bincorrectedroot) | 3 | −Y 250 µm 25-bin; y trim ±92 cm (25×23) |
| [ml_beamshift_pX_100um_100bin.root](#mlbeamshiftpx100um100binroot) | 8 | Particle +X 100 µm, rebinned to 100×100 |
| [ml_beamshift_pX_100um_100bin_corrected.root](#mlbeamshiftpx100um100bincorrectedroot) | 3 | +X 100 µm 100-bin; y trim ±96 cm (100×96) |
| [ml_beamshift_pX_100um_25bin.root](#mlbeamshiftpx100um25binroot) | 8 | Particle +X 100 µm, rebinned to 25×25 |
| [ml_beamshift_pX_100um_25bin_corrected.root](#mlbeamshiftpx100um25bincorrectedroot) | 3 | +X 100 µm 25-bin; y trim ±92 cm (25×23) |
| [ml_beamshift_pY_250um_100bin.root](#mlbeamshiftpy250um100binroot) | 8 | Particle +Y 250 µm, rebinned to 100×100 |
| [ml_beamshift_pY_250um_100bin_corrected.root](#mlbeamshiftpy250um100bincorrectedroot) | 3 | +Y 250 µm 100-bin; y trim ±96 cm (100×96; identical to nY corrected) |
| [ml_beamshift_pY_250um_25bin.root](#mlbeamshiftpy250um25binroot) | 8 | Particle +Y 250 µm, rebinned to 25×25 |
| [ml_beamshift_pY_250um_25bin_corrected.root](#mlbeamshiftpy250um25bincorrectedroot) | 3 | +Y 250 µm 25-bin; y trim ±92 cm (25×23) |
| [ml_nominal_100bin.root](#mlnominal100binroot) | 8 | Particle ML nominal, rebinned to 100×100 (`Shift*` keys) |
| [ml_nominal_100bin_corrected.root](#mlnominal100bincorrectedroot) | 3 | ML nominal 100-bin; y trim ±96 cm (100×96) |
| [ml_nominal_25bin.root](#mlnominal25binroot) | 8 | Particle ML nominal, rebinned to 25×25 |
| [ml_nominal_25bin_corrected.root](#mlnominal25bincorrectedroot) | 3 | ML nominal 25-bin; y trim ±92 cm (25×23) |
| [ml_nominal_mm1_2d_kde.root](#mlnominalmm12dkderoot) | 3 | 2D KDE fit on ML nominal corrected MM1 (100×96) |
| [ml_nominal_mm1_2d_kde_plus_pX_100um_mm1_2d_kde_sA(-0.07625,0)_sB(0.07625,0)_c(0.5,0.5).root](#mlnominalmm12dkdepluspx100ummm12dkdesa0076250sb0076250c0505root) | 3 | combine: 0.5·shifted nominal + 0.5·shifted pX KDEs |
| [ml_tracking_hornb_tilt_down_5sigma_beamshift_nX_100um_trackingon.root](#mltrackinghornbtiltdown5sigmabeamshiftnx100umtrackingonroot) | 8 | ML tracking −X 100 µm (distinct Particle sample) |
| [ml_tracking_hornb_tilt_down_5sigma_beamshift_nY_250um_trackingon.root](#mltrackinghornbtiltdown5sigmabeamshiftny250umtrackingonroot) | 8 | ML tracking −Y 250 µm (distinct Particle sample) |
| [ml_tracking_hornb_tilt_down_5sigma_beamshift_pX_100um_trackingon.root](#mltrackinghornbtiltdown5sigmabeamshiftpx100umtrackingonroot) | 8 | ML tracking +X 100 µm (distinct Particle sample) |
| [ml_tracking_hornb_tilt_down_5sigma_beamshift_pY_250um_trackingon.root](#mltrackinghornbtiltdown5sigmabeamshiftpy250umtrackingonroot) | 8 | ML tracking +Y 250 µm (distinct Particle sample) |
| [ml_tracking_hornb_tilt_down_5sigma_nominal_trackingon.root](#mltrackinghornbtiltdown5sigmanominaltrackingonroot) | 8 | ML tracking nominal (now Particle `Shift*`, ~279 MB) |
| [mz_nominal_100bin_run1.root](#mznominal100binrun1root) | 3 | Mo nominal run1 rebinned to 100×100 |
| [mz_nominal_100bin_run2.root](#mznominal100binrun2root) | 3 | Mo nominal run2 rebinned to 100×100 |
| [mz_nominal_2000bin_run1.root](#mznominal2000binrun1root) | 3 | Mo nominal run1 full-res 2000×2000 (±100 cm) |
| [mz_nominal_2000bin_run1_75x75.root](#mznominal2000binrun175x75root) | 3 | Mo run1 restricted to ±75 cm (1500×1500) |
| [mz_nominal_2000bin_run1_corrected.root](#mznominal2000binrun1correctedroot) | 3 | Mo run1 y-trimmed (±75 cm, 1500×1500) |
| [mz_nominal_2000bin_run2.root](#mznominal2000binrun2root) | 3 | Mo nominal run2 full-res 2000×2000 (±100 cm) |
| [mz_nominal_2000bin_run2_75x75.root](#mznominal2000binrun275x75root) | 3 | Mo run2 restricted to ±75 cm (1500×1500) |
| [mz_nominal_2000bin_run2_corrected.root](#mznominal2000binrun2correctedroot) | 3 | Mo run2 y-trimmed (±75 cm, 1500×1500) |
| [nX_100um_mm1_2d_kde.root](#nx100ummm12dkderoot) | 3 | 2D KDE fit on −X 100 µm corrected MM1 (100×96) |
| [nominal.root](#nominalroot) | 3 | Coarse lab nominal alcove xy (25×25, ±100 cm) |
| [nominal_75x75.root](#nominal75x75root) | 3 | Coarse nominal restricted to ≈±76 cm (19×19) |
| [nominal_75x75_2d_kde.root](#nominal75x752dkderoot) | 3 | 2D KDE fit on `nominal_75x75` MM1 |
| [nominal_corrected.root](#nominalcorrectedroot) | 3 | Coarse nominal with y anomaly trimmed (25×23) |
| [nominal_corrected_2d_kde.root](#nominalcorrected2dkderoot) | 3 | 2D KDE fit on `nominal_corrected` MM1 |
| [pX_100um_mm1_2d_kde.root](#px100ummm12dkderoot) | 3 | 2D KDE fit on +X 100 µm corrected MM1 (100×96) |
| [proton_beam_radius_p_1sigma.root](#protonbeamradiusp1sigmaroot) | 3 | Proton-beam-radius +1σ alcove xy (25×25) |

Files with **no TH2s** (listed for completeness): `kde_chi2_fit.root`, `unweighted_fixed_kde.root`, `weighted_fixed_kde.root`.

Gone from disk since last catalog (do not use): `ml_tracking_nominal_100bin.root`, `ml_tracking_nominal_100bin_corrected.root`, `negX_100um_100bin.root`, `negX_100um_100bin_corrected.root`, `posX_100um_100bin.root`, `posX_100um_100bin_corrected.root`, `jan2026studies_nom_corrected.root`, `jan2026_mm1_2d_kde_shifted_(0.0001,0).root`, `jan2026_2d_kde.root`, `jan2026studies_nominal_corrected.root`, `jan2026studies_nominal_corrected_2d_kde.root`, and the temporary `…trackingon (1).root` Finder copies.

---

## 2d_kde.root

- **Path:** `beam_profile_modeling/root_files/2d_kde.root`
- **Description:** Early `2d_kde.py` fit product trained on `nominalxyposMM1` from `nominal_corrected.root` (25×23 grid). Prefer newer `*_2d_kde.root` files for current campaigns.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `target_hist` | TH2D | 25 × 23 | [-100, 100] | [-100, 84] | 1, 2 | Muon x-y position at Alcove 1 |
| `kde_shape` | TH2D | 25 × 23 | [-100, 100] | [-100, 84] | 1, 2 | mirrored/unmirrored KDE shape |
| `kde_template` | TH2D | 25 × 23 | [-100, 100] | [-100, 84] | 1, 2 | α×KDE(x,y) |

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`). Fit: ρ=3.5; reduced χ²≈1.530.

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

## jan2026_100bin.root

- **Path:** `beam_profile_modeling/root_files/jan2026_100bin.root`
- **Description:** Jan 2026 study nominal alcove xy rebinned to 100×100 (±125 cm) from `jan2026studies_nominal.root`. Input to `remove_y_anomaly.py` → `jan2026_100bin_corrected.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `nominal_xypos_1` | TH2D | 100 × 100 | [-125, 125] | [-125, 125] | 1, 2 | Muon x-y position at Alcove 1 |
| `nominal_xypos_2` | TH2D | 100 × 100 | [-125, 125] | [-125, 125] | 1, 2 | Muon x-y position at Alcove 2 |
| `nominal_xypos_3` | TH2D | 100 × 100 | [-125, 125] | [-125, 125] | 1, 2 | Muon x-y position at Alcove 3 |

---

## jan2026_100bin_corrected.root

- **Path:** `beam_profile_modeling/root_files/jan2026_100bin_corrected.root`
- **Description:** Symmetric y trim (±120 cm) of `jan2026_100bin.root` via `remove_y_anomaly.py` → 100×96. Training target for `jan2026_mm1_2d_kde.root`. Replaces gone `jan2026studies_nom_corrected.root` (200×192).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `nominal_xypos_1` | TH2D | 100 × 96 | [-125, 125] | [-120, 120] | 1 | Muon x-y position at Alcove 1 |
| `nominal_xypos_2` | TH2D | 100 × 96 | [-125, 125] | [-120, 120] | 1 | Muon x-y position at Alcove 2 |
| `nominal_xypos_3` | TH2D | 100 × 96 | [-125, 125] | [-120, 120] | 1 | Muon x-y position at Alcove 3 |

---

## jan2026_mm1_2d_kde.root

- **Path:** `beam_profile_modeling/root_files/jan2026_mm1_2d_kde.root`
- **Description:** `2d_kde.py` fit product on `nominal_xypos_1` from `jan2026_100bin_corrected.root` (Jan 2026 campaign).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `target_hist` | TH2D | 100 × 96 | [-125, 125] | [-120, 120] | 1, 2 | Muon x-y position at Alcove 1 |
| `kde_shape` | TH2D | 100 × 96 | [-125, 125] | [-120, 120] | 1, 2 | mirrored/unmirrored KDE shape |
| `kde_template` | TH2D | 100 × 96 | [-125, 125] | [-120, 120] | 1, 2 | α×KDE(x,y) |

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`). Fit: ρ=3.5; reduced χ²≈1.333.

---

## jan2026studies_nominal.root

- **Path:** `beam_profile_modeling/root_files/jan2026studies_nominal.root`
- **Description:** Jan 2026 beam-studies nominal muon xy at alcoves 1–3. Hist names use underscore form (`nominal_xypos_*`). Source for `jan2026_100bin.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `nominal_xypos_1` | TH2D | 200 × 200 | [-125, 125] | [-125, 125] | 1 | Muon x-y position at Alcove 1 |
| `nominal_xypos_2` | TH2D | 200 × 200 | [-125, 125] | [-125, 125] | 1 | Muon x-y position at Alcove 2 |
| `nominal_xypos_3` | TH2D | 200 × 200 | [-125, 125] | [-125, 125] | 1 | Muon x-y position at Alcove 3 |

---

## ml_beamshift_nX_100um_100bin.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_nX_100um_100bin.root`
- **Description:** All eight TH2s from `…beamshift_nX_100um_trackingon.root`, rebinned 2000→100 via `rebin_TH2.py`. “Particle …” titles; MM1 integral ≈2.33e6.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at HaDeS tracking plane |

---

## ml_beamshift_nX_100um_100bin_corrected.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_nX_100um_100bin_corrected.root`
- **Description:** Symmetric y trim (±96 cm) of alcove xy from `ml_beamshift_nX_100um_100bin.root` via `remove_y_anomaly.py` → 100×96. MM1–MM3 only (no newtracks).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 3 |

---

## ml_beamshift_nX_100um_25bin.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_nX_100um_25bin.root`
- **Description:** All eight TH2s from `…beamshift_nX_100um_trackingon.root`, rebinned 2000→25 via `rebin_TH2.py` (25×25, ±100 cm).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at HaDeS tracking plane |

---

## ml_beamshift_nX_100um_25bin_corrected.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_nX_100um_25bin_corrected.root`
- **Description:** Symmetric y trim (±92 cm) of alcove xy from `ml_beamshift_nX_100um_25bin.root` → 25×23. MM1–MM3 only. Used by `plot_hist.py` defaults.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 3 |

---

## ml_beamshift_nY_250um_100bin.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_nY_250um_100bin.root`
- **Description:** Rebinned 100×100 from `…beamshift_nY_250um_trackingon.root` via `rebin_TH2.py`. Distinct Particle sample; MM1 integral ≈2.33e6.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at HaDeS tracking plane |

---

## ml_beamshift_nY_250um_100bin_corrected.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_nY_250um_100bin_corrected.root`
- **Description:** Symmetric y trim (±96 cm) of alcove xy from `ml_beamshift_nY_250um_100bin.root` → 100×96. MM1–MM3 only. **Bit-identical** to `ml_beamshift_pY_250um_100bin_corrected.root` (suspected copy error).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 3 |

---

## ml_beamshift_nY_250um_25bin.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_nY_250um_25bin.root`
- **Description:** Rebinned 25×25 from `…beamshift_nY_250um_trackingon.root` via `rebin_TH2.py`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at HaDeS tracking plane |

---

## ml_beamshift_nY_250um_25bin_corrected.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_nY_250um_25bin_corrected.root`
- **Description:** Symmetric y trim (±92 cm) of alcove xy from `ml_beamshift_nY_250um_25bin.root` → 25×23. MM1–MM3 only. Distinct from pY 25bin corrected.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 3 |

---

## ml_beamshift_pX_100um_100bin.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_pX_100um_100bin.root`
- **Description:** Rebinned 100×100 from `…beamshift_pX_100um_trackingon.root` via `rebin_TH2.py`. Distinct Particle sample; MM1 integral ≈2.32e6.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at HaDeS tracking plane |

---

## ml_beamshift_pX_100um_100bin_corrected.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_pX_100um_100bin_corrected.root`
- **Description:** Symmetric y trim (±96 cm) of alcove xy from `ml_beamshift_pX_100um_100bin.root` → 100×96. MM1–MM3 only. Training target for `pX_100um_mm1_2d_kde.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 3 |

---

## ml_beamshift_pX_100um_25bin.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_pX_100um_25bin.root`
- **Description:** Rebinned 25×25 from `…beamshift_pX_100um_trackingon.root` via `rebin_TH2.py`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at HaDeS tracking plane |

---

## ml_beamshift_pX_100um_25bin_corrected.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_pX_100um_25bin_corrected.root`
- **Description:** Symmetric y trim (±92 cm) of alcove xy from `ml_beamshift_pX_100um_25bin.root` → 25×23. MM1–MM3 only.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 3 |

---

## ml_beamshift_pY_250um_100bin.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_pY_250um_100bin.root`
- **Description:** Rebinned 100×100 from `…beamshift_pY_250um_trackingon.root` via `rebin_TH2.py`. Distinct Particle sample; MM1 integral ≈2.33e6.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at HaDeS tracking plane |

---

## ml_beamshift_pY_250um_100bin_corrected.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_pY_250um_100bin_corrected.root`
- **Description:** Symmetric y trim (±96 cm) of alcove xy from `ml_beamshift_pY_250um_100bin.root` → 100×96. MM1–MM3 only. **Bit-identical** to `ml_beamshift_nY_250um_100bin_corrected.root` (suspected copy error).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 3 |

---

## ml_beamshift_pY_250um_25bin.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_pY_250um_25bin.root`
- **Description:** Rebinned 25×25 from `…beamshift_pY_250um_trackingon.root` via `rebin_TH2.py`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at HaDeS tracking plane |

---

## ml_beamshift_pY_250um_25bin_corrected.root

- **Path:** `beam_profile_modeling/root_files/ml_beamshift_pY_250um_25bin_corrected.root`
- **Description:** Symmetric y trim (±92 cm) of alcove xy from `ml_beamshift_pY_250um_25bin.root` → 25×23. MM1–MM3 only.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 3 |

---

## ml_nominal_100bin.root

- **Path:** `beam_profile_modeling/root_files/ml_nominal_100bin.root`
- **Description:** All eight TH2s from current `…nominal_trackingon.root` (Particle sample), rebinned 2000→100 via `rebin_TH2.py`. Keys are `Shift*` (not `Nominal*`); MM1 integral ≈2.33e6. Replaces gone `ml_tracking_nominal_100bin.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 100 × 100 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at HaDeS tracking plane |

---

## ml_nominal_100bin_corrected.root

- **Path:** `beam_profile_modeling/root_files/ml_nominal_100bin_corrected.root`
- **Description:** Symmetric y trim (±96 cm) of alcove xy from `ml_nominal_100bin.root` → 100×96. MM1–MM3 only. Current default hist for `2d_kde.py` / `compare_2d_kde_to_hist.py` / `plot_hist_ratio.py` den.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1 | Particle x-y position at Alcove 3 |

---

## ml_nominal_25bin.root

- **Path:** `beam_profile_modeling/root_files/ml_nominal_25bin.root`
- **Description:** All eight TH2s from current `…nominal_trackingon.root`, rebinned 2000→25 via `rebin_TH2.py` (25×25, ±100 cm). Keys `Shift*`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1, 2 | Particle x-y position at HaDeS tracking plane |

---

## ml_nominal_25bin_corrected.root

- **Path:** `beam_profile_modeling/root_files/ml_nominal_25bin_corrected.root`
- **Description:** Symmetric y trim (±92 cm) of alcove xy from `ml_nominal_25bin.root` → 25×23. MM1–MM3 only.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 25 × 23 | [-100, 100] | [-92, 92] | 1 | Particle x-y position at Alcove 3 |

---

## ml_nominal_mm1_2d_kde.root

- **Path:** `beam_profile_modeling/root_files/ml_nominal_mm1_2d_kde.root`
- **Description:** `2d_kde.py` fit product on `ShiftxyposMM1` from `ml_nominal_100bin_corrected.root`. Current ML-campaign KDE; default combine leg A in `shift_kde.py`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `target_hist` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1, 2 | Particle x-y position at Alcove 1 |
| `kde_shape` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1, 2 | mirrored/unmirrored KDE shape |
| `kde_template` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1, 2 | α×KDE(x,y) |

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`). Fit: ρ=3.5; reduced χ²≈1.064.

---

## ml_nominal_mm1_2d_kde_plus_pX_100um_mm1_2d_kde_sA(-0.07625,0)_sB(0.07625,0)_c(0.5,0.5).root

- **Path:** `beam_profile_modeling/root_files/ml_nominal_mm1_2d_kde_plus_pX_100um_mm1_2d_kde_sA(-0.07625,0)_sB(0.07625,0)_c(0.5,0.5).root`
- **Description:** `shift_kde.py` MODE=`combine` product: shift A (−0.07625,0) + B (+0.07625,0), then 0.5·A′+0.5·B′. Inputs `ml_nominal_mm1_2d_kde.root` + `pX_100um_mm1_2d_kde.root`. Same key layout as a 2d_kde fit ROOT.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `target_hist` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1, 2 | Particle x-y position at Alcove 1 |
| `kde_shape` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1, 2 | 0.5*kde_shape + 0.5*kde_shape |
| `kde_template` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1, 2 | 0.5*kde_template + 0.5*kde_template |

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`). Fit: ρ=3.5; reduced χ²≈1.064; combine A=`ml_nominal_mm1_2d_kde.root` B=`pX_100um_mm1_2d_kde.root`; sA=(-0.07625,0.0) sB=(0.07625,0.0); c=(0.5,0.5).

---

## ml_tracking_hornb_tilt_down_5sigma_beamshift_nX_100um_trackingon.root

- **Path:** `beam_profile_modeling/root_files/ml_tracking_hornb_tilt_down_5sigma_beamshift_nX_100um_trackingon.root`
- **Description:** Distinct ML-tracking sample: Horn B tilt −5σ + beam −X 100 µm, tracking on. ~279 MB; MM1 integral ≈2.33e6; “Particle …” titles; keys `Shift*`. Source for `ml_beamshift_nX_100um_{100,25}bin.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at HaDeS tracking plane |

---

## ml_tracking_hornb_tilt_down_5sigma_beamshift_nY_250um_trackingon.root

- **Path:** `beam_profile_modeling/root_files/ml_tracking_hornb_tilt_down_5sigma_beamshift_nY_250um_trackingon.root`
- **Description:** Distinct ML-tracking sample: Horn B tilt −5σ + beam −Y 250 µm, tracking on. ~279 MB; MM1 integral ≈2.33e6; “Particle …” titles. Source for `ml_beamshift_nY_250um_{100,25}bin.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at HaDeS tracking plane |

---

## ml_tracking_hornb_tilt_down_5sigma_beamshift_pX_100um_trackingon.root

- **Path:** `beam_profile_modeling/root_files/ml_tracking_hornb_tilt_down_5sigma_beamshift_pX_100um_trackingon.root`
- **Description:** Distinct ML-tracking sample: Horn B tilt −5σ + beam +X 100 µm, tracking on. ~279 MB; MM1 integral ≈2.32e6; “Particle …” titles. Source for `ml_beamshift_pX_100um_{100,25}bin.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at HaDeS tracking plane |

---

## ml_tracking_hornb_tilt_down_5sigma_beamshift_pY_250um_trackingon.root

- **Path:** `beam_profile_modeling/root_files/ml_tracking_hornb_tilt_down_5sigma_beamshift_pY_250um_trackingon.root`
- **Description:** Distinct ML-tracking sample: Horn B tilt −5σ + beam +Y 250 µm, tracking on. ~279 MB; MM1 integral ≈2.33e6; “Particle …” titles. Source for `ml_beamshift_pY_250um_{100,25}bin.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at HaDeS tracking plane |

---

## ml_tracking_hornb_tilt_down_5sigma_nominal_trackingon.root

- **Path:** `beam_profile_modeling/root_files/ml_tracking_hornb_tilt_down_5sigma_nominal_trackingon.root`
- **Description:** Current ML-tracking nominal beam (Horn B tilt −5σ label, tracking on). **Now** a ~279 MB Particle sample with `Shift*` keys (MM1 ≈2.33e6) — replaces the older ~55 MB “Muon …” / `Nominal*` file. Source for `ml_nominal_{100,25}bin.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `ShiftxyposMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 1 |
| `ShiftxyposMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 2 |
| `ShiftxyposMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Alcove 3 |
| `ShiftnewtracksMM1` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 1 tracking plane |
| `ShiftnewtracksMM2` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 2 tracking plane |
| `ShiftnewtracksMM3` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Horn 3 tracking plane |
| `ShiftnewtracksMM4` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at Decay Pipe tracking plane |
| `ShiftnewtracksMM5` | TH2D | 2000 × 2000 | [-100, 100] | [-100, 100] | 1 | Particle x-y position at HaDeS tracking plane |

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
- **Description:** Mo run1 with y trimmed to ±75 cm (same footprint as `_75x75` here). Older corrected naming for this campaign.

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
- **Description:** Mo run2 restricted to ±75 cm (1500×1500).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `NominalxyposMM1` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 1 |
| `NominalxyposMM2` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 2 |
| `NominalxyposMM3` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 3 |

---

## mz_nominal_2000bin_run2_corrected.root

- **Path:** `beam_profile_modeling/root_files/mz_nominal_2000bin_run2_corrected.root`
- **Description:** Mo run2 with y trimmed to ±75 cm (same footprint as `_75x75` here).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `NominalxyposMM1` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 1 |
| `NominalxyposMM2` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 2 |
| `NominalxyposMM3` | TH2D | 1500 × 1500 | [-75, 75] | [-75, 75] | 1 | Muon x-y position at Alcove 3 |

---

## nX_100um_mm1_2d_kde.root

- **Path:** `beam_profile_modeling/root_files/nX_100um_mm1_2d_kde.root`
- **Description:** `2d_kde.py` fit product on `ShiftxyposMM1` from `ml_beamshift_nX_100um_100bin_corrected.root` (100×96).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `target_hist` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1, 2 | Particle x-y position at Alcove 1 |
| `kde_shape` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1, 2 | mirrored/unmirrored KDE shape |
| `kde_template` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1, 2 | α×KDE(x,y) |

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`). Fit: ρ=3.5; reduced χ²≈1.041.

---

## nominal.root

- **Path:** `beam_profile_modeling/root_files/nominal.root`
- **Description:** Coarse lab nominal alcove MM1–MM3 xy (25×25, ±100 cm). Early / systematics comparison baseline.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `nominalxyposMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 1 |
| `nominalxyposMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 2 |
| `nominalxyposMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 3 |

---

## nominal_75x75.root

- **Path:** `beam_profile_modeling/root_files/nominal_75x75.root`
- **Description:** Coarse nominal with axes restricted to ≈±76 cm (19×19). Source for `nominal_75x75_2d_kde.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `nominalxyposMM1` | TH2D | 19 × 19 | [-76, 76] | [-76, 76] | 1 | Muon x-y position at Alcove 1 |
| `nominalxyposMM2` | TH2D | 19 × 19 | [-76, 76] | [-76, 76] | 1 | Muon x-y position at Alcove 2 |
| `nominalxyposMM3` | TH2D | 19 × 19 | [-76, 76] | [-76, 76] | 1 | Muon x-y position at Alcove 3 |

---

## nominal_75x75_2d_kde.root

- **Path:** `beam_profile_modeling/root_files/nominal_75x75_2d_kde.root`
- **Description:** `2d_kde.py` fit product on `nominalxyposMM1` from `nominal_75x75.root`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `target_hist` | TH2D | 19 × 19 | [-76, 76] | [-76, 76] | 1, 2 | Muon x-y position at Alcove 1 |
| `kde_shape` | TH2D | 19 × 19 | [-76, 76] | [-76, 76] | 1, 2 | mirrored/unmirrored KDE shape |
| `kde_template` | TH2D | 19 × 19 | [-76, 76] | [-76, 76] | 1, 2 | α×KDE(x,y) |

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`). Fit: ρ=3.5; reduced χ²≈0.476.

---

## nominal_corrected.root

- **Path:** `beam_profile_modeling/root_files/nominal_corrected.root`
- **Description:** Coarse nominal with y anomaly trimmed (asymmetric: y to +84 cm) → 25×23. Source for early `2d_kde.root` / `nominal_corrected_2d_kde.root`.

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

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`). Fit: ρ=3.5; reduced χ²≈1.530.

---

## pX_100um_mm1_2d_kde.root

- **Path:** `beam_profile_modeling/root_files/pX_100um_mm1_2d_kde.root`
- **Description:** `2d_kde.py` fit product on `ShiftxyposMM1` from `ml_beamshift_pX_100um_100bin_corrected.root` (100×96). Default combine leg B in `shift_kde.py`.

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `target_hist` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1, 2 | Particle x-y position at Alcove 1 |
| `kde_shape` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1, 2 | mirrored/unmirrored KDE shape |
| `kde_template` | TH2D | 100 × 96 | [-100, 100] | [-96, 96] | 1, 2 | α×KDE(x,y) |

Also contains: `TGraph` (`chi2_vs_rho`, `reduced_chi2_vs_rho`, `alpha_vs_rho`), `TNamed` (`fit_meta`, `stats_meta`). Fit: ρ=3.5; reduced χ²≈1.047.

---

## proton_beam_radius_p_1sigma.root

- **Path:** `beam_profile_modeling/root_files/proton_beam_radius_p_1sigma.root`
- **Description:** Systematics: proton beam radius at +1σ. Coarse alcove MM1–MM3 xy (same binning as `nominal.root`).

| Name | Class | Bins (nx × ny) | x range [cm] | y range [cm] | Cycles | Title |
|------|-------|----------------|---------|---------|--------|-------|
| `proton_beam_radius_p_1sigmaxyposMM1` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 1 |
| `proton_beam_radius_p_1sigmaxyposMM2` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 2 |
| `proton_beam_radius_p_1sigmaxyposMM3` | TH2D | 25 × 25 | [-100, 100] | [-100, 100] | 1 | Muon x-y position at Alcove 3 |
