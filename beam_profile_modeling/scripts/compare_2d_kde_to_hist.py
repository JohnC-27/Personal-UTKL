#!/usr/bin/env python3
"""
Overlay the 2D KDE from a fit ROOT file on a general TH2 histogram.

Display-only: rebins and/or restricts the histogram in memory for plotting and
statistics. Does not write or modify any ROOT files.

Plain fit products: rebuild RooNDKeysPdf from target_hist + fit_meta. Overlay
SURF and projection curves evaluate the PDF on a fine grid (KDE_OVERLAY_POINTS /
KDE_PROJECTION_POINTS). Chi2 uses PDF values at data-bin centers.

Shifted / combined products (shift_kde.py): target_hist is not the model.
Use stored kde_template (bin copy or Interpolate) so means/χ² match
shift_kde; live rebuild would ignore shift_*/combine_*.

Outputs (plot_2d_kde.py style):
  - overlay: data + alpha*KDE surface, KDE alone, region stats (WLS alpha)
  - projections: X/Y with Data/KDE ratio and per-bin chi2; KDE uses
    integral-matched alpha for display (curve, ratio, projection chi2/stats)
  - chi2 contribution lego (standalone; 2D WLS alpha)
"""

import array
import os
import sys
from dataclasses import dataclass

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# ===============================SCRIPT PARAMS===============================
# target hist file and hist name
HIST_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "ml_beamshift_nX_100um_100bin_corrected.root"
)
HIST_NAME = "ShiftxyposMM1"


KDE_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "ml_nominal_mm1_2d_kde_plus_pX_100um_mm1_2d_kde_sA(-0.0912,0)_sB(0.06108,0)_c(0.5,0.5).root"
)

# Rebin option here for hists >100 bins each axis or so
# Rebin by combining adjacent bins (1 = no rebin). Must evenly divide source bins.
# int -> original_nbins / REBIN_FACTOR = new_nbins
REBIN_FACTOR_X = 1
REBIN_FACTOR_Y = 1

# Restriction for getting the "good area" of the KDE (roughly the middle 80% of each axis)
# After rebinning, keep a centered Nx x Ny window for plot + stats (None = all).
# int = bins to keep on that axis after rebin (centered crop so oringinal nbins must be even)
RESTRICT_NBINS_X = 60
RESTRICT_NBINS_Y = 60

OUTPUT_TAG = "2d_kde_vs"
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "plots")

# Plot titles (projections / chi2 append suffixes below)
PLOT_TITLE = "Combined -100 and +100um shifts vs nominal hist"
PLOT_TITLE_OVERLAY = f"{PLOT_TITLE}"
PLOT_TITLE_KDE_TEMPLATE = "KDE"
PLOT_TITLE_X_PROJ = f"{PLOT_TITLE} x projection"
PLOT_TITLE_Y_PROJ = f"{PLOT_TITLE} y projection"
PLOT_TITLE_CHI2 = f"{PLOT_TITLE} #chi^{{2}} contribution"

# Direct RooNDKeysPdf evaluation grid for 2D surfaces.
# kde eval and point plotting scales as this squared
KDE_OVERLAY_POINTS = 200 # PER AXIS -> 100 by 100 = 10,000 points over the plane
# Continuous projection curves: PDF marginalization sample count per axis.
KDE_PROJECTION_POINTS = 200

RATIO_Z_PAD = 1.05
RATIO_Z_MIN_HALF_WIDTH = 0.05

NDKEYS_NO_MIRROR = "a"
NDKEYS_MIRROR_BOTH = "am"
# ===========================================================================


@dataclass
class Th2RegionStats:
  integral: float
  mean_x: float
  mean_y: float
  std_x: float
  std_y: float


@dataclass
class Th1ProjectionStats:
  integral_data: float
  integral_kde: float
  mean_data: float
  mean_kde: float
  chi2: float
  ndf: int
  mean_label: str

  @property
  def reduced_chi2(self) -> float:
    return self.chi2 / max(self.ndf, 1)


@dataclass
class KdeEvalContext:
  x_var: ROOT.RooRealVar
  y_var: ROOT.RooRealVar
  argset: ROOT.RooArgSet
  dataset: ROOT.RooDataSet


@dataclass
class KdeModel:
  ctx: KdeEvalContext
  alpha: float
  mix: float
  use_linear_combo: bool
  pdf_single: ROOT.RooNDKeysPdf | None = None
  pdf_unmirrored: ROOT.RooNDKeysPdf | None = None
  pdf_mirrored: ROOT.RooNDKeysPdf | None = None

  def shape_at(self, x: float, y: float) -> float:
    self.ctx.x_var.setVal(x)
    self.ctx.y_var.setVal(y)
    if self.use_linear_combo:
      u = self.pdf_unmirrored.getVal(self.ctx.argset)
      m = self.pdf_mirrored.getVal(self.ctx.argset)
      return self.mix * u + (1.0 - self.mix) * m
    return self.pdf_single.getVal(self.ctx.argset)

  def scaled_at(self, x: float, y: float) -> float:
    return self.alpha * self.shape_at(x, y)


@dataclass
class StoredShapeModel:
  """Model backed by stored kde_template (for shift/combine fit products)."""

  shape: ROOT.TH2
  alpha: float = 1.0

  def shape_at(self, x: float, y: float) -> float:
    xaxis = self.shape.GetXaxis()
    yaxis = self.shape.GetYaxis()
    if not (xaxis.GetXmin() <= x <= xaxis.GetXmax()):
      return 0.0
    if not (yaxis.GetXmin() <= y <= yaxis.GetXmax()):
      return 0.0
    ix = xaxis.FindFixBin(x)
    iy = yaxis.FindFixBin(y)
    nbx = self.shape.GetNbinsX()
    nby = self.shape.GetNbinsY()
    if ix < 1 or ix > nbx or iy < 1 or iy > nby:
      return 0.0
    # ROOT Interpolate is undefined near the outer bin rim; use bin content there.
    if ix == 1 or ix == nbx or iy == 1 or iy == nby:
      return float(self.shape.GetBinContent(ix, iy))
    return float(self.shape.Interpolate(x, y))

  def scaled_at(self, x: float, y: float) -> float:
    return self.alpha * self.shape_at(x, y)


def fit_has_spatial_ops(meta: dict[str, float | str | int]) -> bool:
  """True if fit_meta records a shift_kde shift and/or combine."""
  keys = (
    "shift_x",
    "shift_y",
    "shift_a_x",
    "shift_a_y",
    "shift_b_x",
    "shift_b_y",
    "combine_coeff_a",
    "combine_coeff_b",
  )
  return any(key in meta for key in keys)


def parse_fit_meta(meta: ROOT.TNamed) -> dict[str, float | str | int]:
  out: dict[str, float | str | int] = {}
  for part in meta.GetTitle().split(";"):
    key, val = part.split("=", 1)
    if key == "pdf":
      out[key] = val
      continue
    try:
      num = float(val)
      out[key] = int(num) if num.is_integer() and "." not in val else num
    except ValueError:
      out[key] = val
  return out


def load_kde_fit(
  filepath: str,
) -> tuple[ROOT.TH2, dict[str, float | str | int], ROOT.TH2, ROOT.TH2]:
  """Load target_hist, fit_meta, kde_shape, and kde_template from a fit product."""
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open {filepath}")

  train_target = tfile.Get("target_hist")
  meta = tfile.Get("fit_meta")
  kde_shape = tfile.Get("kde_shape")
  kde_template = tfile.Get("kde_template")
  if not train_target or not meta or not kde_shape or not kde_template:
    tfile.Close()
    raise KeyError(
      f"missing target_hist, fit_meta, kde_shape, or kde_template in {filepath}"
    )

  train_target = train_target.Clone("kde_train_target")
  train_target.SetDirectory(0)
  kde_shape = kde_shape.Clone("kde_shape_compare")
  kde_shape.SetDirectory(0)
  kde_template = kde_template.Clone("kde_template_compare")
  kde_template.SetDirectory(0)
  meta_dict = parse_fit_meta(meta)
  tfile.Close()
  return train_target, meta_dict, kde_shape, kde_template


def th2_axes_match(a: ROOT.TH2, b: ROOT.TH2) -> bool:
  if a.GetNbinsX() != b.GetNbinsX() or a.GetNbinsY() != b.GetNbinsY():
    return False
  ax, bx = a.GetXaxis(), b.GetXaxis()
  ay, by = a.GetYaxis(), b.GetYaxis()
  edges = (
    (ax.GetXmin(), bx.GetXmin()),
    (ax.GetXmax(), bx.GetXmax()),
    (ay.GetXmin(), by.GetXmin()),
    (ay.GetXmax(), by.GetXmax()),
  )
  for va, vb in edges:
    if abs(va - vb) > 1e-9 * max(1.0, abs(va), abs(vb)):
      return False
  return True


def copy_th2_contents(source: ROOT.TH2, ref_hist: ROOT.TH2, name: str) -> ROOT.TH2D:
  """Clone ref_hist binning and copy source bin contents (axes must match)."""
  out = ref_hist.Clone(name)
  out.SetDirectory(0)
  out.Reset()
  out.SetTitle("KDE shape on data grid")
  for ix in range(1, out.GetNbinsX() + 1):
    for iy in range(1, out.GetNbinsY() + 1):
      out.SetBinContent(ix, iy, source.GetBinContent(ix, iy))
  return out


def histogram_to_weighted_dataset(hist: ROOT.TH2) -> KdeEvalContext:
  xlo = hist.GetXaxis().GetXmin()
  xhi = hist.GetXaxis().GetXmax()
  ylo = hist.GetYaxis().GetXmin()
  yhi = hist.GetYaxis().GetXmax()
  x_var = ROOT.RooRealVar("plot_x", "x [cm]", xlo, xhi)
  y_var = ROOT.RooRealVar("plot_y", "y [cm]", ylo, yhi)
  argset = ROOT.RooArgSet(x_var, y_var)
  w_var = ROOT.RooRealVar("plot_w", "weight", 0.0, 1.0e20)
  dataset = ROOT.RooDataSet(
    "plot_weighted_points",
    "plot_weighted_points",
    argset,
    ROOT.RooFit.WeightVar(w_var),
  )
  for ix in range(1, hist.GetNbinsX() + 1):
    x_var.setVal(hist.GetXaxis().GetBinCenter(ix))
    for iy in range(1, hist.GetNbinsY() + 1):
      content = hist.GetBinContent(ix, iy)
      if content <= 0:
        continue
      y_var.setVal(hist.GetYaxis().GetBinCenter(iy))
      w_var.setVal(content)
      dataset.add(argset, content)
  return KdeEvalContext(x_var=x_var, y_var=y_var, argset=argset, dataset=dataset)


def make_ndkeys_pdf(
  name: str,
  ctx: KdeEvalContext,
  *,
  mirror_options: str,
  rho: float,
) -> ROOT.RooNDKeysPdf:
  return ROOT.RooNDKeysPdf(
    name, name, ctx.argset, ctx.dataset, mirror_options, float(rho)
  )


def build_kde_model(train_target: ROOT.TH2, meta: dict[str, float | str | int]) -> KdeModel:
  ctx = histogram_to_weighted_dataset(train_target)
  rho = float(meta["rho"])
  alpha = float(meta["alpha"])
  use_linear_combo = bool(meta.get("linear_combo", 0))
  mix = float(meta.get("mix", 1.0))
  opt_no = str(meta.get("ndkeys_no_mirror", NDKEYS_NO_MIRROR))
  opt_m = str(meta.get("ndkeys_mirror", NDKEYS_MIRROR_BOTH))
  if use_linear_combo:
    return KdeModel(
      ctx=ctx,
      alpha=alpha,
      mix=mix,
      use_linear_combo=True,
      pdf_unmirrored=make_ndkeys_pdf(
        "compare_kde_unmirrored", ctx, mirror_options=opt_no, rho=rho
      ),
      pdf_mirrored=make_ndkeys_pdf(
        "compare_kde_mirrored", ctx, mirror_options=opt_m, rho=rho
      ),
    )
  return KdeModel(
    ctx=ctx,
    alpha=alpha,
    mix=1.0,
    use_linear_combo=False,
    pdf_single=make_ndkeys_pdf(
      "compare_kde_single", ctx, mirror_options=opt_m, rho=rho
    ),
  )


def evaluate_kde_th2(
  model: KdeModel | StoredShapeModel,
  *,
  n_bins_x: int,
  n_bins_y: int,
  xlo: float,
  xhi: float,
  ylo: float,
  yhi: float,
  name: str,
  title: str,
  scaled: bool = True,
) -> ROOT.TH2D:
  out = ROOT.TH2D(name, title, n_bins_x, xlo, xhi, n_bins_y, ylo, yhi)
  out.SetDirectory(0)
  out.SetStats(0)
  out._hold_model = model
  value_at = model.scaled_at if scaled else model.shape_at
  for ix in range(1, n_bins_x + 1):
    x = out.GetXaxis().GetBinCenter(ix)
    for iy in range(1, n_bins_y + 1):
      y = out.GetYaxis().GetBinCenter(iy)
      out.SetBinContent(ix, iy, value_at(x, y))
  return out


def kde_fine_plot_hist(
  model: KdeModel | StoredShapeModel,
  ref_hist: ROOT.TH2,
  name: str,
) -> ROOT.TH2D:
  """Evaluate model on a fine grid spanning the displayed histogram axes."""
  xtitle = ref_hist.GetXaxis().GetTitle() or "x [cm]"
  ytitle = ref_hist.GetYaxis().GetTitle() or "y [cm]"
  return evaluate_kde_th2(
    model,
    n_bins_x=KDE_OVERLAY_POINTS,
    n_bins_y=KDE_OVERLAY_POINTS,
    xlo=ref_hist.GetXaxis().GetXmin(),
    xhi=ref_hist.GetXaxis().GetXmax(),
    ylo=ref_hist.GetYaxis().GetXmin(),
    yhi=ref_hist.GetYaxis().GetXmax(),
    name=name,
    title=f"#alpha#timesKDE(x,y);{xtitle};{ytitle}",
    scaled=True,
  )


def open_histogram2d(filepath: str, hist_name: str) -> ROOT.TH2:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open ROOT file: {filepath}")

  obj = tfile.Get(hist_name)
  if not obj:
    tfile.Close()
    raise KeyError(f"object {hist_name!r} not found in {filepath}")
  if not obj.InheritsFrom("TH2"):
    tfile.Close()
    raise TypeError(f"{hist_name!r} is not TH2")

  hist = obj.Clone(f"{hist_name}_compare")
  hist.SetDirectory(0)
  tfile.Close()
  return hist


def rebin_histogram2d(hist: ROOT.TH2, factor_x: int, factor_y: int) -> ROOT.TH2:
  if factor_x < 1 or factor_y < 1:
    raise ValueError("REBIN_FACTOR_* must be >= 1")
  if factor_x == 1 and factor_y == 1:
    return hist

  nbx = hist.GetNbinsX()
  nby = hist.GetNbinsY()
  if nbx % factor_x != 0:
    raise ValueError(f"REBIN_FACTOR_X={factor_x} must evenly divide x bins ({nbx})")
  if nby % factor_y != 0:
    raise ValueError(f"REBIN_FACTOR_Y={factor_y} must evenly divide y bins ({nby})")

  rebinned = hist.Rebin2D(factor_x, factor_y, f"{hist.GetName()}_rebin")
  rebinned.SetDirectory(0)
  return rebinned


def _centered_bin_window(n_total: int, n_keep: int | None) -> tuple[int, int]:
  if n_keep is None or n_keep <= 0 or n_keep >= n_total:
    return 1, n_total
  if n_keep > n_total:
    raise ValueError(f"cannot keep {n_keep} bins from axis with {n_total}")
  first = 1 + (n_total - n_keep) // 2
  return first, first + n_keep - 1


def restrict_histogram2d(
  hist: ROOT.TH2,
  n_bins_x: int | None,
  n_bins_y: int | None,
) -> ROOT.TH2:
  """Return an in-memory TH2 containing only the centered restrict window."""
  ix_lo, ix_hi = _centered_bin_window(hist.GetNbinsX(), n_bins_x)
  iy_lo, iy_hi = _centered_bin_window(hist.GetNbinsY(), n_bins_y)
  if ix_lo == 1 and ix_hi == hist.GetNbinsX() and iy_lo == 1 and iy_hi == hist.GetNbinsY():
    return hist

  xaxis = hist.GetXaxis()
  yaxis = hist.GetYaxis()
  x_edges = [xaxis.GetBinLowEdge(i) for i in range(ix_lo, ix_hi + 1)]
  x_edges.append(xaxis.GetBinUpEdge(ix_hi))
  y_edges = [yaxis.GetBinLowEdge(i) for i in range(iy_lo, iy_hi + 1)]
  y_edges.append(yaxis.GetBinUpEdge(iy_hi))

  out = ROOT.TH2D(
    f"{hist.GetName()}_restrict",
    hist.GetTitle(),
    len(x_edges) - 1,
    array.array("d", x_edges),
    len(y_edges) - 1,
    array.array("d", y_edges),
  )
  out.SetDirectory(0)
  out.GetXaxis().SetTitle(xaxis.GetTitle())
  out.GetYaxis().SetTitle(yaxis.GetTitle())
  out.GetZaxis().SetTitle(hist.GetZaxis().GetTitle())
  if hist.GetSumw2N() == 0:
    hist.Sumw2()
  out.Sumw2()

  for ix_new, ix in enumerate(range(ix_lo, ix_hi + 1), start=1):
    for iy_new, iy in enumerate(range(iy_lo, iy_hi + 1), start=1):
      out.SetBinContent(ix_new, iy_new, hist.GetBinContent(ix, iy))
      out.SetBinError(ix_new, iy_new, hist.GetBinError(ix, iy))

  print(
    f"Restricted display/stats region to "
    f"{out.GetNbinsX()}x{out.GetNbinsY()} bins "
    f"(source bins x[{ix_lo},{ix_hi}] y[{iy_lo},{iy_hi}])"
  )
  return out


def th2_region_stats(hist: ROOT.TH2) -> Th2RegionStats:
  total_w = 0.0
  sum_x = 0.0
  sum_y = 0.0
  for ix in range(1, hist.GetNbinsX() + 1):
    x = hist.GetXaxis().GetBinCenter(ix)
    for iy in range(1, hist.GetNbinsY() + 1):
      w = hist.GetBinContent(ix, iy)
      if w <= 0:
        continue
      y = hist.GetYaxis().GetBinCenter(iy)
      total_w += w
      sum_x += w * x
      sum_y += w * y

  if total_w <= 0:
    return Th2RegionStats(0.0, 0.0, 0.0, 0.0, 0.0)

  mean_x = sum_x / total_w
  mean_y = sum_y / total_w
  var_x = 0.0
  var_y = 0.0
  for ix in range(1, hist.GetNbinsX() + 1):
    x = hist.GetXaxis().GetBinCenter(ix)
    dx2 = (x - mean_x) ** 2
    for iy in range(1, hist.GetNbinsY() + 1):
      w = hist.GetBinContent(ix, iy)
      if w <= 0:
        continue
      y = hist.GetYaxis().GetBinCenter(iy)
      var_x += w * dx2
      var_y += w * (y - mean_y) ** 2

  return Th2RegionStats(
    integral=total_w,
    mean_x=mean_x,
    mean_y=mean_y,
    std_x=(var_x / total_w) ** 0.5,
    std_y=(var_y / total_w) ** 0.5,
  )


def optimal_alpha(kde_shape: ROOT.TH2, hist: ROOT.TH2) -> float:
  num, den = 0.0, 0.0
  for ix in range(1, hist.GetNbinsX() + 1):
    for iy in range(1, hist.GetNbinsY() + 1):
      err = hist.GetBinError(ix, iy)
      if err <= 0:
        continue
      observed = hist.GetBinContent(ix, iy)
      shape = kde_shape.GetBinContent(ix, iy)
      w = 1.0 / (err * err)
      num += w * observed * shape
      den += w * shape * shape
  if den <= 0:
    return 1.0
  return num / den


def chi_squared_vs_hist(
  kde_shape: ROOT.TH2,
  hist: ROOT.TH2,
  alpha: float,
) -> tuple[float, int]:
  chi2 = 0.0
  n_used = 0
  for ix in range(1, hist.GetNbinsX() + 1):
    for iy in range(1, hist.GetNbinsY() + 1):
      err = hist.GetBinError(ix, iy)
      if err <= 0:
        continue
      observed = hist.GetBinContent(ix, iy)
      expected = alpha * kde_shape.GetBinContent(ix, iy)
      diff = observed - expected
      chi2 += (diff * diff) / (err * err)
      n_used += 1
  return chi2, n_used


def expand_model_range(model: KdeModel | StoredShapeModel, ref_hist: ROOT.TH2) -> None:
  """Allow live PDF evaluation over the displayed histogram axes."""
  if not isinstance(model, KdeModel):
    return
  xlo = min(model.ctx.x_var.getMin(), ref_hist.GetXaxis().GetXmin())
  xhi = max(model.ctx.x_var.getMax(), ref_hist.GetXaxis().GetXmax())
  ylo = min(model.ctx.y_var.getMin(), ref_hist.GetYaxis().GetXmin())
  yhi = max(model.ctx.y_var.getMax(), ref_hist.GetYaxis().GetXmax())
  model.ctx.x_var.setRange(xlo, xhi)
  model.ctx.y_var.setRange(ylo, yhi)


def fill_kde_shape_on_hist_grid(
  model: KdeModel | StoredShapeModel,
  ref_hist: ROOT.TH2,
  name: str,
) -> ROOT.TH2D:
  """Evaluate the KDE shape at each data-bin center (for chi2 / projections)."""
  if isinstance(model, StoredShapeModel) and th2_axes_match(model.shape, ref_hist):
    out = copy_th2_contents(model.shape, ref_hist, name)
    out._hold_model = model
    return out
  out = ref_hist.Clone(name)
  out.SetDirectory(0)
  out.Reset()
  out.SetTitle("KDE shape on data grid")
  out._hold_model = model
  for ix in range(1, out.GetNbinsX() + 1):
    x = out.GetXaxis().GetBinCenter(ix)
    for iy in range(1, out.GetNbinsY() + 1):
      y = out.GetYaxis().GetBinCenter(iy)
      out.SetBinContent(ix, iy, model.shape_at(x, y))
  return out


def kde_template_histogram(
  kde_shape: ROOT.TH2,
  ref_hist: ROOT.TH2,
  alpha: float,
  name: str,
) -> ROOT.TH2D:
  template = kde_shape.Clone(name)
  template.SetDirectory(0)
  template.Scale(alpha)
  xtitle = ref_hist.GetXaxis().GetTitle() or "x [cm]"
  ytitle = ref_hist.GetYaxis().GetTitle() or "y [cm]"
  template.SetTitle(f"#alpha#timesKDE(x,y);{xtitle};{ytitle}")
  return template


def make_chi2_contrib_th2(
  data: ROOT.TH2,
  template: ROOT.TH2,
  name: str,
) -> ROOT.TH2D:
  contrib = data.Clone(name)
  contrib.SetDirectory(0)
  contrib.Reset()
  contrib.SetStats(0)
  xtitle = data.GetXaxis().GetTitle() or "x [cm]"
  ytitle = data.GetYaxis().GetTitle() or "y [cm]"
  contrib.SetTitle(
    f"{PLOT_TITLE_CHI2};{xtitle};{ytitle};#chi^{{2}}_{{i}}"
  )
  for ix in range(1, data.GetNbinsX() + 1):
    for iy in range(1, data.GetNbinsY() + 1):
      err = data.GetBinError(ix, iy)
      if err <= 0:
        continue
      diff = data.GetBinContent(ix, iy) - template.GetBinContent(ix, iy)
      contrib.SetBinContent(ix, iy, (diff * diff) / (err * err))
  return contrib


def projection_with_errors(hist: ROOT.TH2, axis: str, name: str) -> ROOT.TH1:
  if hist.GetSumw2N() == 0:
    hist.Sumw2()
  if axis == "x":
    proj = hist.ProjectionX(name, 1, hist.GetNbinsY(), "e")
  elif axis == "y":
    proj = hist.ProjectionY(name, 1, hist.GetNbinsX(), "e")
  else:
    raise ValueError(f"axis must be 'x' or 'y', got {axis!r}")
  proj.SetDirectory(0)
  if proj.GetSumw2N() == 0:
    proj.Sumw2()
  return proj


def _hist_axis_centers(axis: ROOT.TAxis) -> list[float]:
  return [axis.GetBinCenter(i) for i in range(1, axis.GetNbins() + 1)]


def _axis_linspace(lo: float, hi: float, n_points: int) -> list[float]:
  if n_points <= 1:
    return [(lo + hi) / 2.0]
  step = (hi - lo) / (n_points - 1)
  pts = [min(max(lo + i * step, lo), hi) for i in range(n_points)]
  pts[0] = lo
  pts[-1] = hi
  return pts


def kde_projection_x_curve(
  model: KdeModel | StoredShapeModel,
  ref_hist: ROOT.TH2,
  name: str,
) -> ROOT.TGraph:
  x_values = _axis_linspace(
    ref_hist.GetXaxis().GetXmin(),
    ref_hist.GetXaxis().GetXmax(),
    KDE_PROJECTION_POINTS,
  )
  y_centers = _hist_axis_centers(ref_hist.GetYaxis())
  graph = ROOT.TGraph(len(x_values))
  graph.SetName(name)
  graph._hold_model = model
  for i, x in enumerate(x_values):
    graph.SetPoint(i, x, sum(model.scaled_at(x, y) for y in y_centers))
  return graph


def kde_projection_y_curve(
  model: KdeModel | StoredShapeModel,
  ref_hist: ROOT.TH2,
  name: str,
) -> ROOT.TGraph:
  y_values = _axis_linspace(
    ref_hist.GetYaxis().GetXmin(),
    ref_hist.GetYaxis().GetXmax(),
    KDE_PROJECTION_POINTS,
  )
  x_centers = _hist_axis_centers(ref_hist.GetXaxis())
  graph = ROOT.TGraph(len(y_values))
  graph.SetName(name)
  graph._hold_model = model
  for i, y in enumerate(y_values):
    graph.SetPoint(i, y, sum(model.scaled_at(x, y) for x in x_centers))
  return graph


def _integral_match_scale(data_integral: float, model_integral: float) -> float:
  """Scale factor so model integral matches data (display α_int / α_WLS)."""
  if model_integral <= 0:
    return 1.0
  return data_integral / model_integral


def _scale_graph_y(graph: ROOT.TGraph, scale: float) -> None:
  if scale == 1.0:
    return
  for i in range(graph.GetN()):
    graph.SetPoint(i, graph.GetPointX(i), graph.GetPointY(i) * scale)


def _print_comparison_stats(
  label: str,
  data: ROOT.TH2,
  alpha: float,
  n_bins: int,
  chi2: float,
  ndf: int,
  reduced: float,
  hist_stats: Th2RegionStats,
  kde_stats: Th2RegionStats,
) -> None:
  print(f"\n{label} (restricted region):")
  print(f"  bins: {data.GetNbinsX()}x{data.GetNbinsY()}")
  print(f"  bins used in chi2: {n_bins}")
  print(f"  profiled alpha: {alpha:.6g}")
  print(
    f"  data  integral={hist_stats.integral:.6g}  "
    f"mean=({hist_stats.mean_x:.6g}, {hist_stats.mean_y:.6g})  "
    f"std=({hist_stats.std_x:.6g}, {hist_stats.std_y:.6g})"
  )
  print(
    f"  model integral={kde_stats.integral:.6g}  "
    f"mean=({kde_stats.mean_x:.6g}, {kde_stats.mean_y:.6g})  "
    f"std=({kde_stats.std_x:.6g}, {kde_stats.std_y:.6g})"
  )
  print(f"  chi2 = {chi2:.6g}, ndf = {ndf}, chi2/ndf = {reduced:.6g}")


def evaluate_comparison(
  label: str,
  data: ROOT.TH2,
  model: KdeModel | StoredShapeModel,
) -> dict:
  if data.GetSumw2N() == 0:
    data.Sumw2()
  expand_model_range(model, data)

  shape_on_grid = fill_kde_shape_on_hist_grid(
    model, data, f"kde_shape_on_grid_{label}"
  )
  alpha = optimal_alpha(shape_on_grid, data)
  model.alpha = alpha
  chi2, n_bins = chi_squared_vs_hist(shape_on_grid, data, alpha)
  ndf = max(n_bins - 1, 1)
  reduced = chi2 / ndf
  template = kde_template_histogram(
    shape_on_grid, data, alpha, f"kde_template_{label}"
  )
  hist_stats = th2_region_stats(data)
  kde_stats = th2_region_stats(template)
  _print_comparison_stats(
    label, data, alpha, n_bins, chi2, ndf, reduced, hist_stats, kde_stats
  )
  return {
    "alpha": alpha,
    "chi2": chi2,
    "ndf": float(ndf),
    "reduced_chi2": reduced,
    "hist_stats": hist_stats,
    "kde_stats": kde_stats,
    "template": template,
    "model": model,
    "proj_x_data": projection_with_errors(data, "x", f"data_px_{label}"),
    "proj_y_data": projection_with_errors(data, "y", f"data_py_{label}"),
    "proj_x_model": projection_with_errors(template, "x", f"model_px_{label}"),
    "proj_y_model": projection_with_errors(template, "y", f"model_py_{label}"),
    "curve_x": kde_projection_x_curve(model, data, f"curve_px_{label}"),
    "curve_y": kde_projection_y_curve(model, data, f"curve_py_{label}"),
  }


def _axis_titles(hist: ROOT.TH2) -> tuple[str, str, str]:
  return (
    hist.GetXaxis().GetTitle() or "x [cm]",
    hist.GetYaxis().GetTitle() or "y [cm]",
    hist.GetZaxis().GetTitle() or "Entries",
  )


def _style_surf3d_axes(hist: ROOT.TH2) -> tuple[str, str, str]:
  xtitle, ytitle, ztitle = _axis_titles(hist)
  for axis in (hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis()):
    axis.SetTitle("")
    axis.SetLabelSize(0.022)
    axis.SetNdivisions(505)
  return xtitle, ytitle, ztitle


def _draw_surf3d_axis_titles(xtitle: str, ytitle: str, ztitle: str) -> None:
  latex = ROOT.TLatex()
  latex.SetNDC(True)
  latex.SetTextFont(42)
  latex.SetTextSize(0.032)
  latex.SetTextAlign(22)
  latex.DrawLatex(0.74, 0.05, xtitle)
  latex.DrawLatex(0.26, 0.05, ytitle)
  latex.SetTextAngle(90)
  latex.DrawLatex(0.97, 0.54, ztitle)


def _configure_surf_pad(pad: ROOT.TPad) -> None:
  pad.SetGrid()
  pad.SetTheta(28)
  pad.SetPhi(60)
  pad.SetLeftMargin(0.07)
  pad.SetRightMargin(0.07)
  pad.SetBottomMargin(0.05)
  pad.SetTopMargin(0.10)


def _style_surf_hist(
  hist: ROOT.TH2,
  *,
  line_color: int,
  line_width: int = 1,
) -> tuple[str, str, str]:
  hist.SetLineColor(line_color)
  hist.SetLineWidth(line_width)
  hist.SetFillStyle(0)
  return _style_surf3d_axes(hist)


def _graph_max_y(graph: ROOT.TGraph) -> float:
  ymax = 0.0
  for i in range(graph.GetN()):
    ymax = max(ymax, graph.GetY()[i])
  return ymax


def _hist_1d_stats(hist: ROOT.TH1) -> tuple[float, float]:
  total = 0.0
  weighted = 0.0
  for i in range(1, hist.GetNbinsX() + 1):
    content = hist.GetBinContent(i)
    if content <= 0:
      continue
    total += content
    weighted += content * hist.GetXaxis().GetBinCenter(i)
  if total <= 0:
    return 0.0, 0.0
  return total, weighted / total


def _compute_chi2_1d(data: ROOT.TH1, model: ROOT.TH1) -> tuple[float, int]:
  chi2 = 0.0
  ndf = 0
  for i in range(1, data.GetNbinsX() + 1):
    err = data.GetBinError(i)
    if err <= 0:
      continue
    diff = data.GetBinContent(i) - model.GetBinContent(i)
    chi2 += (diff * diff) / (err * err)
    ndf += 1
  return chi2, ndf


def _projection_stats(
  data: ROOT.TH1,
  model: ROOT.TH1,
  mean_label: str,
) -> Th1ProjectionStats:
  chi2, ndf = _compute_chi2_1d(data, model)
  int_data, mean_data = _hist_1d_stats(data)
  int_kde, mean_kde = _hist_1d_stats(model)
  return Th1ProjectionStats(
    integral_data=int_data,
    integral_kde=int_kde,
    mean_data=mean_data,
    mean_kde=mean_kde,
    chi2=chi2,
    ndf=ndf,
    mean_label=mean_label,
  )


def _make_ratio_hist(data: ROOT.TH1, model: ROOT.TH1, name: str) -> ROOT.TH1:
  ratio = data.Clone(name)
  ratio.SetDirectory(0)
  if ratio.GetSumw2N() == 0:
    ratio.Sumw2()
  ratio.Divide(model)
  ratio.SetTitle("")
  ratio.SetStats(0)
  ratio.SetMarkerSize(data.GetMarkerSize())
  ratio.SetMarkerStyle(data.GetMarkerStyle())
  ratio.SetMarkerColor(data.GetMarkerColor())
  ratio.SetLineColor(data.GetLineColor())
  ratio.SetLineWidth(data.GetLineWidth())
  ratio.GetYaxis().SetTitle("Data / KDE")
  ratio.GetYaxis().SetNdivisions(505)
  ratio.GetXaxis().SetTitle(data.GetXaxis().GetTitle())
  return ratio


def _make_chi2_contrib_hist(data: ROOT.TH1, model: ROOT.TH1, name: str) -> ROOT.TH1:
  contrib = data.Clone(name)
  contrib.SetDirectory(0)
  contrib.Reset()
  contrib.SetTitle("")
  contrib.SetStats(0)
  contrib.SetMarkerSize(0.0)
  contrib.SetLineColor(ROOT.kRed + 1)
  contrib.SetFillColor(ROOT.kRed + 1)
  contrib.SetFillStyle(1001)
  contrib.GetYaxis().SetTitle("#chi^{2} contrib")
  contrib.GetXaxis().SetTitle(data.GetXaxis().GetTitle())
  for i in range(1, data.GetNbinsX() + 1):
    err = data.GetBinError(i)
    if err <= 0:
      continue
    diff = data.GetBinContent(i) - model.GetBinContent(i)
    contrib.SetBinContent(i, (diff * diff) / (err * err))
  return contrib


def _ratio_y_range(ratio: ROOT.TH1, pad_frac: float = 0.08) -> None:
  ymax = 1.0
  ymin = 1.0
  for i in range(1, ratio.GetNbinsX() + 1):
    val = ratio.GetBinContent(i)
    err = ratio.GetBinError(i)
    if err <= 0 and val == 0:
      continue
    ymax = max(ymax, val + err)
    ymin = min(ymin, val - err)
  span = max(ymax - ymin, 0.2)
  margin = pad_frac * span
  ratio.GetYaxis().SetRangeUser(ymin - margin, ymax + margin)


def _chi2_y_range(chi2_hist: ROOT.TH1, pad_frac: float = 0.12) -> None:
  ymax = chi2_hist.GetMaximum()
  if ymax <= 0:
    ymax = 1.0
  chi2_hist.GetYaxis().SetRangeUser(0.0, ymax * (1.0 + pad_frac))


def _draw_unity_line(hist: ROOT.TH1) -> ROOT.TLine:
  line = ROOT.TLine(
    hist.GetXaxis().GetXmin(),
    1.0,
    hist.GetXaxis().GetXmax(),
    1.0,
  )
  line.SetLineStyle(2)
  line.SetLineColor(ROOT.kBlack)
  line.Draw()
  hist._unity_line = line
  return line


def _info_pad_box(x1: float, y1: float, x2: float, y2: float) -> ROOT.TPave:
  box = ROOT.TPave(x1, y1, x2, y2, 1, "NDC")
  box.SetFillColor(ROOT.kWhite)
  box.SetFillStyle(1001)
  box.SetBorderSize(1)
  box.SetLineColor(ROOT.kBlack)
  return box


def _draw_projection_stats_box(pad: ROOT.TPad, stats: Th1ProjectionStats) -> None:
  pad.cd()
  box = _info_pad_box(0.28, 0.16, 0.72, 0.55)
  box.SetName("proj_stats_box")
  box.Draw()
  pad._stats_box = box

  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextFont(42)
  latex.SetTextSize(0.034)

  x_label = 0.31
  x_data = 0.46
  x_kde = 0.60
  y_header = 0.50
  y0 = 0.42
  dy = 0.06

  latex.SetTextAlign(23)
  latex.DrawLatex(x_data, y_header, "data")
  latex.DrawLatex(x_kde, y_header, "KDE")

  rows = (
    ("integral", f"{stats.integral_data:.4g}", f"{stats.integral_kde:.4g}"),
    (stats.mean_label, f"{stats.mean_data:.5g}", f"{stats.mean_kde:.5g}"),
  )
  for i, (label, data_val, kde_val) in enumerate(rows):
    y = y0 - i * dy
    latex.SetTextAlign(13)
    latex.DrawLatex(x_label, y, label)
    latex.SetTextAlign(23)
    latex.DrawLatex(x_data, y, data_val)
    latex.DrawLatex(x_kde, y, kde_val)

  y_chi = y0 - 2.3 * dy
  latex.SetTextAlign(12)
  latex.DrawLatex(
    x_label,
    y_chi,
    f"#chi^{{2}} = {stats.chi2:.4g},  ndf = {stats.ndf}",
  )
  latex.DrawLatex(
    x_label,
    y_chi - dy,
    f"#chi^{{2}}/ndf = {stats.reduced_chi2:.4g}",
  )
  pad._stats_latex = latex


def _style_projection_curve(curve: ROOT.TGraph) -> None:
  curve.SetLineColor(ROOT.kBlue + 1)
  curve.SetLineWidth(1)


def _stats_delta(hist_stats: Th2RegionStats, kde_stats: Th2RegionStats) -> Th2RegionStats:
  return Th2RegionStats(
    integral=kde_stats.integral - hist_stats.integral,
    mean_x=kde_stats.mean_x - hist_stats.mean_x,
    mean_y=kde_stats.mean_y - hist_stats.mean_y,
    std_x=kde_stats.std_x - hist_stats.std_x,
    std_y=kde_stats.std_y - hist_stats.std_y,
  )


def _param_lines(fit: dict, meta: dict[str, float | str | int]) -> list[str]:
  lines = [
    f"#rho = {float(meta.get('rho', float('nan'))):.5g},  #alpha = {fit['alpha']:.5g}",
    (
      f"#chi^{{2}} = {fit['chi2']:.4g},  #chi^{{2}}/ndf = {fit['reduced_chi2']:.4g},  "
      f"ndf = {fit['ndf']:.0f}"
    ),
  ]
  if meta.get("linear_combo", 0):
    mix = float(meta.get("mix", float("nan")))
    lines.append(f"mix = {mix:.4g} (unmirrored),  {1.0 - mix:.4g} (mirrored)")
  return lines


def _draw_overlay_stats_table(
  latex: ROOT.TLatex,
  hist_stats: Th2RegionStats,
  kde_stats: Th2RegionStats,
) -> None:
  delta = _stats_delta(hist_stats, kde_stats)
  x_label, x_int, x_mx, x_my = 0.10, 0.24, 0.34, 0.44
  y_header, y_hist = 0.82, 0.58
  latex.SetTextSize(0.15)
  latex.SetTextAlign(23)
  latex.DrawLatex(x_int, y_header, "integral")
  latex.DrawLatex(x_mx, y_header, "mean x")
  latex.DrawLatex(x_my, y_header, "mean y")
  latex.SetTextAlign(13)
  latex.DrawLatex(x_label, y_hist, "histogram")
  latex.DrawLatex(x_label, y_hist - 0.2, "KDE")
  latex.DrawLatex(x_label, y_hist - 0.4, "difference")
  latex.SetTextAlign(23)
  for row, stats in enumerate((hist_stats, kde_stats, delta)):
    y = y_hist - 0.2 * row
    latex.DrawLatex(x_int, y, f"{stats.integral:.4g}")
    latex.DrawLatex(x_mx, y, f"{stats.mean_x:.5g}")
    latex.DrawLatex(x_my, y, f"{stats.mean_y:.5g}")


def _draw_overlay_params(latex: ROOT.TLatex, param_lines: list[str]) -> None:
  y_param = 0.54
  latex.SetTextAlign(12)
  latex.SetTextSize(0.15)
  latex.DrawLatex(0.68, 0.82, "KDE Parameters")
  for line in param_lines:
    latex.DrawLatex(0.58, y_param, line)
    y_param -= 0.20


def _draw_stats_and_params(
  pad: ROOT.TPad,
  fit: dict,
  meta: dict[str, float | str | int],
) -> None:
  pad.cd()
  pad.SetFillStyle(0)

  stats_box = _info_pad_box(0.07, 0.04, 0.52, 0.94)
  params_box = _info_pad_box(0.54, 0.04, 0.96, 0.94)
  stats_box.Draw()
  params_box.Draw()
  pad._stats_table_box = stats_box
  pad._params_box = params_box

  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextFont(42)
  _draw_overlay_stats_table(latex, fit["hist_stats"], fit["kde_stats"])
  _draw_overlay_params(latex, _param_lines(fit, meta))


def plot_overlay(
  data: ROOT.TH2,
  fit: dict,
  meta: dict[str, float | str | int],
  outfile: str,
) -> None:
  model: KdeModel | StoredShapeModel = fit["model"]
  data_plot = data.Clone("overlay_data")
  data_plot.SetDirectory(0)
  data_plot.SetStats(0)
  if isinstance(model, StoredShapeModel):
    print(
      f"Evaluating fine KDE surface ({KDE_OVERLAY_POINTS}x{KDE_OVERLAY_POINTS}) "
      "via stored kde_template..."
    )
  else:
    print(
      f"Evaluating fine KDE surface ({KDE_OVERLAY_POINTS}x{KDE_OVERLAY_POINTS}) "
      "via RooNDKeysPdf..."
    )
  kde_surf = kde_fine_plot_hist(model, data, "overlay_kde_surf")
  kde_only = kde_surf.Clone("overlay_kde_only")
  kde_only.SetDirectory(0)

  canvas = ROOT.TCanvas("c_overlay", "2D KDE overlay", 3000, 1640)
  pad_info = ROOT.TPad("pad_info", "", 0.0, 0.0, 1.0, 0.22)
  pad_info.SetFillStyle(0)
  pad_info.Draw()
  pad_left = ROOT.TPad("pad_left", "", 0.0, 0.22, 0.5, 1.0)
  pad_left.Draw()
  pad_right = ROOT.TPad("pad_right", "", 0.5, 0.22, 1.0, 1.0)
  pad_right.Draw()

  pad_left.cd()
  _configure_surf_pad(pad_left)
  data_plot.SetTitle(PLOT_TITLE_OVERLAY)
  axis_titles = _style_surf_hist(data_plot, line_color=ROOT.kBlue + 1)
  kde_axes = _style_surf_hist(kde_surf, line_color=ROOT.kRed + 1, line_width=1)
  data_plot.Draw("LEGO")
  kde_surf.SetLineColorAlpha(ROOT.kRed + 1, 0.2)
  kde_surf.Draw("SURF SAME")
  _draw_surf3d_axis_titles(*axis_titles)

  leg = ROOT.TLegend(0.12, 0.82, 0.42, 0.92)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.04)
  leg.AddEntry(data_plot, "Data", "l")
  leg.AddEntry(kde_surf, "#alpha#timesKDE(x,y)", "l")
  leg.Draw()

  pad_right.cd()
  _configure_surf_pad(pad_right)
  kde_only.SetTitle(PLOT_TITLE_KDE_TEMPLATE)
  _style_surf_hist(kde_only, line_color=ROOT.kBlue + 1)
  kde_only.Draw("SURF")
  _draw_surf3d_axis_titles(*kde_axes)

  _draw_stats_and_params(pad_info, fit, meta)
  canvas._keepalive = [data_plot, kde_surf, kde_only, leg]
  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def _make_projection_pads(name_prefix: str) -> tuple[ROOT.TPad, ROOT.TPad, ROOT.TPad]:
  pad_chi2 = ROOT.TPad(f"{name_prefix}_chi2_pad", "", 0.0, 0.0, 1.0, 0.17)
  pad_chi2.SetTopMargin(0.04)
  pad_chi2.SetBottomMargin(0.32)
  pad_chi2.SetLeftMargin(0.12)
  pad_chi2.SetRightMargin(0.04)
  pad_chi2.SetGridy()
  pad_ratio = ROOT.TPad(f"{name_prefix}_ratio_pad", "", 0.0, 0.17, 1.0, 0.34)
  pad_ratio.SetTopMargin(0.02)
  pad_ratio.SetBottomMargin(0.32)
  pad_ratio.SetLeftMargin(0.12)
  pad_ratio.SetRightMargin(0.04)
  pad_ratio.SetGridy()
  pad_main = ROOT.TPad(f"{name_prefix}_main", "", 0.0, 0.34, 1.0, 1.0)
  pad_main.SetBottomMargin(0.12)
  pad_main.SetTopMargin(0.06)
  pad_main.SetLeftMargin(0.12)
  pad_main.SetRightMargin(0.04)
  pad_main.SetGridy()
  return pad_main, pad_ratio, pad_chi2


def _draw_projection_main(
  pad: ROOT.TPad,
  data: ROOT.TH1,
  model_curve: ROOT.TGraph,
  title: str,
) -> ROOT.TLegend:
  pad.cd()
  data.SetStats(0)
  data.SetTitle(title)
  data.SetMarkerSize(0.8)
  data.SetMarkerColor(ROOT.kBlack)
  data.SetLineColor(ROOT.kBlack)
  data.SetLineWidth(1)
  _style_projection_curve(model_curve)
  data.GetXaxis().SetLabelSize(0.045)
  data.GetXaxis().SetTitleSize(0.045)
  data.GetXaxis().SetTitleOffset(1.0)
  data.GetXaxis().SetNdivisions(505)
  ymax = max(data.GetMaximum(), _graph_max_y(model_curve))
  ymin = min(0.0, data.GetMinimum())
  span = max(ymax - ymin, 1.0)
  data.GetYaxis().SetRangeUser(ymin - 0.05 * span, ymax + 0.15 * span)
  data.GetYaxis().SetLabelSize(0.045)
  data.GetYaxis().SetTitleSize(0.045)
  data.GetYaxis().SetTitleOffset(1.2)
  data.Draw("E1 HIST")
  model_curve.Draw("L SAME")
  leg = ROOT.TLegend(0.82, 0.72, 0.98, 0.88)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.035)
  leg.AddEntry(data, "Data", "lep")
  leg.AddEntry(model_curve, "#alpha#timesKDE", "l")
  leg.Draw()
  return leg


def _draw_projection_column(
  parent_pad: ROOT.TPad,
  data: ROOT.TH1,
  model_curve: ROOT.TGraph,
  model_bins: ROOT.TH1,
  title: str,
  stats: Th1ProjectionStats,
  name_prefix: str,
) -> None:
  ratio = _make_ratio_hist(data, model_bins, f"{name_prefix}_ratio")
  chi2_hist = _make_chi2_contrib_hist(data, model_bins, f"{name_prefix}_chi2")
  parent_pad.cd()
  parent_pad.SetFillStyle(0)
  pad_main, pad_ratio, pad_chi2 = _make_projection_pads(name_prefix)
  pad_chi2.Draw()
  pad_ratio.Draw()
  pad_main.Draw()

  pad_main._legend = _draw_projection_main(pad_main, data, model_curve, title)
  pad_ratio.cd()
  ratio.GetXaxis().SetLabelSize(0.0)
  ratio.GetXaxis().SetTitleSize(0.0)
  ratio.GetYaxis().SetLabelSize(0.11)
  ratio.GetYaxis().SetTitleSize(0.11)
  ratio.GetYaxis().SetTitleOffset(0.42)
  _ratio_y_range(ratio)
  ratio.Draw("E1 HIST")
  pad_ratio._unity_line = _draw_unity_line(ratio)

  pad_chi2.cd()
  chi2_hist.GetXaxis().SetLabelSize(0.0)
  chi2_hist.GetXaxis().SetTitleSize(0.0)
  chi2_hist.GetYaxis().SetLabelSize(0.11)
  chi2_hist.GetYaxis().SetTitleSize(0.11)
  chi2_hist.GetYaxis().SetTitleOffset(0.48)
  _chi2_y_range(chi2_hist)
  chi2_hist.Draw("HIST")
  _draw_projection_stats_box(pad_main, stats)
  parent_pad._projection_pads = (pad_main, pad_ratio, pad_chi2)
  parent_pad._projection_ratio = ratio
  parent_pad._projection_chi2 = chi2_hist


def plot_projections(data: ROOT.TH2, fit: dict, outfile: str) -> None:
  hx_data = fit["proj_x_data"]
  hy_data = fit["proj_y_data"]
  hx_kde = fit["proj_x_model"]
  hy_kde = fit["proj_y_model"]
  gx = fit["curve_x"]
  gy = fit["curve_y"]
  # Display: rescale WLS-α model to integral-matched α (height matches data).
  scale = _integral_match_scale(
    fit["hist_stats"].integral, fit["kde_stats"].integral
  )
  hx_kde.Scale(scale)
  hy_kde.Scale(scale)
  _scale_graph_y(gx, scale)
  _scale_graph_y(gy, scale)
  print(
    f"Projection display: integral-match scale={scale:.6g} "
    f"(α_int / α_WLS; overlay/2D χ² still use WLS α)"
  )
  xtitle, ytitle, ztitle = _axis_titles(data)
  hx_data.GetXaxis().SetTitle(xtitle)
  hy_data.GetXaxis().SetTitle(ytitle)
  hx_data.GetYaxis().SetTitle(ztitle)
  hy_data.GetYaxis().SetTitle(ztitle)

  x_stats = _projection_stats(hx_data, hx_kde, "mean x")
  y_stats = _projection_stats(hy_data, hy_kde, "mean y")

  canvas = ROOT.TCanvas("c_proj", "2D KDE projections", 2800, 1240)
  canvas.Divide(2, 1)
  pad_x = canvas.GetPad(1)
  pad_y = canvas.GetPad(2)
  _draw_projection_column(
    pad_x, hx_data, gx, hx_kde, PLOT_TITLE_X_PROJ, x_stats, f"{data.GetName()}_px"
  )
  _draw_projection_column(
    pad_y, hy_data, gy, hy_kde, PLOT_TITLE_Y_PROJ, y_stats, f"{data.GetName()}_py"
  )
  canvas._projection_data = (hx_data, hy_data, gx, gy, hx_kde, hy_kde)
  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def plot_chi2_contrib_lego(
  data: ROOT.TH2,
  fit: dict,
  region_label: str,
  outfile: str,
) -> None:
  contrib = make_chi2_contrib_th2(data, fit["template"], "chi2_contrib_compare")
  zmax = contrib.GetMaximum()
  if zmax <= 0:
    zmax = 1.0
  contrib.GetZaxis().SetRangeUser(0.0, zmax * 1.05)
  contrib.GetZaxis().SetTitle("#chi^{2}_{i}")
  contrib.GetXaxis().SetTitle(data.GetXaxis().GetTitle() or "x [cm]")
  contrib.GetXaxis().SetTitleOffset(2)
  contrib.GetYaxis().SetTitle(data.GetYaxis().GetTitle() or "y [cm]")

  canvas = ROOT.TCanvas("c_chi2_contrib", "#chi^{2} contributions", 1200, 900)
  canvas.SetRightMargin(0.14)
  canvas.SetLeftMargin(0.12)
  canvas.SetBottomMargin(0.12)
  canvas.SetTopMargin(0.10)
  canvas.SetTheta(28)
  canvas.SetPhi(60)
  contrib.Draw("LEGO")

  latex = ROOT.TLatex()
  latex.SetNDC(True)
  latex.SetTextFont(42)
  latex.SetTextSize(0.035)
  latex.SetTextAlign(12)
  latex.DrawLatex(0.08, 0.92, region_label)
  latex.DrawLatex(
    0.08,
    0.87,
    f"#chi^{{2}}={fit['chi2']:.2f}, #chi^{{2}}/ndf={fit['reduced_chi2']:.3f}",
  )
  canvas._keepalive = [contrib, latex]
  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def region_label(hist: ROOT.TH2) -> str:
  return f"{hist.GetNbinsX()}#times{hist.GetNbinsY()} bins"


def output_paths(hist: ROOT.TH2) -> tuple[str, str, str]:
  tag = (
    f"{OUTPUT_TAG}_{HIST_NAME}_"
    f"{hist.GetNbinsX()}x{hist.GetNbinsY()}"
  )
  return (
    os.path.join(OUTPUT_DIR, f"{tag}_overlay.pdf"),
    os.path.join(OUTPUT_DIR, f"{tag}_projections.pdf"),
    os.path.join(OUTPUT_DIR, f"{tag}_chi2_contrib.pdf"),
  )


def prepare_histogram(hist: ROOT.TH2) -> ROOT.TH2:
  print(
    f"Loaded {HIST_NAME!r}: {hist.GetNbinsX()}x{hist.GetNbinsY()} bins "
    f"from {HIST_ROOT_FILE}"
  )
  if hist.GetSumw2N() == 0:
    hist.Sumw2()
  hist = rebin_histogram2d(hist, REBIN_FACTOR_X, REBIN_FACTOR_Y)
  if REBIN_FACTOR_X != 1 or REBIN_FACTOR_Y != 1:
    print(
      f"Rebinned in memory by {REBIN_FACTOR_X}x{REBIN_FACTOR_Y} -> "
      f"{hist.GetNbinsX()}x{hist.GetNbinsY()} bins"
    )
  hist = restrict_histogram2d(hist, RESTRICT_NBINS_X, RESTRICT_NBINS_Y)
  return hist


def build_compare_model(
  train_target: ROOT.TH2,
  meta: dict[str, float | str | int],
  kde_template: ROOT.TH2,
) -> KdeModel | StoredShapeModel:
  """Live PDF for plain fits; stored template when shift/combine ops are present."""
  if fit_has_spatial_ops(meta):
    # shift_kde bake-in: per-leg alphas live in kde_template; overall scale is
    # re-profiled against the comparison hist. Do not rebuild from target_hist.
    print(
      "fit_meta has shift/combine keys: using stored kde_template "
      f"({kde_template.GetNbinsX()}x{kde_template.GetNbinsY()})"
    )
    return StoredShapeModel(shape=kde_template, alpha=1.0)
  print(
    f"Rebuilding RooNDKeysPdf from training hist "
    f"{train_target.GetNbinsX()}x{train_target.GetNbinsY()} "
    f"(rho={meta.get('rho')}, linear_combo={bool(meta.get('linear_combo', 0))})"
  )
  return build_kde_model(train_target, meta)


def main() -> int:
  data = prepare_histogram(open_histogram2d(HIST_ROOT_FILE, HIST_NAME))
  print(f"Loading KDE fit from {KDE_ROOT_FILE}")
  train_target, meta, _kde_shape, kde_template = load_kde_fit(KDE_ROOT_FILE)
  model = build_compare_model(train_target, meta, kde_template)
  fit = evaluate_comparison("compare", data, model)

  os.makedirs(OUTPUT_DIR, exist_ok=True)
  overlay_path, proj_path, chi2_path = output_paths(data)
  label = region_label(data)
  plot_overlay(data, fit, meta, overlay_path)
  plot_projections(data, fit, proj_path)
  plot_chi2_contrib_lego(data, fit, label, chi2_path)
  return 0


if __name__ == "__main__":
  sys.exit(main())
