#!/usr/bin/env python3
"""
Shift and/or linearly combine 2D KDEs from 2d_kde fit ROOT files.

Modes (MODE):
  "shift"   — translate training bin centers by (SHIFT_X, SHIFT_Y), rebuild
              RooNDKeysPdf with the original fit_meta params, then evaluate the
              translated density back onto the original target_hist axes.
  "combine" — do the same for KDE_A and KDE_B independently, then form
              c1*A' + c2*B' on A's original grid.

No TH2 Interpolate of an existing KDE. Training uses an axis-shifted copy of
the target (same bin heights); the saved kde_shape / kde_template use the
original lab binning. Mass that leaves the original window leaves a vacated
strip (near-empty bins); Δmean on that fixed window need not equal the shift.
"""

import os
import sys
from dataclasses import dataclass

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


# ===============================SCRIPT PARAMS===============================
# "shift" | "combine"
MODE = "shift"

INPUT_FILE_NAME = "pX_100um_mm1_2d_kde"  # .root ADDED BELOW
INPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", f"{INPUT_FILE_NAME}.root"
)

# Enter as mm. *0.1 converts to cm (nominal axis units).
SHIFT_X = 1.54976 * 0.1
SHIFT_Y = 0.0 * 0.1

# Combine mode: load two fit products, shift+rebuild each, then c1*A' + c2*B'.
# Primary (A) supplies scan graphs, base fit_meta, and the original eval grid.
COMBINE_FILE_A_NAME = "ml_nominal_mm1_2d_kde"  # .root added below
COMBINE_FILE_B_NAME = "pX_100um_mm1_2d_kde"
COMBINE_ROOT_A = os.path.join(
  os.path.dirname(__file__), "..", "root_files", f"{COMBINE_FILE_A_NAME}.root"
)
COMBINE_ROOT_B = os.path.join(
  os.path.dirname(__file__), "..", "root_files", f"{COMBINE_FILE_B_NAME}.root"
)

SHIFT_A_X = -0.0912
SHIFT_A_Y = 0.0
SHIFT_B_X = 0.06108
SHIFT_B_Y = 0.0
COEFF_A = 0.5
COEFF_B = 0.5

OUTPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__),
  "..",
  "root_files",
  f"{INPUT_FILE_NAME}_shifted_({SHIFT_X:g},{SHIFT_Y:g}).root",
)
COMBINE_OUTPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__),
  "..",
  "root_files",
  f"{COMBINE_FILE_A_NAME}_plus_{COMBINE_FILE_B_NAME}"
  f"_sA({SHIFT_A_X:g},{SHIFT_A_Y:g})"
  f"_sB({SHIFT_B_X:g},{SHIFT_B_Y:g})"
  f"_c({COEFF_A:g},{COEFF_B:g}).root",
)

NDKEYS_NO_MIRROR = "a"
NDKEYS_MIRROR_BOTH = "am"
# ===========================================================================

REQUIRED_KEYS = (
  "target_hist",
  "chi2_vs_rho",
  "reduced_chi2_vs_rho",
  "alpha_vs_rho",
  "kde_shape",
  "kde_template",
  "fit_meta",
  "stats_meta",
)


@dataclass
class Th2Stats:
  integral: float
  mean_x: float
  mean_y: float


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


def parse_fit_meta(meta: ROOT.TNamed) -> dict[str, float | str | int]:
  out: dict[str, float | str | int] = {}
  for part in meta.GetTitle().split(";"):
    if "=" not in part:
      continue
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


def shift_th2_axes(
  hist: ROOT.TH2,
  shift_x: float,
  shift_y: float,
  out_name: str,
) -> ROOT.TH2D:
  """Same bin heights; axis limits (and thus bin centers) translated by (dx, dy)."""
  xaxis = hist.GetXaxis()
  yaxis = hist.GetYaxis()
  out = ROOT.TH2D(
    out_name,
    f"{hist.GetTitle()} shifted ({shift_x:g},{shift_y:g})",
    hist.GetNbinsX(),
    xaxis.GetXmin() + shift_x,
    xaxis.GetXmax() + shift_x,
    hist.GetNbinsY(),
    yaxis.GetXmin() + shift_y,
    yaxis.GetXmax() + shift_y,
  )
  out.SetDirectory(0)
  for ix in range(1, hist.GetNbinsX() + 1):
    for iy in range(1, hist.GetNbinsY() + 1):
      out.SetBinContent(ix, iy, hist.GetBinContent(ix, iy))
      out.SetBinError(ix, iy, hist.GetBinError(ix, iy))
  out.ResetStats()
  return out


def histogram_to_weighted_dataset(hist: ROOT.TH2) -> KdeEvalContext:
  xlo = hist.GetXaxis().GetXmin()
  xhi = hist.GetXaxis().GetXmax()
  ylo = hist.GetYaxis().GetXmin()
  yhi = hist.GetYaxis().GetXmax()
  x_var = ROOT.RooRealVar("x", "x [cm]", xlo, xhi)
  y_var = ROOT.RooRealVar("y", "y [cm]", ylo, yhi)
  argset = ROOT.RooArgSet(x_var, y_var)
  w_var = ROOT.RooRealVar("w", "weight", 0.0, 1.0e20)
  dataset = ROOT.RooDataSet(
    "weighted_points",
    "weighted_points",
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
  if dataset.numEntries() == 0:
    raise ValueError("histogram has no positive bin content; cannot build KDE")
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


def build_kde_model(target: ROOT.TH2, meta: dict[str, float | str | int]) -> KdeModel:
  ctx = histogram_to_weighted_dataset(target)
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
        "shift_kde_unmirrored", ctx, mirror_options=opt_no, rho=rho
      ),
      pdf_mirrored=make_ndkeys_pdf(
        "shift_kde_mirrored", ctx, mirror_options=opt_m, rho=rho
      ),
    )
  return KdeModel(
    ctx=ctx,
    alpha=alpha,
    mix=1.0,
    use_linear_combo=False,
    pdf_single=make_ndkeys_pdf(
      "shift_kde_single", ctx, mirror_options=opt_m, rho=rho
    ),
  )


def expand_model_range(model: KdeModel, ref_hist: ROOT.TH2) -> None:
  """Allow PDF evaluation over ref_hist axes (union with training range)."""
  xlo = min(model.ctx.x_var.getMin(), ref_hist.GetXaxis().GetXmin())
  xhi = max(model.ctx.x_var.getMax(), ref_hist.GetXaxis().GetXmax())
  ylo = min(model.ctx.y_var.getMin(), ref_hist.GetYaxis().GetXmin())
  yhi = max(model.ctx.y_var.getMax(), ref_hist.GetYaxis().GetXmax())
  model.ctx.x_var.setRange(xlo, xhi)
  model.ctx.y_var.setRange(ylo, yhi)


def evaluate_shape_on_hist(
  model: KdeModel,
  ref_hist: ROOT.TH2,
  name: str,
  title: str,
) -> ROOT.TH2D:
  """Evaluate unnormalized KDE shape at ref_hist bin centers."""
  out = ROOT.TH2D(
    name,
    title,
    ref_hist.GetNbinsX(),
    ref_hist.GetXaxis().GetXmin(),
    ref_hist.GetXaxis().GetXmax(),
    ref_hist.GetNbinsY(),
    ref_hist.GetYaxis().GetXmin(),
    ref_hist.GetYaxis().GetXmax(),
  )
  out.SetDirectory(0)
  out._hold_model = model
  for ix in range(1, out.GetNbinsX() + 1):
    x = out.GetXaxis().GetBinCenter(ix)
    for iy in range(1, out.GetNbinsY() + 1):
      y = out.GetYaxis().GetBinCenter(iy)
      out.SetBinContent(ix, iy, model.shape_at(x, y))
  return out


def scaled_template(shape: ROOT.TH2, alpha: float, name: str) -> ROOT.TH2:
  template = shape.Clone(name)
  template.SetDirectory(0)
  template.Scale(alpha)
  template.SetTitle("#alpha#timesKDE(x,y);x [cm];y [cm]")
  return template


def rebuild_shifted_kde(
  target: ROOT.TH2,
  meta: dict[str, float | str | int],
  shift_x: float,
  shift_y: float,
) -> tuple[ROOT.TH2D, ROOT.TH2, ROOT.TH2, KdeModel]:
  """Train on axis-shifted target; evaluate shape/template on original axes.

  Returns (shifted_train, shape, template, model). Saved products should use
  the original target as target_hist; shape/template match its binning.
  """
  shifted_train = shift_th2_axes(target, shift_x, shift_y, "shifted_train")
  model = build_kde_model(shifted_train, meta)
  expand_model_range(model, target)
  title = (
    "mirrored/unmirrored KDE shape" if model.use_linear_combo else "KDE shape"
  )
  shape = evaluate_shape_on_hist(model, target, "kde_shape", title)
  template = scaled_template(shape, model.alpha, "kde_template")
  # Keep training hist alive for the model dataset.
  model._hold_train = shifted_train
  return shifted_train, shape, template, model


def add_th2s(
  hist_a: ROOT.TH2,
  hist_b: ROOT.TH2,
  coeff_a: float,
  coeff_b: float,
  out_name: str,
) -> ROOT.TH2:
  """Return coeff_a * hist_a + coeff_b * hist_b; axes must match hist_a."""
  if hist_a.GetNbinsX() != hist_b.GetNbinsX() or hist_a.GetNbinsY() != hist_b.GetNbinsY():
    raise ValueError(
      f"bin count mismatch: {hist_a.GetName()} "
      f"({hist_a.GetNbinsX()}x{hist_a.GetNbinsY()}) vs {hist_b.GetName()} "
      f"({hist_b.GetNbinsX()}x{hist_b.GetNbinsY()})"
    )
  out = hist_a.Clone(out_name)
  out.SetDirectory(0)
  out.SetTitle(f"{coeff_a:g}*{hist_a.GetName()} + {coeff_b:g}*{hist_b.GetName()}")
  for ix in range(1, out.GetNbinsX() + 1):
    for iy in range(1, out.GetNbinsY() + 1):
      out.SetBinContent(
        ix,
        iy,
        coeff_a * hist_a.GetBinContent(ix, iy)
        + coeff_b * hist_b.GetBinContent(ix, iy),
      )
  out.ResetStats()
  return out


def evaluate_model_on_hist(
  model: KdeModel,
  ref_hist: ROOT.TH2,
  name: str,
  title: str,
  *,
  scaled: bool,
) -> ROOT.TH2D:
  shape = evaluate_shape_on_hist(model, ref_hist, f"{name}_shape", title)
  if scaled:
    return scaled_template(shape, model.alpha, name)
  shape.SetName(name)
  return shape


def th2_distribution_stats(hist: ROOT.TH2) -> Th2Stats:
  integral = 0.0
  mean_x_num = 0.0
  mean_y_num = 0.0
  for ix in range(1, hist.GetNbinsX() + 1):
    x = hist.GetXaxis().GetBinCenter(ix)
    for iy in range(1, hist.GetNbinsY() + 1):
      weight = hist.GetBinContent(ix, iy)
      if weight <= 0:
        continue
      y = hist.GetYaxis().GetBinCenter(iy)
      integral += weight
      mean_x_num += weight * x
      mean_y_num += weight * y
  if integral <= 0:
    return Th2Stats(0.0, 0.0, 0.0)
  return Th2Stats(
    integral=integral,
    mean_x=mean_x_num / integral,
    mean_y=mean_y_num / integral,
  )


def stats_meta_string(hist_stats: Th2Stats, kde_stats: Th2Stats) -> str:
  return (
    f"hist_integral={hist_stats.integral};"
    f"hist_mean_x={hist_stats.mean_x};hist_mean_y={hist_stats.mean_y};"
    f"kde_integral={kde_stats.integral};"
    f"kde_mean_x={kde_stats.mean_x};kde_mean_y={kde_stats.mean_y}"
  )


def load_fit_objects(filepath: str) -> dict[str, ROOT.TObject]:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open ROOT file: {filepath}")
  objects: dict[str, ROOT.TObject] = {}
  for name in REQUIRED_KEYS:
    obj = tfile.Get(name)
    if not obj:
      tfile.Close()
      raise KeyError(f"missing {name!r} in {filepath}")
    clone = obj.Clone(name)
    if hasattr(clone, "SetDirectory"):
      clone.SetDirectory(0)
    objects[name] = clone
  tfile.Close()
  return objects


def augment_fit_meta_shift(
  meta: ROOT.TNamed,
  shift_x: float,
  shift_y: float,
) -> ROOT.TNamed:
  title = meta.GetTitle()
  for key in ("shift_x=", "shift_y="):
    if key in title:
      raise ValueError(f"fit_meta already contains {key!r}; refusing to double-shift")
  return ROOT.TNamed(
    "fit_meta",
    f"{title};shift_x={shift_x};shift_y={shift_y}",
  )


def augment_fit_meta_combine(
  meta: ROOT.TNamed,
  file_a: str,
  file_b: str,
  shift_a_x: float,
  shift_a_y: float,
  shift_b_x: float,
  shift_b_y: float,
  coeff_a: float,
  coeff_b: float,
) -> ROOT.TNamed:
  title = meta.GetTitle()
  for key in ("combine_coeff_a=", "combine_coeff_b="):
    if key in title:
      raise ValueError(f"fit_meta already contains {key!r}; refusing to re-combine")
  base_a = os.path.basename(file_a)
  base_b = os.path.basename(file_b)
  return ROOT.TNamed(
    "fit_meta",
    f"{title};combine_file_a={base_a};combine_file_b={base_b};"
    f"shift_a_x={shift_a_x};shift_a_y={shift_a_y};"
    f"shift_b_x={shift_b_x};shift_b_y={shift_b_y};"
    f"combine_coeff_a={coeff_a};combine_coeff_b={coeff_b}",
  )


def write_fit_product(
  outfile: str,
  objects: dict[str, ROOT.TObject],
  target_hist: ROOT.TH2,
  kde_shape: ROOT.TH2,
  kde_template: ROOT.TH2,
  fit_meta: ROOT.TNamed,
  stats_meta: ROOT.TNamed,
) -> None:
  out_dir = os.path.dirname(os.path.abspath(outfile))
  os.makedirs(out_dir, exist_ok=True)
  fout = ROOT.TFile.Open(outfile, "RECREATE")
  if not fout or fout.IsZombie():
    raise OSError(f"cannot create output file: {outfile}")

  fout.cd()
  target_out = target_hist.Clone("target_hist")
  target_out.SetDirectory(fout)
  target_out.Write()

  for name in ("chi2_vs_rho", "reduced_chi2_vs_rho", "alpha_vs_rho"):
    obj = objects[name].Clone(name)
    if hasattr(obj, "SetDirectory"):
      obj.SetDirectory(fout)
    obj.Write()

  shape = kde_shape.Clone("kde_shape")
  shape.SetDirectory(fout)
  shape.Write()

  template = kde_template.Clone("kde_template")
  template.SetDirectory(fout)
  template.Write()

  fit_meta.Write()
  stats_meta.Write()
  fout.Close()


def print_mean_shift_check(
  label: str,
  before: Th2Stats,
  after: Th2Stats,
  shift_x: float,
  shift_y: float,
) -> None:
  dx = after.mean_x - before.mean_x
  dy = after.mean_y - before.mean_y
  print(
    f"  {label}: mean before=({before.mean_x:.6g}, {before.mean_y:.6g}) "
    f"after=({after.mean_x:.6g}, {after.mean_y:.6g}) "
    f"Δ=({dx:.6g}, {dy:.6g}) shift=({shift_x:g}, {shift_y:g}) "
    f"(fixed-window Δ need not equal shift)"
  )


def run_shift() -> int:
  print(
    f"Mode: shift (rebuild on shifted train; eval on original axes)\n"
    f"Reading fit product from {INPUT_ROOT_FILE}\n"
    f"Shift (x, y) = ({SHIFT_X:g}, {SHIFT_Y:g}) cm\n"
    f"Writing {OUTPUT_ROOT_FILE}"
  )

  objects = load_fit_objects(INPUT_ROOT_FILE)
  target = objects["target_hist"]
  if not target.InheritsFrom("TH2"):
    raise TypeError("target_hist is not TH2")
  meta = parse_fit_meta(objects["fit_meta"])

  orig_hist_stats = th2_distribution_stats(target)
  orig_kde_stats = th2_distribution_stats(objects["kde_template"])

  _shifted_train, shape, template, model = rebuild_shifted_kde(
    target, meta, SHIFT_X, SHIFT_Y
  )

  hist_stats = orig_hist_stats
  kde_stats = th2_distribution_stats(template)
  print(
    f"  eval grid: x=[{target.GetXaxis().GetXmin():g},{target.GetXaxis().GetXmax():g}] "
    f"y=[{target.GetYaxis().GetXmin():g},{target.GetYaxis().GetXmax():g}] "
    f"(original axes)"
  )
  print_mean_shift_check(
    "kde_template", orig_kde_stats, kde_stats, SHIFT_X, SHIFT_Y
  )

  fit_meta = augment_fit_meta_shift(objects["fit_meta"], SHIFT_X, SHIFT_Y)
  stats_meta = ROOT.TNamed("stats_meta", stats_meta_string(hist_stats, kde_stats))
  write_fit_product(
    OUTPUT_ROOT_FILE,
    objects,
    target,
    shape,
    template,
    fit_meta,
    stats_meta,
  )

  print(
    f"Wrote {OUTPUT_ROOT_FILE}\n"
    f"  kde_template integral={kde_stats.integral:.6g} "
    f"mean=({kde_stats.mean_x:.6g}, {kde_stats.mean_y:.6g})"
  )
  _ = model
  return 0


def run_combine() -> int:
  print(
    f"Mode: combine (shift+rebuild each, eval on A's original grid, then sum)\n"
    f"A: {COMBINE_ROOT_A}\n"
    f"   shift=({SHIFT_A_X:g}, {SHIFT_A_Y:g}) cm  coeff={COEFF_A:g}\n"
    f"B: {COMBINE_ROOT_B}\n"
    f"   shift=({SHIFT_B_X:g}, {SHIFT_B_Y:g}) cm  coeff={COEFF_B:g}\n"
    f"Writing {COMBINE_OUTPUT_ROOT_FILE}"
  )

  objects_a = load_fit_objects(COMBINE_ROOT_A)
  objects_b = load_fit_objects(COMBINE_ROOT_B)
  target_a = objects_a["target_hist"]
  target_b = objects_b["target_hist"]
  if not target_a.InheritsFrom("TH2") or not target_b.InheritsFrom("TH2"):
    raise TypeError("target_hist is not TH2")

  meta_a = parse_fit_meta(objects_a["fit_meta"])
  meta_b = parse_fit_meta(objects_b["fit_meta"])

  _train_a, shape_a, template_a, model_a = rebuild_shifted_kde(
    target_a, meta_a, SHIFT_A_X, SHIFT_A_Y
  )
  _train_b, _shape_b, _template_b, model_b = rebuild_shifted_kde(
    target_b, meta_b, SHIFT_B_X, SHIFT_B_Y
  )

  # Both PDFs on A's original grid (PDF getVal; no TH2 Interpolate).
  expand_model_range(model_b, target_a)
  shape_b_on_a = evaluate_model_on_hist(
    model_b,
    target_a,
    "kde_shape_b_on_a",
    "KDE B on A grid",
    scaled=False,
  )
  template_b_on_a = evaluate_model_on_hist(
    model_b,
    target_a,
    "kde_template_b_on_a",
    "KDE B template on A grid",
    scaled=True,
  )

  shape = add_th2s(shape_a, shape_b_on_a, COEFF_A, COEFF_B, "kde_shape")
  template = add_th2s(template_a, template_b_on_a, COEFF_A, COEFF_B, "kde_template")
  template.SetTitle("#alpha#timesKDE(x,y);x [cm];y [cm]")

  hist_stats = th2_distribution_stats(target_a)
  kde_stats = th2_distribution_stats(template)
  fit_meta = augment_fit_meta_combine(
    objects_a["fit_meta"],
    COMBINE_ROOT_A,
    COMBINE_ROOT_B,
    SHIFT_A_X,
    SHIFT_A_Y,
    SHIFT_B_X,
    SHIFT_B_Y,
    COEFF_A,
    COEFF_B,
  )
  stats_meta = ROOT.TNamed("stats_meta", stats_meta_string(hist_stats, kde_stats))

  write_fit_product(
    COMBINE_OUTPUT_ROOT_FILE,
    objects_a,
    target_a,
    shape,
    template,
    fit_meta,
    stats_meta,
  )

  print(
    f"Wrote {COMBINE_OUTPUT_ROOT_FILE}\n"
    f"  kde_template integral={kde_stats.integral:.6g} "
    f"mean=({kde_stats.mean_x:.6g}, {kde_stats.mean_y:.6g})"
  )
  _ = (model_a, model_b)
  return 0


def main() -> int:
  if MODE == "shift":
    return run_shift()
  if MODE == "combine":
    return run_combine()
  raise ValueError(f"MODE must be 'shift' or 'combine', got {MODE!r}")


if __name__ == "__main__":
  sys.exit(main())
