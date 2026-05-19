import ROOT
import array
import random
ROOT.gErrorIgnoreLevel = ROOT.kWarning  # or kError to hide warnings too


nominal_file = ROOT.TFile.Open("/Users/johnculbertson/Documents/Personal-UTKL/root-test-files/nominal.root", "READ")
## three TH2D named nominalxyposMM1, 2, 3

mm1 = nominal_file.Get("nominalxyposMM1")
mm2 = nominal_file.Get("nominalxyposMM2")
mm3 = nominal_file.Get("nominalxyposMM3")

mm1x = mm1.ProjectionX()
mm1y = mm1.ProjectionY()
mm2x = mm2.ProjectionX()
mm2y = mm2.ProjectionY()
mm3x = mm3.ProjectionX()
mm3y = mm3.ProjectionY()


DEFAULT_KDE_KERNEL = "gaussian"
DEFAULT_KDE_ADAPTIVE = True

KDE_KERNELS = ("gaussian", "epanechnikov", "biweight", "cosinearch")


def build_kde_options(kernel, adaptive=True):
  k = kernel.lower()
  if k not in KDE_KERNELS:
    raise ValueError(f"unknown kernel {kernel!r}; choose from {sorted(KDE_KERNELS)}")
  iteration = "Adaptive" if adaptive else "Fixed"
  kernel_tag = {"gaussian": "Gaussian", "epanechnikov": "Epanechnikov",
                "biweight": "Biweight", "cosinearch": "CosineArch"}[k]
  return (
    f"KernelType:{kernel_tag};Iteration:{iteration};"
    "Mirror:noMirror;Binning:RelaxedBinning"
  )


def hist_bin_centers_weights(hist):
  xs, ws = [], []
  for i in range(1, hist.GetNbinsX() + 1):
    c = hist.GetBinContent(i)
    if c <= 0:
      continue
    xs.append(hist.GetBinCenter(i))
    ws.append(c)
  return xs, ws


# TKDE from histogram bin centers/contents. bandwidth is global tuning factor rho.
# adaptive=True: local bandwidth per point (ROOT Iteration:Adaptive).
def weighted_kde(hist, bandwidth, kernel=DEFAULT_KDE_KERNEL, adaptive=DEFAULT_KDE_ADAPTIVE):
  xs, ws = hist_bin_centers_weights(hist)
  if not xs:
    return None
  xlo = hist.GetXaxis().GetXmin()
  xhi = hist.GetXaxis().GetXmax()
  da = array.array("d", xs)
  wa = array.array("d", ws)
  options = build_kde_options(kernel, adaptive)
  kde = ROOT.TKDE(len(da), da, wa, xlo, xhi, options, bandwidth)
  kde._kde_adaptive = adaptive
  kde._kde_kernel = kernel
  return kde


def local_bandwidths(kde, hist):
  # Per-point bandwidths (adaptive) or constant (fixed). Returns (x[], h[]).
  xs, _ = hist_bin_centers_weights(hist)
  if not xs:
    return [], []
  kde.GetFunction()
  if getattr(kde, "_kde_adaptive", True):
    weights = kde.GetAdaptiveWeights()
    if weights is None:
      return xs, []
    hs = [weights[i] for i in range(len(xs))]
  else:
    h0 = kde.GetFixedWeight()
    hs = [h0] * len(xs)
  return xs, hs


def adaptive_bandwidth_graph(kde, hist):
  xs, hs = local_bandwidths(kde, hist)
  if not xs or not hs:
    return None
  g = ROOT.TGraph(len(xs))
  for i, (x, h) in enumerate(zip(xs, hs)):
    g.SetPoint(i, x, h)
  g.SetTitle("local bandwidth h(x);x;h(x)")
  g.SetLineColor(ROOT.kRed + 1)
  g.SetLineWidth(2)
  return g


# TF1::Integral with tight epsrel often fails on TKDE (GSL status 18 over [-100,100]).
def integrate_tf1(func, xlo, xhi, n=2000):
  if n < 2:
    n = 2
  dx = (xhi - xlo) / (n - 1)
  total = 0.5 * (func.Eval(xlo) + func.Eval(xhi))
  x = xlo
  for _ in range(1, n - 1):
    x += dx
    total += func.Eval(x)
  return total * dx


# TKDE::GetFunction() is not guaranteed to integrate to 1 on [xlo, xhi]. Histogram
# density from Scale(1/Integral("width")) does. Renormalize so ∫ kde dx = 1.
def renorm_kde_tf1(raw_func, xlo, xhi, npx=2000):
  if raw_func is None:
    return None
  raw_func.SetNpx(npx)
  den = integrate_tf1(raw_func, xlo, xhi, npx)
  if den <= 0:
    return raw_func

  def scaled(x, p):
    return p[0] * raw_func.Eval(x[0])

  name = "kde_pdf_%08x" % random.getrandbits(32)
  out = ROOT.TF1(name, scaled, xlo, xhi, 1)
  out.SetParameter(0, 1.0 / den)
  out.SetNpx(npx)
  out._hold_raw = raw_func
  return out


# does the MSE for a function and a 1D histogram at the CENTER of each bin
# (function, 1D histogram)
# returns the MSE
def MSE(func, hist):
  nb = hist.GetNbinsX()
  if nb == 0:
    return 0.0
  sse = 0.0
  for i in range(1, nb + 1):
    hist_value = hist.GetBinContent(i)
    func_value = func.Eval(hist.GetBinCenter(i))
    sse += (hist_value - func_value) ** 2
  return sse / nb


# copy past of MSE() but weighted by 1/the value of the histagram at each point
def chi_sq(func, hist):
  nb = hist.GetNbinsX()
  if nb == 0:
    return 0.0
  s = 0.0
  for i in range(1, nb + 1):
    hist_value = hist.GetBinContent(i)
    func_value = func.Eval(hist.GetBinCenter(i))
    err = hist.GetBinError(i)
    if err <= 0:
      continue
    s += ((hist_value - func_value) ** 2) / (err * err)
  return s


# does a quasi integral to compute ISE for a function and a 1D histogram
# (function, 1D histogram, divisions/bin)
# returns the ISE
def quasi_ISE(func, hist, divisions):
  ise = 0.0
  for i in range(1, hist.GetNbinsX() + 1):
    a = hist.GetBinLowEdge(i)
    b = a + hist.GetBinWidth(i)
    dens = hist.GetBinContent(i)
    dx = (b - a) / divisions
    for j in range(divisions):
      x = a + (j + 0.5) * dx
      d = func.Eval(x) - dens
      ise += d * d * dx
  return ise


# Build normalized pdf KDE for one bandwidth; returns TF1 or None.
def kde_pdf_for_bandwidth(hist, bandwidth, kernel=DEFAULT_KDE_KERNEL,
                          adaptive=DEFAULT_KDE_ADAPTIVE, npx=2000):
  kde = weighted_kde(hist, bandwidth, kernel=kernel, adaptive=adaptive)
  if kde is None:
    return None
  xlo = hist.GetXaxis().GetXmin()
  xhi = hist.GetXaxis().GetXmax()
  pdf = renorm_kde_tf1(kde.GetFunction(), xlo, xhi, npx)
  if pdf is not None:
    pdf._hold_kde = kde
    pdf._kde_adaptive = adaptive
    pdf._kde_kernel = kernel
  return pdf


# Scan bandwidths; returns list of {bandwidth, mse, ise, chi2} (failed fits omitted).
def scan_kde_bandwidths(hist, bw_low, bw_high, step, kernel=DEFAULT_KDE_KERNEL,
                        adaptive=DEFAULT_KDE_ADAPTIVE, ise_divisions=100, npx=2000):
  results = []
  bw = bw_low
  while bw <= bw_high + 1e-12 * max(abs(bw_high), 1.0):
    kde_fn = kde_pdf_for_bandwidth(hist, bw, kernel=kernel, adaptive=adaptive, npx=npx)
    if kde_fn is not None:
      results.append({
        "bandwidth": bw,
        "mse": MSE(kde_fn, hist),
        "ise": quasi_ISE(kde_fn, hist, ise_divisions),
        "chi2": chi_sq(kde_fn, hist),
      })
    bw += step
  return results


# Best bandwidth from scan_kde_bandwidths; metric is "mse", "ise", or "chi2".
def find_optimal_bandwidth(hist, bw_low, bw_high, step, kernel=DEFAULT_KDE_KERNEL,
                           adaptive=DEFAULT_KDE_ADAPTIVE, metric="mse",
                           ise_divisions=100, npx=2000):
  rows = scan_kde_bandwidths(hist, bw_low, bw_high, step, kernel, adaptive,
                              ise_divisions, npx)
  if not rows:
    return None, rows
  best = min(rows, key=lambda r: r[metric])
  return best["bandwidth"], rows


def print_bandwidth_summary(kde, hist, rho):
  xs, hs = local_bandwidths(kde, hist)
  if not hs:
    return
  hmin, hmax = min(hs), max(hs)
  hmean = sum(hs) / len(hs)
  mode = "adaptive" if getattr(kde, "_kde_adaptive", True) else "fixed"
  print(f"KDE bandwidth mode: {mode} (rho={rho})")
  print(f"  local h(x): min={hmin:.4g}, max={hmax:.4g}, mean={hmean:.4g} ({len(hs)} points)")


def optimized_kde(hist, min_bw, kernel=DEFAULT_KDE_KERNEL, adaptive=DEFAULT_KDE_ADAPTIVE):
  hist_normalized = hist.Clone(hist.GetTitle())
  hist_normalized.SetDirectory(0)
  den_int = hist_normalized.Integral("width")
  if den_int > 0:
    hist_normalized.Scale(1.0 / den_int)

  print("KDE kernel:", kernel)
  optimal_bw, bw_scan = find_optimal_bandwidth(
    hist_normalized, min_bw, 1.0, 0.005, kernel=kernel, adaptive=adaptive, metric="mse")
  print("optimal rho", optimal_bw)

  kde_fn = None
  if optimal_bw is not None:
    kde_fn = kde_pdf_for_bandwidth(
      hist_normalized, optimal_bw, kernel=kernel, adaptive=adaptive)
    if kde_fn and kde_fn._hold_kde:
      print_bandwidth_summary(kde_fn._hold_kde, hist_normalized, optimal_bw)

  if kde_fn:
    print("MSE (bin centers, normalized mm1x vs KDE):", MSE(kde_fn, hist_normalized))
    print("chi sq sum:", chi_sq(kde_fn, hist_normalized))

  return kde_fn


def pol2_fn(hist):
  hist_normalized = hist.Clone(hist.GetTitle())
  hist_normalized.SetDirectory(0)
  den_int = hist_normalized.Integral("width")
  if den_int > 0:
    hist_normalized.Scale(1.0 / den_int)


  pol2 = hist_normalized.Clone("pol4")
  pol2.SetDirectory(0)
  pol2_fit = pol2.Fit("pol4", "SQ0")
  pol2_fn = pol2.GetFunction("pol4") if pol2_fit and pol2_fit.Status() == 0 else None

  if pol2_fn:
    pol2_fn._hold_pol2 = pol2
    print("\npol2 coefficients (p0 + p1*x + p2*x^2):")
    for i in range(pol2_fn.GetNpar()):
      print(f"  p{i} = {pol2_fn.GetParameter(i):.6g} ± {pol2_fn.GetParError(i):.6g}")
    print("MSE (normalized mm1x vs pol2):", MSE(pol2_fn, hist_normalized))
    print("REDUCED chi sq sum:", chi_sq(pol2_fn, hist_normalized) / 4)
  return pol2_fn



def format_output(hist, kde, pol2, adaptive=DEFAULT_KDE_ADAPTIVE):
  # Drawing: KDE + local bandwidth (left) vs polynomial fit (right)
  c1 = ROOT.TCanvas("c1", "KDE vs pol4", 1400, 600)
  c1.Divide(2, 1)

  hist_normalized = hist.Clone(hist.GetTitle())
  hist_normalized.SetDirectory(0)
  den_int = hist_normalized.Integral("width")
  if den_int > 0:
    hist_normalized.Scale(1.0 / den_int)

  c1.cd(1)
  pad_pdf = ROOT.TPad("pdf", "pdf", 0, 0.32, 1, 1)
  pad_bw = ROOT.TPad("bw", "bw", 0, 0, 1, 0.32)
  pad_pdf.Draw()
  pad_bw.Draw()

  pad_pdf.cd()
  mode = "adaptive" if adaptive else "fixed"
  hist_normalized.SetMarkerSize(0.8)
  hist_normalized.Draw("E1 HIST")
  hist_normalized.SetTitle(f"KDE fit ({mode} bandwidth)")
  tkde = getattr(kde, "_hold_kde", None) if kde else None
  if kde:
    kde.SetLineColor(ROOT.kBlue)
    kde.SetLineWidth(2)
    kde.Draw("SAME C")

  pad_bw.cd()
  if tkde:
    bw_graph = adaptive_bandwidth_graph(tkde, hist_normalized)
    if bw_graph:
      bw_graph.Draw("AL")
      pad_bw.SetGridy()

  c1.cd(2)
  hist_pol2 = hist_normalized.Clone()
  hist_pol2.SetDirectory(0)
  hist_pol2.SetMarkerSize(0.8)
  hist_pol2.Draw("E1 HIST")
  hist_pol2.SetTitle("pol4 fit")
  if pol2:
    pol2.SetLineColor(ROOT.kRed)
    pol2.SetLineWidth(2)
    pol2.Draw("SAME")

  c1.Update()
  input()




hist = mm1y

# seems like 0.18-0.22 are best for gaussian, not sure for others
# 0.4 is about the max before it becomes nonsensical. The max is currently 1 in the def of optimized_kde
min_bw = 0.2
format_output(
  hist,
  optimized_kde(hist, min_bw, adaptive=DEFAULT_KDE_ADAPTIVE),
  pol2_fn(hist),
  adaptive=DEFAULT_KDE_ADAPTIVE,
)
