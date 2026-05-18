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


# Weighted TKDE: bin centers + contents. Unweighted TKDE(...) expects x samples, not bin counts.
# (1D histogram, bandwidth (float))
# returns the kde object
def weighted_kde(hist, bandwidth):
  xs, ws = [], []
  for i in range(1, hist.GetNbinsX() + 1):
    c = hist.GetBinContent(i)
    if c <= 0:
      continue
    xs.append(hist.GetBinCenter(i))
    ws.append(c)
  xlo = hist.GetXaxis().GetXmin()
  xhi = hist.GetXaxis().GetXmax()
  if len(xs) == 0:
    return None
  da = array.array("d", xs)
  wa = array.array("d", ws)
  kde = ROOT.TKDE(len(da), da, wa, xlo, xhi, "", bandwidth)
  kde.SetKernelType(ROOT.TKDE.kEpanechnikov)
  return kde


# TKDE::GetFunction() is not guaranteed to integrate to 1 on [xlo, xhi]. Histogram
# density from Scale(1/Integral("width")) does. Renormalize so ∫ kde dx = 1.
def renorm_kde_tf1(raw_func, xlo, xhi, npx=2000, epsrel=1e-9):
  if raw_func is None:
    return None
  raw_func.SetNpx(npx)
  den = raw_func.Integral(xlo, xhi, epsrel)
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
def kde_pdf_for_bandwidth(hist, bandwidth, npx=2000, epsrel=1e-9):
  kde = weighted_kde(hist, bandwidth)
  if kde is None:
    return None
  xlo = hist.GetXaxis().GetXmin()
  xhi = hist.GetXaxis().GetXmax()
  pdf = renorm_kde_tf1(kde.GetFunction(), xlo, xhi, npx, epsrel)
  if pdf is not None:
    pdf._hold_kde = kde
  return pdf


# Scan bandwidths; returns list of {bandwidth, mse, ise, chi2} (failed fits omitted).
def scan_kde_bandwidths(hist, bw_low, bw_high, step, ise_divisions=100, npx=2000):
  results = []
  bw = bw_low
  while bw <= bw_high + 1e-12 * max(abs(bw_high), 1.0):
    kde_fn = kde_pdf_for_bandwidth(hist, bw, npx=npx)
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
def find_optimal_bandwidth(hist, bw_low, bw_high, step, metric="mse", ise_divisions=100, npx=2000):
  rows = scan_kde_bandwidths(hist, bw_low, bw_high, step, ise_divisions, npx)
  if not rows:
    return None, rows
  best = min(rows, key=lambda r: r[metric])
  return best["bandwidth"], rows


def optimized_kde(hist):
  hist_normalized = hist.Clone(hist.GetTitle())
  hist_normalized.SetDirectory(0)
  den_int = hist_normalized.Integral("width")
  if den_int > 0:
    hist_normalized.Scale(1.0 / den_int)

  optimal_bw, bw_scan = find_optimal_bandwidth(hist_normalized, 0.21, 1.0, 0.005, metric="mse")
  print("optimal bandwidth", optimal_bw)

  kde_fn = kde_pdf_for_bandwidth(hist_normalized, optimal_bw) if optimal_bw is not None else None

  if kde_fn:
    print("MSE (bin centers, normalized mm1x vs KDE):", MSE(kde_fn, hist_normalized))
    print("chi sq sum:", chi_sq(kde_fn, hist_normalized))

  return kde_fn


def pol9_fn(hist):
  hist_normalized = hist.Clone(hist.GetTitle())
  hist_normalized.SetDirectory(0)
  den_int = hist_normalized.Integral("width")
  if den_int > 0:
    hist_normalized.Scale(1.0 / den_int)

  pol9 = hist_normalized.Clone("pol9")
  pol9.SetDirectory(0)
  pol9_fit = pol9.Fit("pol9", "SQ0")
  pol9_fn = pol9.GetFunction("pol9") if pol9_fit and pol9_fit.Status() == 0 else None

  if pol9_fn:
    pol9_fn._hold_pol9 = pol9
    print("\nMSE (normalized mm1x vs pol9):", MSE(pol9_fn, hist_normalized))
    print("REDUCED chi sq sum: ", (chi_sq(pol9_fn, hist_normalized) / 9))
  return pol9_fn



def format_output(hist, kde, pol9):
  # Drawing: KDE (left) vs degree-9 polynomial (right)
  c1 = ROOT.TCanvas("c1", "KDE vs pol9", 1400, 600)
  c1.Divide(2, 1)

  hist_normalized = hist.Clone(hist.GetTitle())
  hist_normalized.SetDirectory(0)
  den_int = hist_normalized.Integral("width")
  if den_int > 0:
    hist_normalized.Scale(1.0 / den_int)

  c1.cd(1)
  hist_normalized.Draw("HIST")
  hist_normalized.SetTitle("KDE fit")
  if kde:
    kde.SetLineColor(ROOT.kBlue)
    kde.SetLineWidth(2)
    kde.Draw("SAME C")

  c1.cd(2)
  hist_pol9 = hist_normalized.Clone()
  hist_pol9.SetDirectory(0)
  hist_pol9.Draw("HIST")
  hist_pol9.SetTitle("pol9 fit")
  if pol9:
    pol9.SetLineColor(ROOT.kRed)
    pol9.SetLineWidth(2)
    pol9.Draw("SAME")

  c1.Update()
  input()




hist = mm1x
format_output(hist, optimized_kde(hist), pol9_fn(hist))
