import ROOT
import array

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
  return kde


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
    sse += (func_value - hist_value) ** 2
  return sse / nb


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


mm1x_kde_obj = weighted_kde(mm1x, 1.9)
mm1x_kde = mm1x_kde_obj.GetFunction()
mm1x_kde = mm1x_kde / mm1x_kde.Integral(mm1x_kde.GetXaxis().GetXmin(), mm1x_kde.GetXaxis().GetXmax())



mm1x_normalized = mm1x.Clone("mm1x_normalized")
mm1x_normalized.SetDirectory(0)
den_int = mm1x_normalized.Integral("width")
if den_int > 0:
  mm1x_normalized.Scale(1.0 / den_int)


if mm1x_kde:
  print("MSE (bin centers, normalized mm1x vs KDE):", MSE(mm1x_kde, mm1x_normalized))
  print("ISE:", quasi_ISE(mm1x_kde, mm1x_normalized, 100))


# Drawing
c1 = ROOT.TCanvas("c1", "TITLE", 1200, 800)
mm1x_normalized.Draw("HIST")
if mm1x_kde:
  mm1x_kde.SetLineColor(ROOT.kBlue)
  mm1x_kde.Draw("SAME C")

c1.Update()
input()
