import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kWarning  # or kError to hide warnings too


nominal_file = ROOT.TFile.Open("/Users/johnculbertson/Documents/Personal-UTKL/root-test-files/nominal.root", "READ")
## three TH2D named nominalxyposMM1, 2, 3

mm1 = nominal_file.Get("nominalxyposMM1")
mm2 = nominal_file.Get("nominalxyposMM2")
mm3 = nominal_file.Get("nominalxyposMM3")

c1 = ROOT.TCanvas("c1", "nominal xy pos MM1,2,3", 1500, 500)
c1.Divide(3, 1)

for i, h in enumerate([mm1, mm2, mm3], 1):
  c1.cd(i)
  h.SetStats(0)
  h.DrawCopy("COLZ")

c1.Update()
input()
