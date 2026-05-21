import ROOT

file = ROOT.TFile.Open("/Users/johnculbertson/Documents/Personal-UTKL/root-test-files/horn_current_p_1sigma.root", "READ")
## horn_current_p_1sigmaxyposMM1, 2, 3
## 


h1 = file.Get("horn_current_p_1sigmaxyposMM1")
h2 = file.Get("horn_current_p_1sigmaxyposMM2")
h3 = file.Get("horn_current_p_1sigmaxyposMM3")

h1x = h1.ProjectionX()
h2x = h2.ProjectionX()
h3x = h3.ProjectionX()


c1 = ROOT.TCanvas('c1', "TITLE")
pad1 = ROOT.TPad("pad1", "1", 0, 0, 0.5, 0.5)
pad2 = ROOT.TPad("pad2", "2", 0.5, 0.5, 1, 1)
pad3 = ROOT.TPad("pad3", "3", 0, 0.5, .45, 1)
pad1.Draw()
pad2.Draw()
pad3.Draw()

pad1.cd()
h1x.DrawCopy()

pad2.cd()
h2x.DrawCopy()

pad3.cd()
h3x.DrawCopy()

c1.Update()
input()