import ROOT

file = ROOT.TFile.Open("/Users/johnculbertson/Documents/Personal-UTKL/root-test-files/proton_beam_radius_p_1sigma.root", "READ")

# proton_beam_raduis_p_1sigmaxyposMM1, 2, 3


file.ls()

h1 = file.Get("proton_beam_radius_p_1sigmaxyposMM1")
h2 = file.Get("proton_beam_radius_p_1sigmaxyposMM2")
h3 = file.Get("proton_beam_radius_p_1sigmaxyposMM3")

c1 = ROOT.TCanvas('c1', "TITLE")
pad1 = ROOT.TPad("pad1", "1", 0, 0, 0.5, 0.5)
pad2 = ROOT.TPad("pad2", "2", 0.5, 0.5, 1, 1)
pad3 = ROOT.TPad("pad3", "3", 0, 0.5, .45, 1)
pad1.Draw()
pad2.Draw()
pad3.Draw()

pad1.cd()
h1.DrawCopy()

pad2.cd()
h2.DrawCopy()

pad3.cd()
h3.DrawCopy()

c1.Update()
input()