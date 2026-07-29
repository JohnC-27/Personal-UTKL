import ROOT

horn_file = ROOT.TFile.Open("/Users/johnculbertson/Documents/Personal-UTKL/root-test-files/horn_current_p_1sigma.root", "READ")
nominal_file = ROOT.TFile.Open("/Users/johnculbertson/Documents/Personal-UTKL/root-test-files/nominal.root", "READ")

h1 = nominal_file.Get("nominalxyposMM1")
h2 = horn_file.Get("horn_current_p_1sigmaxyposMM1")


c1 = ROOT.TCanvas('c1', "TITLE")
pad1 = ROOT.TPad("pad1", "Alcove1", 0, 0, 1, 0.5)
pad2 = ROOT.TPad("pad2", "Alcove2", 0, 0.5, 1, 1)
pad1.Draw()
pad2.Draw()

pad1.cd()
h1.DrawCopy()

pad2.cd()
h2.DrawCopy()

c1.Update()
input()