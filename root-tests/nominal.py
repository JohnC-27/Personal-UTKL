import ROOT

nominal_file = ROOT.TFile.Open("/Users/johnculbertson/Documents/Personal-UTKL/root-test-files/nominal.root", "READ")
## nominalxyposMM1, 2, 3

h1 = nominal_file.Get("nominalxyposMM1")
h2 = nominal_file.Get("nominalxyposMM2")
h3 = nominal_file.Get("nominalxyposMM3")

h1x = h1.ProjectionX()
h1y = h1.ProjectionY()
h2x = h2.ProjectionX()
h2y = h2.ProjectionY()
h3x = h3.ProjectionX()
h3y = h3.ProjectionY()


h1x_fit = h1x.Fit("pol3", "S")



c1 = ROOT.TCanvas('c1', "TITLE", 1200, 600)
c1.Divide(3,1)

c1.cd(1)
h1.Draw()

c1.cd(2)
h1x.Draw()

c1.cd(3)
h1y.Draw()

c1.Update()
input()
