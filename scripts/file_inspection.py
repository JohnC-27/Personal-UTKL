import ROOT
import os

INPUT_ROOT_FILE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "kde_chi2_fit.root"
)
f = ROOT.TFile.Open(INPUT_ROOT_FILE, "READ")
f.ls()