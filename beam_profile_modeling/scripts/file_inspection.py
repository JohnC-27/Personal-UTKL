import os
import sys

import ROOT

INPUT_ROOT_FILE = os.path.abspath(
  os.path.join(
    os.path.dirname(__file__), "..", "root_files", "mz_nominal_100bin_run1.root"
  )
)

if not os.path.isfile(INPUT_ROOT_FILE):
  sys.exit(f"ROOT file not found: {INPUT_ROOT_FILE}")

f = ROOT.TFile.Open(INPUT_ROOT_FILE, "READ")
if not f or f.IsZombie():
  sys.exit(f"cannot open {INPUT_ROOT_FILE}")

f.ls()

print("\nTH2 sizes:")
for key in f.GetListOfKeys():
  obj = key.ReadObj()
  if not obj or not obj.InheritsFrom("TH2"):
    continue
  nx = obj.GetNbinsX()
  ny = obj.GetNbinsY()
  xaxis = obj.GetXaxis()
  yaxis = obj.GetYaxis()
  print(
    f"  {obj.GetName()} [{obj.ClassName()}]: "
    f"{nx} x {ny} bins "
    f"(x: [{xaxis.GetXmin()}, {xaxis.GetXmax()}], "
    f"y: [{yaxis.GetXmin()}, {yaxis.GetXmax()}])"
  )

f.Close()
