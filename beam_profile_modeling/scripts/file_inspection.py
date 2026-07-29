import os
import sys

import ROOT

INPUT_ROOT_FILE = os.path.abspath(
  os.path.join(
    os.path.dirname(__file__), "..", "root_files", "Jan2026studies_nominal.root"
  )
)

if not os.path.isfile(INPUT_ROOT_FILE):
  sys.exit(f"ROOT file not found: {INPUT_ROOT_FILE}")

f = ROOT.TFile.Open(INPUT_ROOT_FILE, "READ")
if not f or f.IsZombie():
  sys.exit(f"cannot open {INPUT_ROOT_FILE}")

f.ls()
