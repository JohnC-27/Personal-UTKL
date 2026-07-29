import pandas as pd
import os

INPUT_TXT_FILE = os.path.abspath(
  os.path.join(
    os.path.dirname(__file__), "..", "root_files", "Jan2026studies_nominal.txt"
  )
)
OUTPUT_CSV_FILE = os.path.abspath(
  os.path.join(
    os.path.dirname(__file__), "..", "root_files", "raw_xypos_nominal.csv"
  )
)

# Load the tab-separated file
df = pd.read_csv(INPUT_TXT_FILE, sep='\t')

# Save it as a comma-separated file
df.to_csv(OUTPUT_CSV_FILE, index=False)