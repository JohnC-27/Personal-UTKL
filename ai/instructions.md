# Instructions / rules (agent ↔ human)

Standing constraints for work in **Personal-UTKL**. Humans: edit this file when expectations change. Agents: follow it unless the human overrides in the current message.

## Orientation

1. Read `ai/active_task.md` before coding.
2. Read `docs/overview.md`, `docs/structure.md`, and `docs/goals.md` for project shape.
3. Prefer facts in `ai/context_log.md` over guessing filenames or hist names.
4. After meaningful work, append `ai/work_log.md`. After debug work, append `ai/debugging_log.md`. Keep `active_task.md` current.

## Scope discipline

- Change only files needed for the stated task. Do not “clean up” unrelated scripts, plots, or ROOT files.
- Do not delete files without explicit human permission.
- Do not commit or push unless the human asks.
- Do not edit outside the repo (or outside the prompt’s scope) without permission.
- Prefer READ-ONLY inspection when writing is not required.

## Code conventions

- **Indent with two spaces** (never four).
- New functions must be under 60 lines.
- 
- Match existing script style: top-of-file path/config constants, dataclasses where already used, ROOT batch mode for plot/fit scripts.
- Regenerate plots via scripts over hand-editing `plots/` binaries.
- Do not introduce new absolute home-directory paths; use paths relative to the script/`repo root` like neighboring files.
- Do not add unsolicited markdown docs beyond what was requested; when updating `ai/` or `docs/`, keep entries factual and scannable.

## Physics / analysis hygiene

- Treat ROOT hist names and file names as fragile: verify with `file_inspection.py` or `TFile` listing before assuming.
- Always reference ROOT documentation to verify function and class operation.
- When changing KDE or projection math, check **normalization** (integrals, bin widths, marginal scales) — scale mismatches have been a recurring failure mode.
- Distinguish **fit scripts** (write ROOT) from **plot scripts** (read fit ROOT, write `plots/`).
- Summarize durable debug findings from `.cursor/debug-*.log` into `ai/debugging_log.md`; do not rely on NDJSON dumps as the only record.

## Communication

Write log entries so a **new agent or a tired human** can resume:

- Lead with outcome / status, then evidence.
- Use concrete paths, hist names, parameter values, and commands.
- Separate **observation** vs **hypothesis** vs **decision**.
- Mark uncertainty explicitly (`Unknown`, `Suspected`, `Confirmed`).
- Prefer short tables and bullets over long prose.

## Entry templates

Use the section templates in each log file. Newest entries **on top** (reverse chronological) unless a file says otherwise.

## Environment notes

- Python venv: `.venv/` (NumPy/SciPy/Matplotlib/Pandas available).
- Most `scripts/` require a working **ROOT / PyROOT** install in the shell that runs them.
- `swic/` analysis is primarily NumPy/SciPy/Matplotlib and may not need ROOT.
