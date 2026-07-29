# Debugging log

Hypothesis-driven record for bugs and anomalies. Goal: another agent can continue without replaying the whole investigation.

## How to append

Newest cases on top. One case per incident; update the same case until closed.

```
### YYYY-MM-DD — <short symptom title>
- **Status:** open | monitoring | resolved | wontfix
- **Symptom:** what was observed (plot, exception, bad number)
- **Repro:** command / script / inputs
- **Hypotheses:**
  | ID | Claim | Result |
  |----|-------|--------|
  | H1 | ... | confirmed / rejected / inconclusive |
- **Evidence:** key measurements, log lines, paths (not full NDJSON dumps)
- **Root cause:** if known
- **Fix:** what changed
- **Verification:** how you know it is fixed
- **Residual risk:** what might still be wrong
```

---

### Template only — no open case yet

_No formal debugging case opened at docs bootstrap. Prior raw probes exist in `.cursor/debug-8940bd.log` and `.cursor/debug-8dd097.log` (projection / KDE diagnostics in `plot_2d_kde.py`). If that work resumes, promote a real case entry here with H1/H2 results._

### Seed: projection scale checks (`plot_2d_kde.py`)

- **Status:** monitoring (historical; not fully promoted to a closed case)
- **Symptom:** KDE projection integrals (`sum_kde_px` / `sum_kde_py`) can disagree with data projection sums by a large factor unless scale factors align; pointwise `kde_px` vs `data_px` at a probe x can look close while sums do not.
- **Repro:** run `scripts/plot_2d_kde.py` with debug logging to `.cursor/debug-8940bd.log` (`hypothesisId` H1-H2, location `plot_2d_kde.py:plot_projections`).
- **Hypotheses:**
  | ID | Claim | Result |
  |----|-------|--------|
  | H1 | Projection scale related to bin width / grid factor (e.g. 0.125) | inconclusive in this log alone — needs human confirmation of intended convention |
  | H2 | Marginalization grid (`n_integrate`) or asymmetric `x_scale`/`y_scale` causes sum mismatch | inconclusive — logs show trials with `x_scale`/`y_scale` near 0.125 / 0.115 |
- **Evidence:** NDJSON entries in `.cursor/debug-8940bd.log` with fields `data_px`, `kde_px`, `sum_data_px`, `sum_kde_px`, `x_scale`, `y_scale`, `n_integrate`
- **Root cause:** not recorded as confirmed in `ai/` yet
- **Fix:** none logged here
- **Verification:** n/a
- **Residual risk:** regenerating projection plots without checking integral ratios may hide a global scale error
