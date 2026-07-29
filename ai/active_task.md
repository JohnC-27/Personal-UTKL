# Active task

Single source of truth for **what we are doing right now**. Replace the body when the task changes; move completed work into `work_log.md`.

---

## Status

`idle` — no active analysis task after docs bootstrap.

| Field | Value |
|-------|--------|
| Status | `idle` \| `in_progress` \| `blocked` \| `done` |
| Owner | human / agent / shared |
| Started | — |
| Target done | — |

## Objective

_(One sentence: what “done” means.)_

Bootstrap complete: `docs/` and `ai/` initial markdown exist for agent ↔ human coordination. Awaiting next analysis or engineering task from the human.

## Acceptance criteria

- [x] `docs/` describes overview, structure, and goals
- [x] `ai/` has instructions, active task, work log, context log, debugging log
- [ ] _(Next task criteria go here)_

## In scope

- Documentation and AI working files only (this bootstrap)

## Out of scope

- Script/ROOT/plot changes unless a new active task says so

## Context pointers

- `docs/overview.md`, `docs/structure.md`, `docs/goals.md`
- Recent KDE plot work lived in `scripts/plot_2d_kde.py` and `.cursor/debug-8940bd.log` (projection scale checks)

## Blockers

None.

## Next actions

1. Human: set the next objective in this file (or tell the agent to).
2. Agent: on new task, flip Status to `in_progress` and fill Objective + Acceptance criteria before editing code.

---

## How to update (agents)

When starting work:

```
Status: in_progress
Objective: <clear done-state>
Acceptance criteria: checklist
Next actions: ordered steps
```

When blocked: set Status `blocked`, name the blocker, and what input is needed.
When finished: set Status `done`, then append a short entry to `work_log.md` and reset this file toward `idle` or the next task.
