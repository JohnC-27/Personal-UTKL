# AI folder index

## Project purpose:
- The purpose of this project is to analyze and model 2d histograms.

## ai/ logs and documentation purpose:
- Working memory for **agents and humans**. Update these files over burying decisions in chat only.

| File | Purpose | Update when |
|------|---------|-------------|
| [instructions.md](instructions.md) | Standing rules and conventions | Rules change |
| [active_task.md](active_task.md) | Current objective, acceptance criteria, blockers | Task starts / changes / completes |
| [work_log.md](work_log.md) | Chronological what-was-done | End of meaningful session or PR-sized chunk |
| [context_log.md](context_log.md) | Durable facts, decisions, “gotchas” | Something new must survive the next session |
| [debugging_log.md](debugging_log.md) | Hypothesis-driven debug record | Investigating or resolving a bug |

## Boot sequence for agents: 
`ai/instructions.md` → `ai/active_task.md` → `docs/overview.md` + `docs/structure.md` → skim latest entries in the three logs.
