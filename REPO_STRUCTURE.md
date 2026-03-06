# Repository Structure and Reorganization Plan

This document defines a lightweight structure for maintaining this research monorepo without breaking historical experiments.

## 1) Domain ownership map

| Domain | Primary path(s) | Notes |
|---|---|---|
| Free-energy methods | `Free_Energy/`, `free_energy/` | Two trees currently coexist; `Free_Energy/` should be preferred for active updates. |
| Parallel programming | `Parallel/` | Contains source, tests, docs, and third-party code. |
| Molecular interaction studies | `C12E2_work/` | Includes Hamaker and archived calculations. |
| Build-system learning | `CMake_practice/` | Educational/reference examples. |
| Experimental drafts | `WIP/` | Unstable and exploratory. |

## 2) Reorganization policy

To avoid disrupting old workflows, use **soft reorganization** first:

- Improve discoverability with docs and consistent folder READMEs.
- Add ignore rules to reduce noise from generated output.
- Introduce naming conventions for new work.
- Delay hard moves/renames unless actively migrating scripts.

## 3) Naming conventions for new additions

- Use lowercase with underscores for new top-level folders if needed.
- Put runnable source under `src/` and tests under `test/` or `tests/` inside a project folder.
- Co-locate per-project README with build/run commands.
- Avoid committing generated HTML docs, binaries, and local build dirs unless they are explicit release artifacts.

## 4) Incremental cleanup checklist

- [ ] Add/refresh `README.md` in each major top-level folder.
- [ ] Standardize build directories to `<project>/build/` and ignore them in git.
- [ ] Identify duplicate projects across `Free_Energy/` and `free_energy/` for potential consolidation.
- [ ] Tag legacy-only folders with a short `README.legacy.md`.
- [ ] Add CI checks selectively for actively maintained modules.

## 5) Suggested next migration milestone

**Milestone A (non-breaking):**

1. Document entry points for `Parallel/` and `Free_Energy/`.
2. Add smoke tests for at least one representative module in each area.
3. Capture known-good compiler/toolchain versions in per-module docs.

