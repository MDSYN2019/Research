# Research Monorepo (MD, Free-Energy, and Parallel Computing Work)

This repository combines several years of experimental and reference work related to:

- molecular dynamics and free-energy estimation (including Jarzynski/Bennett workflows),
- MPI/OpenMP/C++ experimentation,
- CMake learning and cookbook-style examples,
- prototype and work-in-progress research code.

The previous root README was minimal and outdated (last updated 2021). This version provides a current navigation map and maintenance conventions so the repo is easier to understand and evolve.

## Repository layout

Top-level directories and their main purpose:

- `Free_Energy/` — primary free-energy analysis work (Jarzynski, Bennett, and related computation).
- `free_energy/` — legacy/alternate free-energy tree (kept for compatibility and historical reference).
- `Parallel/` — MPI/OpenMP and modern C++ experiments, with docs/tests/third-party code.
- `C12E2_work/` — molecular interaction and Hamaker-related calculations plus archived material.
- `CMake_practice/` — cookbook-style CMake examples and notes.
- `Graph/` — graphing/plot-related materials.
- `WIP/` — active prototypes and unfinished studies.
- `workshop/` — workshop examples and headers.

## Current state and maintenance notes

This is a research-first monorepo with mixed maturity levels. Some folders are polished; others are snapshots.

To keep the repository manageable:

1. **Do not remove historical folders lightly** — many are referenced in notes/scripts.
2. **Prefer adding clear README files in subfolders** when introducing or reviving work.
3. **Keep generated artifacts out of version control** (see root `.gitignore`).
4. **Treat `WIP/` as unstable** and keep runnable/documented code in domain folders when ready.

## Quick navigation tasks

- Find CMake entry points:
  - `rg --files -g 'CMakeLists.txt'`
- Find tests:
  - `rg --files -g '*test*' Parallel/test Free_Energy free_energy`
- Identify documentation roots:
  - `rg --files -g 'README*'`

## Legacy references

Earlier project notes referenced these classic free-energy papers:

1. *Nonequilibrium equality for free energy differences* (Phys. Rev. Lett. 78, 2690).
2. *Free energy reconstruction from nonequilibrium single-molecule pulling experiments* (PNAS 98, 7, 3658–3661).
3. *Free energy calculation from steered molecular dynamics simulations using Jarzynski's equality* (JCP 119, 6, 2003).
4. *Bias and error in estimates of equilibrium free-energy differences from nonequilibrium measurements* (PNAS 100, 22, 12564–12569).

---

If you are onboarding to this repository, start with [`REPO_STRUCTURE.md`](REPO_STRUCTURE.md) for a higher-level organization and cleanup roadmap.
