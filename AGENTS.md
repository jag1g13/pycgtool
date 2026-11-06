# AGENTS

## General

- Favour simplicity in implementation where possible.
- Don't add additional dependencies / tools unless they're useful.
- Ask before adding any new dependencies / tools.
- When a change affects features, setup steps, project structure, tech stack, or data model, flag that README.md / docs/ need updating and propose the wording, but ask before editing prose — AI-written prose tends to read poorly.

## Dependencies & tooling

- Manage dependencies with `uv` only (`pyproject.toml` + `uv.lock`); never reintroduce poetry or edit the lock by hand.
- Dependency updates must respect the `exclude-newer` cooldown in `[tool.uv]`; use `uv lock --upgrade` and run the full test suite afterwards.
- The minimum supported Python version is set by `requires-python`; bumping it (or the CI matrix) requires explicit confirmation, and `xdrlib`-style stdlib backport shims must be re-checked when it changes.
- Reference data (e.g. `src/pycgtool/data/`) is scientific input: never regenerate or silently reformat it; document the provenance of any change.

## Scientific correctness

- PyCGTOOL produces physical models: never change numerical behaviour (parameters, units, tolerances, averaging, mapping rules) casually. State the intent and verify against reference results before and after.
- Keep units consistent and explicit in code and docs (mdtraj uses nm/ps/K; Gromacs `.gro`/`.itp` files have their own conventions). Convert at boundaries, never implicitly.
- Use fixed random seeds in tests; no test may depend on wall-clock time, network access, or machine ordering.
- Prefer vectorised NumPy idioms, but profile before optimising; correctness and readability beat cleverness in research code.
- Guard numerical edge cases (division by zero, empty selections, singular fits) with clear error messages, not silent fallbacks.

## Testing & quality gates

- Before any commit: `uv run pytest tests/` (all pass) and `uv run ruff check` (clean). CI enforces coverage `--fail-under=80`; new features need tests and must not drop coverage.
- Bug fixes need a regression test that fails without the fix.
- Don't commit during manual/interactive testing until given confirmation.

## Git workflow

- Commit after each completed task.
- Use Conventional Commits format for commit messages.
- When committing, credit AI assistance using the Linux kernel convention: an "Assisted-by: <harness> <model>" trailer.
- Commits must be signed - if the key isn't loaded, ask for it.
- Don't force push commits without asking.
- Never commit generated or local-only artefacts: `dist/`, `.venv/`, `*.egg-info/`, `.cauldron/`, local scratch data.
