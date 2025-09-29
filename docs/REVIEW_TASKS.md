# Suggested Follow-up Tasks

## Typographical Fix
- Correct the typo "Menadżer pakietów Python" to "Menedżer pakietów Pythona" in `INSTALLATION.md` (line 12).
  - Location: INSTALLATION.md line 12

## Bug Fix
- Update `Grid.add_particle` to respect the configured `max_particles` instead of the hard-coded limit of 10,000 so simulations with different capacities don't reject particles prematurely.
  - Location: backend/sim/core/grid.py lines 100-108

## Documentation Alignment
- Adjust the backend startup instructions (e.g., `README.md` Quick Start) to use `python api/server.py` or add an `__init__.py` so that the documented `python -m api.server` command works.
  - Location: README.md lines 41-45, backend/api/ (missing __init__.py)

## Test Enhancement
- Extend `test_save_snapshot` to verify that the snapshot file is actually written to disk (and clean it up afterwards) instead of only checking the HTTP response payload.
  - Location: backend/tests/test_api.py lines 190-203

