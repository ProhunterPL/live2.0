#!/usr/bin/env python3
"""Check which runs have thermodynamics failures"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.validation import TruthFilter, ValidationLevel
import json

filter = TruthFilter(ValidationLevel.MEDIUM)

failures = []
for scenario_dir in Path('results/phase2b_additional').iterdir():
    if not scenario_dir.is_dir():
        continue
    
    for run_dir in scenario_dir.iterdir():
        if not run_dir.is_dir() or not run_dir.name.startswith('run_'):
            continue
        
        result = filter.filter_run(str(run_dir))
        thermo_status = result.filters['thermodynamics'].status
        
        if thermo_status.value == 'FAIL':
            failures.append({
                'run': str(run_dir),
                'scenario': scenario_dir.name,
                'errors': result.filters['thermodynamics'].errors,
                'details': result.filters['thermodynamics'].details
            })

print(f"Thermodynamics failures: {len(failures)}")
print("\nDetails:")
for f in failures[:10]:
    print(f"\n{f['run']}:")
    print(f"  Errors: {f['errors']}")
    print(f"  Details: {f['details']}")

# Save to JSON
with open('validation_reports/phase2b/thermodynamics_failures.json', 'w') as out:
    json.dump(failures, out, indent=2)

print(f"\n✅ Saved to validation_reports/phase2b/thermodynamics_failures.json")

