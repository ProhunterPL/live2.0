#!/usr/bin/env python3
"""
Fill Results Section 3.3 (Autocatalytic Cycles) with analysis data

Usage:
    python scripts/fill_autocatalysis_section.py \
        --data paper/results_data \
        --template paper/manuscript_draft.tex \
        --output paper/manuscript_draft_filled.tex
"""

import json
import argparse
import re
from pathlib import Path
from typing import Dict


def load_analysis_data(data_dir: Path) -> Dict:
    """Load all analysis data files"""
    data = {}
    
    # Load scenario comparison
    comparison_file = data_dir / "scenario_comparison.json"
    if comparison_file.exists():
        with open(comparison_file) as f:
            data['comparison'] = json.load(f)
    
    # Load individual scenario analyses
    scenarios = ['miller_urey_extended', 'hydrothermal_extended', 'formamide_extended']
    for scenario in scenarios:
        scenario_file = data_dir / f"{scenario}_analysis.json"
        if scenario_file.exists():
            with open(scenario_file) as f:
                data[scenario] = json.load(f)
    
    return data


def extract_autocatalysis_values(data: Dict) -> Dict:
    """Extract all values needed for Section 3.3"""
    values = {}
    
    comp = data.get('comparison', {})
    auto_comp = comp.get('autocatalysis_comparison', {})
    
    # Total unique cycles (sum across all scenarios)
    miller = auto_comp.get('miller_urey_extended', {})
    hydro = auto_comp.get('hydrothermal_extended', {})
    form = auto_comp.get('formamide_extended', {})
    
    total_cycles = (
        miller.get('total_cycles', 0) +
        hydro.get('total_cycles', 0) +
        form.get('total_cycles', 0)
    )
    values['total_unique_cycles'] = f"{total_cycles:,}"
    
    # Cycles per run by scenario
    miller_cycles_per_run = miller.get('cycles_per_run', 0)
    hydro_cycles_per_run = hydro.get('cycles_per_run', 0)
    form_cycles_per_run = form.get('cycles_per_run', 0)
    
    # Get std from individual analyses
    miller_data = data.get('miller_urey_extended', {}).get('autocatalysis', {})
    hydro_data = data.get('hydrothermal_extended', {}).get('autocatalysis', {})
    form_data = data.get('formamide_extended', {}).get('autocatalysis', {})
    
    miller_std = miller_data.get('cycles_per_run_std', 0)
    hydro_std = hydro_data.get('cycles_per_run_std', 0)
    form_std = form_data.get('cycles_per_run_std', 0)
    
    values['miller_cycles'] = f"{miller_cycles_per_run:.0f} ± {miller_std:.0f}"
    values['hydro_cycles'] = f"{hydro_cycles_per_run:.0f} ± {hydro_std:.0f}"
    values['form_cycles'] = f"{form_cycles_per_run:.0f} ± {form_std:.0f}"
    
    # Cycle types (aggregate across all scenarios)
    all_cycle_types = {
        'direct': 0,
        'indirect': 0,
        'hypercycle': 0
    }
    
    for scenario_data in [miller, hydro, form]:
        cycle_types = scenario_data.get('cycle_types', {})
        for cycle_type, count in cycle_types.items():
            all_cycle_types[cycle_type] = all_cycle_types.get(cycle_type, 0) + count
    
    total_by_type = sum(all_cycle_types.values())
    if total_by_type > 0:
        values['simple_loops_pct'] = f"{100 * (all_cycle_types.get('direct', 0) / total_by_type):.1f}"
        # Medium loops = indirect cycles
        values['medium_loops_pct'] = f"{100 * (all_cycle_types.get('indirect', 0) / total_by_type):.1f}"
        # Complex = hypercycles
        values['complex_loops_pct'] = f"{100 * (all_cycle_types.get('hypercycle', 0) / total_by_type):.1f}"
        values['direct_count'] = f"{all_cycle_types.get('direct', 0):,}"
    else:
        values['simple_loops_pct'] = "0.0"
        values['medium_loops_pct'] = "0.0"
        values['complex_loops_pct'] = "0.0"
        values['direct_count'] = "0"
    
    # Amplification factors (aggregate across all scenarios)
    all_amplifications = []
    for scenario_data in [miller_data, hydro_data, form_data]:
        amp_stats = scenario_data.get('amplification_stats', {})
        # We'll use min/max/median from aggregated stats
        if amp_stats:
            all_amplifications.append(amp_stats)
    
    # Aggregate amplification stats
    if all_amplifications:
        all_mins = [s.get('min', 0) for s in all_amplifications if s.get('min')]
        all_maxs = [s.get('max', 0) for s in all_amplifications if s.get('max')]
        all_medians = [s.get('median', 0) for s in all_amplifications if s.get('median')]
        all_q25s = [s.get('quartiles', [0, 0, 0])[0] for s in all_amplifications if s.get('quartiles')]
        all_q75s = [s.get('quartiles', [0, 0, 0])[2] for s in all_amplifications if s.get('quartiles')]
        
        if all_mins:
            values['amp_min'] = f"{min(all_mins):.2f}"
        else:
            values['amp_min'] = "1.0"
        
        if all_maxs:
            values['amp_max'] = f"{max(all_maxs):.0f}"
        else:
            values['amp_max'] = "1"
        
        if all_medians:
            values['amp_median'] = f"{sum(all_medians) / len(all_medians):.2f}"
        else:
            values['amp_median'] = "1.0"
        
        if all_q25s and all_q75s:
            values['amp_iqr'] = f"{sum(all_q25s) / len(all_q25s):.2f}-{sum(all_q75s) / len(all_q75s):.2f}"
        else:
            values['amp_iqr'] = "1.0-1.0"
    else:
        values['amp_min'] = "1.0"
        values['amp_max'] = "1"
        values['amp_median'] = "1.0"
        values['amp_iqr'] = "1.0-1.0"
    
    # Placeholders for strongest amplifiers (to be filled manually or from detailed analysis)
    values['strongest_molecules'] = "[molecule names]"
    values['strongest_scenario'] = "[scenario]"
    values['max_amplification'] = f"{values['amp_max']}"
    
    # Formose-like cycles (placeholder - would need specific detection)
    values['formose_runs'] = "[N]"
    values['glycolaldehyde_amp'] = "[XX]"
    values['glycolaldehyde_time'] = "[time]"
    
    return values


def fill_section_3_3(template: str, values: Dict) -> str:
    """Fill Section 3.3 placeholders"""
    
    # Map placeholders to values
    replacements = {
        r'\[XX\]\s*unique autocatalytic cycles': f"{values['total_unique_cycles']} unique autocatalytic cycles",
        r'\[X ± Y\]\s*cycles/run\)': f"{values['miller_cycles']} cycles/run)",
        r'\[A ± B\]': values['hydro_cycles'],
        r'\[C ± D\]': values['form_cycles'],
        r'\[XX\]%': f"{values['simple_loops_pct']}%",
        r'\[YY\]%': f"{values['medium_loops_pct']}%",
        r'\[ZZ\]%': f"{values['complex_loops_pct']}%",
        r'\[N\]\s*instances': f"{values['direct_count']} instances",
        r'\[min\]': values['amp_min'],
        r'\[max\]': values['amp_max'],
        r'\[X\.X\]\s*\(IQR:': f"{values['amp_median']} (IQR:",
        r'\[Y\.Y\]-\[Z\.Z\]\)': f"{values['amp_iqr']})",
        r'\[molecule names\]': values['strongest_molecules'],
        r'\[scenario\]': values['strongest_scenario'],
        r'\[XX\]-fold amplification': f"{values['max_amplification']}-fold amplification",
        r'\[N\]\s*formamide runs': f"{values['formose_runs']} formamide runs",
        r'\[XX\]-fold over \[time\]': f"{values['glycolaldehyde_amp']}-fold over {values['glycolaldehyde_time']}",
    }
    
    result = template
    for pattern, replacement in replacements.items():
        result = re.sub(pattern, replacement, result)
    
    return result


def main():
    parser = argparse.ArgumentParser(description="Fill Section 3.3 with autocatalysis data")
    parser.add_argument('--data', required=True, help='Data directory (paper/results_data)')
    parser.add_argument('--template', required=True, help='Template file (paper/manuscript_draft.tex)')
    parser.add_argument('--output', help='Output file (default: template with _filled suffix)')
    
    args = parser.parse_args()
    
    data_dir = Path(args.data)
    template_file = Path(args.template)
    
    if not template_file.exists():
        print(f"Error: Template file not found: {template_file}")
        return 1
    
    # Load data
    print("Loading analysis data...")
    data = load_analysis_data(data_dir)
    
    # Extract values
    print("Extracting autocatalysis values...")
    values = extract_autocatalysis_values(data)
    
    print("\nExtracted values:")
    for key, val in values.items():
        print(f"  {key}: {val}")
    
    # Load template
    print(f"\nLoading template: {template_file}")
    with open(template_file) as f:
        template = f.read()
    
    # Fill Section 3.3
    print("Filling Section 3.3...")
    filled = fill_section_3_3(template, values)
    
    # Save output
    if args.output:
        output_file = Path(args.output)
    else:
        output_file = template_file.parent / f"{template_file.stem}_filled.tex"
    
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(filled)
    
    print(f"\n✅ Filled manuscript saved to: {output_file}")
    print("\n⚠️  Note: Some placeholders may need manual review:")
    print("   - [molecule names] (strongest amplifiers)")
    print("   - [scenario] (where strongest found)")
    print("   - [N] formamide runs (formose-like cycles)")
    print("   - [XX]-fold over [time] (glycolaldehyde amplification)")
    
    return 0


if __name__ == "__main__":
    exit(main())

