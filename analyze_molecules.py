#!/usr/bin/env python3
"""
Simple molecule analysis script
"""
import json

def analyze_molecules():
    # Load results
    with open('results/overnight_test_2025-10-14_16-59-45/results.json', 'r') as f:
        data = json.load(f)
    
    molecules = data.get('molecules_detected', [])
    
    print("=== MOLECULES DETECTED ===")
    print(f"Total molecules: {len(molecules)}")
    print()
    
    for i, m in enumerate(molecules):
        print(f"{i+1}. ID: {m['id']}")
        print(f"   Formula: {m['formula']}")
        print(f"   Mass: {m['mass']:.3f}")
        print(f"   Complexity: {m['complexity']:.1f}")
        print(f"   Count: {m['count']}")
        print(f"   First detected: {m['first_detected']:.1f}")
        print(f"   Last seen: {m['last_seen']:.1f}")
        print()
    
    # Statistics
    total_instances = sum(m['count'] for m in molecules)
    avg_mass = sum(m['mass'] for m in molecules) / len(molecules)
    avg_complexity = sum(m['complexity'] for m in molecules) / len(molecules)
    
    print("=== STATISTICS ===")
    print(f"Total instances: {total_instances:,}")
    print(f"Average mass: {avg_mass:.3f}")
    print(f"Average complexity: {avg_complexity:.1f}")
    print()
    
    # Most common
    sorted_molecules = sorted(molecules, key=lambda x: x['count'], reverse=True)
    print("=== TOP 3 MOST COMMON ===")
    for i, m in enumerate(sorted_molecules[:3]):
        print(f"{i+1}. {m['formula'][:20]}... (count: {m['count']:,})")

if __name__ == "__main__":
    analyze_molecules()
