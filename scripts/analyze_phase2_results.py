"""
Phase 2 Results Analyzer
========================

Analyzes simulation results from Phase 2:
- Extracts novel molecules
- Uses MatcherV2 for PubChem identification
- Generates molecule catalog
- Identifies reaction networks
- Finds autocatalytic cycles

Usage:
    python scripts/analyze_phase2_results.py --scenario miller_urey
    python scripts/analyze_phase2_results.py --all
"""

import sys
import json
import argparse
from pathlib import Path
from typing import List, Dict, Set
from collections import Counter, defaultdict
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from matcher.matcher_v2 import MatcherV2
from matcher.confidence import Reliability


class Phase2ResultsAnalyzer:
    """Analyzes Phase 2 simulation results"""
    
    def __init__(self, output_dir: str = "results/phase2"):
        self.output_dir = Path(output_dir)
        
        # Initialize MatcherV2
        print("[+] Initializing PubChem MatcherV2...")
        self.matcher = MatcherV2(
            classifier_model='data/atom_classifier.pkl',
            use_ml_classifier=True  # Use ML for best results
        )
        print("    ✓ Matcher ready")
        
        self.molecule_catalog = []
        self.reaction_networks = {}
    
    def analyze_scenario(self, scenario: str) -> Dict:
        """
        Analyze all runs for a scenario
        
        Returns:
            analysis_results dict
        """
        scenario_dir = self.output_dir / scenario
        
        if not scenario_dir.exists():
            print(f"❌ Scenario directory not found: {scenario_dir}")
            return {}
        
        print("\n" + "=" * 70)
        print(f"ANALYZING: {scenario.upper()}")
        print("=" * 70)
        
        # Find all run directories
        run_dirs = sorted([d for d in scenario_dir.iterdir() if d.is_dir()])
        print(f"Found {len(run_dirs)} runs")
        
        # Collect molecules from all runs
        all_molecules = []
        all_reactions = []
        
        for i, run_dir in enumerate(run_dirs):
            print(f"\n[{i+1}/{len(run_dirs)}] Processing {run_dir.name}...")
            
            # Load molecules
            molecules_file = run_dir / "molecules.json"
            if molecules_file.exists():
                with open(molecules_file, 'r') as f:
                    molecules = json.load(f)
                all_molecules.extend(molecules)
                print(f"    Found {len(molecules)} molecules")
            
            # Load reactions
            reactions_file = run_dir / "reactions.json"
            if reactions_file.exists():
                with open(reactions_file, 'r') as f:
                    reactions = json.load(f)
                all_reactions.extend(reactions)
                print(f"    Found {len(reactions)} reactions")
        
        print(f"\n{'='*70}")
        print(f"Total molecules collected: {len(all_molecules)}")
        print(f"Total reactions collected: {len(all_reactions)}")
        
        # Filter unique molecules
        unique_molecules = self._filter_unique_molecules(all_molecules)
        print(f"Unique molecules: {len(unique_molecules)}")
        
        # Match molecules to PubChem
        print(f"\n[+] Matching molecules to PubChem...")
        matched_molecules = self._match_molecules_to_pubchem(unique_molecules)
        
        # Analyze reactions
        print(f"\n[+] Analyzing reaction networks...")
        reaction_analysis = self._analyze_reactions(all_reactions)
        
        # Identify autocatalytic cycles
        print(f"\n[+] Searching for autocatalytic cycles...")
        autocatalytic_cycles = self._find_autocatalytic_cycles(all_reactions)
        
        # Generate statistics
        stats = self._generate_statistics(
            unique_molecules,
            matched_molecules,
            reaction_analysis,
            autocatalytic_cycles
        )
        
        # Save results
        results = {
            'scenario': scenario,
            'timestamp': datetime.now().isoformat(),
            'num_runs': len(run_dirs),
            'total_molecules': len(all_molecules),
            'unique_molecules': len(unique_molecules),
            'matched_molecules': matched_molecules,
            'reactions': reaction_analysis,
            'autocatalytic_cycles': autocatalytic_cycles,
            'statistics': stats
        }
        
        # Save to file
        output_file = self.output_dir / f"{scenario}_analysis.json"
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\n✅ Analysis saved to: {output_file}")
        
        return results
    
    def _filter_unique_molecules(self, molecules: List[Dict]) -> List[Dict]:
        """Filter unique molecules by formula"""
        unique = {}
        
        for mol in molecules:
            formula = mol.get('formula', '')
            if formula not in unique:
                unique[formula] = mol
            else:
                # Keep the one with more information
                if len(mol.get('bonds', [])) > len(unique[formula].get('bonds', [])):
                    unique[formula] = mol
        
        return list(unique.values())
    
    def _match_molecules_to_pubchem(self, molecules: List[Dict]) -> List[Dict]:
        """Match molecules to PubChem using MatcherV2"""
        matched = []
        
        print(f"  Matching {len(molecules)} unique molecules...")
        
        for i, mol in enumerate(molecules):
            if i % 10 == 0:
                print(f"    Progress: {i}/{len(molecules)}")
            
            try:
                result = self.matcher.match_cluster(
                    cluster=mol,
                    top_n=3,
                    min_similarity=0.5
                )
                
                if result.success:
                    matched_result = {
                        'formula': mol.get('formula', ''),
                        'pubchem_cid': result.pubchem_cid,
                        'pubchem_name': result.pubchem_name,
                        'pubchem_formula': result.pubchem_formula,
                        'pubchem_smiles': result.pubchem_smiles,
                        'similarity': result.similarity_score.overall,
                        'confidence': result.confidence.confidence_score,
                        'reliability': result.confidence.reliability.value,
                        'is_novel': result.confidence.reliability == Reliability.LOW
                    }
                    matched.append(matched_result)
            
            except Exception as e:
                print(f"    Warning: Failed to match {mol.get('formula', 'unknown')}: {e}")
        
        print(f"  ✓ Successfully matched: {len(matched)}/{len(molecules)}")
        
        # Sort by reliability
        matched.sort(key=lambda x: (
            {'high': 3, 'medium': 2, 'low': 1, 'invalid': 0}.get(x['reliability'], 0),
            x['similarity']
        ), reverse=True)
        
        return matched
    
    def _analyze_reactions(self, reactions: List[Dict]) -> Dict:
        """Analyze reaction networks"""
        if not reactions:
            return {'count': 0, 'reaction_types': {}, 'most_common': []}
        
        # Count reaction types
        reaction_types = Counter()
        reactant_product_pairs = []
        
        for reaction in reactions:
            rtype = reaction.get('type', 'unknown')
            reaction_types[rtype] += 1
            
            reactants = reaction.get('reactants', [])
            products = reaction.get('products', [])
            if reactants and products:
                reactant_product_pairs.append((tuple(reactants), tuple(products)))
        
        # Most common reactions
        most_common = Counter(reactant_product_pairs).most_common(10)
        
        return {
            'count': len(reactions),
            'reaction_types': dict(reaction_types),
            'most_common_reactions': [
                {'reactants': list(r), 'products': list(p), 'count': c}
                for (r, p), c in most_common
            ]
        }
    
    def _find_autocatalytic_cycles(self, reactions: List[Dict]) -> List[Dict]:
        """Find autocatalytic cycles in reaction network"""
        if not reactions:
            return []
        
        # Build reaction graph
        graph = defaultdict(set)  # product -> reactants that produce it
        reverse_graph = defaultdict(set)  # reactant -> products it makes
        
        for reaction in reactions:
            reactants = frozenset(reaction.get('reactants', []))
            products = frozenset(reaction.get('products', []))
            
            for product in products:
                graph[product].add(reactants)
                for reactant in reactants:
                    reverse_graph[reactant].add(product)
        
        # Search for cycles (simplified DFS)
        cycles = []
        visited = set()
        
        def dfs_cycle(node, path):
            if node in path:
                # Found cycle
                cycle_start = path.index(node)
                cycle = path[cycle_start:]
                if len(cycle) >= 2:
                    cycles.append(cycle)
                return
            
            if node in visited:
                return
            
            visited.add(node)
            path.append(node)
            
            # Follow products
            for product in reverse_graph.get(node, []):
                dfs_cycle(product, path.copy())
        
        # Start DFS from each molecule
        for molecule in list(graph.keys())[:100]:  # Limit search
            dfs_cycle(molecule, [])
        
        # Remove duplicates
        unique_cycles = []
        seen = set()
        for cycle in cycles:
            cycle_key = tuple(sorted(cycle))
            if cycle_key not in seen:
                seen.add(cycle_key)
                unique_cycles.append({
                    'molecules': cycle,
                    'length': len(cycle)
                })
        
        print(f"    Found {len(unique_cycles)} potential autocatalytic cycles")
        
        return unique_cycles[:20]  # Return top 20
    
    def _generate_statistics(self,
                            unique_molecules: List[Dict],
                            matched_molecules: List[Dict],
                            reaction_analysis: Dict,
                            autocatalytic_cycles: List[Dict]) -> Dict:
        """Generate summary statistics"""
        # Formula distribution
        formulas = [mol.get('formula', 'unknown') for mol in unique_molecules]
        formula_counts = Counter(formulas).most_common(20)
        
        # Atom counts
        atom_counts = Counter()
        for mol in unique_molecules:
            atoms = mol.get('atoms', [])
            for atom in atoms:
                if isinstance(atom, str):
                    atom_counts[atom] += 1
                elif isinstance(atom, dict):
                    atom_counts[atom.get('element', 'unknown')] += 1
        
        # Size distribution
        sizes = [len(mol.get('atoms', [])) for mol in unique_molecules]
        size_dist = Counter(sizes)
        
        # Match quality
        high_conf = sum(1 for m in matched_molecules if m['reliability'] == 'high')
        medium_conf = sum(1 for m in matched_molecules if m['reliability'] == 'medium')
        low_conf = sum(1 for m in matched_molecules if m['reliability'] == 'low')
        
        return {
            'total_unique_molecules': len(unique_molecules),
            'matched_molecules': len(matched_molecules),
            'match_quality': {
                'high_confidence': high_conf,
                'medium_confidence': medium_conf,
                'low_confidence': low_conf
            },
            'top_20_formulas': [{'formula': f, 'count': c} for f, c in formula_counts],
            'atom_distribution': dict(atom_counts),
            'size_distribution': dict(size_dist),
            'total_reactions': reaction_analysis.get('count', 0),
            'autocatalytic_cycles_found': len(autocatalytic_cycles)
        }
    
    def analyze_all_scenarios(self) -> Dict:
        """Analyze all scenarios"""
        scenarios = ['miller_urey', 'hydrothermal', 'formamide']
        
        all_results = {}
        
        for scenario in scenarios:
            results = self.analyze_scenario(scenario)
            all_results[scenario] = results
        
        # Generate comparison
        print("\n\n" + "=" * 70)
        print("PHASE 2 ANALYSIS - COMPARISON")
        print("=" * 70)
        
        for scenario, results in all_results.items():
            stats = results.get('statistics', {})
            print(f"\n{scenario.upper()}:")
            print(f"  Unique molecules: {stats.get('total_unique_molecules', 0)}")
            print(f"  Matched to PubChem: {stats.get('matched_molecules', 0)}")
            print(f"  High confidence: {stats.get('match_quality', {}).get('high_confidence', 0)}")
            print(f"  Reactions: {stats.get('total_reactions', 0)}")
            print(f"  Autocatalytic cycles: {stats.get('autocatalytic_cycles_found', 0)}")
        
        print("\n" + "=" * 70)
        
        # Save comparison
        comparison_file = self.output_dir / "phase2_comparison.json"
        with open(comparison_file, 'w') as f:
            json.dump(all_results, f, indent=2)
        
        print(f"\n✅ Comparison saved to: {comparison_file}")
        
        return all_results


def main():
    parser = argparse.ArgumentParser(description="Phase 2 Results Analyzer")
    parser.add_argument('--scenario', choices=['miller_urey', 'hydrothermal', 'formamide'],
                       help='Scenario to analyze')
    parser.add_argument('--all', action='store_true',
                       help='Analyze all scenarios')
    parser.add_argument('--output', default='results/phase2',
                       help='Results directory (default: results/phase2)')
    
    args = parser.parse_args()
    
    if not args.scenario and not args.all:
        parser.error("Must specify --scenario or --all")
    
    # Initialize analyzer
    analyzer = Phase2ResultsAnalyzer(output_dir=args.output)
    
    # Analyze
    if args.all:
        analyzer.analyze_all_scenarios()
    else:
        analyzer.analyze_scenario(args.scenario)
    
    print("\n✅ Analysis complete!")


if __name__ == "__main__":
    main()

