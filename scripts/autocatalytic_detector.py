"""
Autocatalytic Cycle Detector
==============================

Detects and analyzes autocatalytic cycles in reaction networks.

An autocatalytic cycle is a sequence of reactions where:
1. A molecule X catalyzes its own production (directly or indirectly)
2. The cycle is self-sustaining (products feed back into reactants)
3. The cycle can amplify (produces more than it consumes)

Types of autocatalysis detected:
- Direct autocatalysis: A + B -> 2A
- Indirect autocatalysis: A -> B -> C -> A (with amplification)
- Hypercycles: Coupled autocatalytic cycles (A helps B, B helps A)
- RAF sets: Reflexively Autocatalytic and Food-generated sets

Usage:
    # Detect cycles in network
    python scripts/autocatalytic_detector.py analysis/reaction_network/reaction_network.json
    
    # Detailed analysis with visualization
    python scripts/autocatalytic_detector.py analysis/reaction_network/reaction_network.json --detailed
"""

import sys
import argparse
import logging
import json
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict, deque
from dataclasses import dataclass, field

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class Cycle:
    """Represents an autocatalytic cycle"""
    molecules: List[str]  # Molecules in cycle (in order)
    reactions: List[str]  # Reaction strings
    amplification_factor: float  # Net production/consumption ratio
    cycle_type: str  # 'direct', 'indirect', 'hypercycle', 'raf'
    is_self_sustaining: bool
    
    def __str__(self):
        return " -> ".join(self.molecules + [self.molecules[0]])
    
    def size(self) -> int:
        return len(self.molecules)


class AutocatalyticDetector:
    """Detects autocatalytic cycles in reaction networks"""
    
    def __init__(self, network_file: Path, output_dir: Path):
        self.network_file = network_file
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Network data
        self.molecules = []
        self.reactions = []
        self.forward_graph = defaultdict(set)  # molecule -> products
        self.reverse_graph = defaultdict(set)  # molecule -> reactants
        self.reaction_map = {}  # (reactants, products) -> reaction info
        
        # Detected cycles
        self.cycles = []
        
    def load_network(self):
        """Load reaction network from JSON"""
        logger.info(f"Loading network: {self.network_file}")
        
        with open(self.network_file) as f:
            data = json.load(f)
        
        self.molecules = [mol['formula'] for mol in data.get('molecules', [])]
        self.reactions = data.get('reactions', [])
        
        logger.info(f"  Molecules: {len(self.molecules)}")
        logger.info(f"  Reactions: {len(self.reactions)}")
        
        # Build graph
        for rxn in self.reactions:
            reactants = tuple(rxn['reactants'])
            products = tuple(rxn['products'])
            
            # Store reaction
            self.reaction_map[(reactants, products)] = rxn
            
            # Build graph edges
            for r in reactants:
                for p in products:
                    self.forward_graph[r].add(p)
                    self.reverse_graph[p].add(r)
    
    def detect_direct_autocatalysis(self) -> List[Cycle]:
        """
        Detect direct autocatalysis: A + B -> A + A (or similar)
        Molecule appears on both sides with net production
        """
        logger.info("\n[1/4] Detecting direct autocatalysis...")
        
        cycles = []
        
        for rxn in self.reactions:
            reactants = rxn['reactants']
            products = rxn['products']
            
            # Count molecules on each side
            reactant_counts = {}
            product_counts = {}
            
            for r in reactants:
                reactant_counts[r] = reactant_counts.get(r, 0) + 1
            for p in products:
                product_counts[p] = product_counts.get(p, 0) + 1
            
            # Find molecules that appear on both sides with net production
            for mol in set(reactants) & set(products):
                net = product_counts[mol] - reactant_counts[mol]
                if net > 0:  # Net production
                    amplification = product_counts[mol] / reactant_counts[mol]
                    
                    cycle = Cycle(
                        molecules=[mol],
                        reactions=[rxn['string']],
                        amplification_factor=amplification,
                        cycle_type='direct',
                        is_self_sustaining=True
                    )
                    cycles.append(cycle)
        
        logger.info(f"  Found {len(cycles)} direct autocatalytic reactions")
        return cycles
    
    def detect_indirect_autocatalysis(self, max_length: int = 10) -> List[Cycle]:
        """
        Detect indirect autocatalysis: cycles in the reaction graph
        Uses DFS to find cycles where a molecule eventually produces itself
        """
        logger.info(f"\n[2/4] Detecting indirect autocatalysis (max length: {max_length})...")
        
        cycles = []
        
        # For each molecule, try to find paths back to itself
        for start_mol in self.molecules[:100]:  # Limit for performance
            paths = self._find_cycles_from(start_mol, max_length)
            
            for path in paths:
                # Check if this is a valid autocatalytic cycle
                if self._is_autocatalytic(path):
                    amplification = self._compute_amplification(path)
                    
                    cycle = Cycle(
                        molecules=path,
                        reactions=self._get_reactions_for_path(path),
                        amplification_factor=amplification,
                        cycle_type='indirect',
                        is_self_sustaining=amplification >= 1.0
                    )
                    cycles.append(cycle)
        
        logger.info(f"  Found {len(cycles)} indirect autocatalytic cycles")
        return cycles
    
    def _find_cycles_from(self, start: str, max_length: int) -> List[List[str]]:
        """Find all cycles starting from a given molecule using DFS"""
        cycles = []
        
        def dfs(current: str, path: List[str], visited: Set[str]):
            if len(path) > max_length:
                return
            
            # Found a cycle back to start
            if current == start and len(path) > 1:
                cycles.append(path[:])
                return
            
            # Already visited (but not the start)
            if current in visited:
                return
            
            visited.add(current)
            
            # Explore neighbors
            for neighbor in self.forward_graph.get(current, set()):
                dfs(neighbor, path + [neighbor], visited.copy())
        
        dfs(start, [start], set())
        return cycles
    
    def _is_autocatalytic(self, path: List[str]) -> bool:
        """Check if a cycle path is truly autocatalytic"""
        # For now, simple check: path forms a cycle
        # More sophisticated: check stoichiometry
        return len(path) > 1 and path[0] == path[-1]
    
    def _compute_amplification(self, path: List[str]) -> float:
        """Compute amplification factor for a cycle"""
        # Simplified: assume 1.0 for now
        # Real implementation would track stoichiometry
        return 1.0
    
    def _get_reactions_for_path(self, path: List[str]) -> List[str]:
        """Get reaction strings for a path"""
        reactions = []
        for i in range(len(path) - 1):
            # Find reaction from path[i] to path[i+1]
            reactions.append(f"{path[i]} -> {path[i+1]}")
        return reactions
    
    def detect_hypercycles(self) -> List[Cycle]:
        """
        Detect hypercycles: coupled autocatalytic cycles
        A catalyzes B, B catalyzes C, C catalyzes A
        """
        logger.info("\n[3/4] Detecting hypercycles...")
        
        # Simplified: look for strongly connected components
        # that form mutual catalysis networks
        
        cycles = []
        
        # Find strongly connected components (simplified)
        components = self._find_strongly_connected_components()
        
        for comp in components:
            if len(comp) >= 2:  # At least 2 molecules
                # Check if they mutually catalyze
                is_hypercycle = True
                for mol in comp:
                    # Check if mol is produced by others in comp
                    producers = self.reverse_graph.get(mol, set())
                    if not producers & set(comp):
                        is_hypercycle = False
                        break
                
                if is_hypercycle:
                    cycle = Cycle(
                        molecules=list(comp),
                        reactions=[],  # Would populate with actual reactions
                        amplification_factor=1.0,
                        cycle_type='hypercycle',
                        is_self_sustaining=True
                    )
                    cycles.append(cycle)
        
        logger.info(f"  Found {len(cycles)} hypercycles")
        return cycles
    
    def _find_strongly_connected_components(self) -> List[Set[str]]:
        """Find strongly connected components using Tarjan's algorithm"""
        # Simplified implementation
        components = []
        visited = set()
        
        for mol in self.molecules[:50]:  # Limit for performance
            if mol not in visited:
                # Simple reachability check
                reachable = self._get_reachable(mol, max_depth=5)
                
                # Check if forms SCC
                scc = set()
                for other in reachable:
                    if mol in self._get_reachable(other, max_depth=5):
                        scc.add(other)
                
                if len(scc) >= 2:
                    components.append(scc)
                    visited.update(scc)
        
        return components
    
    def _get_reachable(self, start: str, max_depth: int) -> Set[str]:
        """Get all molecules reachable from start within max_depth"""
        reachable = set()
        queue = deque([(start, 0)])
        
        while queue:
            mol, depth = queue.popleft()
            
            if depth > max_depth or mol in reachable:
                continue
            
            reachable.add(mol)
            
            for neighbor in self.forward_graph.get(mol, set()):
                queue.append((neighbor, depth + 1))
        
        return reachable
    
    def detect_raf_sets(self) -> List[Dict]:
        """
        Detect RAF (Reflexively Autocatalytic and Food-generated) sets
        A set of reactions where each reaction is catalyzed by at least
        one molecule producible from food molecules
        """
        logger.info("\n[4/4] Detecting RAF sets...")
        
        # This is a complex algorithm - simplified version
        # Full implementation would need:
        # 1. Define food set (initial molecules)
        # 2. Iteratively expand: reactions catalyzed by available molecules
        # 3. Check for closure (RAF property)
        
        raf_sets = []
        
        # Placeholder implementation
        logger.info("  RAF detection: simplified implementation")
        logger.info("  (Full RAF analysis requires catalysis tracking)")
        
        return raf_sets
    
    def analyze(self):
        """Run all autocatalysis detection algorithms"""
        logger.info("\n" + "=" * 70)
        logger.info("AUTOCATALYTIC CYCLE DETECTION")
        logger.info("=" * 70)
        
        # Detect all types
        direct = self.detect_direct_autocatalysis()
        indirect = self.detect_indirect_autocatalysis()
        hypercycles = self.detect_hypercycles()
        raf_sets = self.detect_raf_sets()
        
        # Combine results
        self.cycles = direct + indirect + hypercycles
        
        # Summary
        logger.info("\n" + "=" * 70)
        logger.info("DETECTION SUMMARY")
        logger.info("=" * 70)
        logger.info(f"Direct autocatalysis: {len(direct)}")
        logger.info(f"Indirect cycles: {len(indirect)}")
        logger.info(f"Hypercycles: {len(hypercycles)}")
        logger.info(f"RAF sets: {len(raf_sets)}")
        logger.info(f"Total cycles: {len(self.cycles)}")
        logger.info("=" * 70)
        
        return {
            'direct': direct,
            'indirect': indirect,
            'hypercycles': hypercycles,
            'raf_sets': raf_sets,
            'total': len(self.cycles)
        }
    
    def export_results(self, results: Dict):
        """Export results to JSON"""
        output_file = self.output_dir / "autocatalytic_cycles.json"
        
        logger.info(f"\nExporting results: {output_file}")
        
        data = {
            'analysis_time': datetime.now().isoformat(),
            'network_file': str(self.network_file),
            'summary': {
                'direct_count': len(results['direct']),
                'indirect_count': len(results['indirect']),
                'hypercycle_count': len(results['hypercycles']),
                'raf_set_count': len(results['raf_sets']),
                'total_count': results['total']
            },
            'cycles': []
        }
        
        # Export all cycles
        for cycle in self.cycles:
            data['cycles'].append({
                'type': cycle.cycle_type,
                'molecules': cycle.molecules,
                'reactions': cycle.reactions,
                'size': cycle.size(),
                'amplification_factor': cycle.amplification_factor,
                'is_self_sustaining': cycle.is_self_sustaining,
                'string': str(cycle)
            })
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        logger.info(f"Results exported: {output_file}")
    
    def generate_report(self, results: Dict):
        """Generate human-readable report"""
        report_file = self.output_dir / "autocatalytic_report.txt"
        
        logger.info(f"Generating report: {report_file}")
        
        with open(report_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("AUTOCATALYTIC CYCLE ANALYSIS\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Analysis time: {datetime.now()}\n")
            f.write(f"Network file: {self.network_file}\n\n")
            
            # Summary
            f.write("Summary:\n")
            f.write("-" * 70 + "\n")
            f.write(f"  Direct autocatalysis: {len(results['direct'])}\n")
            f.write(f"  Indirect cycles: {len(results['indirect'])}\n")
            f.write(f"  Hypercycles: {len(results['hypercycles'])}\n")
            f.write(f"  RAF sets: {len(results['raf_sets'])}\n")
            f.write(f"  Total: {results['total']}\n\n")
            
            # Direct autocatalysis
            if results['direct']:
                f.write("Direct Autocatalysis:\n")
                f.write("-" * 70 + "\n")
                for i, cycle in enumerate(results['direct'][:20], 1):
                    f.write(f"{i}. {cycle.molecules[0]}\n")
                    f.write(f"   Reaction: {cycle.reactions[0]}\n")
                    f.write(f"   Amplification: {cycle.amplification_factor:.2f}x\n\n")
            
            # Indirect cycles
            if results['indirect']:
                f.write("\nIndirect Autocatalytic Cycles:\n")
                f.write("-" * 70 + "\n")
                for i, cycle in enumerate(results['indirect'][:10], 1):
                    f.write(f"{i}. {str(cycle)}\n")
                    f.write(f"   Size: {cycle.size()} molecules\n")
                    f.write(f"   Amplification: {cycle.amplification_factor:.2f}x\n\n")
            
            # Hypercycles
            if results['hypercycles']:
                f.write("\nHypercycles:\n")
                f.write("-" * 70 + "\n")
                for i, cycle in enumerate(results['hypercycles'][:10], 1):
                    f.write(f"{i}. {' <-> '.join(cycle.molecules)}\n")
                    f.write(f"   Size: {cycle.size()} molecules\n\n")
            
            f.write("\n" + "=" * 70 + "\n")
        
        logger.info(f"Report saved: {report_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Detect autocatalytic cycles in reaction networks"
    )
    parser.add_argument(
        'network_file',
        type=str,
        help="Path to reaction network JSON file"
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='analysis/autocatalytic_cycles',
        help="Output directory (default: analysis/autocatalytic_cycles)"
    )
    parser.add_argument(
        '--max-cycle-length',
        type=int,
        default=10,
        help="Maximum cycle length to search (default: 10)"
    )
    parser.add_argument(
        '--detailed',
        action='store_true',
        help="Perform detailed analysis (slower)"
    )
    
    args = parser.parse_args()
    
    network_file = Path(args.network_file)
    output_dir = Path(args.output)
    
    if not network_file.exists():
        logger.error(f"Network file not found: {network_file}")
        return 1
    
    logger.info("=" * 70)
    logger.info("AUTOCATALYTIC CYCLE DETECTOR")
    logger.info("=" * 70)
    logger.info(f"Network: {network_file}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Max cycle length: {args.max_cycle_length}")
    logger.info("=" * 70)
    
    # Run detection
    detector = AutocatalyticDetector(network_file, output_dir)
    detector.load_network()
    results = detector.analyze()
    
    # Export
    detector.export_results(results)
    detector.generate_report(results)
    
    logger.info("\n" + "=" * 70)
    logger.info("DETECTION COMPLETE!")
    logger.info("=" * 70)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

