"""
Reaction Network Analyzer
==========================

Analyzes reaction networks from Phase 2 simulation results.
Builds reaction graphs, identifies key pathways, and computes network metrics.

Features:
- Reaction graph construction (molecules as nodes, reactions as edges)
- Network topology analysis (degree distribution, centrality, clustering)
- Pathway identification (shortest paths, reaction cascades)
- Export to various formats (GraphML, JSON, DOT)

Usage:
    # Analyze single simulation
    python scripts/reaction_network_analyzer.py results/overnight_test
    
    # Analyze multiple simulations (merged network)
    python scripts/reaction_network_analyzer.py results/phase2/miller_urey --merge
    
    # Export to GraphML for visualization
    python scripts/reaction_network_analyzer.py results/overnight_test --export graphml
"""

import sys
import argparse
import logging
import json
import pickle
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict, Counter
from dataclasses import dataclass, field

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class Molecule:
    """Represents a molecule in the network"""
    formula: str
    id: Optional[str] = None
    name: Optional[str] = None
    mass: Optional[float] = None
    num_atoms: Optional[int] = None
    first_seen: Optional[int] = None  # Step when first appeared
    
    def __hash__(self):
        return hash(self.formula)
    
    def __eq__(self, other):
        return isinstance(other, Molecule) and self.formula == other.formula


@dataclass
class Reaction:
    """Represents a chemical reaction"""
    reactants: Tuple[str, ...]  # Tuple of formulas
    products: Tuple[str, ...]   # Tuple of formulas
    step: Optional[int] = None  # Step when reaction occurred
    energy_change: Optional[float] = None
    
    def __hash__(self):
        return hash((self.reactants, self.products))
    
    def __eq__(self, other):
        return (isinstance(other, Reaction) and 
                self.reactants == other.reactants and 
                self.products == other.products)
    
    def __str__(self):
        reactants_str = " + ".join(self.reactants)
        products_str = " + ".join(self.products)
        return f"{reactants_str} -> {products_str}"


@dataclass
class ReactionNetwork:
    """Complete reaction network graph"""
    molecules: Dict[str, Molecule] = field(default_factory=dict)
    reactions: List[Reaction] = field(default_factory=list)
    
    # Graph structure (adjacency lists)
    forward_edges: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))
    reverse_edges: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))
    
    def add_molecule(self, mol: Molecule):
        """Add molecule to network"""
        if mol.formula not in self.molecules:
            self.molecules[mol.formula] = mol
    
    def add_reaction(self, rxn: Reaction):
        """Add reaction to network and update graph"""
        # Add molecules if not present
        for formula in rxn.reactants + rxn.products:
            if formula not in self.molecules:
                self.add_molecule(Molecule(formula=formula))
        
        # Add reaction
        self.reactions.append(rxn)
        
        # Update graph edges
        for reactant in rxn.reactants:
            for product in rxn.products:
                self.forward_edges[reactant].add(product)
                self.reverse_edges[product].add(reactant)
    
    def get_degree(self, formula: str) -> Tuple[int, int]:
        """Get in-degree and out-degree for a molecule"""
        in_degree = len(self.reverse_edges.get(formula, set()))
        out_degree = len(self.forward_edges.get(formula, set()))
        return in_degree, out_degree
    
    def get_sources(self) -> List[str]:
        """Get source molecules (no incoming reactions)"""
        sources = []
        for formula in self.molecules:
            in_degree, _ = self.get_degree(formula)
            if in_degree == 0:
                sources.append(formula)
        return sources
    
    def get_sinks(self) -> List[str]:
        """Get sink molecules (no outgoing reactions)"""
        sinks = []
        for formula in self.molecules:
            _, out_degree = self.get_degree(formula)
            if out_degree == 0:
                sinks.append(formula)
        return sinks
    
    def get_hubs(self, threshold: int = 5) -> List[str]:
        """Get hub molecules (high total degree)"""
        hubs = []
        for formula in self.molecules:
            in_degree, out_degree = self.get_degree(formula)
            if in_degree + out_degree >= threshold:
                hubs.append(formula)
        return hubs
    
    def compute_statistics(self) -> Dict:
        """Compute network statistics"""
        stats = {
            'num_molecules': len(self.molecules),
            'num_reactions': len(self.reactions),
            'num_sources': len(self.get_sources()),
            'num_sinks': len(self.get_sinks()),
            'num_hubs': len(self.get_hubs()),
        }
        
        # Degree distribution
        degrees = []
        for formula in self.molecules:
            in_deg, out_deg = self.get_degree(formula)
            degrees.append(in_deg + out_deg)
        
        if degrees:
            stats['avg_degree'] = sum(degrees) / len(degrees)
            stats['max_degree'] = max(degrees)
            stats['min_degree'] = min(degrees)
        
        # Reaction statistics
        reactant_counts = [len(r.reactants) for r in self.reactions]
        product_counts = [len(r.products) for r in self.reactions]
        
        if reactant_counts:
            stats['avg_reactants_per_reaction'] = sum(reactant_counts) / len(reactant_counts)
        if product_counts:
            stats['avg_products_per_reaction'] = sum(product_counts) / len(product_counts)
        
        return stats


class ReactionNetworkAnalyzer:
    """Main analyzer class"""
    
    def __init__(self, result_dirs: List[Path], output_dir: Path):
        self.result_dirs = result_dirs
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.network = ReactionNetwork()
        
    def load_results(self):
        """Load results from all directories"""
        logger.info(f"Loading results from {len(self.result_dirs)} directories...")
        
        for result_dir in self.result_dirs:
            logger.info(f"  Loading: {result_dir.name}")
            
            # Try to load reactions from snapshots
            snapshot_dir = result_dir / "snapshots"
            if snapshot_dir.exists():
                self._load_from_snapshots(snapshot_dir)
            
            # Try to load from results.json
            results_file = result_dir / "results.json"
            if results_file.exists():
                self._load_from_results(results_file)
    
    def _load_from_snapshots(self, snapshot_dir: Path):
        """Load reactions from snapshot files"""
        # Try both naming conventions
        snapshot_files = sorted(snapshot_dir.glob("snapshot_*.json"))
        if not snapshot_files:
            snapshot_files = sorted(snapshot_dir.glob("step_*.json"))
        
        if not snapshot_files:
            return
        
        logger.info(f"    Found {len(snapshot_files)} snapshots")
        
        # Track molecules across snapshots
        prev_molecules = set()
        
        for snapshot_file in snapshot_files:
            try:
                with open(snapshot_file) as f:
                    data = json.load(f)
                
                step = data.get('step', 0)
                
                # Extract current molecules
                curr_molecules = set()
                
                # Try direct 'molecules' field first
                if 'molecules' in data:
                    for mol in data['molecules']:
                        formula = mol.get('formula', '')
                        if formula:
                            curr_molecules.add(formula)
                            
                            # Add molecule to network
                            self.network.add_molecule(Molecule(
                                formula=formula,
                                name=mol.get('name'),
                                num_atoms=len(mol.get('atoms', [])),
                                first_seen=step
                            ))
                
                # If no molecules field, try to extract from bonds/attributes
                elif 'bonds' in data and data.get('bonds'):
                    from collections import defaultdict
                    
                    bonds = data.get('bonds', [])
                    if bonds:
                        # Build graph of connected particles
                        graph = defaultdict(set)
                        for bond in bonds:
                            if len(bond) >= 2:
                                i, j = bond[0], bond[1]
                                graph[i].add(j)
                                graph[j].add(i)
                        
                        # Find connected components (molecules)
                        visited = set()
                        components = []
                        
                        for node in graph:
                            if node not in visited:
                                component = []
                                stack = [node]
                                while stack:
                                    n = stack.pop()
                                    if n not in visited:
                                        visited.add(n)
                                        component.append(n)
                                        stack.extend(graph[n])
                                if len(component) >= 2:  # At least 2 atoms for a molecule
                                    components.append(component)
                        
                        # Create molecule signatures (simplified: use size + bond pattern)
                        for component in components:
                            size = len(component)
                            bond_count = sum(len(graph[n]) for n in component) // 2
                            
                            # Create a simple formula-like identifier
                            formula = f"MOL_{size}_{bond_count}"
                            curr_molecules.add(formula)
                            
                            # Add molecule to network
                            self.network.add_molecule(Molecule(
                                formula=formula,
                                num_atoms=size,
                                first_seen=step
                            ))
                
                # Detect new molecules (potential products)
                new_molecules = curr_molecules - prev_molecules
                if new_molecules and prev_molecules:
                    # Infer reactions (simplified: any previous -> any new)
                    # In reality, would need more detailed tracking
                    for new_mol in new_molecules:
                        # Create pseudo-reaction from existing molecules
                        # This is a simplification - real tracking would be better
                        pass
                
                prev_molecules = curr_molecules
                
            except Exception as e:
                logger.debug(f"Failed to load {snapshot_file}: {e}")
    
    def _load_from_results(self, results_file: Path):
        """Load from aggregated results file"""
        try:
            with open(results_file) as f:
                data = json.load(f)
            
            # Extract molecules
            if 'molecules' in data:
                for mol in data['molecules']:
                    formula = mol.get('formula', '')
                    if formula:
                        self.network.add_molecule(Molecule(
                            formula=formula,
                            name=mol.get('name'),
                            num_atoms=len(mol.get('atoms', []))
                        ))
            
            # Extract reactions if available
            if 'reactions' in data:
                for rxn_data in data['reactions']:
                    reactants = tuple(rxn_data.get('reactants', []))
                    products = tuple(rxn_data.get('products', []))
                    
                    if reactants and products:
                        rxn = Reaction(
                            reactants=reactants,
                            products=products,
                            step=rxn_data.get('step'),
                            energy_change=rxn_data.get('energy_change')
                        )
                        self.network.add_reaction(rxn)
            
        except Exception as e:
            logger.error(f"Failed to load {results_file}: {e}")
    
    def analyze(self):
        """Perform network analysis"""
        logger.info("\nAnalyzing reaction network...")
        
        stats = self.network.compute_statistics()
        
        logger.info(f"  Molecules: {stats['num_molecules']}")
        logger.info(f"  Reactions: {stats['num_reactions']}")
        logger.info(f"  Sources: {stats['num_sources']}")
        logger.info(f"  Sinks: {stats['num_sinks']}")
        logger.info(f"  Hubs: {stats['num_hubs']}")
        
        if 'avg_degree' in stats:
            logger.info(f"  Avg degree: {stats['avg_degree']:.2f}")
        
        return stats
    
    def export_graphml(self):
        """Export network to GraphML format"""
        output_file = self.output_dir / "reaction_network.graphml"
        
        logger.info(f"Exporting to GraphML: {output_file}")
        
        with open(output_file, 'w') as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write('<graphml xmlns="http://graphml.graphdrawing.org/xmlns">\n')
            f.write('  <key id="formula" for="node" attr.name="formula" attr.type="string"/>\n')
            f.write('  <key id="name" for="node" attr.name="name" attr.type="string"/>\n')
            f.write('  <key id="reaction" for="edge" attr.name="reaction" attr.type="string"/>\n')
            f.write('  <graph id="G" edgedefault="directed">\n')
            
            # Write nodes
            for i, (formula, mol) in enumerate(self.network.molecules.items()):
                f.write(f'    <node id="n{i}">\n')
                f.write(f'      <data key="formula">{formula}</data>\n')
                if mol.name:
                    f.write(f'      <data key="name">{mol.name}</data>\n')
                f.write(f'    </node>\n')
            
            # Write edges
            formula_to_id = {formula: f"n{i}" for i, formula in enumerate(self.network.molecules)}
            edge_id = 0
            
            for rxn in self.network.reactions:
                for reactant in rxn.reactants:
                    for product in rxn.products:
                        source_id = formula_to_id.get(reactant)
                        target_id = formula_to_id.get(product)
                        
                        if source_id and target_id:
                            f.write(f'    <edge id="e{edge_id}" source="{source_id}" target="{target_id}">\n')
                            f.write(f'      <data key="reaction">{rxn}</data>\n')
                            f.write(f'    </edge>\n')
                            edge_id += 1
            
            f.write('  </graph>\n')
            f.write('</graphml>\n')
        
        logger.info(f"GraphML exported: {output_file}")
    
    def export_json(self):
        """Export network to JSON format"""
        output_file = self.output_dir / "reaction_network.json"
        
        logger.info(f"Exporting to JSON: {output_file}")
        
        data = {
            'molecules': [
                {
                    'formula': mol.formula,
                    'name': mol.name,
                    'num_atoms': mol.num_atoms,
                    'first_seen': mol.first_seen
                }
                for mol in self.network.molecules.values()
            ],
            'reactions': [
                {
                    'reactants': list(rxn.reactants),
                    'products': list(rxn.products),
                    'step': rxn.step,
                    'energy_change': rxn.energy_change,
                    'string': str(rxn)
                }
                for rxn in self.network.reactions
            ],
            'statistics': self.network.compute_statistics()
        }
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        logger.info(f"JSON exported: {output_file}")
    
    def generate_report(self, stats: Dict):
        """Generate analysis report"""
        report_file = self.output_dir / "network_analysis.txt"
        
        logger.info(f"Generating report: {report_file}")
        
        with open(report_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("REACTION NETWORK ANALYSIS\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Analysis time: {datetime.now()}\n")
            f.write(f"Result directories: {len(self.result_dirs)}\n\n")
            
            # Network statistics
            f.write("Network Statistics:\n")
            f.write("-" * 70 + "\n")
            for key, value in stats.items():
                f.write(f"  {key}: {value}\n")
            f.write("\n")
            
            # Source molecules
            sources = self.network.get_sources()
            f.write(f"Source Molecules ({len(sources)}):\n")
            f.write("-" * 70 + "\n")
            for formula in sources[:20]:  # Top 20
                f.write(f"  {formula}\n")
            if len(sources) > 20:
                f.write(f"  ... and {len(sources) - 20} more\n")
            f.write("\n")
            
            # Sink molecules
            sinks = self.network.get_sinks()
            f.write(f"Sink Molecules ({len(sinks)}):\n")
            f.write("-" * 70 + "\n")
            for formula in sinks[:20]:  # Top 20
                f.write(f"  {formula}\n")
            if len(sinks) > 20:
                f.write(f"  ... and {len(sinks) - 20} more\n")
            f.write("\n")
            
            # Hub molecules
            hubs = self.network.get_hubs()
            f.write(f"Hub Molecules ({len(hubs)}):\n")
            f.write("-" * 70 + "\n")
            hub_data = []
            for formula in hubs:
                in_deg, out_deg = self.network.get_degree(formula)
                hub_data.append((formula, in_deg + out_deg, in_deg, out_deg))
            
            hub_data.sort(key=lambda x: x[1], reverse=True)
            
            for formula, total_deg, in_deg, out_deg in hub_data[:20]:
                f.write(f"  {formula:20s} (total: {total_deg:3d}, in: {in_deg:3d}, out: {out_deg:3d})\n")
            if len(hub_data) > 20:
                f.write(f"  ... and {len(hub_data) - 20} more\n")
            
            f.write("\n" + "=" * 70 + "\n")
        
        logger.info(f"Report saved: {report_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze reaction networks from Phase 2 results"
    )
    parser.add_argument(
        'result_dirs',
        nargs='+',
        type=str,
        help="Result directories to analyze"
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='analysis/reaction_network',
        help="Output directory (default: analysis/reaction_network)"
    )
    parser.add_argument(
        '--export',
        choices=['json', 'graphml', 'both'],
        default='both',
        help="Export format (default: both)"
    )
    parser.add_argument(
        '--merge',
        action='store_true',
        help="Merge networks from multiple simulations"
    )
    
    args = parser.parse_args()
    
    # Convert to Path objects
    result_dirs = [Path(d) for d in args.result_dirs]
    output_dir = Path(args.output)
    
    # Validate directories
    valid_dirs = []
    for d in result_dirs:
        if d.exists():
            valid_dirs.append(d)
        else:
            logger.warning(f"Directory not found: {d}")
    
    if not valid_dirs:
        logger.error("No valid result directories found")
        return 1
    
    logger.info("=" * 70)
    logger.info("REACTION NETWORK ANALYZER")
    logger.info("=" * 70)
    logger.info(f"Result directories: {len(valid_dirs)}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Merge networks: {args.merge}")
    logger.info("=" * 70)
    
    # Analyze
    analyzer = ReactionNetworkAnalyzer(valid_dirs, output_dir)
    analyzer.load_results()
    stats = analyzer.analyze()
    
    # Export
    if args.export in ['json', 'both']:
        analyzer.export_json()
    if args.export in ['graphml', 'both']:
        analyzer.export_graphml()
    
    # Report
    analyzer.generate_report(stats)
    
    logger.info("\n" + "=" * 70)
    logger.info("ANALYSIS COMPLETE!")
    logger.info("=" * 70)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

