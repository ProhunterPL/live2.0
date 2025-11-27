"""Debug script to test autocatalysis detection with full pipeline"""
import json
import networkx as nx
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.sim.analysis.autocatalysis_detector import AutocatalysisDetector

# Load run_1 data
run_dir = Path("results/phase2b_additional/miller_urey_extended/run_1")
network_file = run_dir / "reaction_network.json"

with open(network_file) as f:
    r = json.load(f)

# Build graph
G = nx.DiGraph()
for edge in r.get('edges', []):
    G.add_edge(edge['source'], edge['target'])

# Get abundance history
ah = r.get('abundance_history', {})
mol_names = {m['id']: m.get('formula', '') for m in r.get('molecules', [])}

print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
print(f"Abundance history: {len(ah)} molecules")

# Detect cycles with full pipeline
det = AutocatalysisDetector(min_amplification=1.1, cycle_timeout=60)  # Shorter timeout for testing
cycles = det.detect_cycles_in_network(G, ah, mol_names)

print(f"\nDetected {len(cycles)} autocatalytic cycles")
for i, cycle in enumerate(cycles[:5]):
    print(f"\nCycle {i+1}:")
    print(f"  Nodes: {cycle.nodes}")
    print(f"  Amplification: {cycle.amplification_factor:.3f}")
    print(f"  Type: {cycle.cycle_type}")

