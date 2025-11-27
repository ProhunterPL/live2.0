"""Debug script to test autocatalysis detection on single run"""
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
if ah:
    sample_key = list(ah.keys())[0]
    print(f"Sample history length: {len(ah[sample_key])}")
    print(f"Sample history: {ah[sample_key]}")

# Find cycles
det = AutocatalysisDetector(min_amplification=1.1)
cycles = det._find_all_cycles(G)
print(f"\nTotal cycles found: {len(cycles)}")

if cycles:
    # Check first few cycles
    for i, cycle in enumerate(cycles[:5]):
        print(f"\nCycle {i+1}: {cycle}")
        # Check which nodes are in abundance_history
        in_history = [n in ah for n in cycle]
        print(f"  Nodes in abundance_history: {in_history}")
        
        # Check amplification for each node
        for node in cycle:
            if node in ah:
                history = ah[node]
                print(f"  Node {node}: history_len={len(history)}, history={history}")
                if len(history) >= 3:
                    if len(history) < 6:
                        early = history[0] if history[0] > 0 else sum(history[:max(1, len(history)//2)]) / max(1, len(history)//2) + 1e-10
                        late = sum(history[-max(1, len(history)//2):]) / max(1, len(history)//2)
                    else:
                        early = sum(history[:len(history)//3]) / (len(history)//3) + 1e-10
                        late = sum(history[2*len(history)//3:]) / (len(history) - 2*len(history)//3)
                    amp = late / early if early > 0 else 0.0
                    print(f"    Amplification: {amp:.3f} (threshold: {det.min_amplification})")
                    print(f"    Passes: {amp >= det.min_amplification}")
        
        # Check if cycle is autocatalytic
        is_auto = det._is_autocatalytic(cycle, G, ah)
        print(f"  Is autocatalytic: {is_auto}")

