"""
Network Visualizer
==================

Creates visualizations of reaction networks and autocatalytic cycles.

Generates:
- Network topology plots (nodes = molecules, edges = reactions)
- Degree distribution histograms
- Autocatalytic cycle diagrams
- Heatmaps of connectivity
- Interactive HTML visualizations

Requirements:
    pip install matplotlib networkx

Usage:
    # Visualize network from JSON
    python scripts/network_visualizer.py analysis/reaction_network/reaction_network.json
    
    # Include autocatalytic cycles
    python scripts/network_visualizer.py \
        analysis/reaction_network/reaction_network.json \
        --cycles analysis/autocatalytic_cycles/autocatalytic_cycles.json
    
    # Generate interactive HTML
    python scripts/network_visualizer.py \
        analysis/reaction_network/reaction_network.json \
        --interactive
"""

import sys
import argparse
import logging
import json
from pathlib import Path
from typing import List, Dict, Set, Tuple, Optional
from collections import Counter

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class NetworkVisualizer:
    """Visualizes reaction networks"""
    
    def __init__(self, network_file: Path, output_dir: Path, cycles_file: Optional[Path] = None):
        self.network_file = network_file
        self.cycles_file = cycles_file
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.network_data = None
        self.cycles_data = None
        
    def load_data(self):
        """Load network and cycle data"""
        logger.info(f"Loading network: {self.network_file}")
        
        with open(self.network_file) as f:
            self.network_data = json.load(f)
        
        if self.cycles_file and self.cycles_file.exists():
            logger.info(f"Loading cycles: {self.cycles_file}")
            with open(self.cycles_file) as f:
                self.cycles_data = json.load(f)
    
    def create_degree_distribution(self):
        """Create degree distribution histogram"""
        logger.info("Creating degree distribution plot...")
        
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.warning("matplotlib not installed, skipping plot")
            return
        
        # Extract molecules and reactions
        molecules = {mol['formula']: mol for mol in self.network_data.get('molecules', [])}
        reactions = self.network_data.get('reactions', [])
        
        # Calculate degrees
        in_degree = Counter()
        out_degree = Counter()
        
        for rxn in reactions:
            for reactant in rxn['reactants']:
                out_degree[reactant] += 1
            for product in rxn['products']:
                in_degree[product] += 1
        
        # Total degree
        total_degree = Counter()
        for mol in molecules:
            total_degree[mol] = in_degree[mol] + out_degree[mol]
        
        # Plot
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        
        # In-degree distribution
        if in_degree:
            degrees = list(in_degree.values())
            axes[0].hist(degrees, bins=20, edgecolor='black')
            axes[0].set_xlabel('In-degree')
            axes[0].set_ylabel('Frequency')
            axes[0].set_title('In-degree Distribution')
            axes[0].grid(True, alpha=0.3)
        
        # Out-degree distribution
        if out_degree:
            degrees = list(out_degree.values())
            axes[1].hist(degrees, bins=20, edgecolor='black')
            axes[1].set_xlabel('Out-degree')
            axes[1].set_ylabel('Frequency')
            axes[1].set_title('Out-degree Distribution')
            axes[1].grid(True, alpha=0.3)
        
        # Total degree distribution
        if total_degree:
            degrees = list(total_degree.values())
            axes[2].hist(degrees, bins=20, edgecolor='black')
            axes[2].set_xlabel('Total degree')
            axes[2].set_ylabel('Frequency')
            axes[2].set_title('Total Degree Distribution')
            axes[2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        output_file = self.output_dir / "degree_distribution.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Saved: {output_file}")
    
    def create_network_topology(self, max_nodes: int = 100):
        """Create network topology visualization"""
        logger.info("Creating network topology plot...")
        
        try:
            import matplotlib.pyplot as plt
            import networkx as nx
        except ImportError:
            logger.warning("matplotlib or networkx not installed, skipping plot")
            return
        
        # Build NetworkX graph
        G = nx.DiGraph()
        
        molecules = self.network_data.get('molecules', [])[:max_nodes]
        reactions = self.network_data.get('reactions', [])
        
        # Add nodes
        for mol in molecules:
            G.add_node(mol['formula'], name=mol.get('name', ''))
        
        # Add edges
        edge_count = {}
        for rxn in reactions:
            for reactant in rxn['reactants']:
                for product in rxn['products']:
                    if reactant in G.nodes and product in G.nodes:
                        edge = (reactant, product)
                        edge_count[edge] = edge_count.get(edge, 0) + 1
        
        for (source, target), count in edge_count.items():
            G.add_edge(source, target, weight=count)
        
        # Layout
        logger.info(f"  Computing layout for {G.number_of_nodes()} nodes...")
        
        if G.number_of_nodes() > 0:
            # Use spring layout for small graphs, circular for large
            if G.number_of_nodes() <= 30:
                pos = nx.spring_layout(G, k=2, iterations=50)
            else:
                pos = nx.circular_layout(G)
            
            # Plot
            fig, ax = plt.subplots(figsize=(16, 12))
            
            # Draw nodes
            node_sizes = [G.degree(node) * 100 + 100 for node in G.nodes()]
            nx.draw_networkx_nodes(
                G, pos, node_size=node_sizes, 
                node_color='lightblue', alpha=0.7, ax=ax
            )
            
            # Draw edges
            nx.draw_networkx_edges(
                G, pos, edge_color='gray', 
                alpha=0.3, arrows=True, ax=ax,
                arrowsize=10, arrowstyle='->'
            )
            
            # Draw labels (only for small graphs)
            if G.number_of_nodes() <= 50:
                labels = {node: node for node in G.nodes()}
                nx.draw_networkx_labels(
                    G, pos, labels, font_size=8, ax=ax
                )
            
            ax.set_title(f'Reaction Network Topology ({G.number_of_nodes()} molecules, {G.number_of_edges()} connections)', 
                        fontsize=14, fontweight='bold')
            ax.axis('off')
            
            plt.tight_layout()
            output_file = self.output_dir / "network_topology.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Saved: {output_file}")
        else:
            logger.warning("No nodes in graph, skipping topology plot")
    
    def create_cycle_diagrams(self):
        """Create diagrams of autocatalytic cycles"""
        if not self.cycles_data:
            logger.info("No cycle data available, skipping cycle diagrams")
            return
        
        logger.info("Creating autocatalytic cycle diagrams...")
        
        try:
            import matplotlib.pyplot as plt
            import networkx as nx
        except ImportError:
            logger.warning("matplotlib or networkx not installed, skipping plot")
            return
        
        cycles = self.cycles_data.get('cycles', [])[:12]  # Top 12 cycles
        
        if not cycles:
            logger.warning("No cycles to visualize")
            return
        
        # Create subplots
        n_cycles = len(cycles)
        n_cols = min(4, n_cycles)
        n_rows = (n_cycles + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
        if n_rows == 1 and n_cols == 1:
            axes = [[axes]]
        elif n_rows == 1:
            axes = [axes]
        
        for idx, cycle in enumerate(cycles):
            row = idx // n_cols
            col = idx % n_cols
            ax = axes[row][col] if n_rows > 1 else axes[col]
            
            # Build cycle graph
            G = nx.DiGraph()
            molecules = cycle['molecules']
            
            for i, mol in enumerate(molecules):
                G.add_node(mol)
                if i < len(molecules) - 1:
                    G.add_edge(mol, molecules[i+1])
                elif len(molecules) > 1:
                    # Close the cycle
                    G.add_edge(mol, molecules[0])
            
            # Layout
            if len(molecules) <= 3:
                pos = nx.circular_layout(G)
            else:
                pos = nx.spring_layout(G, k=1, iterations=30)
            
            # Draw
            nx.draw_networkx_nodes(G, pos, node_color='lightcoral', node_size=800, ax=ax)
            nx.draw_networkx_edges(G, pos, edge_color='red', arrows=True, 
                                  arrowsize=15, arrowstyle='->', ax=ax, width=2)
            
            labels = {node: node[:8] for node in G.nodes()}  # Truncate long formulas
            nx.draw_networkx_labels(G, pos, labels, font_size=9, ax=ax)
            
            ax.set_title(f"{cycle['type'].upper()} (size {cycle['size']})", 
                        fontsize=10, fontweight='bold')
            ax.axis('off')
        
        # Remove empty subplots
        for idx in range(n_cycles, n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            ax = axes[row][col] if n_rows > 1 else axes[col]
            ax.axis('off')
        
        plt.tight_layout()
        output_file = self.output_dir / "autocatalytic_cycles.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Saved: {output_file}")
    
    def create_statistics_summary(self):
        """Create visual summary of network statistics"""
        logger.info("Creating statistics summary...")
        
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.warning("matplotlib not installed, skipping plot")
            return
        
        stats = self.network_data.get('statistics', {})
        
        if not stats:
            logger.warning("No statistics available")
            return
        
        # Create figure
        fig = plt.figure(figsize=(12, 6))
        
        # Text summary
        ax = fig.add_subplot(111)
        ax.axis('off')
        
        summary_text = "REACTION NETWORK STATISTICS\n" + "="*50 + "\n\n"
        
        for key, value in stats.items():
            label = key.replace('_', ' ').title()
            if isinstance(value, float):
                summary_text += f"{label:40s}: {value:.2f}\n"
            else:
                summary_text += f"{label:40s}: {value}\n"
        
        # Add cycle statistics if available
        if self.cycles_data:
            summary = self.cycles_data.get('summary', {})
            summary_text += "\n" + "="*50 + "\n"
            summary_text += "AUTOCATALYTIC CYCLES\n" + "="*50 + "\n\n"
            
            for key, value in summary.items():
                label = key.replace('_', ' ').title()
                summary_text += f"{label:40s}: {value}\n"
        
        ax.text(0.1, 0.9, summary_text, 
               transform=ax.transAxes,
               fontsize=11,
               verticalalignment='top',
               fontfamily='monospace')
        
        output_file = self.output_dir / "statistics_summary.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Saved: {output_file}")
    
    def create_interactive_html(self):
        """Create interactive HTML visualization"""
        logger.info("Creating interactive HTML visualization...")
        
        # Simple HTML with embedded data
        output_file = self.output_dir / "network_interactive.html"
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Reaction Network Visualization</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #333;
            border-bottom: 3px solid #4CAF50;
            padding-bottom: 10px;
        }}
        .stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .stat-card {{
            background-color: #f9f9f9;
            padding: 15px;
            border-radius: 5px;
            border-left: 4px solid #4CAF50;
        }}
        .stat-value {{
            font-size: 24px;
            font-weight: bold;
            color: #4CAF50;
        }}
        .stat-label {{
            color: #666;
            font-size: 14px;
        }}
        .section {{
            margin: 30px 0;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
        }}
        th, td {{
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #4CAF50;
            color: white;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🔬 Reaction Network Analysis</h1>
        
        <div class="section">
            <h2>Network Statistics</h2>
            <div class="stats">
"""
        
        # Add statistics
        stats = self.network_data.get('statistics', {})
        for key, value in list(stats.items())[:6]:
            label = key.replace('_', ' ').title()
            html_content += f"""
                <div class="stat-card">
                    <div class="stat-value">{value if isinstance(value, int) else f'{value:.2f}'}</div>
                    <div class="stat-label">{label}</div>
                </div>
"""
        
        html_content += """
            </div>
        </div>
        
        <div class="section">
            <h2>Top Molecules</h2>
            <table>
                <tr>
                    <th>#</th>
                    <th>Formula</th>
                    <th>Name</th>
                    <th>Atoms</th>
                </tr>
"""
        
        # Add top molecules
        molecules = self.network_data.get('molecules', [])[:20]
        for i, mol in enumerate(molecules, 1):
            html_content += f"""
                <tr>
                    <td>{i}</td>
                    <td><strong>{mol['formula']}</strong></td>
                    <td>{mol.get('name', 'Unknown')}</td>
                    <td>{mol.get('num_atoms', 'N/A')}</td>
                </tr>
"""
        
        html_content += """
            </table>
        </div>
"""
        
        # Add cycle information if available
        if self.cycles_data:
            cycles = self.cycles_data.get('cycles', [])[:10]
            html_content += """
        <div class="section">
            <h2>Autocatalytic Cycles</h2>
            <table>
                <tr>
                    <th>#</th>
                    <th>Type</th>
                    <th>Size</th>
                    <th>Cycle</th>
                </tr>
"""
            for i, cycle in enumerate(cycles, 1):
                cycle_str = cycle['string'][:100] + ('...' if len(cycle['string']) > 100 else '')
                html_content += f"""
                <tr>
                    <td>{i}</td>
                    <td><strong>{cycle['type']}</strong></td>
                    <td>{cycle['size']}</td>
                    <td>{cycle_str}</td>
                </tr>
"""
            html_content += """
            </table>
        </div>
"""
        
        html_content += """
    </div>
</body>
</html>
"""
        
        with open(output_file, 'w') as f:
            f.write(html_content)
        
        logger.info(f"Saved: {output_file}")
    
    def visualize(self, interactive: bool = False):
        """Run all visualizations"""
        logger.info("\n" + "=" * 70)
        logger.info("NETWORK VISUALIZATION")
        logger.info("=" * 70)
        
        self.load_data()
        
        # Create visualizations
        self.create_degree_distribution()
        self.create_network_topology()
        self.create_cycle_diagrams()
        self.create_statistics_summary()
        
        if interactive:
            self.create_interactive_html()
        
        logger.info("\n" + "=" * 70)
        logger.info("VISUALIZATION COMPLETE!")
        logger.info(f"Output directory: {self.output_dir}")
        logger.info("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description="Visualize reaction networks"
    )
    parser.add_argument(
        'network_file',
        type=str,
        help="Path to reaction network JSON file"
    )
    parser.add_argument(
        '-c', '--cycles',
        type=str,
        help="Path to autocatalytic cycles JSON file"
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='analysis/visualizations',
        help="Output directory (default: analysis/visualizations)"
    )
    parser.add_argument(
        '--interactive',
        action='store_true',
        help="Generate interactive HTML visualization"
    )
    parser.add_argument(
        '--max-nodes',
        type=int,
        default=100,
        help="Maximum nodes in topology plot (default: 100)"
    )
    
    args = parser.parse_args()
    
    network_file = Path(args.network_file)
    cycles_file = Path(args.cycles) if args.cycles else None
    output_dir = Path(args.output)
    
    if not network_file.exists():
        logger.error(f"Network file not found: {network_file}")
        return 1
    
    logger.info("=" * 70)
    logger.info("NETWORK VISUALIZER")
    logger.info("=" * 70)
    logger.info(f"Network: {network_file}")
    if cycles_file:
        logger.info(f"Cycles: {cycles_file}")
    logger.info(f"Output: {output_dir}")
    logger.info("=" * 70)
    
    # Run visualization
    visualizer = NetworkVisualizer(network_file, output_dir, cycles_file)
    visualizer.visualize(interactive=args.interactive)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

