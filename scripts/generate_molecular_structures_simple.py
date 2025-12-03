#!/usr/bin/env python3
"""
Simple Molecular Structures Generator (No RDKit Dependency)
============================================================

Generates molecular structures panel without requiring RDKit/NumPy 2.x compatibility.
Uses simple text-based visualization or basic molecular info.

Usage:
    python scripts/generate_molecular_structures_simple.py \
        --molecules-file results/phase2b_additional/miller_urey_extended/run_1/molecules.json \
        --output paper/figures/molecular_structures_panel.png
"""

import sys
import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Set publication-quality defaults
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def query_pubchem_simple(formula: str) -> Optional[Dict]:
    """
    Simple PubChem query using REST API (no RDKit needed)
    Returns basic info: name, CID, formula
    """
    try:
        import requests
        import urllib.parse
        
        # Query by formula
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula/{formula}/JSON"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            cids = data.get('PC_Compounds', [])
            if cids:
                cid = cids[0].get('id', {}).get('id', {}).get('cid')
                if cid:
                    # Get name
                    name_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
                    name_response = requests.get(name_url, timeout=10)
                    if name_response.status_code == 200:
                        name_data = name_response.json()
                        name = name_data.get('PropertyTable', {}).get('Properties', [{}])[0].get('IUPACName', 'Unknown')
                        return {
                            'cid': cid,
                            'name': name,
                            'formula': formula
                        }
    except Exception as e:
        logger.debug(f"PubChem query failed for {formula}: {e}")
    
    return None


def create_molecular_structures_panel(molecules: List[Dict], output_path: Path):
    """Create visualization panel with molecular structures (text-based)"""
    n_molecules = len(molecules)
    cols = min(3, n_molecules)
    rows = (n_molecules + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(15, 5*rows))
    fig.suptitle('Example Molecular Structures Detected in Simulations', fontsize=14, fontweight='bold')
    
    if rows == 1:
        axes = axes.reshape(1, -1) if cols > 1 else [axes]
    elif cols == 1:
        axes = axes.reshape(-1, 1)
    
    for idx, mol_data in enumerate(molecules):
        row = idx // cols
        col = idx % cols
        ax = axes[row, col] if rows > 1 else axes[col]
        
        formula = mol_data.get('formula', 'Unknown')
        pubchem_info = mol_data.get('pubchem', {})
        name = pubchem_info.get('name', 'Unknown compound')
        cid = pubchem_info.get('cid', 'N/A')
        abundance = mol_data.get('abundance', 0)
        
        # Display molecule info
        ax.text(0.5, 0.9, f"{formula}", transform=ax.transAxes,
                ha='center', fontsize=16, fontweight='bold')
        ax.text(0.5, 0.8, f"{name}", transform=ax.transAxes,
                ha='center', fontsize=11, wrap=True)
        ax.text(0.5, 0.7, f"PubChem CID: {cid}", transform=ax.transAxes,
                ha='center', fontsize=9)
        ax.text(0.5, 0.6, f"Abundance: {abundance}", transform=ax.transAxes,
                ha='center', fontsize=9)
        
        # Add simple molecular structure representation (ASCII art or formula)
        ax.text(0.5, 0.4, f"Structure:\n{formula}", transform=ax.transAxes,
                ha='center', fontsize=10, family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
    
    # Hide unused subplots
    for idx in range(n_molecules, rows * cols):
        row = idx // cols
        col = idx % cols
        ax = axes[row, col] if rows > 1 else axes[col]
        ax.axis('off')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate molecular structures panel (simple version, no RDKit)"
    )
    parser.add_argument(
        '--molecules-file',
        type=str,
        required=True,
        help='Path to molecules.json file'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='paper/figures/molecular_structures_panel.png',
        help='Output path for panel'
    )
    parser.add_argument(
        '--top-n',
        type=int,
        default=5,
        help='Number of top molecules to visualize (default: 5)'
    )
    
    args = parser.parse_args()
    
    molecules_file = Path(args.molecules_file)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if not molecules_file.exists():
        logger.error(f"Molecules file not found: {molecules_file}")
        return 1
    
    logger.info("="*70)
    logger.info("GENERATING MOLECULAR STRUCTURES PANEL (Simple Version)")
    logger.info("="*70)
    logger.info(f"Input: {molecules_file}")
    logger.info(f"Output: {output_path}")
    logger.info("="*70)
    
    # Load molecules
    with open(molecules_file, 'r') as f:
        molecules_data = json.load(f)
    
    # Handle different JSON formats
    if isinstance(molecules_data, list):
        # Direct list of molecules
        molecules = molecules_data
    elif isinstance(molecules_data, dict):
        # Dict with 'molecules' key
        molecules = molecules_data.get('molecules', [])
        if isinstance(molecules, dict):
            molecules = list(molecules.values())
    else:
        logger.error(f"Unexpected molecules.json format: {type(molecules_data)}")
        return 1
    
    # If molecules.json is empty, try to extract from snapshots
    if len(molecules) == 0:
        logger.warning("molecules.json is empty - trying to extract from snapshots...")
        try:
            # Add project root to path
            project_root = Path(__file__).parent.parent
            sys.path.insert(0, str(project_root))
            
            from backend.sim.molecule_extractor import extract_molecules_from_results
            
            result = extract_molecules_from_results(
                str(molecules_file.parent),
                output_dir=str(molecules_file.parent / "analysis"),
                export_for_matcher=False
            )
            molecules = result.get('molecules', [])
            logger.info(f"  Extracted {len(molecules)} molecules from snapshots")
        except Exception as e:
            logger.error(f"  Failed to extract from snapshots: {e}")
            logger.error("  💡 Try running: python scripts/extract_hydrothermal_molecules.py")
            logger.error("  💡 Or use existing analysis data if available")
            logger.error("  Cannot generate molecular structures panel without molecules")
            return 1
    
    if len(molecules) == 0:
        logger.error("No molecules found - cannot generate panel")
        return 1
    
    # Sort by abundance or count
    if len(molecules) > 0:
        # Try different abundance fields
        if 'abundance' in molecules[0]:
            top_molecules = sorted(molecules, key=lambda x: x.get('abundance', 0), reverse=True)[:args.top_n]
        elif 'count' in molecules[0]:
            top_molecules = sorted(molecules, key=lambda x: x.get('count', 0), reverse=True)[:args.top_n]
        else:
            top_molecules = molecules[:args.top_n]
    else:
        top_molecules = molecules[:args.top_n]
    
    logger.info(f"Selected {len(top_molecules)} molecules")
    
    # Query PubChem for each molecule
    enriched_molecules = []
    for i, mol in enumerate(top_molecules):
        formula = mol.get('formula', '')
        logger.info(f"  [{i+1}/{len(top_molecules)}] Querying PubChem for {formula}...")
        
        pubchem_info = query_pubchem_simple(formula)
        if pubchem_info:
            logger.info(f"    ✅ Found: {pubchem_info['name']} (CID: {pubchem_info['cid']})")
        else:
            logger.warning(f"    ⚠️  No PubChem match found")
            pubchem_info = {'name': 'Unknown compound', 'cid': 'N/A', 'formula': formula}
        
        enriched_molecules.append({
            'formula': formula,
            'abundance': mol.get('abundance', 0),
            'pubchem': pubchem_info
        })
    
    # Generate panel
    create_molecular_structures_panel(enriched_molecules, output_path)
    
    logger.info("\n" + "="*70)
    logger.info("✅ MOLECULAR STRUCTURES PANEL GENERATED!")
    logger.info("="*70)
    logger.info(f"Output: {output_path}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

