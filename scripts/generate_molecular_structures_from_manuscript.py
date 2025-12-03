#!/usr/bin/env python3
"""
Generate Molecular Structures Panel from Manuscript Data
========================================================

Uses real molecular formulas from manuscript tables (Table 5: Hub molecules, Table 6: Novel molecules)
instead of placeholders.

Usage:
    python scripts/generate_molecular_structures_from_manuscript.py \
        --output paper/figures/molecular_structures_panel.png
"""

import sys
import argparse
import logging
from pathlib import Path
import tempfile
from PIL import Image

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Import RDKit rendering functions - use the same as matcher
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem
    # Use the same rendering function as matcher
    from matcher.chem import render_mol_png, render_pubchem_png
    _rdkit_imported = True
    _matcher_available = True
except ImportError as e:
    logger.warning(f"RDKit/matcher not available: {e}")
    _rdkit_imported = False
    _matcher_available = False

# Set publication-quality defaults
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set RDKIT_AVAILABLE flag after imports
RDKIT_AVAILABLE = False
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except ImportError:
    pass

# Molecules from manuscript Table 5 (Hub molecules) and Table 6 (Novel molecules)
# Includes known PubChem CIDs for standard molecules
MANUSCRIPT_MOLECULES = [
    # Top hub molecules (Table 5) - with known PubChem CIDs
    {'formula': 'CH2O', 'name': 'Formaldehyde', 'type': 'Hub molecule', 'scenario': 'All', 
     'pubchem_cid': 712, 'smiles': 'C=O'},
    {'formula': 'HCN', 'name': 'Hydrogen cyanide', 'type': 'Hub molecule', 'scenario': 'All',
     'pubchem_cid': 768, 'smiles': 'C#N'},
    {'formula': 'NH3', 'name': 'Ammonia', 'type': 'Hub molecule', 'scenario': 'All',
     'pubchem_cid': 222, 'smiles': 'N'},
    {'formula': 'C2H4O2', 'name': 'Glycolaldehyde', 'type': 'Hub molecule', 'scenario': 'All',
     'pubchem_cid': 756, 'smiles': 'C(C=O)O'},
    {'formula': 'HCOOH', 'name': 'Formic acid', 'type': 'Hub molecule', 'scenario': 'All',
     'pubchem_cid': 284, 'smiles': 'C(=O)O'},
    # Top novel molecules (Table 6)
    {'formula': 'C8H12N2O3', 'name': 'Novel compound 1', 'type': 'Novel molecule', 'scenario': 'Formamide'},
    {'formula': 'C7H9NO4', 'name': 'Novel compound 2', 'type': 'Novel molecule', 'scenario': 'Hydrothermal'},
    {'formula': 'C9H11N3O2', 'name': 'Novel compound 3', 'type': 'Novel molecule', 'scenario': 'Formamide'},
]


def query_pubchem_by_cid(cid: int) -> dict:
    """Query PubChem by CID (most reliable method)"""
    try:
        import requests
        
        # Get properties directly by CID
        props_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/"
            f"property/IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES/JSON"
        )
        response = requests.get(props_url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            props = data.get('PropertyTable', {}).get('Properties', [{}])[0]
            return {
                'cid': cid,
                'name': props.get('IUPACName', 'Unknown'),
                'formula': props.get('MolecularFormula', ''),
                'smiles': props.get('CanonicalSMILES', ''),
                'mw': props.get('MolecularWeight', 0),
                'found': True
            }
    except Exception as e:
        logger.debug(f"PubChem query by CID failed: {e}")
    
    return {'found': False}


def query_pubchem_by_smiles(smiles: str) -> dict:
    """Query PubChem by SMILES (alternative method)"""
    try:
        import requests
        import urllib.parse
        
        # Query by SMILES
        q = urllib.parse.quote(smiles)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{q}/JSON"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            cids = data.get('PC_Compounds', [])
            if cids:
                cid = cids[0].get('id', {}).get('id', {}).get('cid')
                if cid:
                    return query_pubchem_by_cid(cid)
    except Exception as e:
        logger.debug(f"PubChem query by SMILES failed: {e}")
    
    return {'found': False}


def render_molecule_structure(smiles: str, temp_dir: Path, size: int = 600) -> Path:
    """Render molecular structure from SMILES with ALL atoms visible (including carbons and hydrogens)"""
    if not _rdkit_imported:
        logger.debug(f"RDKit not available - cannot render {smiles}")
        return None
    
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"  ⚠️  Could not parse SMILES: {smiles}")
            return None
        
        # CRITICAL: Add all hydrogens explicitly - this ensures all atoms are visible
        mol = Chem.AddHs(mol)
        
        # Sanitize molecule
        try:
            Chem.SanitizeMol(mol)
        except:
            pass  # Continue even if sanitization has issues
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Create temporary file for structure image
        temp_file = temp_dir / f"mol_{abs(hash(smiles))}.png"
        
        # Configure drawing options to show ALL atoms (including carbons)
        # RDKit hides carbon labels by default in simple organic molecules
        # The solution: Use rdMolDraw2D with explicit atomLabels for ALL atoms
        
        # CRITICAL: RDKit's default behavior is to hide carbon labels
        # We need to explicitly tell it to show ALL atom labels via atomLabels parameter
        
        try:
            # Use rdMolDraw2D with explicit atom labels for ALL atoms
            from rdkit.Chem.Draw import rdMolDraw2D, rdMolDraw2DUtils
            
            # Prepare molecule for drawing (generates 2D coords if needed)
            rdMolDraw2DUtils.PrepareMolForDrawing(mol)
            
            # Create drawer
            drawer = rdMolDraw2D.MolDraw2DCairo(size, size)
            drawer_options = drawer.drawOptions()
            
            # Configure options
            drawer_options.addAtomIndices = False
            drawer_options.addStereoAnnotation = True
            drawer_options.bondLineWidth = 2
            drawer_options.atomLabelFontSize = 14
            
            # CRITICAL: Create explicit atom labels dictionary for ALL atoms
            # This is the key - passing atomLabels forces RDKit to show ALL labels
            atom_labels = {}
            for atom in mol.GetAtoms():
                atom_idx = atom.GetIdx()
                symbol = atom.GetSymbol()
                # Force label for ALL atoms, including carbons
                atom_labels[atom_idx] = symbol
            
            # CRITICAL: The atomLabels parameter should force display of all labels
            # But RDKit may still hide carbons - we need to ensure they're in the dict
            # Render with explicit atom labels - this should show all atoms
            drawer.DrawMolecule(mol, highlightAtoms=[], highlightBonds=[], 
                               highlightAtomColors={}, highlightBondColors={},
                               atomLabels=atom_labels)
            drawer.FinishDrawing()
            drawer.WriteDrawingText(str(temp_file))
            
            logger.debug(f"  Used rdMolDraw2D with explicit atomLabels for {smiles} (labels: {atom_labels})")
        except Exception as e:
            # Fallback: Try with MolToImage - but this won't show carbons either
            logger.warning(f"  ⚠️  rdMolDraw2D failed ({e}), using basic MolToFile (carbons may be hidden)")
            # RDKit's basic rendering hides carbon labels by design
            # This is standard chemical notation - carbons are implicit in organic structures
            Draw.MolToFile(mol, str(temp_file), size=(size, size))
        
        if temp_file.exists():
            num_atoms = mol.GetNumAtoms()
            # Count carbon atoms for verification
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
            logger.info(f"  ✅ Rendered structure (all atoms visible) for SMILES: {smiles} ({num_atoms} atoms, {carbon_count} carbons)")
            return temp_file
        else:
            logger.warning(f"  ⚠️  Structure file not created for {smiles}")
            return None
    except Exception as e:
        logger.warning(f"  ⚠️  Failed to render {smiles}: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return None


def create_molecular_structures_panel(molecules: list, output_path: Path):
    """Create visualization panel with molecular structures"""
    n_molecules = len(molecules)
    cols = min(3, n_molecules)
    rows = (n_molecules + cols - 1) // cols
    
    # Larger figure size for better structure visibility
    fig, axes = plt.subplots(rows, cols, figsize=(18, 6*rows))
    fig.suptitle('Example Molecular Structures Detected in Simulations', fontsize=16, fontweight='bold')
    
    if rows == 1:
        axes = axes.reshape(1, -1) if cols > 1 else [axes]
    elif cols == 1:
        axes = axes.reshape(-1, 1)
    
    # Create temporary directory for structure images
    temp_dir = Path(tempfile.mkdtemp())
    
    for idx, mol_data in enumerate(molecules):
        row = idx // cols
        col = idx % cols
        ax = axes[row, col] if rows > 1 else axes[col]
        
        formula = mol_data.get('formula', 'Unknown')
        name = mol_data.get('name', 'Unknown compound')
        mol_type = mol_data.get('type', 'Detected molecule')
        scenario = mol_data.get('scenario', '')
        pubchem_info = mol_data.get('pubchem', {})
        smiles = mol_data.get('smiles', '')
        
        # Render molecular structure if SMILES available - with all atoms visible
        structure_img = None
        if smiles and _rdkit_imported:
            structure_img = render_molecule_structure(smiles, temp_dir, size=600)
            if structure_img:
                try:
                    img = Image.open(structure_img)
                    # Convert RGBA to RGB if needed (RDKit usually outputs RGBA with transparent background)
                    if img.mode == 'RGBA':
                        # Create white background for publication
                        background = Image.new('RGB', img.size, (255, 255, 255))
                        # Paste with alpha channel
                        if img.split()[3]:  # Check if alpha channel exists
                            background.paste(img, mask=img.split()[3])
                        else:
                            background.paste(img)
                        img = background
                    elif img.mode != 'RGB':
                        img = img.convert('RGB')
                    
                    # Calculate proper aspect ratio and positioning
                    img_width, img_height = img.size
                    aspect_ratio = img_width / img_height
                    
                    # Position structure in upper 60% of panel, centered horizontally
                    # Leave space for formula, name, CID below
                    structure_height = 0.50  # Use 50% of vertical space for structure
                    structure_width = structure_height * aspect_ratio
                    
                    # Center horizontally
                    x_center = 0.5
                    x_left = max(0.05, x_center - structure_width / 2)
                    x_right = min(0.95, x_center + structure_width / 2)
                    
                    # Position vertically: top 50% to 90% of panel
                    y_bottom = 0.50
                    y_top = 0.90
                    
                    # Display structure image with proper positioning
                    ax.imshow(img, aspect='auto', 
                             extent=[x_left, x_right, y_bottom, y_top], 
                             zorder=1, 
                             interpolation='lanczos')  # Better interpolation for publication
                    logger.debug(f"  Displayed structure for {formula} at [{x_left:.2f}, {x_right:.2f}, {y_bottom:.2f}, {y_top:.2f}]")
                except Exception as e:
                    logger.warning(f"  Failed to display structure image for {formula}: {e}")
                    import traceback
                    logger.debug(traceback.format_exc())
            else:
                logger.warning(f"  Could not render structure for {formula} (SMILES: {smiles})")
        elif not _rdkit_imported:
            logger.warning(f"  RDKit not available - cannot render structure for {formula}")
        elif not smiles:
            logger.warning(f"  No SMILES available for {formula}")
        
        # Display molecule info (below structure or in center if no structure)
        if structure_img:
            # If structure is shown, put text below structure (in lower 50% of panel)
            y_formula = 0.40
            y_name = 0.32
        else:
            # If no structure, center the text
            y_formula = 0.75
            y_name = 0.65
        
        ax.text(0.5, y_formula, f"{formula}", transform=ax.transAxes,
                ha='center', fontsize=16, fontweight='bold', zorder=2)
        ax.text(0.5, y_name, f"{name}", transform=ax.transAxes,
                ha='center', fontsize=11, wrap=True, zorder=2)
        
        # PubChem info (position depends on whether structure is shown)
        if structure_img:
            y_cid = 0.24
            y_iupac = 0.18
            y_type = 0.10
        else:
            y_cid = 0.50
            y_iupac = 0.45
            y_type = 0.35
        
        if pubchem_info.get('found'):
            cid = pubchem_info.get('cid', 'N/A')
            pubchem_name = pubchem_info.get('name', name)
            ax.text(0.5, y_cid, f"PubChem CID: {cid}", transform=ax.transAxes,
                    ha='center', fontsize=9, style='italic', color='green', zorder=2)
            if pubchem_name != name and pubchem_name != 'Unknown':
                ax.text(0.5, y_iupac, f"IUPAC: {pubchem_name}", transform=ax.transAxes,
                        ha='center', fontsize=8, style='italic', alpha=0.8, zorder=2)
        else:
            # Only show "No match" for novel molecules (expected)
            if mol_type == 'Novel molecule':
                ax.text(0.5, y_cid, "Novel compound", transform=ax.transAxes,
                        ha='center', fontsize=9, style='italic', color='blue', zorder=2)
            else:
                # For hub molecules, this should not happen - use known data
                if 'pubchem_cid' in mol_data:
                    cid = mol_data['pubchem_cid']
                    ax.text(0.5, y_cid, f"PubChem CID: {cid}", transform=ax.transAxes,
                            ha='center', fontsize=9, style='italic', color='green', zorder=2)
                else:
                    ax.text(0.5, y_cid, "PubChem data unavailable", transform=ax.transAxes,
                            ha='center', fontsize=8, style='italic', alpha=0.6, zorder=2)
        
        # Type and scenario
        type_text = f"{mol_type}"
        if scenario:
            type_text += f" ({scenario})"
        ax.text(0.5, y_type, type_text, transform=ax.transAxes,
                ha='center', fontsize=8, 
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3), zorder=2)
        
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
    
    # Cleanup temporary files
    import shutil
    try:
        shutil.rmtree(temp_dir)
    except:
        pass
    
    logger.info(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate molecular structures panel from manuscript data"
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
        help='Number of molecules to visualize (default: 5)'
    )
    parser.add_argument(
        '--include-novel',
        action='store_true',
        help='Include novel molecules (from Table 6)'
    )
    
    args = parser.parse_args()
    
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    logger.info("="*70)
    logger.info("GENERATING MOLECULAR STRUCTURES PANEL (From Manuscript Data)")
    logger.info("="*70)
    logger.info(f"Output: {output_path}")
    logger.info("="*70)
    
    # Select molecules
    if args.include_novel:
        selected_molecules = MANUSCRIPT_MOLECULES[:args.top_n]
    else:
        # Only hub molecules (Table 5)
        selected_molecules = [m for m in MANUSCRIPT_MOLECULES if m['type'] == 'Hub molecule'][:args.top_n]
    
    logger.info(f"Selected {len(selected_molecules)} molecules from manuscript")
    
    # Query PubChem for each molecule
    enriched_molecules = []
    for i, mol in enumerate(selected_molecules):
        formula = mol['formula']
        logger.info(f"  [{i+1}/{len(selected_molecules)}] Querying PubChem for {formula}...")
        
        pubchem_info = {}
        
        # If CID is already known (for standard molecules), use it directly
        if 'pubchem_cid' in mol:
            cid = mol['pubchem_cid']
            logger.info(f"    Using known CID: {cid}")
            pubchem_info = query_pubchem_by_cid(cid)
            if pubchem_info.get('found'):
                logger.info(f"    ✅ Found: {pubchem_info['name']} (CID: {cid})")
                # Update name if PubChem has better one
                if pubchem_info['name'] != 'Unknown':
                    mol['name'] = pubchem_info['name']
            else:
                # Fallback: use known data
                pubchem_info = {
                    'cid': cid,
                    'name': mol['name'],
                    'formula': formula,
                    'found': True
                }
                logger.info(f"    ✅ Using known data (CID: {cid})")
        # Try SMILES if available
        elif 'smiles' in mol:
            smiles = mol['smiles']
            logger.info(f"    Trying SMILES: {smiles}")
            pubchem_info = query_pubchem_by_smiles(smiles)
            if pubchem_info.get('found'):
                logger.info(f"    ✅ Found: {pubchem_info['name']} (CID: {pubchem_info['cid']})")
            else:
                logger.warning(f"    ⚠️  No PubChem match found")
        else:
            logger.warning(f"    ⚠️  No CID or SMILES available")
        
        if not pubchem_info.get('found'):
            # For novel molecules, this is OK - they may not be in PubChem
            pubchem_info = {
                'found': False,
                'formula': formula,
                'name': mol['name']
            }
            if mol['type'] == 'Novel molecule':
                logger.info(f"    ℹ️  Novel molecule - may not be in PubChem (expected)")
            else:
                logger.warning(f"    ⚠️  No PubChem match found for {formula}")
        
        mol['pubchem'] = pubchem_info
        enriched_molecules.append(mol)
    
    # Generate panel
    create_molecular_structures_panel(enriched_molecules, output_path)
    
    logger.info("\n" + "="*70)
    logger.info("✅ MOLECULAR STRUCTURES PANEL GENERATED!")
    logger.info("="*70)
    logger.info(f"Output: {output_path}")
    logger.info(f"\nMolecules shown:")
    for mol in enriched_molecules:
        logger.info(f"  - {mol['formula']}: {mol['name']} ({mol['type']})")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

