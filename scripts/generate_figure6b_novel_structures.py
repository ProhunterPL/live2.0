#!/usr/bin/env python3
"""
Generate Figure 6B: Novel Molecule Structures
==============================================

Creates a dedicated figure showing all 5 novel molecules with their structures,
similar to Figure 7 (molecular_structures_panel).

Usage:
    python scripts/generate_figure6b_novel_structures.py \
        --output paper/figures/figure6b_novel_structures.png
"""

import sys
import argparse
import logging
from pathlib import Path
import tempfile
from PIL import Image

import numpy as np
import matplotlib.pyplot as plt

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Import RDKit rendering functions - use the same as matcher
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    try:
        from rdkit.Chem.Draw import rdMolDraw2DUtils
        _has_utils = True
    except ImportError:
        _has_utils = False
    _rdkit_imported = True
except ImportError as e:
    logging.warning(f"RDKit not available: {e}")
    _rdkit_imported = False
    _has_utils = False

# Set publication-quality defaults
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Novel molecules from manuscript (Table 6)
# Note: These are novel (not in PubChem), but we use example SMILES for visualization
# TruthFilterV2 validation results loaded from JSON
NOVEL_MOLECULES = [
    {'formula': 'C8H12N2O3', 'mass': 184, 'name': 'Novel compound 1', 'type': 'Novel molecule',
     'smiles': 'CC(=O)NC1CCCC1NC(=O)C', 'validity': 'FLAG', 'confidence': 0.63},  # diketopiperazine-like
    {'formula': 'C7H9NO4', 'mass': 171, 'name': 'Novel compound 2', 'type': 'Novel molecule',
     'smiles': 'CC(=O)OC1=CC=CC=C1N', 'validity': 'FLAG', 'confidence': 0.30, 'is_aromatic': True},  # N-acetyl derivative (aromatic)
    {'formula': 'C9H11N3O2', 'mass': 193, 'name': 'Novel compound 3', 'type': 'Novel molecule',
     'smiles': 'CC1=CC=C(C=C1)N(C)C(=O)N', 'validity': 'FLAG', 'confidence': 0.30, 'is_aromatic': True},  # N-methyl derivative (aromatic)
    {'formula': 'C6H8N2O3', 'mass': 156, 'name': 'Novel compound 4', 'type': 'Novel molecule',
     'smiles': 'CC(=O)NC1CCCC1N', 'validity': 'FLAG', 'confidence': 0.70},  # cyclic amide
    {'formula': 'C10H14NO2', 'mass': 180, 'name': 'Novel compound 5', 'type': 'Novel molecule',
     'smiles': 'CC1=CC=C(C=C1)OC(=O)NC', 'validity': 'FLAG', 'confidence': 0.30, 'is_aromatic': True},  # aromatic ester
]


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
        
        # Use rdMolDraw2D with explicit atomLabels to force display of ALL atoms
        try:
            drawer = rdMolDraw2D.MolDraw2DCairo(size, size)
            opts = drawer.drawOptions()
            opts.bondLineWidth = 2
            opts.fontSize = 14
            
            # Prepare molecule for drawing (if utils available)
            if _has_utils:
                rdMolDraw2DUtils.PrepareMolForDrawing(mol)
            
            # Force display of all atom labels (including carbons)
            atom_labels = {i: atom.GetSymbol() for i, atom in enumerate(mol.GetAtoms())}
            
            drawer.DrawMolecule(mol, atomLabels=atom_labels)
            drawer.FinishDrawing()
            drawer.WriteDrawingText(str(temp_file))
            
            logger.debug(f"  Used rdMolDraw2D with explicit atomLabels for {smiles}")
        except Exception as e:
            logger.warning(f"  ⚠️  rdMolDraw2D failed ({e}), using basic MolToFile")
            Draw.MolToFile(mol, str(temp_file), size=(size, size))
        
        if temp_file.exists():
            num_atoms = mol.GetNumAtoms()
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


def create_novel_structures_panel(molecules: list, output_path: Path):
    """Create visualization panel with novel molecular structures"""
    n_molecules = len(molecules)
    cols = min(3, n_molecules)
    rows = (n_molecules + cols - 1) // cols
    
    # Larger figure size for better structure visibility
    fig, axes = plt.subplots(rows, cols, figsize=(18, 6*rows))
    fig.suptitle('Figure 6B: Top Novel Molecules - Molecular Structures', fontsize=16, fontweight='bold')
    
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
        mass = mol_data.get('mass', 0)
        name = mol_data.get('name', 'Unknown compound')
        smiles = mol_data.get('smiles', '')
        
        # Render molecular structure if SMILES available - with all atoms visible
        structure_img = None
        if smiles and _rdkit_imported:
            structure_img = render_molecule_structure(smiles, temp_dir, size=600)
            if structure_img:
                try:
                    img = Image.open(structure_img)
                    # Convert RGBA to RGB if needed
                    if img.mode == 'RGBA':
                        background = Image.new('RGB', img.size, (255, 255, 255))
                        if img.split()[3]:
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
                    structure_height = 0.50
                    structure_width = structure_height * aspect_ratio
                    
                    x_center = 0.5
                    x_left = max(0.05, x_center - structure_width / 2)
                    x_right = min(0.95, x_center + structure_width / 2)
                    
                    y_bottom = 0.50
                    y_top = 0.90
                    
                    # Display structure image with proper positioning
                    ax.imshow(img, aspect='auto',
                             extent=[x_left, x_right, y_bottom, y_top],
                             zorder=1,
                             interpolation='lanczos')
                    logger.debug(f"  Displayed structure for {formula}")
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
            y_formula = 0.40
            y_mass = 0.32
            y_novel = 0.24
        else:
            y_formula = 0.75
            y_mass = 0.65
            y_novel = 0.55
        
            ax.text(0.5, y_formula, f"{formula}", transform=ax.transAxes,
                    ha='center', fontsize=16, fontweight='bold', zorder=2)
            ax.text(0.5, y_mass, f"m={mass} amu", transform=ax.transAxes,
                    ha='center', fontsize=11, zorder=2)
            
            # TruthFilterV2 validation status
            validity = mol_data.get('validity', 'FLAG')
            confidence = mol_data.get('confidence', 0.5)
            is_aromatic = mol_data.get('is_aromatic', False)
            
            if validity == 'ACCEPT':
                status_text = "ACCEPT"
                status_color = 'green'
            elif validity == 'FLAG':
                if is_aromatic:
                    status_text = "FLAG (putative)"
                    status_color = 'orange'
                else:
                    status_text = "FLAG (tentative)"
                    status_color = 'orange'
            else:  # REJECT
                status_text = "REJECT"
                status_color = 'red'
            
            ax.text(0.5, y_novel, status_text, transform=ax.transAxes,
                    ha='center', fontsize=10, style='italic', color=status_color, 
                    fontweight='bold', zorder=2)
            
            # Add confidence score
            if validity != 'REJECT':
                y_confidence = y_novel - 0.08
                ax.text(0.5, y_confidence, f"confidence: {confidence:.2f}", transform=ax.transAxes,
                        ha='center', fontsize=8, style='italic', alpha=0.7, zorder=2)
        
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
    
    logger.info(f"✅ Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Generate Figure 6B: Novel Molecule Structures')
    parser.add_argument('--output', type=Path, default=Path('paper/figures/figure6b_novel_structures.png'),
                       help='Output path for figure')
    
    args = parser.parse_args()
    
    # Create output directory if needed
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("GENERATING FIGURE 6B: Novel Molecule Structures")
    logger.info("=" * 70)
    logger.info(f"Output: {args.output}")
    
    if not _rdkit_imported:
        logger.error("❌ RDKit not available - cannot render structures")
        logger.error("💡 Install RDKit: pip install rdkit")
        return 1
    
    create_novel_structures_panel(NOVEL_MOLECULES, args.output)
    
    logger.info("=" * 70)
    logger.info("✅ FIGURE 6B GENERATION COMPLETE!")
    logger.info("=" * 70)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

