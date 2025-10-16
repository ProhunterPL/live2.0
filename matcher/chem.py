"""
RDKit and PubChem utilities for LIVE 2.0 cluster-to-molecule matching.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import requests
import urllib.parse
from typing import Optional, Dict, List


def choose_symbol(deg: int, given_label: Optional[str] = None, mass: float = 1.0, charge_magnitude: float = 0.0) -> str:
    """
    Heuristic mapping from node degree, mass, and charge to atomic symbol.
    
    CRITICAL: Hydrogen can ONLY have degree 1 (one bond). Any atom with degree > 1
    cannot be hydrogen, regardless of mass.
    
    If given_label is a valid element symbol, use it.
    Otherwise, use degree as PRIMARY criterion:
    - degree 1 → could be H (if terminal) or terminal C
    - degree 2 → O, N, or C in chain (NEVER H)
    - degree 3 → N or C (NEVER H)
    - degree ≥4 → C (NEVER H)
    
    Mass is used as secondary criterion:
    - mass ≈ 12 (10-14) → C (carbon)
    - mass ≈ 14 (13-15) → N (nitrogen)
    - mass ≈ 16 (15-17) → O (oxygen)
    - mass ≈ 1 (0.5-2) → H only if degree == 1
    
    This is a proxy mapping for clusters without real chemical symbols.
    """
    valid_elements = {"H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"}
    
    if given_label and given_label in valid_elements:
        return given_label
    
    # CRITICAL: Hydrogen can ONLY have degree 1
    # Any atom with degree > 1 CANNOT be hydrogen
    if deg > 1:
        # Use mass to choose between C, N, O for multi-bonded atoms
        if 10.0 <= mass <= 14.0:
            return "C"
        elif 13.0 <= mass <= 15.0:
            return "N"
        elif 15.0 <= mass <= 17.0:
            return "O"
        
        # Fallback based on degree for multi-bonded atoms
        if deg == 2:
            # In organic chemistry, degree 2 is often O in chains (e.g., ethers, carbonyls)
            # But could also be C in chains or N in rings
            return "O" if charge_magnitude < 0.3 else "N"
        elif deg == 3:
            # Degree 3 is typically N (amines) or C (sp2/sp3)
            return "N" if charge_magnitude > 0.3 else "C"
        else:  # deg >= 4
            # Only carbon can have 4 or more bonds in organic chemistry
            return "C"
    
    # Degree == 1: could be H (terminal hydrogen) or terminal C
    # Use mass to decide
    if 0.5 <= mass <= 2.0:
        # Low mass, degree 1 → likely hydrogen
        return "H"
    elif 10.0 <= mass <= 14.0:
        # Carbon mass, degree 1 → terminal carbon (e.g., methyl group)
        return "C"
    elif 15.0 <= mass <= 17.0:
        # Oxygen mass, degree 1 → hydroxyl or similar
        return "O"
    
    # Default fallback for degree 1
    # In organic molecules, terminal atoms are usually C or O
    return "C"


def json_to_mol(cluster_json: dict) -> Chem.Mol:
    """
    Convert cluster JSON (nodes + bonds) to RDKit Mol object.
    
    Expected JSON structure:
    {
        "nodes": [{"id": 0, "label": "A", ...}, ...],
        "bonds": [{"a": 0, "b": 2, "order": 1}, ...]
    }
    """
    nodes = cluster_json["nodes"]
    bonds = cluster_json["bonds"]
    
    # Calculate node degrees
    deg = {n["id"]: 0 for n in nodes}
    for b in bonds:
        deg[b["a"]] += 1
        deg[b["b"]] += 1
    
    rw = Chem.RWMol()
    idx_map = {}
    
    # Add atoms
    for n in nodes:
        # Extract mass and charge for better element selection
        mass = n.get("mass", 1.0)
        charge = n.get("charge", [0.0, 0.0, 0.0])
        charge_magnitude = (charge[0]**2 + charge[1]**2 + charge[2]**2)**0.5
        
        sym = choose_symbol(deg[n["id"]], n.get("label"), mass, charge_magnitude)
        idx_map[n["id"]] = rw.AddAtom(Chem.Atom(sym))
    
    # Add bonds
    order_map = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE
    }
    for b in bonds:
        bond_order = b.get("order", 1)
        bond_type = order_map.get(bond_order, Chem.BondType.SINGLE)
        rw.AddBond(idx_map[b["a"]], idx_map[b["b"]], bond_type)
    
    mol = rw.GetMol()
    
    # Try to sanitize (kekulize, etc.)
    try:
        AllChem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    except Exception as e:
        print(f"Warning: Sanitization failed: {e}")
        # Continue anyway - some molecules might not be perfect
    
    # CRITICAL: Add implicit hydrogens to satisfy valence requirements
    # Without this, N3 becomes N3 instead of N3H3, C2 becomes C2 instead of C2H6, etc.
    try:
        mol = Chem.AddHs(mol)
    except Exception as e:
        print(f"Warning: Adding hydrogens failed: {e}")
    
    return mol


def mol_to_smiles(mol: Chem.Mol) -> str:
    """
    Convert RDKit Mol to canonical SMILES string.
    
    CRITICAL: Remove explicit hydrogens before generating SMILES for PubChem search.
    PubChem uses implicit hydrogens:
    - [H]N([H])[H] → N (ammonia)
    - [H][H] → [HH] or just don't search
    - [H]N1N([H])N1[H] → N1NN1
    """
    try:
        # Remove explicit hydrogens for PubChem-compatible SMILES
        mol_no_h = Chem.RemoveHs(mol)
        return Chem.MolToSmiles(mol_no_h, canonical=True)
    except Exception as e:
        print(f"Warning: SMILES generation failed: {e}")
        # Fallback: try with explicit H
        return Chem.MolToSmiles(mol, canonical=True)


def pubchem_similar_top(smiles: str, threshold: int = 75) -> Optional[Dict]:
    """
    Search PubChem for the most similar compound to the given SMILES.
    
    Returns a dict with:
    - cid: PubChem Compound ID
    - name: IUPAC name
    - formula: Molecular formula
    - mw: Molecular weight
    - smiles: Canonical SMILES from PubChem
    - inchikey: InChIKey
    
    Returns None if no match found.
    """
    try:
        q = urllib.parse.quote(smiles)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/{q}/JSON?Threshold={threshold}"
        
        r = requests.get(url, timeout=25)
        r.raise_for_status()
        
        data = r.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        
        if not cids:
            return None
        
        cid = cids[0]
        
        # Get properties for the top result
        props_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/"
            f"property/IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES,InChIKey/JSON"
        )
        p_response = requests.get(props_url, timeout=25)
        p_response.raise_for_status()
        
        props = p_response.json()["PropertyTable"]["Properties"][0]
        
        return {
            "cid": cid,
            "name": props.get("IUPACName", "unknown"),
            "formula": props.get("MolecularFormula", ""),
            "mw": props.get("MolecularWeight", 0),
            "smiles": props.get("CanonicalSMILES", ""),
            "inchikey": props.get("InChIKey", "")
        }
    
    except requests.exceptions.RequestException as e:
        print(f"PubChem API error: {e}")
        return None
    except Exception as e:
        print(f"Error querying PubChem: {e}")
        return None


def render_pubchem_png(smiles: str, out_png: str, size: int = 512):
    """
    Render a PNG image of a molecule from its SMILES string.
    
    Args:
        smiles: SMILES string
        out_png: Output PNG file path
        size: Image size (width and height)
    """
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    
    Draw.MolToFile(mol, out_png, size=(size, size))


def render_mol_png(mol: Chem.Mol, out_png: str, size: int = 512):
    """
    Render a PNG image of an RDKit Mol object.
    
    Args:
        mol: RDKit Mol object
        out_png: Output PNG file path
        size: Image size (width and height)
    """
    Draw.MolToFile(mol, out_png, size=(size, size))


def export_to_mol_file(mol: Chem.Mol, out_mol: str):
    """
    Export molecule to MDL Molfile format (.mol).
    
    This format is widely used in chemistry software and databases.
    
    Args:
        mol: RDKit Mol object
        out_mol: Output .mol file path
    """
    try:
        Chem.MolToMolFile(mol, out_mol)
        print(f"✓ Exported to MOL format: {out_mol}")
    except Exception as e:
        print(f"⚠️  Warning: Could not export to MOL format: {e}")


def export_to_xyz_file(mol: Chem.Mol, out_xyz: str, stamp: str = ""):
    """
    Export molecule to XYZ format with 3D coordinates.
    
    The XYZ format is a simple Cartesian coordinate format widely used
    in computational chemistry. This function:
    1. Generates 3D coordinates using ETKDG method
    2. Writes atom symbols and positions
    
    Args:
        mol: RDKit Mol object
        out_xyz: Output .xyz file path
        stamp: Optional timestamp/identifier for file header
    """
    try:
        # Generate 3D coordinates using ETKDG (Experimental Torsion Knowledge Distance Geometry)
        # This is a fast, rule-based 3D coordinate generation method
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        if result == -1:
            # Embedding failed, try without ETKDG parameters
            print("⚠️  ETKDG embedding failed, trying basic embedding...")
            result = AllChem.EmbedMolecule(mol)
            
        if result == -1:
            print("⚠️  Warning: Could not generate 3D coordinates, skipping XYZ export")
            return
        
        # Optimize the geometry (optional but recommended)
        AllChem.UFFOptimizeMolecule(mol)
        
        # Write XYZ file
        with open(out_xyz, 'w', encoding='utf-8') as f:
            conf = mol.GetConformer()
            num_atoms = mol.GetNumAtoms()
            
            # XYZ format:
            # Line 1: Number of atoms
            # Line 2: Comment line
            # Following lines: Element X Y Z
            f.write(f"{num_atoms}\n")
            f.write(f"LIVE 2.0 export {stamp}\n")
            
            for atom_idx in range(num_atoms):
                atom = mol.GetAtomWithIdx(atom_idx)
                symbol = atom.GetSymbol()
                pos = conf.GetAtomPosition(atom_idx)
                
                f.write(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")
        
        print(f"✓ Exported to XYZ format: {out_xyz}")
        
    except Exception as e:
        print(f"⚠️  Warning: Could not export to XYZ format: {e}")


def export_all_formats(mol: Chem.Mol, base_path: str, stamp: str = ""):
    """
    Export molecule to all supported formats.
    
    This is a convenience function that exports to:
    - .mol (MDL Molfile)
    - .xyz (XYZ Cartesian coordinates)
    
    Args:
        mol: RDKit Mol object
        base_path: Base path without extension (e.g., "matches/cluster_123")
        stamp: Optional timestamp/identifier
    """
    from pathlib import Path
    
    base = Path(base_path)
    
    # Export to MOL format
    mol_path = base.with_suffix('.mol')
    export_to_mol_file(mol, str(mol_path))
    
    # Export to XYZ format
    xyz_path = base.with_suffix('.xyz')
    export_to_xyz_file(mol, str(xyz_path), stamp=stamp or base.stem)

