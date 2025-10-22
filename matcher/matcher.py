"""
LIVE 2.0 Cluster to PubChem Matcher - CLI tool.

Usage:
    python matcher/matcher.py <path_to_cluster_image.png>

Requirements:
    - A JSON file with the same name must exist (e.g., cluster.json for cluster.png)
    - JSON must contain 'nodes' and 'bonds' data

Output:
    - matches/<name>_match.png: Comparison panel
    - matches/<name>_match.json: Match metadata
"""
import json
import sys
from pathlib import Path
import traceback

# Fix encoding for Windows console
if sys.platform == 'win32':
    import io
    if sys.stdout.encoding != 'utf-8':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    if sys.stderr.encoding != 'utf-8':
        sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

from chem import (
    json_to_mol, mol_to_smiles, pubchem_similar_top, 
    render_mol_png, render_pubchem_png, export_all_formats
)
from compose import compose_panel_with_metadata


def run_for_pair(img_path: str):
    """
    Process a cluster image and its metadata JSON.
    
    Steps:
    1. Load cluster metadata from JSON
    2. Convert cluster to RDKit Mol → SMILES
    3. Query PubChem for similar compounds
    4. Render both molecules
    5. Compose comparison panel
    6. Save results
    """
    p = Path(img_path)
    
    # Check for metadata JSON
    meta_path = p.with_suffix(".json")
    if not meta_path.exists():
        print(f"❌ Error: Missing metadata JSON next to image: {meta_path}")
        sys.exit(1)
    
    print(f"📂 Processing: {p.name}")
    print(f"📄 Metadata: {meta_path.name}")
    
    # Load metadata
    try:
        with open(meta_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except Exception as e:
        print(f"❌ Error loading JSON: {e}")
        sys.exit(1)
    
    # Validate JSON structure
    if "nodes" not in data or "bonds" not in data:
        print("❌ Error: JSON must contain 'nodes' and 'bonds' fields")
        sys.exit(1)
    
    print(f"   Nodes: {len(data['nodes'])}, Bonds: {len(data['bonds'])}")
    
    # Convert to RDKit Mol
    try:
        mol = json_to_mol(data)
        smiles = mol_to_smiles(mol)
        print(f"🧪 SMILES: {smiles}")
    except Exception as e:
        print(f"❌ Error converting to SMILES: {e}")
        traceback.print_exc()
        sys.exit(1)
    
    # Export to standard chemical formats (.mol and .xyz)
    out_dir = Path("matches")
    out_dir.mkdir(exist_ok=True)
    base_path = out_dir / p.stem
    print("📦 Exporting to chemical file formats...")
    export_all_formats(mol, str(base_path), stamp=data.get("id", p.stem))
    
    # Render LIVE cluster as RDKit molecule
    tmp_left = p.parent / f"{p.stem}_live_rdkit.png"
    try:
        render_mol_png(mol, str(tmp_left), size=512)
        print(f"✓ Rendered LIVE cluster: {tmp_left.name}")
    except Exception as e:
        print(f"❌ Error rendering LIVE cluster: {e}")
        sys.exit(1)
    
    # Query PubChem
    print("🔍 Searching PubChem...")
    try:
        top = pubchem_similar_top(smiles, threshold=80)
    except Exception as e:
        print(f"❌ Error querying PubChem: {e}")
        top = None
    
    # Handle PubChem results
    if top and top.get("smiles"):
        print(f"✓ Found match: CID {top['cid']} - {top.get('name', 'unknown')}")
        pubchem_png = p.parent / f"{p.stem}_pubchem.png"
        try:
            render_pubchem_png(top["smiles"], str(pubchem_png), size=512)
            print(f"✓ Rendered PubChem molecule: {pubchem_png.name}")
            
            # Export PubChem molecule to chemical formats as well
            from rdkit import Chem
            pubchem_mol = Chem.MolFromSmiles(top["smiles"])
            if pubchem_mol:
                pubchem_base = out_dir / f"{p.stem}_pubchem"
                print("📦 Exporting PubChem molecule to chemical formats...")
                export_all_formats(pubchem_mol, str(pubchem_base), stamp=f"PubChem_CID_{top['cid']}")
        except Exception as e:
            print(f"⚠️  Warning: Could not render PubChem molecule: {e}")
            pubchem_png = tmp_left  # Fallback to live cluster
            top = {}
    else:
        print("⚠️  No PubChem match found (threshold=80)")
        pubchem_png = tmp_left  # Fallback
        top = {}
    
    # Compose comparison panel
    panel_path = out_dir / f"{p.stem}_match.png"
    
    try:
        if top:
            compose_panel_with_metadata(
                str(tmp_left),
                str(pubchem_png),
                data.get("metadata", {}),
                top,
                str(panel_path)
            )
        else:
            # No match case
            from compose import compose_panel
            compose_panel(
                str(tmp_left),
                str(tmp_left),
                "LIVE 2.0 Cluster",
                "No PubChem Match",
                str(panel_path),
                f"Size: {data.get('metadata', {}).get('size', '?')}",
                "Try adjusting threshold"
            )
    except Exception as e:
        print(f"❌ Error composing panel: {e}")
        traceback.print_exc()
        sys.exit(1)
    
    # Save result metadata
    result = {
        "cluster_id": data.get("id", p.stem),
        "smiles_live": smiles,
        "pubchem_top": top if top else None,
        "panel_path": str(panel_path),
        "cluster_metadata": data.get("metadata", {})
    }
    
    result_json_path = out_dir / f"{p.stem}_match.json"
    with open(result_json_path, 'w', encoding='utf-8') as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    
    print(f"✓ Saved results: {result_json_path.name}")
    print(f"✅ Done! Panel: {panel_path}")
    
    # Print result as JSON for programmatic use
    print("\n" + "="*60)
    print(json.dumps(result, ensure_ascii=False, indent=2))


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python matcher/matcher.py <path_to_cluster_image.png>")
        print("\nExample:")
        print("  python matcher/matcher.py clustershots/2025-10-11_123456.png")
        sys.exit(1)
    
    img_path = sys.argv[1]
    
    if not Path(img_path).exists():
        print(f"❌ Error: Image file not found: {img_path}")
        sys.exit(1)
    
    run_for_pair(img_path)


if __name__ == "__main__":
    main()

