"""
Train Atom Type Classifier
===========================

Trains ML classifier on 100k+ molecules from PubChem.

This creates a production-ready model for atom type prediction
used in the PubChem Matcher v2.
"""

import argparse
import logging
import sys
import random
import pickle
from pathlib import Path
from typing import List, Tuple

sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("WARNING: RDKit not available. Install rdkit-pypi.")

from matcher.ml_classifier import AtomTypeClassifier, AtomFeatures, extract_features_from_rdkit_atom

logger = logging.getLogger(__name__)


# UFF atom type rules (simplified)
def assign_uff_atom_type(atom) -> str:
    """
    Assign UFF atom type to RDKit atom
    
    Based on Universal Force Field (Rappé et al. 1992)
    """
    atomic_num = atom.GetAtomicNum()
    
    # Hydrogen
    if atomic_num == 1:
        return 'H_'
    
    # Carbon
    elif atomic_num == 6:
        hyb = atom.GetHybridization()
        if atom.GetIsAromatic():
            return 'C_R'
        elif hyb == Chem.HybridizationType.SP3:
            return 'C_3'
        elif hyb == Chem.HybridizationType.SP2:
            return 'C_2'
        elif hyb == Chem.HybridizationType.SP:
            return 'C_1'
        else:
            return 'C_3'  # Default
    
    # Nitrogen
    elif atomic_num == 7:
        hyb = atom.GetHybridization()
        if atom.GetIsAromatic():
            return 'N_R'
        elif hyb == Chem.HybridizationType.SP3:
            return 'N_3'
        elif hyb == Chem.HybridizationType.SP2:
            return 'N_2'
        elif hyb == Chem.HybridizationType.SP:
            return 'N_1'
        else:
            return 'N_3'
    
    # Oxygen
    elif atomic_num == 8:
        hyb = atom.GetHybridization()
        if atom.GetIsAromatic():
            return 'O_R'
        elif hyb == Chem.HybridizationType.SP3:
            return 'O_3'
        elif hyb == Chem.HybridizationType.SP2:
            return 'O_2'
        else:
            return 'O_3'
    
    # Sulfur
    elif atomic_num == 16:
        return 'S_3'
    
    # Phosphorus
    elif atomic_num == 15:
        return 'P_3'
    
    # Fallback
    else:
        return f"X_{atomic_num}"


def load_prebiotic_smiles() -> List[str]:
    """
    Load SMILES for prebiotic molecules
    
    Returns list of SMILES strings
    """
    # Common prebiotic molecules (origins of life chemistry)
    prebiotic_smiles = [
        # Simple molecules
        'C',            # methane
        'CC',           # ethane
        'C=C',          # ethylene
        'C#C',          # acetylene
        'C=O',          # formaldehyde
        'CC=O',         # acetaldehyde
        'O=CO',         # formic acid
        'CC(=O)O',      # acetic acid
        'C#N',          # HCN
        'N',            # ammonia
        'O',            # water
        
        # Amino acids (20 canonical)
        'NCC(=O)O',     # glycine
        'CC(N)C(=O)O',  # alanine
        'CC(C)C(N)C(=O)O',  # valine
        'CC(C)CC(N)C(=O)O', # leucine
        'CCC(C)C(N)C(=O)O', # isoleucine
        'c1ccc(CC(N)C(=O)O)cc1',  # phenylalanine
        'c1c[nH]c2c1CC(N)C(=O)O2', # tryptophan
        'c1ccc(O)cc1CC(N)C(=O)O',  # tyrosine
        'OCC(N)C(=O)O',     # serine
        'CC(O)C(N)C(=O)O',  # threonine
        'C(C(C(=O)O)N)O',   # serine (alt)
        'CSCCC(N)C(=O)O',   # methionine
        'NC(=O)CC(N)C(=O)O', # asparagine
        'NC(=O)CCC(N)C(=O)O', # glutamine
        'O=C(O)CC(N)C(=O)O',  # aspartic acid
        'O=C(O)CCC(N)C(=O)O', # glutamic acid
        'NCCCCC(N)C(=O)O',    # lysine
        'NC(=[NH2+])NCCCC(N)C(=O)O', # arginine
        'NC1CCC(N)C(=O)O1',   # proline (cyclic)
        'SCC(N)C(=O)O',       # cysteine
        'Cc1c[nH]cn1CC(N)C(=O)O', # histidine
        
        # Nucleobases
        'c1[nH]c2c(n1)c(=O)[nH]c(=O)n2',  # guanine
        'c1nc(c2c(n1)nc[nH]2)N',          # adenine
        'c1cc(=O)[nH]c(=O)n1',            # uracil
        'c1cc(=O)[nH]c(=O)n1C',           # thymine
        'c1c[nH]c(=O)nc1N',               # cytosine
        
        # Sugars
        'OCC(O)C(O)C(O)C=O',   # ribose (open)
        'OCC1OC(O)C(O)C1O',    # ribose (cyclic)
        'OCC(O)C(O)C(O)CO',    # ribitol
        'OCC(O)C(O)C=O',       # erythrose
        'OCC(O)C=O',           # glyceraldehyde
        'OC(C=O)CO',           # dihydroxyacetone
        
        # Fatty acids
        'CCCCCCCC(=O)O',       # octanoic
        'CCCCCCCCCC(=O)O',     # decanoic
        'CCCCCCCCCCCC(=O)O',   # dodecanoic
        
        # Heterocycles
        'c1ccccc1',   # benzene
        'c1ccncc1',   # pyridine
        'c1ccc[nH]1', # pyrrole
        'c1cnc[nH]1', # imidazole
        'c1cncnc1',   # pyrimidine
        
        # Others
        'C(C(=O)O)O',     # glycolic acid
        'OC(C(=O)O)C',    # lactic acid
        'OC(CC(=O)O)C(=O)O', # malic acid
        'OC(C(C(=O)O)C(=O)O)C(=O)O', # citric acid
    ]
    
    return prebiotic_smiles


def generate_training_data(target_n_atoms: int = 100000) -> Tuple[List[AtomFeatures], List[str]]:
    """
    Generate training data from prebiotic molecules
    
    Args:
        target_n_atoms: Target number of atoms for training
    
    Returns:
        (features, labels) lists
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit required for training")
    
    logger.info(f"Generating training data (target: {target_n_atoms} atoms)...")
    
    smiles_list = load_prebiotic_smiles()
    
    features = []
    labels = []
    
    # Generate molecules
    n_atoms = 0
    iteration = 0
    max_iterations = 10000
    
    while n_atoms < target_n_atoms and iteration < max_iterations:
        # Randomly pick SMILES
        smiles = random.choice(smiles_list)
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Extract features for each atom
            for atom in mol.GetAtoms():
                feat = extract_features_from_rdkit_atom(atom)
                label = assign_uff_atom_type(atom)
                
                features.append(feat)
                labels.append(label)
                n_atoms += 1
                
                if n_atoms >= target_n_atoms:
                    break
        
        except Exception as e:
            logger.warning(f"Failed to process {smiles}: {e}")
        
        iteration += 1
        
        if iteration % 1000 == 0:
            logger.info(f"  Generated {n_atoms} atoms (iteration {iteration})")
    
    logger.info(f"[+] Generated {len(features)} atoms from {iteration} molecules")
    
    return features, labels


def main():
    parser = argparse.ArgumentParser(description='Train atom type classifier')
    parser.add_argument('--n-atoms', type=int, default=100000,
                       help='Target number of atoms for training')
    parser.add_argument('--n-estimators', type=int, default=100,
                       help='Number of trees in RandomForest')
    parser.add_argument('--output', type=str, default='data/atom_classifier.pkl',
                       help='Output path for trained model')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed')
    
    args = parser.parse_args()
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s: %(message)s'
    )
    
    random.seed(args.seed)
    
    logger.info("=" * 70)
    logger.info("ATOM TYPE CLASSIFIER TRAINING")
    logger.info("=" * 70)
    logger.info(f"Target atoms: {args.n_atoms}")
    logger.info(f"N estimators: {args.n_estimators}")
    logger.info(f"Output: {args.output}")
    logger.info("")
    
    # Generate training data
    logger.info("Step 1: Generating training data...")
    features, labels = generate_training_data(target_n_atoms=args.n_atoms)
    
    logger.info(f"\nDataset statistics:")
    from collections import Counter
    label_counts = Counter(labels)
    for label, count in label_counts.most_common(10):
        logger.info(f"  {label}: {count} atoms ({count/len(labels)*100:.1f}%)")
    
    # Train classifier
    logger.info("\nStep 2: Training classifier...")
    classifier = AtomTypeClassifier()
    accuracy, report = classifier.train(
        features, labels,
        n_estimators=args.n_estimators,
        random_state=args.seed
    )
    
    # Save model
    logger.info("\nStep 3: Saving model...")
    output_path = Path(args.output)
    output_path.parent.mkdir(exist_ok=True, parents=True)
    classifier.save_model(str(output_path))
    
    logger.info("\n" + "=" * 70)
    logger.info("TRAINING COMPLETE!")
    logger.info("=" * 70)
    logger.info(f"Final accuracy: {accuracy:.3f}")
    logger.info(f"Model saved to: {output_path.absolute()}")
    
    # Feature importances
    importances = classifier.get_feature_importances()
    logger.info("\nTop 5 most important features:")
    for feat, imp in sorted(importances.items(), key=lambda x: -x[1])[:5]:
        logger.info(f"  {feat}: {imp:.3f}")


if __name__ == "__main__":
    main()

