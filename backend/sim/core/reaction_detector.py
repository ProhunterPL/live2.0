"""
Reaction Detector Module
=========================

Detects chemical reactions in simulation trajectories by analyzing:
- Bond formation/breaking events
- Molecular connectivity changes
- Product formation

This is crucial for validating benchmark reactions against literature data.
"""

import numpy as np
import logging
from typing import List, Dict, Set, Tuple, Optional
from dataclasses import dataclass, field
from collections import defaultdict

logger = logging.getLogger(__name__)


@dataclass
class Molecule:
    """Represents a detected molecule"""
    atom_indices: Set[int]
    formula: str
    mass: float
    bonds: Set[Tuple[int, int]]
    
    def __hash__(self):
        return hash((self.formula, frozenset(self.atom_indices)))
    
    def __eq__(self, other):
        return self.formula == other.formula and self.atom_indices == other.atom_indices


@dataclass
class ReactionEvent:
    """Represents a detected chemical reaction"""
    step: int
    time: float
    reactants: List[Molecule]
    products: List[Molecule]
    reaction_type: str = "unknown"
    delta_G: Optional[float] = None
    
    def __repr__(self):
        reactant_formulas = " + ".join([m.formula for m in self.reactants])
        product_formulas = " + ".join([m.formula for m in self.products])
        return f"t={self.time:.2f}: {reactant_formulas} -> {product_formulas}"


class BondGraph:
    """
    Graph representation of molecular bonds
    
    Used for detecting molecules via connected components
    """
    
    def __init__(self, n_atoms: int):
        self.n_atoms = n_atoms
        self.adjacency: Dict[int, Set[int]] = defaultdict(set)
    
    def add_bond(self, atom_i: int, atom_j: int):
        """Add a bond between atoms i and j"""
        self.adjacency[atom_i].add(atom_j)
        self.adjacency[atom_j].add(atom_i)
    
    def get_connected_components(self) -> List[Set[int]]:
        """
        Find all molecules (connected components)
        
        Returns list of sets, each set contains atom indices in one molecule
        """
        visited = set()
        components = []
        
        for start_atom in range(self.n_atoms):
            if start_atom in visited:
                continue
            
            # BFS to find all atoms in this molecule
            component = set()
            queue = [start_atom]
            
            while queue:
                atom = queue.pop(0)
                if atom in visited:
                    continue
                
                visited.add(atom)
                component.add(atom)
                
                # Add neighbors
                for neighbor in self.adjacency.get(atom, []):
                    if neighbor not in visited:
                        queue.append(neighbor)
            
            if component:
                components.append(component)
        
        return components
    
    def get_bonds_in_molecule(self, atom_set: Set[int]) -> Set[Tuple[int, int]]:
        """Get all bonds within a molecule"""
        bonds = set()
        for atom_i in atom_set:
            for atom_j in self.adjacency.get(atom_i, []):
                if atom_j in atom_set:
                    bond = tuple(sorted([atom_i, atom_j]))
                    bonds.add(bond)
        return bonds


class ReactionDetector:
    """
    Detects chemical reactions in simulation trajectories
    
    Usage:
        detector = ReactionDetector()
        
        for step in simulation:
            detector.update(step, positions, bonds, attributes)
        
        reactions = detector.get_detected_reactions()
    """
    
    def __init__(self, 
                 bond_threshold: float = 1.8,  # Angstrom
                 min_reaction_time: float = 0.1):  # ps
        """
        Initialize reaction detector
        
        Args:
            bond_threshold: Distance threshold for bond detection (Angstrom)
            min_reaction_time: Minimum time between reactions (ps)
        """
        self.bond_threshold = bond_threshold
        self.min_reaction_time = min_reaction_time
        
        self.reactions: List[ReactionEvent] = []
        self.molecule_history: List[List[Molecule]] = []
        self.previous_molecules: Optional[List[Molecule]] = None
        self.last_reaction_time = -np.inf
        
        # Statistics
        self.total_reactions_detected = 0
        self.reaction_types_count: Dict[str, int] = defaultdict(int)
    
    def detect_molecules(self, 
                        n_atoms: int,
                        positions: np.ndarray,
                        bonds: List[Tuple[int, int]],
                        attributes: np.ndarray) -> List[Molecule]:
        """
        Detect all molecules in current state
        
        Args:
            n_atoms: Number of atoms
            positions: Atom positions (n_atoms, 3)
            bonds: List of (i, j) bond pairs
            attributes: Atom attributes (n_atoms, 5) [type, mass, charge, ...]
        
        Returns:
            List of detected molecules
        """
        # Build bond graph
        graph = BondGraph(n_atoms)
        for atom_i, atom_j in bonds:
            graph.add_bond(atom_i, atom_j)
        
        # Find connected components (molecules)
        components = graph.get_connected_components()
        
        molecules = []
        for component in components:
            if not component:
                continue
            
            # Get molecule properties
            atom_indices = component
            bonds_in_mol = graph.get_bonds_in_molecule(atom_indices)
            
            # Calculate formula (simplified - counts atoms)
            atom_types = []
            total_mass = 0.0
            for idx in atom_indices:
                atom_type = int(attributes[idx, 0])
                mass = attributes[idx, 1]
                atom_types.append(atom_type)
                total_mass += mass
            
            # Generate formula string (e.g. "C2H4O2")
            formula = self._generate_formula(atom_types, attributes[list(atom_indices)])
            
            molecule = Molecule(
                atom_indices=atom_indices,
                formula=formula,
                mass=total_mass,
                bonds=bonds_in_mol
            )
            molecules.append(molecule)
        
        return molecules
    
    def _generate_formula(self, atom_types: List[int], attributes: np.ndarray) -> str:
        """
        Generate molecular formula from atom types
        
        Note: This is simplified. Full implementation would map atom types
        to element symbols (C, H, O, N, etc.)
        """
        # Count each atom type
        type_counts = defaultdict(int)
        for atom_type in atom_types:
            type_counts[atom_type] += 1
        
        # Sort by type (conventionally: C, H, then others)
        # For now, just sort by type number
        formula_parts = []
        for atom_type in sorted(type_counts.keys()):
            count = type_counts[atom_type]
            # Map type to symbol (simplified - would use full periodic table)
            symbol = self._atom_type_to_symbol(atom_type)
            if count == 1:
                formula_parts.append(symbol)
            else:
                formula_parts.append(f"{symbol}{count}")
        
        return "".join(formula_parts)
    
    def _atom_type_to_symbol(self, atom_type: int) -> str:
        """
        Map atom type to element symbol
        
        Note: Placeholder - full implementation would use proper element mapping
        """
        # Simplified mapping (extend as needed)
        type_to_symbol = {
            0: "C",
            1: "H",
            2: "O",
            3: "N",
            4: "S",
            5: "P"
        }
        return type_to_symbol.get(atom_type, f"X{atom_type}")
    
    def update(self, 
               step: int,
               time: float,
               n_atoms: int,
               positions: np.ndarray,
               bonds: List[Tuple[int, int]],
               attributes: np.ndarray):
        """
        Update detector with new simulation state
        
        Detects reactions by comparing current molecules with previous state
        """
        # Detect molecules in current state
        current_molecules = self.detect_molecules(n_atoms, positions, bonds, attributes)
        
        # Store molecule history
        self.molecule_history.append(current_molecules)
        
        # Check for reactions (if we have previous state)
        if self.previous_molecules is not None:
            self._detect_reactions(step, time, self.previous_molecules, current_molecules)
        
        self.previous_molecules = current_molecules
    
    def _detect_reactions(self,
                         step: int,
                         time: float,
                         previous: List[Molecule],
                         current: List[Molecule]):
        """
        Detect reactions by comparing previous and current molecular states
        """
        # Only detect if enough time has passed
        if time - self.last_reaction_time < self.min_reaction_time:
            return
        
        # Convert to sets for comparison
        prev_set = set(previous)
        curr_set = set(current)
        
        # Check for changes
        disappeared = prev_set - curr_set  # Reactants
        appeared = curr_set - prev_set      # Products
        
        if disappeared or appeared:
            # Reaction detected!
            reaction = ReactionEvent(
                step=step,
                time=time,
                reactants=list(disappeared),
                products=list(appeared),
                reaction_type=self._classify_reaction(disappeared, appeared)
            )
            
            self.reactions.append(reaction)
            self.total_reactions_detected += 1
            self.reaction_types_count[reaction.reaction_type] += 1
            self.last_reaction_time = time
            
            logger.debug(f"Reaction detected at step {step}: {reaction}")
    
    def _classify_reaction(self, reactants: Set[Molecule], products: Set[Molecule]) -> str:
        """
        Classify reaction type based on reactants and products
        
        Types:
        - condensation: A + B -> C
        - dissociation: A -> B + C
        - isomerization: A -> B (same formula)
        - complex: other
        """
        n_reactants = len(reactants)
        n_products = len(products)
        
        if n_reactants == 2 and n_products == 1:
            return "condensation"
        elif n_reactants == 1 and n_products == 2:
            return "dissociation"
        elif n_reactants == 1 and n_products == 1:
            # Check if same formula (isomerization)
            r = list(reactants)[0]
            p = list(products)[0]
            if r.formula == p.formula:
                return "isomerization"
            else:
                return "rearrangement"
        else:
            return "complex"
    
    def get_detected_reactions(self) -> List[ReactionEvent]:
        """Get all detected reactions"""
        return self.reactions
    
    def get_reactions_by_type(self, reaction_type: str) -> List[ReactionEvent]:
        """Get reactions of specific type"""
        return [r for r in self.reactions if r.reaction_type == reaction_type]
    
    def get_product_yields(self, target_formula: str) -> Dict[str, float]:
        """
        Calculate yields for a specific product
        
        Args:
            target_formula: Chemical formula to search for (e.g. "C2H4O2")
        
        Returns:
            Dict with statistics: count, first_appearance, etc.
        """
        count = 0
        first_appearance = None
        last_appearance = None
        
        for reaction in self.reactions:
            for product in reaction.products:
                if product.formula == target_formula:
                    count += 1
                    if first_appearance is None:
                        first_appearance = reaction.time
                    last_appearance = reaction.time
        
        return {
            'count': count,
            'first_appearance': first_appearance,
            'last_appearance': last_appearance,
            'formula': target_formula
        }
    
    def count_molecules(self, formula: str, at_step: Optional[int] = None) -> int:
        """
        Count molecules with given formula
        
        Args:
            formula: Chemical formula
            at_step: Specific step to count at (default: final step)
        
        Returns:
            Number of molecules
        """
        if at_step is None:
            at_step = len(self.molecule_history) - 1
        
        if at_step >= len(self.molecule_history):
            return 0
        
        molecules = self.molecule_history[at_step]
        return sum(1 for m in molecules if m.formula == formula)
    
    def get_statistics(self) -> Dict:
        """Get detector statistics"""
        return {
            'total_reactions': self.total_reactions_detected,
            'reaction_types': dict(self.reaction_types_count),
            'snapshots': len(self.molecule_history),
            'unique_molecules': len(set(
                m.formula for molecules in self.molecule_history for m in molecules
            ))
        }
    
    def summary(self) -> str:
        """Generate human-readable summary"""
        stats = self.get_statistics()
        
        lines = []
        lines.append("=" * 70)
        lines.append("REACTION DETECTOR SUMMARY")
        lines.append("=" * 70)
        lines.append(f"Total reactions detected: {stats['total_reactions']}")
        lines.append(f"Snapshots analyzed: {stats['snapshots']}")
        lines.append(f"Unique molecules observed: {stats['unique_molecules']}")
        lines.append("")
        lines.append("Reaction types:")
        for rtype, count in stats['reaction_types'].items():
            lines.append(f"  {rtype}: {count}")
        lines.append("=" * 70)
        
        return "\n".join(lines)


# Example usage
if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    
    # Create detector
    detector = ReactionDetector()
    
    # Simulate some data
    print("Testing ReactionDetector with mock data...")
    
    # Step 1: Two separate molecules (C + H)
    n_atoms_1 = 2
    positions_1 = np.array([[0, 0, 0], [3, 0, 0]])  # Far apart
    bonds_1 = []  # No bonds
    attributes_1 = np.array([
        [0, 12.0, 0, 0, 0],  # C
        [1, 1.0, 0, 0, 0]    # H
    ])
    
    detector.update(0, 0.0, n_atoms_1, positions_1, bonds_1, attributes_1)
    
    # Step 2: Bonded (C-H formed)
    bonds_2 = [(0, 1)]  # C-H bond
    detector.update(1, 0.5, n_atoms_1, positions_1, bonds_2, attributes_1)
    
    # Print summary
    print("\n" + detector.summary())
    
    # Print reactions
    print("\nDetected reactions:")
    for reaction in detector.get_detected_reactions():
        print(f"  {reaction}")

