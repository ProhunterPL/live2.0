"""
Snapshot Processor for dataset export.

Processes simulation snapshots to extract molecules and states.
Used by ReactionTrajectoryExporter to compare consecutive snapshots.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Set, Tuple
from collections import defaultdict

logger = logging.getLogger(__name__)


class SnapshotProcessor:
    """
    Processes simulation snapshots to extract molecules and states.
    
    Used by ReactionTrajectoryExporter to compare consecutive snapshots.
    """
    
    def load_snapshot(self, snapshot_file: Path) -> Dict:
        """
        Load snapshot and extract molecular state.
        
        Args:
            snapshot_file: Path to snapshot JSON file
        
        Returns:
            Dict with:
            - step: int
            - time: float
            - particles: List (positions, types, etc.)
            - bonds: List[List] (bond data)
            - molecules: List[Dict] (detected molecules)
        """
        try:
            with open(snapshot_file, 'r') as f:
                data = json.load(f)
        except Exception as e:
            logger.error(f"Failed to load snapshot {snapshot_file.name}: {e}")
            return {
                "step": 0,
                "time": 0.0,
                "particles": [],
                "bonds": [],
                "molecules": []
            }
        
        # Extract molecules from bonds (connected components)
        bonds = data.get("bonds", [])
        attributes = data.get("attributes", [])
        molecules = self._extract_molecules_from_bonds(bonds, attributes)
        
        return {
            "step": data.get("step", 0),
            "time": data.get("time", 0.0),
            "particles": data.get("particles", []),
            "bonds": bonds,
            "molecules": molecules
        }
    
    def _extract_molecules_from_bonds(
        self, 
        bonds: List[List], 
        attributes: List
    ) -> List[Dict]:
        """
        Extract molecules from bond graph (connected components).
        
        Similar to build_reaction_network_from_snapshots.extract_molecules_from_snapshot()
        
        Args:
            bonds: List of bond data (each bond is [i, j, ...] or similar)
            attributes: List of particle attributes
        
        Returns:
            List of molecule dicts with formula, atoms, etc.
        """
        if not bonds:
            return []
        
        # Build graph
        graph = defaultdict(set)
        for bond in bonds:
            if len(bond) >= 2:
                i, j = int(bond[0]), int(bond[1])
                graph[i].add(j)
                graph[j].add(i)
        
        if not graph:
            return []
        
        # Find connected components
        visited = set()
        molecules = []
        
        for node in graph:
            if node in visited:
                continue
            
            # BFS to find all atoms in this molecule
            component = []
            stack = [node]
            while stack:
                n = stack.pop()
                if n not in visited:
                    visited.add(n)
                    component.append(n)
                    stack.extend(graph[n])
            
            if len(component) >= 2:  # At least 2 atoms for a molecule
                # Create molecule dict
                molecule = self._create_molecule(component, graph, attributes)
                molecules.append(molecule)
        
        return molecules
    
    def _create_molecule(
        self, 
        component: List[int], 
        graph: Dict[int, Set[int]], 
        attributes: List
    ) -> Dict:
        """
        Create molecule dict from connected component.
        
        Args:
            component: List of atom indices in molecule
            graph: Bond graph
            attributes: Particle attributes
        
        Returns:
            Molecule dict with formula, atoms, bonds, etc.
        """
        # Get bonds within this molecule
        molecule_bonds = []
        for atom_i in component:
            for atom_j in graph.get(atom_i, []):
                if atom_j in component:
                    bond = tuple(sorted([atom_i, atom_j]))
                    if bond not in molecule_bonds:
                        molecule_bonds.append(bond)
        
        # Calculate formula (simplified - count atom types)
        atom_types = []
        total_mass = 0.0
        
        for idx in component:
            if idx < len(attributes):
                attr = attributes[idx]
                if isinstance(attr, list) and len(attr) > 0:
                    atom_type = int(attr[0]) if isinstance(attr[0], (int, float)) else 0
                    mass = float(attr[1]) if len(attr) > 1 else 1.0
                else:
                    atom_type = 0
                    mass = 1.0
            else:
                atom_type = 0
                mass = 1.0
            
            atom_types.append(atom_type)
            total_mass += mass
        
        # Generate formula string
        formula = self._generate_formula(atom_types)
        
        return {
            "formula": formula,
            "atoms": component,
            "bonds": molecule_bonds,
            "mass": total_mass,
            "size": len(component)
        }
    
    def _generate_formula(self, atom_types: List[int]) -> str:
        """
        Generate molecular formula from atom types.
        
        Simplified implementation - maps atom types to element symbols.
        
        Args:
            atom_types: List of atom type numbers
        
        Returns:
            Formula string (e.g., "C2H4O2")
        """
        # Count each atom type
        from collections import Counter
        type_counts = Counter(atom_types)
        
        # Map type to symbol (simplified - would use full periodic table)
        type_to_symbol = {
            0: "C",
            1: "H",
            2: "O",
            3: "N",
            4: "S",
            5: "P"
        }
        
        # Build formula (conventionally: C, H, then others)
        formula_parts = []
        
        # Add C first if present
        if 0 in type_counts:
            count = type_counts[0]
            if count == 1:
                formula_parts.append("C")
            else:
                formula_parts.append(f"C{count}")
        
        # Add H second if present
        if 1 in type_counts:
            count = type_counts[1]
            if count == 1:
                formula_parts.append("H")
            else:
                formula_parts.append(f"H{count}")
        
        # Add other elements
        for atom_type in sorted(type_counts.keys()):
            if atom_type in (0, 1):  # Already added
                continue
            count = type_counts[atom_type]
            symbol = type_to_symbol.get(atom_type, f"X{atom_type}")
            if count == 1:
                formula_parts.append(symbol)
            else:
                formula_parts.append(f"{symbol}{count}")
        
        return "".join(formula_parts) if formula_parts else "UNKNOWN"

