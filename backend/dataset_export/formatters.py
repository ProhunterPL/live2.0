"""
Format converters for dataset export.

Supports:
- Parquet (streaming for large datasets)
- JSON (for smaller datasets)
- GraphML (for network graphs)
- CSV (for molecules)
"""

import json
import logging
from pathlib import Path
from typing import List, Dict, Optional, Any

logger = logging.getLogger(__name__)

# Try to import optional dependencies
try:
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
    PARQUET_AVAILABLE = True
except ImportError:
    PARQUET_AVAILABLE = False
    logger.debug("pandas/pyarrow not available, Parquet export disabled (optional)")

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    logger.warning("networkx not available, GraphML export disabled")


class FormatterFactory:
    """Factory for creating formatters"""
    
    @staticmethod
    def create(format_type: str):
        """
        Create formatter for given format type.
        
        Args:
            format_type: "parquet", "json", "graphml", or "csv"
        
        Returns:
            Formatter instance
        """
        if format_type == "parquet":
            if not PARQUET_AVAILABLE:
                raise ImportError("pandas/pyarrow required for Parquet export")
            return ParquetFormatter()
        elif format_type == "json":
            return JSONFormatter()
        elif format_type == "graphml":
            if not NETWORKX_AVAILABLE:
                raise ImportError("networkx required for GraphML export")
            return GraphMLFormatter()
        elif format_type == "csv":
            if not PARQUET_AVAILABLE:
                raise ImportError("pandas required for CSV export")
            return CSVFormatter()
        else:
            raise ValueError(f"Unknown format: {format_type}")


class ParquetFormatter:
    """Streaming Parquet writer for large datasets"""
    
    def open_stream(self, output_path: str):
        """
        Open Parquet file for streaming writes.
        
        Args:
            output_path: Path to output Parquet file
        
        Returns:
            Context manager for streaming writes
        """
        return ParquetStreamWriter(output_path)
    
    def export_graph(self, graph, output_path: str):
        """
        Export NetworkX graph to Parquet (convert to edge list).
        
        Args:
            graph: NetworkX graph
            output_path: Path to output file
        """
        if not NETWORKX_AVAILABLE:
            raise ImportError("networkx required for graph export")
        
        # Convert graph to edge list DataFrame
        edges = []
        for u, v, data in graph.edges(data=True):
            edge_dict = {
                "source": str(u),
                "target": str(v),
            }
            edge_dict.update(data)
            edges.append(edge_dict)
        
        if not edges:
            logger.warning("Empty graph, creating empty Parquet file")
            df = pd.DataFrame(columns=["source", "target"])
        else:
            df = pd.DataFrame(edges)
        
        df.to_parquet(output_path, index=False)
        logger.info(f"Exported graph to Parquet: {output_path} ({len(edges)} edges)")
    
    def export_molecules(self, molecules: List[Dict], output_path: str, include_graphs: bool):
        """
        Export molecules to Parquet.
        
        Args:
            molecules: List of molecule dicts
            output_path: Path to output file
            include_graphs: Whether to include graph data (flattened)
        """
        if not molecules:
            logger.warning("No molecules to export")
            df = pd.DataFrame()
        else:
            # Flatten nested structures
            rows = []
            for mol in molecules:
                row = {
                    "id": mol.get("id", ""),
                    "formula": mol.get("formula", ""),
                    "smiles": mol.get("smiles", ""),
                    "novelty_score": mol.get("novelty_score", 0.0),
                    "complexity": mol.get("complexity", 0),
                    "mass": mol.get("mass", 0.0),
                    "count": mol.get("count", 0),
                }
                
                # Flatten discovery_conditions
                conditions = mol.get("discovery_conditions", {})
                if conditions:
                    row["discovery_scenario"] = conditions.get("scenario", "")
                    row["discovery_run_id"] = conditions.get("run_id", "")
                    row["discovery_step"] = conditions.get("step", 0)
                    row["discovery_temperature"] = conditions.get("temperature", 0.0)
                
                # Flatten properties
                props = mol.get("properties", {})
                if props:
                    row["properties_mass"] = props.get("mass", 0.0)
                    row["properties_count"] = props.get("count", 0)
                    row["properties_first_detected"] = props.get("first_detected", 0.0)
                
                rows.append(row)
            
            df = pd.DataFrame(rows)
        
        df.to_parquet(output_path, index=False)
        logger.info(f"Exported {len(molecules)} molecules to Parquet: {output_path}")


class ParquetStreamWriter:
    """Context manager for streaming Parquet writes"""
    
    def __init__(self, output_path: str):
        self.output_path = Path(output_path)
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        self.records = []
        self.schema = None
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        # Write all records at once (could be optimized for true streaming)
        if self.records:
            df = pd.DataFrame(self.records)
            df.to_parquet(self.output_path, index=False)
            logger.info(f"Wrote {len(self.records)} records to {self.output_path}")
        else:
            logger.warning(f"No records to write to {self.output_path}")
        return False
    
    def write(self, record: Dict):
        """Write a single record (buffered)"""
        self.records.append(record)


class JSONFormatter:
    """JSON formatter (non-streaming for smaller datasets)"""
    
    def open_stream(self, output_path: str):
        """Open JSON file for streaming writes"""
        return JSONStreamWriter(output_path)
    
    def export_graph(self, graph, output_path: str):
        """
        Export NetworkX graph to JSON.
        
        Args:
            graph: NetworkX graph
            output_path: Path to output file
        """
        if not NETWORKX_AVAILABLE:
            raise ImportError("networkx required for graph export")
        
        # Use NetworkX node_link_data format
        data = nx.node_link_data(graph)
        
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        logger.info(f"Exported graph to JSON: {output_path} ({graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges)")
    
    def export_molecules(self, molecules: List[Dict], output_path: str, include_graphs: bool):
        """
        Export molecules to JSON.
        
        Args:
            molecules: List of molecule dicts
            output_path: Path to output file
            include_graphs: Whether to include graph data
        """
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Optionally remove graph data if not needed
        export_molecules = molecules
        if not include_graphs:
            export_molecules = []
            for mol in molecules:
                mol_copy = mol.copy()
                mol_copy.pop("graph", None)  # Safe - won't raise if key doesn't exist
                export_molecules.append(mol_copy)
        
        with open(output_file, 'w') as f:
            json.dump(export_molecules, f, indent=2)
        
        logger.info(f"Exported {len(molecules)} molecules to JSON: {output_path}")


class JSONStreamWriter:
    """Context manager for streaming JSON writes"""
    
    def __init__(self, output_path: str):
        self.output_path = Path(output_path)
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        self.records = []
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        # Write all records at once
        with open(self.output_path, 'w') as f:
            json.dump(self.records, f, indent=2)
        logger.info(f"Wrote {len(self.records)} records to {self.output_path}")
        return False
    
    def write(self, record: Dict):
        """Write a single record (buffered)"""
        self.records.append(record)


class GraphMLFormatter:
    """GraphML formatter for network graphs"""
    
    def export_graph(self, graph, output_path: str):
        """
        Export NetworkX graph to GraphML.
        
        Args:
            graph: NetworkX graph
            output_path: Path to output file
        """
        if not NETWORKX_AVAILABLE:
            raise ImportError("networkx required for GraphML export")
        
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        nx.write_graphml(graph, output_file)
        logger.info(f"Exported graph to GraphML: {output_path} ({graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges)")
    
    def open_stream(self, output_path: str):
        """GraphML doesn't support streaming, use export_graph instead"""
        raise NotImplementedError("GraphML doesn't support streaming, use export_graph()")
    
    def export_molecules(self, molecules: List[Dict], output_path: str, include_graphs: bool):
        """GraphML is for graphs only"""
        raise NotImplementedError("GraphML is for graphs only, use export_graph()")


class CSVFormatter:
    """CSV formatter for molecules"""
    
    def export_molecules(self, molecules: List[Dict], output_path: str, include_graphs: bool):
        """
        Export molecules to CSV.
        
        Args:
            molecules: List of molecule dicts
            output_path: Path to output file
            include_graphs: Whether to include graph data (flattened)
        """
        if not molecules:
            logger.warning("No molecules to export")
            df = pd.DataFrame()
        else:
            # Flatten nested structures
            rows = []
            for mol in molecules:
                row = {
                    "id": mol.get("id", ""),
                    "formula": mol.get("formula", ""),
                    "smiles": mol.get("smiles", ""),
                    "novelty_score": mol.get("novelty_score", 0.0),
                    "complexity": mol.get("complexity", 0),
                    "mass": mol.get("mass", 0.0),
                    "count": mol.get("count", 0),
                }
                
                # Flatten discovery_conditions
                conditions = mol.get("discovery_conditions", {})
                if conditions:
                    row["discovery_scenario"] = conditions.get("scenario", "")
                    row["discovery_run_id"] = conditions.get("run_id", "")
                    row["discovery_step"] = conditions.get("step", 0)
                    row["discovery_temperature"] = conditions.get("temperature", 0.0)
                
                # Flatten properties
                props = mol.get("properties", {})
                if props:
                    row["properties_mass"] = props.get("mass", 0.0)
                    row["properties_count"] = props.get("count", 0)
                    row["properties_first_detected"] = props.get("first_detected", 0.0)
                
                rows.append(row)
            
            df = pd.DataFrame(rows)
        
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        df.to_csv(output_file, index=False)
        logger.info(f"Exported {len(molecules)} molecules to CSV: {output_path}")
    
    def open_stream(self, output_path: str):
        """CSV doesn't support streaming in this implementation"""
        raise NotImplementedError("CSV doesn't support streaming, use export_molecules()")
    
    def export_graph(self, graph, output_path: str):
        """CSV is for molecules only"""
        raise NotImplementedError("CSV is for molecules only, use export_molecules()")

