"""
Catalog system for Live 2.0 simulation
Tracks novel substances, their properties, and emergence patterns
"""

import taichi as ti
import numpy as np
import json
import time
from typing import Dict, List, Tuple, Optional, Any
from collections import defaultdict
from .graphs import MolecularGraph, GraphCatalog

class SubstanceRecord:
    """Record of a discovered substance"""
    
    def __init__(self, graph: MolecularGraph, timestamp: float, 
                 properties: Dict[str, Any] = None):
        self.graph = graph
        self.timestamp = timestamp
        self.properties = properties or {}
        
        # Generate unique ID
        self.id = self._generate_id()
        
        # Track occurrences
        self.occurrence_count = 1
        self.first_seen = timestamp
        self.last_seen = timestamp
        
        # Stability metrics
        self.lifetime_sum = 0.0
        self.max_lifetime = 0.0
        self.min_lifetime = float('inf')
    
    def _generate_id(self) -> str:
        """Generate unique ID for this substance"""
        canonical_form = self.graph.get_canonical_form()
        timestamp_str = str(int(self.timestamp * 1000))  # Millisecond precision
        return f"SUB_{canonical_form[:8]}_{timestamp_str}"
    
    def update_occurrence(self, timestamp: float, lifetime: float = 0.0):
        """Update occurrence statistics"""
        self.occurrence_count += 1
        self.last_seen = timestamp
        
        # Update lifetime statistics
        self.lifetime_sum += lifetime
        self.max_lifetime = max(self.max_lifetime, lifetime)
        self.min_lifetime = min(self.min_lifetime, lifetime)
    
    def get_average_lifetime(self) -> float:
        """Get average lifetime of this substance"""
        if self.occurrence_count == 0:
            return 0.0
        return self.lifetime_sum / self.occurrence_count
    
    def to_dict(self) -> Dict:
        """Convert to dictionary representation"""
        return {
            'id': self.id,
            'graph': self.graph.to_dict(),
            'timestamp': self.timestamp,
            'properties': self.properties,
            'occurrence_count': self.occurrence_count,
            'first_seen': self.first_seen,
            'last_seen': self.last_seen,
            'lifetime_stats': {
                'sum': self.lifetime_sum,
                'max': self.max_lifetime,
                'min': self.min_lifetime,
                'average': self.get_average_lifetime()
            }
        }
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'SubstanceRecord':
        """Create from dictionary representation"""
        graph = MolecularGraph.from_dict(data['graph'])
        record = cls(graph, data['timestamp'], data['properties'])
        record.id = data['id']
        record.occurrence_count = data['occurrence_count']
        record.first_seen = data['first_seen']
        record.last_seen = data['last_seen']
        
        lifetime_stats = data['lifetime_stats']
        record.lifetime_sum = lifetime_stats['sum']
        record.max_lifetime = lifetime_stats['max']
        record.min_lifetime = lifetime_stats['min']
        
        return record

class SubstanceCatalog:
    """Catalog of all discovered substances"""
    
    def __init__(self):
        self.substances: Dict[str, SubstanceRecord] = {}
        self.graph_catalog = GraphCatalog()
        
        # Statistics
        self.total_discoveries = 0
        self.novel_discoveries = 0
        self.start_time = time.time()
        
        # Temporal tracking with size limits to prevent memory leaks
        self.discovery_timeline: List[Tuple[float, str]] = []
        self.novelty_rate_history: List[Tuple[float, float]] = []
        
        # Maximum number of timeline entries to keep (prevents memory leaks)
        self.max_timeline_entries = 1000  # Reduced from 10000
        self.max_history_entries = 100    # Reduced from 1000
    
    def add_substance(self, graph: MolecularGraph, timestamp: float = None,
                     properties: Dict[str, Any] = None) -> Tuple[bool, str]:
        """Add a substance to the catalog. Returns (is_novel, substance_id)"""
        if timestamp is None:
            timestamp = time.time()
        
        # Check if this is a novel graph
        is_novel_graph = self.graph_catalog.add_graph(graph, timestamp)
        
        # Create or update substance record
        canonical_form = graph.get_canonical_form()
        
        if canonical_form in self.substances:
            # Update existing substance
            record = self.substances[canonical_form]
            record.update_occurrence(timestamp)
            is_novel = False
            substance_id = record.id
        else:
            # Create new substance record
            record = SubstanceRecord(graph, timestamp, properties)
            self.substances[canonical_form] = record
            is_novel = True
            substance_id = record.id
        
        # Update statistics
        self.total_discoveries += 1
        if is_novel:
            self.novel_discoveries += 1
        
        # Update timeline with memory management
        self.discovery_timeline.append((timestamp, substance_id))
        if len(self.discovery_timeline) > self.max_timeline_entries:
            # Remove oldest entries to prevent memory leaks
            self.discovery_timeline = self.discovery_timeline[-self.max_timeline_entries:]
        
        # Update novelty rate history with memory management
        novelty_rate = self.get_novelty_rate()
        self.novelty_rate_history.append((timestamp, novelty_rate))
        if len(self.novelty_rate_history) > self.max_history_entries:
            # Remove oldest entries to prevent memory leaks
            self.novelty_rate_history = self.novelty_rate_history[-self.max_history_entries:]
        
        return is_novel, substance_id
    
    def get_novelty_rate(self, window_size: int = 100) -> float:
        """Get current novelty rate"""
        if self.total_discoveries < window_size:
            return self.novel_discoveries / max(self.total_discoveries, 1)
        
        # Count novel discoveries in recent window
        recent_discoveries = self.discovery_timeline[-window_size:]
        recent_novel = sum(1 for _, substance_id in recent_discoveries 
                         if self._is_recent_novel(substance_id))
        
        return recent_novel / window_size
    
    def _is_recent_novel(self, substance_id: str) -> bool:
        """Check if a substance was novel when first discovered"""
        # This is a simplified check - in practice, you'd track this more precisely
        for canonical_form, record in self.substances.items():
            if record.id == substance_id:
                return record.occurrence_count == 1
        return False
    
    def get_substance_by_id(self, substance_id: str) -> Optional[SubstanceRecord]:
        """Get substance record by ID"""
        for record in self.substances.values():
            if record.id == substance_id:
                return record
        return None
    
    def get_recent_substances(self, count: int = 10) -> List[SubstanceRecord]:
        """Get most recently discovered substances"""
        recent_discoveries = self.discovery_timeline[-count:]
        substances = []
        
        for _, substance_id in recent_discoveries:
            substance = self.get_substance_by_id(substance_id)
            if substance:
                substances.append(substance)
        
        return substances
    
    def get_substances_by_size(self, min_size: int = None, max_size: int = None) -> List[SubstanceRecord]:
        """Get substances filtered by size"""
        filtered = []
        
        for record in self.substances.values():
            size = record.graph.num_nodes
            
            if min_size is not None and size < min_size:
                continue
            if max_size is not None and size > max_size:
                continue
            
            filtered.append(record)
        
        return filtered
    
    def get_most_common_substances(self, count: int = 10) -> List[SubstanceRecord]:
        """Get most frequently occurring substances"""
        sorted_substances = sorted(self.substances.values(), 
                                 key=lambda x: x.occurrence_count, 
                                 reverse=True)
        return sorted_substances[:count]
    
    def get_catalog_stats(self) -> Dict:
        """Get comprehensive catalog statistics"""
        if not self.substances:
            return {
                'total_substances': 0,
                'novel_substances': 0,
                'total_discoveries': 0,
                'novelty_rate': 0.0,
                'average_substance_size': 0.0,
                'discovery_rate': 0.0,
                'runtime': 0.0
            }
        
        # Basic statistics
        total_substances = len(self.substances)
        novel_substances = self.novel_discoveries
        total_discoveries = self.total_discoveries
        novelty_rate = self.get_novelty_rate()
        
        # Size statistics
        sizes = [record.graph.num_nodes for record in self.substances.values()]
        average_size = sum(sizes) / len(sizes)
        
        # Temporal statistics
        runtime = time.time() - self.start_time
        discovery_rate = total_discoveries / max(runtime, 1)
        
        # Complexity statistics
        complexities = [record.graph.density for record in self.substances.values()]
        average_complexity = sum(complexities) / len(complexities)
        
        # Lifetime statistics
        lifetimes = [record.get_average_lifetime() for record in self.substances.values()]
        average_lifetime = sum(lifetimes) / len(lifetimes)
        
        return {
            'total_substances': int(total_substances),
            'novel_substances': int(novel_substances),
            'total_discoveries': int(total_discoveries),
            'novelty_rate': float(novelty_rate),
            'average_substance_size': float(average_size),
            'discovery_rate': float(discovery_rate),
            'runtime': float(runtime),
            'average_complexity': float(average_complexity),
            'average_lifetime': float(average_lifetime),
            'size_distribution': self._get_size_distribution(),
            'complexity_distribution': self._get_complexity_distribution()
        }
    
    def _get_size_distribution(self) -> Dict[int, int]:
        """Get distribution of substance sizes"""
        size_counts = defaultdict(int)
        for record in self.substances.values():
            size_counts[record.graph.num_nodes] += 1
        return dict(size_counts)
    
    def _get_complexity_distribution(self) -> Dict[str, int]:
        """Get distribution of substance complexities"""
        complexity_bins = {
            'low': 0,      # density < 0.3
            'medium': 0,   # 0.3 <= density < 0.7
            'high': 0      # density >= 0.7
        }
        
        for record in self.substances.values():
            density = record.graph.density
            if density < 0.3:
                complexity_bins['low'] += 1
            elif density < 0.7:
                complexity_bins['medium'] += 1
            else:
                complexity_bins['high'] += 1
        
        return complexity_bins
    
    def export_catalog(self, filename: str):
        """Export catalog to JSON file"""
        data = {
            'substances': {canonical_form: record.to_dict() 
                          for canonical_form, record in self.substances.items()},
            'statistics': self.get_catalog_stats(),
            'discovery_timeline': self.discovery_timeline,
            'novelty_rate_history': self.novelty_rate_history
        }
        
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)
    
    def import_catalog(self, filename: str):
        """Import catalog from JSON file"""
        with open(filename, 'r') as f:
            data = json.load(f)
        
        # Clear existing data
        self.substances.clear()
        self.graph_catalog.clear()
        self.discovery_timeline.clear()
        self.novelty_rate_history.clear()
        
        # Import substances
        for canonical_form, substance_data in data['substances'].items():
            record = SubstanceRecord.from_dict(substance_data)
            self.substances[canonical_form] = record
            
            # Rebuild graph catalog
            self.graph_catalog.add_graph(record.graph, record.timestamp)
        
        # Import timeline and history
        self.discovery_timeline = data.get('discovery_timeline', [])
        self.novelty_rate_history = data.get('novelty_rate_history', [])
        
        # Update statistics
        self.total_discoveries = len(self.discovery_timeline)
        self.novel_discoveries = sum(1 for _, substance_id in self.discovery_timeline 
                                   if self._is_recent_novel(substance_id))
    
    def get_novel_substances(self, count: int = 10) -> List[SubstanceRecord]:
        """Get recent novel substances"""
        return self.get_recent_substances(count)
    
    def cleanup_old_data(self, max_age_hours: float = 24.0):
        """Clean up old data to prevent memory leaks"""
        current_time = time.time()
        cutoff_time = current_time - (max_age_hours * 3600)
        
        # Clean up old timeline entries
        self.discovery_timeline = [
            (timestamp, substance_id) 
            for timestamp, substance_id in self.discovery_timeline 
            if timestamp > cutoff_time
        ]
        
        # Clean up old novelty rate history
        self.novelty_rate_history = [
            (timestamp, rate) 
            for timestamp, rate in self.novelty_rate_history 
            if timestamp > cutoff_time
        ]
        
        # Clean up old substances that haven't been seen recently
        substances_to_remove = []
        for canonical_form, record in self.substances.items():
            if record.last_seen < cutoff_time and record.occurrence_count == 1:
                substances_to_remove.append(canonical_form)
        
        for canonical_form in substances_to_remove:
            del self.substances[canonical_form]
    
    def get_memory_usage_stats(self) -> Dict:
        """Get memory usage statistics for monitoring"""
        return {
            'timeline_entries': len(self.discovery_timeline),
            'history_entries': len(self.novelty_rate_history),
            'substances_count': len(self.substances),
            'max_timeline_entries': self.max_timeline_entries,
            'max_history_entries': self.max_history_entries,
            'memory_pressure': len(self.discovery_timeline) / self.max_timeline_entries
        }
    
    def clear(self):
        """Clear all catalog data"""
        self.substances.clear()
        self.graph_catalog.clear()
        self.discovery_timeline.clear()
        self.novelty_rate_history.clear()
        self.total_discoveries = 0
        self.novel_discoveries = 0
        self.start_time = time.time()
    
    def get_catalog_stats(self) -> Dict:
        """Get comprehensive catalog statistics"""
        current_time = time.time()
        runtime_hours = (current_time - self.start_time) / 3600
        
        return {
            'total_substances': len(self.substances),
            'total_novel': self.novel_discoveries,
            'total_discovered': self.total_discoveries,
            'novelty_rate': self.get_novelty_rate(),
            'novelty_rate_100': self.get_novelty_rate(100),
            'novelty_rate_1000': self.get_novelty_rate(1000),
            'runtime_hours': runtime_hours,
            'discoveries_per_hour': self.total_discoveries / max(runtime_hours, 0.01),
            'novel_per_hour': self.novel_discoveries / max(runtime_hours, 0.01)
        }
    
    def clear(self):
        """Clear the catalog"""
        self.substances.clear()
        self.graph_catalog.clear()
        self.total_discoveries = 0
        self.novel_discoveries = 0
        self.discovery_timeline.clear()
        self.novelty_rate_history.clear()
        self.start_time = time.time()
