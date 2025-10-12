"""
Memory management utilities for Live 2.0
Optimizes Taichi memory usage and prevents memory leaks
"""

import taichi as ti
import gc
import psutil
import os
import logging
from typing import Dict, Any, Optional

logger = logging.getLogger(__name__)

class TaichiMemoryManager:
    """Manages Taichi memory allocation and cleanup"""
    
    def __init__(self):
        self.field_registry: Dict[str, Any] = {}
        self.memory_stats = {
            'peak_memory_mb': 0.0,
            'current_memory_mb': 0.0,
            'taichi_memory_mb': 0.0,
            'gc_collections': 0
        }
        self.cleanup_threshold_mb = 1000  # Cleanup when memory exceeds 1GB
        
    def register_field(self, name: str, field: Any) -> None:
        """Register a Taichi field for tracking"""
        self.field_registry[name] = field
        logger.debug(f"Registered Taichi field: {name}")
    
    def unregister_field(self, name: str) -> None:
        """Unregister a Taichi field"""
        if name in self.field_registry:
            del self.field_registry[name]
            logger.debug(f"Unregistered Taichi field: {name}")
    
    def cleanup_unused_fields(self) -> None:
        """Clean up unused Taichi fields"""
        try:
            # Force garbage collection
            collected = gc.collect()
            self.memory_stats['gc_collections'] += collected
            
            # Clear Taichi cache if available
            if hasattr(ti, 'clear_all_gradients'):
                ti.clear_all_gradients()
            
            logger.debug(f"Memory cleanup: collected {collected} objects")
        except Exception as e:
            logger.warning(f"Memory cleanup failed: {e}")
    
    def get_memory_stats(self) -> Dict[str, float]:
        """Get current memory statistics"""
        try:
            process = psutil.Process(os.getpid())
            memory_info = process.memory_info()
            
            self.memory_stats['current_memory_mb'] = memory_info.rss / 1024 / 1024
            self.memory_stats['peak_memory_mb'] = max(
                self.memory_stats['peak_memory_mb'],
                self.memory_stats['current_memory_mb']
            )
            
            # Estimate Taichi memory usage
            self.memory_stats['taichi_memory_mb'] = self._estimate_taichi_memory()
            
        except Exception as e:
            logger.warning(f"Failed to get memory stats: {e}")
        
        return self.memory_stats.copy()
    
    def _estimate_taichi_memory(self) -> float:
        """Estimate Taichi memory usage"""
        try:
            total_memory = 0.0
            
            for name, field in self.field_registry.items():
                if hasattr(field, 'shape') and hasattr(field, 'dtype'):
                    # Estimate memory usage based on shape and dtype
                    size = 1
                    for dim in field.shape:
                        size *= dim
                    
                    # Estimate bytes per element
                    if field.dtype == ti.f32:
                        bytes_per_element = 4
                    elif field.dtype == ti.f64:
                        bytes_per_element = 8
                    elif field.dtype == ti.i32:
                        bytes_per_element = 4
                    else:
                        bytes_per_element = 4  # Default
                    
                    field_memory = size * bytes_per_element / 1024 / 1024  # MB
                    total_memory += field_memory
            
            return total_memory
        except Exception as e:
            logger.warning(f"Failed to estimate Taichi memory: {e}")
            return 0.0
    
    def should_cleanup(self) -> bool:
        """Check if memory cleanup is needed"""
        stats = self.get_memory_stats()
        return stats['current_memory_mb'] > self.cleanup_threshold_mb
    
    def optimize_memory(self) -> Dict[str, Any]:
        """Perform memory optimization"""
        before_stats = self.get_memory_stats()
        
        # Cleanup unused fields
        self.cleanup_unused_fields()
        
        # Force garbage collection
        gc.collect()
        
        after_stats = self.get_memory_stats()
        
        optimization_result = {
            'before_memory_mb': before_stats['current_memory_mb'],
            'after_memory_mb': after_stats['current_memory_mb'],
            'memory_freed_mb': before_stats['current_memory_mb'] - after_stats['current_memory_mb'],
            'fields_registered': len(self.field_registry),
            'gc_collections': after_stats['gc_collections']
        }
        
        logger.info(f"Memory optimization: freed {optimization_result['memory_freed_mb']:.1f} MB")
        return optimization_result

# Global memory manager instance
memory_manager = TaichiMemoryManager()

def register_taichi_field(name: str, field: Any) -> None:
    """Register a Taichi field with the memory manager"""
    memory_manager.register_field(name, field)

def unregister_taichi_field(name: str) -> None:
    """Unregister a Taichi field from the memory manager"""
    memory_manager.unregister_field(name)

def get_memory_stats() -> Dict[str, float]:
    """Get current memory statistics"""
    return memory_manager.get_memory_stats()

def optimize_memory() -> Dict[str, Any]:
    """Perform memory optimization"""
    return memory_manager.optimize_memory()

def should_cleanup_memory() -> bool:
    """Check if memory cleanup is needed"""
    return memory_manager.should_cleanup()
