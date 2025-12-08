"""
Helper functions for dataset export module.

Provides utilities for:
- Resolving run patterns to directories
- Validating run directories
- Progress tracking
"""

import logging
from pathlib import Path
from typing import List, Tuple, Callable, Optional, Dict

logger = logging.getLogger(__name__)


def resolve_run_patterns(patterns: List[str], base_dir: Path) -> List[Path]:
    """
    Resolve run patterns to actual directories.
    
    Examples:
        "miller_urey_extended/run_*" -> all run_* dirs in miller_urey_extended
        "hydrothermal_extended/run_1" -> single run
        "miller_urey_extended/run_1" -> single run
    
    Args:
        patterns: List of pattern strings (can contain wildcards)
        base_dir: Base directory containing results
    
    Returns:
        List of resolved run directories (sorted)
    """
    run_dirs = []
    
    for pattern in patterns:
        if "*" in pattern:
            # Glob pattern
            matches = list(base_dir.glob(pattern))
            if not matches:
                logger.warning(f"No matches for pattern: {pattern}")
            run_dirs.extend(matches)
        else:
            # Direct path
            run_dir = base_dir / pattern
            if run_dir.exists() and run_dir.is_dir():
                run_dirs.append(run_dir)
            else:
                logger.warning(f"Run directory not found: {run_dir}")
    
    # Remove duplicates and sort
    run_dirs = sorted(set(run_dirs))
    logger.info(f"Resolved {len(patterns)} pattern(s) to {len(run_dirs)} run directory(ies)")
    
    return run_dirs


def validate_run(run_dir: Path) -> Tuple[bool, str]:
    """
    Validate that a run directory has required files.
    
    Args:
        run_dir: Path to run directory
    
    Returns:
        (is_valid, error_message)
    """
    if not run_dir.exists():
        return False, f"Directory does not exist: {run_dir}"
    
    if not run_dir.is_dir():
        return False, f"Not a directory: {run_dir}"
    
    # Check for required files
    required = ["results.json"]
    for req in required:
        req_path = run_dir / req
        if not req_path.exists():
            return False, f"Missing required file: {req}"
    
    # Check for optional but important files
    optional = ["snapshots/", "molecules.json", "reaction_network.json"]
    missing_optional = []
    for opt in optional:
        opt_path = run_dir / opt
        if not opt_path.exists():
            missing_optional.append(opt)
    
    if missing_optional:
        logger.debug(f"Run {run_dir.name} missing optional files: {missing_optional}")
    
    return True, ""


def apply_filters(run_dirs: List[Path], filters: Optional[Dict]) -> List[Path]:
    """
    Apply filters to run directories.
    
    Supported filters:
        - scenario: str - filter by scenario name
        - min_steps: int - minimum steps in results.json
        - temperature_range: Tuple[float, float] - temperature range
    
    Args:
        run_dirs: List of run directories
        filters: Optional dict of filters
    
    Returns:
        Filtered list of run directories
    """
    if not filters:
        return run_dirs
    
    filtered = []
    
    for run_dir in run_dirs:
        # Check scenario filter
        if "scenario" in filters:
            scenario_name = run_dir.parent.name
            if scenario_name != filters["scenario"]:
                continue
        
        # Check min_steps filter
        if "min_steps" in filters:
            results_file = run_dir / "results.json"
            if results_file.exists():
                try:
                    import json
                    with open(results_file) as f:
                        data = json.load(f)
                    config = data.get("configuration", {})
                    max_steps = config.get("max_steps", 0)
                    if max_steps < filters["min_steps"]:
                        continue
                except Exception as e:
                    logger.warning(f"Failed to check min_steps for {run_dir.name}: {e}")
                    continue
        
        # Check temperature range filter
        if "temperature_range" in filters:
            temp_range = filters["temperature_range"]
            
            # Validate temperature_range format
            if not isinstance(temp_range, (tuple, list)) or len(temp_range) != 2:
                logger.error(f"Invalid temperature_range format: {temp_range}. Must be tuple/list of 2 elements.")
                continue
            
            if temp_range[0] > temp_range[1]:
                logger.error(f"Invalid temperature_range: {temp_range[0]} > {temp_range[1]}. First element must be <= second.")
                continue
            
            results_file = run_dir / "results.json"
            if results_file.exists():
                try:
                    import json
                    with open(results_file) as f:
                        data = json.load(f)
                    config = data.get("configuration", {})
                    temperature = config.get("temperature", 298.0)
                    if not (temp_range[0] <= temperature <= temp_range[1]):
                        continue
                except Exception as e:
                    logger.warning(f"Failed to check temperature for {run_dir.name}: {e}")
                    continue
        
        filtered.append(run_dir)
    
    logger.info(f"Applied filters: {len(run_dirs)} -> {len(filtered)} runs")
    return filtered


def create_progress_callback(total: int) -> Callable:
    """
    Create a progress callback that logs to console.
    
    Args:
        total: Total number of items to process
    
    Returns:
        Callback function(current, total, message)
    """
    def callback(current: int, total_items: int, message: str):
        if total_items > 0:
            pct = 100 * current / total_items
            logger.info(f"[{pct:.1f}%] {message} ({current}/{total_items})")
        else:
            logger.info(f"[?] {message}")
    
    return callback


def load_run_config(run_dir: Path) -> Dict:
    """
    Load configuration from results.json.
    
    Args:
        run_dir: Path to run directory
    
    Returns:
        Dict with configuration (empty dict if not found)
    """
    results_file = run_dir / "results.json"
    if not results_file.exists():
        logger.warning(f"No results.json in {run_dir.name}")
        return {}
    
    try:
        import json
        with open(results_file) as f:
            data = json.load(f)
        return data.get("configuration", {})
    except Exception as e:
        logger.error(f"Failed to load config from {run_dir.name}: {e}")
        return {}

