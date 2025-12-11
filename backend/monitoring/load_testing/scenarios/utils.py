"""
Utilities for load testing.
"""

import json
from typing import Dict, Optional


def load_test_config(config_file: Optional[str] = None) -> Dict:
    """
    Load test configuration from file.
    
    Args:
        config_file: Path to config file (JSON)
    
    Returns:
        Config dict
    """
    if config_file and os.path.exists(config_file):
        with open(config_file, "r") as f:
            return json.load(f)
    
    # Default config
    return {
        "legally": {
            "users": [100, 500, 1000, 5000],
            "spawn_rate": 10,
            "ramp_up_time": 60,
            "duration": 300,  # 5 minutes
            "target_rps": 1000
        },
        "live2": {
            "users": [50, 200, 500],
            "spawn_rate": 5,
            "ramp_up_time": 120,
            "duration": 600,  # 10 minutes
            "target_rps": 200
        }
    }


def generate_load_test_report(results_file: str, output_file: str):
    """
    Generate load test report from results.
    
    Args:
        results_file: Path to Locust results file
        output_file: Path to output report file
    """
    # TODO: Implement report generation
    # This would parse Locust results and generate HTML/PDF report
    pass
