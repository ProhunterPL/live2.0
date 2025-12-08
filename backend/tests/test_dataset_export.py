"""
Tests for dataset export module.

Tests basic functionality of dataset export components.
"""

import pytest
import json
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

from backend.dataset_export import DatasetExporter
from backend.dataset_export.utils import (
    resolve_run_patterns,
    validate_run,
    apply_filters,
    load_run_config
)
from backend.dataset_export.snapshot_processor import SnapshotProcessor
from backend.dataset_export.formatters import FormatterFactory


@pytest.fixture
def temp_results_dir():
    """Create temporary results directory structure"""
    with tempfile.TemporaryDirectory() as tmpdir:
        base_dir = Path(tmpdir) / "results" / "phase2b_additional"
        base_dir.mkdir(parents=True)
        
        # Create a sample run
        run_dir = base_dir / "miller_urey_extended" / "run_1"
        run_dir.mkdir(parents=True)
        
        # Create results.json
        results_file = run_dir / "results.json"
        results_file.write_text(json.dumps({
            "scenario": "miller_urey_extended",
            "configuration": {
                "max_steps": 500000,
                "temperature": 298.0
            },
            "molecules_detected": []
        }))
        
        # Create snapshots directory
        snapshots_dir = run_dir / "snapshots"
        snapshots_dir.mkdir()
        
        # Create a sample snapshot
        snapshot_file = snapshots_dir / "step_0.json"
        snapshot_file.write_text(json.dumps({
            "step": 0,
            "time": 0.0,
            "particles": [[0, 0, 0], [1, 1, 1]],
            "bonds": [[0, 1, 1.0]],
            "attributes": [[0, 12.0], [1, 1.0]]
        }))
        
        yield base_dir


def test_resolve_run_patterns(temp_results_dir):
    """Test resolving run patterns"""
    # Test glob pattern
    patterns = ["miller_urey_extended/run_*"]
    runs = resolve_run_patterns(patterns, temp_results_dir)
    assert len(runs) == 1
    assert runs[0].name == "run_1"
    
    # Test direct path
    patterns = ["miller_urey_extended/run_1"]
    runs = resolve_run_patterns(patterns, temp_results_dir)
    assert len(runs) == 1


def test_validate_run(temp_results_dir):
    """Test run validation"""
    run_dir = temp_results_dir / "miller_urey_extended" / "run_1"
    is_valid, error = validate_run(run_dir)
    assert is_valid
    assert error == ""
    
    # Test invalid run
    invalid_dir = temp_results_dir / "invalid_run"
    invalid_dir.mkdir()
    is_valid, error = validate_run(invalid_dir)
    assert not is_valid
    assert "results.json" in error


def test_apply_filters(temp_results_dir):
    """Test applying filters"""
    run_dir = temp_results_dir / "miller_urey_extended" / "run_1"
    runs = [run_dir]
    
    # Test scenario filter
    filters = {"scenario": "miller_urey_extended"}
    filtered = apply_filters(runs, filters)
    assert len(filtered) == 1
    
    filters = {"scenario": "other_scenario"}
    filtered = apply_filters(runs, filters)
    assert len(filtered) == 0
    
    # Test min_steps filter
    filters = {"min_steps": 100000}
    filtered = apply_filters(runs, filters)
    assert len(filtered) == 1
    
    filters = {"min_steps": 1000000}
    filtered = apply_filters(runs, filters)
    assert len(filtered) == 0


def test_load_run_config(temp_results_dir):
    """Test loading run configuration"""
    run_dir = temp_results_dir / "miller_urey_extended" / "run_1"
    config = load_run_config(run_dir)
    assert config["max_steps"] == 500000
    assert config["temperature"] == 298.0


def test_snapshot_processor(temp_results_dir):
    """Test snapshot processor"""
    processor = SnapshotProcessor()
    
    snapshot_file = temp_results_dir / "miller_urey_extended" / "run_1" / "snapshots" / "step_0.json"
    state = processor.load_snapshot(snapshot_file)
    
    assert state["step"] == 0
    assert state["time"] == 0.0
    assert len(state["bonds"]) == 1
    assert len(state["molecules"]) > 0  # Should extract at least one molecule


def test_formatter_factory():
    """Test formatter factory"""
    # Test JSON formatter
    formatter = FormatterFactory.create("json")
    assert formatter is not None
    
    # Test invalid format
    with pytest.raises(ValueError):
        FormatterFactory.create("invalid_format")


def test_dataset_exporter_init(temp_results_dir):
    """Test DatasetExporter initialization"""
    with tempfile.TemporaryDirectory() as tmpdir:
        exporter = DatasetExporter(
            str(temp_results_dir),
            output_dir=str(Path(tmpdir) / "datasets")
        )
        
        assert exporter.base_results_dir == temp_results_dir
        assert exporter.output_dir.exists()


def test_dataset_exporter_export_reaction_trajectories(temp_results_dir):
    """Test exporting reaction trajectories"""
    with tempfile.TemporaryDirectory() as tmpdir:
        exporter = DatasetExporter(
            str(temp_results_dir),
            output_dir=str(Path(tmpdir) / "datasets")
        )
        
        output_path = str(Path(tmpdir) / "test_trajectories.json")
        
        # This will likely produce empty results, but should not crash
        result_path = exporter.export_reaction_trajectories(
            runs=["miller_urey_extended/run_1"],
            output_format="json",
            output_path=output_path
        )
        
        assert result_path == output_path
        assert Path(output_path).exists()


def test_dataset_exporter_export_novel_molecules(temp_results_dir):
    """Test exporting novel molecules"""
    with tempfile.TemporaryDirectory() as tmpdir:
        exporter = DatasetExporter(
            str(temp_results_dir),
            output_dir=str(Path(tmpdir) / "datasets")
        )
        
        output_path = str(Path(tmpdir) / "test_molecules.json")
        
        # This will likely produce empty results, but should not crash
        result_path = exporter.export_novel_molecules(
            runs=["miller_urey_extended/run_1"],
            output_format="json",
            output_path=output_path
        )
        
        assert result_path == output_path
        assert Path(output_path).exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

