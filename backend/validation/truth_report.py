"""
Truth Report Generator: Generowanie raportów walidacji
"""
import json
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
from datetime import datetime
from backend.validation.types import TruthResult, FilterStatus

logger = logging.getLogger(__name__)


class TruthReportGenerator:
    """
    Generate validation reports in JSON and Markdown formats
    """
    
    def __init__(self, output_dir: Optional[str] = None):
        """
        Initialize report generator
        
        Args:
            output_dir: Output directory for reports (optional)
        """
        self.output_dir = Path(output_dir) if output_dir else None
    
    def generate_report(self, result: TruthResult) -> Dict[str, Any]:
        """
        Generate report for single run
        
        Args:
            result: TruthResult from validation
            
        Returns:
            Dict with report data
        """
        report = {
            'run_id': result.run_id,
            'run_path': result.run_path,
            'validation_level': result.validation_level.value,
            'overall_status': result.overall_status.value,
            'timestamp': result.timestamp or datetime.now().timestamp(),
            'filters': {},
            'filtered_results': result.filtered_results,
            'summary': result.summary
        }
        
        # Add filter results
        for filter_name, filter_result in result.filters.items():
            report['filters'][filter_name] = {
                'status': filter_result.status.value,
                'details': filter_result.details,
                'warnings': filter_result.warnings,
                'errors': filter_result.errors
            }
        
        return report
    
    def generate_summary(self, results: List[TruthResult]) -> Dict[str, Any]:
        """
        Generate summary statistics for multiple runs
        
        Args:
            results: List of TruthResult objects
            
        Returns:
            Dict with summary statistics
        """
        if not results:
            return {'error': 'No results to summarize'}
        
        total_runs = len(results)
        passed_runs = sum(1 for r in results if r.is_pass())
        warning_runs = sum(1 for r in results if r.overall_status == FilterStatus.WARNING)
        failed_runs = sum(1 for r in results if r.overall_status == FilterStatus.FAIL)
        
        # Aggregate filter statistics
        filter_stats = {}
        for result in results:
            for filter_name, filter_result in result.filters.items():
                if filter_name not in filter_stats:
                    filter_stats[filter_name] = {
                        'pass': 0,
                        'warning': 0,
                        'fail': 0
                    }
                filter_stats[filter_name][filter_result.status.value.lower()] += 1
        
        # Aggregate molecule statistics
        total_molecules = 0
        total_filtered = 0
        for result in results:
            mol_filter = result.filters.get('molecule_filter')
            if mol_filter:
                total_molecules += mol_filter.details.get('total_original', 0)
                total_filtered += mol_filter.details.get('total_filtered', 0)
        
        retention_rate = total_filtered / total_molecules if total_molecules > 0 else 0.0
        
        summary = {
            'total_runs': total_runs,
            'passed_runs': passed_runs,
            'warning_runs': warning_runs,
            'failed_runs': failed_runs,
            'pass_rate': passed_runs / total_runs if total_runs > 0 else 0.0,
            'filter_statistics': filter_stats,
            'molecule_statistics': {
                'total_original': total_molecules,
                'total_filtered': total_filtered,
                'retention_rate': retention_rate
            },
            'timestamp': datetime.now().timestamp()
        }
        
        return summary
    
    def save_json(self, report: Dict[str, Any], filename: str) -> Path:
        """
        Save report as JSON
        
        Args:
            report: Report dict
            filename: Output filename
            
        Returns:
            Path to saved file
        """
        if self.output_dir:
            output_path = self.output_dir / filename
            self.output_dir.mkdir(parents=True, exist_ok=True)
        else:
            output_path = Path(filename)
        
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        logger.info(f"Saved JSON report to {output_path}")
        return output_path
    
    def save_markdown(self, result: TruthResult, filename: str) -> Path:
        """
        Save report as Markdown
        
        Args:
            result: TruthResult
            filename: Output filename
            
        Returns:
            Path to saved file
        """
        if self.output_dir:
            output_path = self.output_dir / filename
            self.output_dir.mkdir(parents=True, exist_ok=True)
        else:
            output_path = Path(filename)
        
        lines = []
        lines.append(f"# Truth-Filter Validation Report: {result.run_id}")
        lines.append("")
        lines.append(f"**Validation Level**: {result.validation_level.value}")
        lines.append(f"**Overall Status**: {result.overall_status.value}")
        lines.append(f"**Run Path**: `{result.run_path}`")
        lines.append("")
        
        # Status icon
        status_icon = {
            FilterStatus.PASS: "✅",
            FilterStatus.WARNING: "⚠️",
            FilterStatus.FAIL: "❌"
        }
        lines.append(f"## {status_icon.get(result.overall_status, '?')} Overall: {result.overall_status.value}")
        lines.append("")
        
        # Filter results
        lines.append("## Filter Results")
        lines.append("")
        for filter_name, filter_result in result.filters.items():
            icon = status_icon.get(filter_result.status, '?')
            lines.append(f"### {icon} {filter_name.replace('_', ' ').title()}")
            lines.append("")
            lines.append(f"**Status**: {filter_result.status.value}")
            lines.append("")
            
            # Details
            if filter_result.details:
                lines.append("**Details**:")
                for key, value in filter_result.details.items():
                    if isinstance(value, float):
                        lines.append(f"- {key}: {value:.3f}")
                    else:
                        lines.append(f"- {key}: {value}")
                lines.append("")
            
            # Warnings
            if filter_result.warnings:
                lines.append("**Warnings**:")
                for warning in filter_result.warnings:
                    lines.append(f"- ⚠️ {warning}")
                lines.append("")
            
            # Errors
            if filter_result.errors:
                lines.append("**Errors**:")
                for error in filter_result.errors:
                    lines.append(f"- ❌ {error}")
                lines.append("")
        
        # Summary
        if result.summary:
            lines.append("## Summary")
            lines.append("")
            for key, value in result.summary.items():
                if isinstance(value, float):
                    lines.append(f"- **{key}**: {value:.3f}")
                else:
                    lines.append(f"- **{key}**: {value}")
            lines.append("")
        
        # Write file
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        
        logger.info(f"Saved Markdown report to {output_path}")
        return output_path
    
    def save_summary_markdown(self, summary: Dict[str, Any], filename: str) -> Path:
        """
        Save summary as Markdown
        
        Args:
            summary: Summary dict
            filename: Output filename
            
        Returns:
            Path to saved file
        """
        if self.output_dir:
            output_path = self.output_dir / filename
            self.output_dir.mkdir(parents=True, exist_ok=True)
        else:
            output_path = Path(filename)
        
        lines = []
        lines.append("# Truth-Filter Validation Summary")
        lines.append("")
        lines.append(f"**Total Runs**: {summary['total_runs']}")
        lines.append(f"**Passed**: {summary['passed_runs']} ({summary['pass_rate']:.1%})")
        lines.append(f"**Warnings**: {summary['warning_runs']}")
        lines.append(f"**Failed**: {summary['failed_runs']}")
        lines.append("")
        
        # Filter statistics
        lines.append("## Filter Statistics")
        lines.append("")
        for filter_name, stats in summary.get('filter_statistics', {}).items():
            lines.append(f"### {filter_name.replace('_', ' ').title()}")
            lines.append(f"- Pass: {stats['pass']}")
            lines.append(f"- Warning: {stats['warning']}")
            lines.append(f"- Fail: {stats['fail']}")
            lines.append("")
        
        # Molecule statistics
        mol_stats = summary.get('molecule_statistics', {})
        if mol_stats:
            lines.append("## Molecule Statistics")
            lines.append("")
            lines.append(f"- Total Original: {mol_stats['total_original']}")
            lines.append(f"- Total Filtered: {mol_stats['total_filtered']}")
            lines.append(f"- Retention Rate: {mol_stats['retention_rate']:.1%}")
            lines.append("")
        
        # Write file
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        
        logger.info(f"Saved summary Markdown to {output_path}")
        return output_path

