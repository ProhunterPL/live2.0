#!/usr/bin/env python3
"""
Manuscript Validation Script
Checks LaTeX manuscript for common issues:
- Cross-references (refs vs labels)
- Figure/table citations
- File existence
- Placeholders
- Basic consistency checks
"""

import re
import os
from pathlib import Path
from typing import List, Dict, Set, Tuple

# Paths
PAPER_DIR = Path(__file__).parent
MANUSCRIPT = PAPER_DIR / "manuscript_draft.tex"
FIGURES_DIR = PAPER_DIR / "figures"
TABLES_DIR = PAPER_DIR / "tables"

class ManuscriptChecker:
    def __init__(self, manuscript_path: Path):
        self.manuscript_path = manuscript_path
        self.content = manuscript_path.read_text(encoding='utf-8')
        self.errors = []
        self.warnings = []
        self.info = []
    
    def check_cross_references(self):
        """Check that all \ref{} have corresponding \label{}"""
        print("\n=== Checking Cross-References ===")
        
        # Find all references
        refs = set(re.findall(r'\\ref\{([^}]+)\}', self.content))
        
        # Find all labels in main file
        labels = set(re.findall(r'\\label\{([^}]+)\}', self.content))
        
        # Also check labels in imported table files
        table_includes = re.findall(r'\\input\{tables/([^}]+)\}', self.content)
        for table_name in table_includes:
            table_file = TABLES_DIR / table_name
            if table_file.exists():
                table_content = table_file.read_text(encoding='utf-8')
                table_labels = set(re.findall(r'\\label\{([^}]+)\}', table_content))
                labels.update(table_labels)
        
        # Check for missing labels
        missing_labels = refs - labels
        if missing_labels:
            self.errors.append(f"Missing labels for references: {missing_labels}")
            print(f"❌ ERROR: Missing labels: {missing_labels}")
        else:
            print(f"✅ All {len(refs)} references have corresponding labels")
        
        # Check for unused labels
        unused_labels = labels - refs
        if unused_labels:
            self.warnings.append(f"Unused labels: {unused_labels}")
            print(f"⚠️  WARNING: Unused labels: {unused_labels}")
        
        return len(missing_labels) == 0
    
    def check_figures(self):
        """Check that all figure files exist and are referenced"""
        print("\n=== Checking Figures ===")
        
        # Find all figure includes
        figure_includes = re.findall(r'\\includegraphics.*?\{([^}]+)\}', self.content)
        figure_refs = set(re.findall(r'\\ref\{fig:([^}]+)\}', self.content))
        figure_labels = set(re.findall(r'\\label\{fig:([^}]+)\}', self.content))
        
        # Check file existence
        missing_files = []
        for fig_path in figure_includes:
            # Remove path prefix if present
            fig_name = Path(fig_path).name
            fig_file = FIGURES_DIR / fig_name
            if not fig_file.exists():
                missing_files.append(fig_path)
                self.errors.append(f"Missing figure file: {fig_path}")
                print(f"❌ ERROR: Missing figure: {fig_path}")
        
        if not missing_files:
            print(f"✅ All {len(figure_includes)} figure files exist")
        
        # Check that all figures are referenced
        missing_refs = figure_labels - figure_refs
        if missing_refs:
            self.warnings.append(f"Figures not referenced in text: {missing_refs}")
            print(f"⚠️  WARNING: Unreferenced figures: {missing_refs}")
        
        return len(missing_files) == 0
    
    def check_tables(self):
        """Check that all table files exist and are referenced"""
        print("\n=== Checking Tables ===")
        
        # Find all table includes
        table_includes = re.findall(r'\\input\{tables/([^}]+)\}', self.content)
        table_refs = set(re.findall(r'\\ref\{tab:([^}]+)\}', self.content))
        table_labels = set(re.findall(r'\\label\{tab:([^}]+)\}', self.content))
        
        # Check file existence
        missing_files = []
        for table_name in table_includes:
            table_file = TABLES_DIR / table_name
            if not table_file.exists():
                missing_files.append(table_name)
                self.errors.append(f"Missing table file: {table_name}")
                print(f"❌ ERROR: Missing table: {table_name}")
        
        if not missing_files:
            print(f"✅ All {len(table_includes)} table files exist")
        
        # Check that all tables are referenced
        missing_refs = table_labels - table_refs
        if missing_refs:
            self.warnings.append(f"Tables not referenced in text: {missing_refs}")
            print(f"⚠️  WARNING: Unreferenced tables: {missing_refs}")
        
        return len(missing_files) == 0
    
    def check_placeholders(self):
        """Check for remaining placeholders"""
        print("\n=== Checking Placeholders ===")
        
        # Common placeholder patterns
        placeholder_patterns = [
            r'\[Author Name\]',
            r'\[Institution\]',
            r'\[email\]',
            r'\[collaborators\]',
            r'\[funding sources\]',
            r'\[AWS/local cluster\]',
            r'\[GitHub repository\]',
            r'\[username\]',
            r'\[Zenodo DOI.*?\]',
            r'\[data DOI.*?\]',
            r'\[XX\]',
            r'\[YY\]',
            r'\[XXX,XXX\]',
        ]
        
        found_placeholders = []
        for pattern in placeholder_patterns:
            matches = re.findall(pattern, self.content, re.IGNORECASE)
            if matches:
                found_placeholders.extend(matches)
        
        if found_placeholders:
            unique_placeholders = set(found_placeholders)
            self.warnings.append(f"Remaining placeholders: {unique_placeholders}")
            print(f"⚠️  WARNING: Found {len(unique_placeholders)} placeholder types:")
            for ph in sorted(unique_placeholders):
                count = found_placeholders.count(ph)
                print(f"   - {ph} ({count} occurrences)")
        else:
            print("✅ No placeholders found")
        
        return len(found_placeholders) == 0
    
    def check_citations(self):
        """Check for citation issues"""
        print("\n=== Checking Citations ===")
        
        # Find all citations
        citations = set(re.findall(r'\\citep\{([^}]+)\}', self.content))
        citations.update(re.findall(r'\\citet\{([^}]+)\}', self.content))
        
        # Check if references.bib exists
        bib_file = PAPER_DIR / "references.bib"
        if not bib_file.exists():
            self.errors.append("Missing references.bib file")
            print("❌ ERROR: references.bib not found")
            return False
        
        # Read bibliography
        bib_content = bib_file.read_text(encoding='utf-8')
        
        # Check for missing references
        missing_refs = []
        for cite in citations:
            # Check if citation exists in .bib (simplified check)
            cite_key = cite.split(',')[0].strip()
            if cite_key not in bib_content:
                missing_refs.append(cite_key)
        
        if missing_refs:
            self.warnings.append(f"Possible missing references: {missing_refs}")
            print(f"⚠️  WARNING: {len(missing_refs)} citations not found in references.bib")
        else:
            print(f"✅ All {len(citations)} citations found in references.bib")
        
        return len(missing_refs) == 0
    
    def check_sections(self):
        """Check section structure"""
        print("\n=== Checking Section Structure ===")
        
        # Expected sections
        expected_sections = [
            r'\\section\{Introduction\}',
            r'\\section\{Methods\}',
            r'\\section\{Results\}',
            r'\\section\{Discussion\}',
            r'\\section\{Conclusions\}',
        ]
        
        found_sections = []
        for pattern in expected_sections:
            if re.search(pattern, self.content):
                found_sections.append(pattern)
        
        if len(found_sections) == len(expected_sections):
            print(f"✅ All {len(expected_sections)} main sections present")
        else:
            missing = set(expected_sections) - set(found_sections)
            self.warnings.append(f"Missing sections: {missing}")
            print(f"⚠️  WARNING: Missing sections: {missing}")
        
        return len(found_sections) == len(expected_sections)
    
    def check_word_count(self):
        """Estimate word count (rough)"""
        print("\n=== Checking Word Count ===")
        
        # Remove LaTeX commands and comments
        text = self.content
        text = re.sub(r'%.*', '', text)  # Remove comments
        text = re.sub(r'\\[a-zA-Z]+\{?[^}]*\}?', '', text)  # Remove LaTeX commands
        text = re.sub(r'\{[^}]*\}', '', text)  # Remove braces
        words = len(text.split())
        
        print(f"📊 Estimated word count: ~{words:,} words")
        print(f"   (Target: ~6,000 words)")
        
        if words < 5000:
            self.warnings.append(f"Word count ({words}) below target (6000)")
        elif words > 7000:
            self.warnings.append(f"Word count ({words}) above target (6000)")
        else:
            print("✅ Word count within acceptable range")
        
        return True
    
    def run_all_checks(self):
        """Run all validation checks"""
        print("=" * 60)
        print("Manuscript Validation Report")
        print("=" * 60)
        
        checks = [
            ("Cross-References", self.check_cross_references),
            ("Figures", self.check_figures),
            ("Tables", self.check_tables),
            ("Placeholders", self.check_placeholders),
            ("Citations", self.check_citations),
            ("Sections", self.check_sections),
            ("Word Count", self.check_word_count),
        ]
        
        results = {}
        for name, check_func in checks:
            try:
                results[name] = check_func()
            except Exception as e:
                self.errors.append(f"Error in {name} check: {e}")
                print(f"❌ ERROR in {name}: {e}")
                results[name] = False
        
        # Summary
        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        print(f"Errors: {len(self.errors)}")
        print(f"Warnings: {len(self.warnings)}")
        print(f"Info: {len(self.info)}")
        
        if self.errors:
            print("\n❌ ERRORS:")
            for i, error in enumerate(self.errors, 1):
                print(f"  {i}. {error}")
        
        if self.warnings:
            print("\n⚠️  WARNINGS:")
            for i, warning in enumerate(self.warnings, 1):
                print(f"  {i}. {warning}")
        
        # Return success status
        return len(self.errors) == 0

def main():
    checker = ManuscriptChecker(MANUSCRIPT)
    success = checker.run_all_checks()
    
    if success:
        print("\n✅ Manuscript validation completed with no errors!")
        return 0
    else:
        print("\n❌ Manuscript validation found errors. Please fix before submission.")
        return 1

if __name__ == "__main__":
    exit(main())

