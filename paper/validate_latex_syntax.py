#!/usr/bin/env python3
"""
Basic LaTeX Syntax Validator
Checks for common LaTeX errors without requiring full compilation
"""

import re
from pathlib import Path

PAPER_DIR = Path(__file__).parent
MANUSCRIPT = PAPER_DIR / "manuscript_draft.tex"

def check_balanced_braces(content):
    """Check for balanced braces"""
    errors = []
    lines = content.split('\n')
    
    for i, line in enumerate(lines, 1):
        # Count braces (ignore escaped braces)
        line_clean = re.sub(r'\\{', '', line)
        line_clean = re.sub(r'\\}', '', line_clean)
        
        open_braces = line_clean.count('{')
        close_braces = line_clean.count('}')
        
        if open_braces != close_braces:
            errors.append(f"Line {i}: Unbalanced braces ({{: {open_braces}, }}: {close_braces})")
    
    return errors

def check_balanced_commands(content):
    """Check for balanced LaTeX commands"""
    errors = []
    
    # Check \begin/\end pairs
    begins = re.findall(r'\\begin\{([^}]+)\}', content)
    ends = re.findall(r'\\end\{([^}]+)\}', content)
    
    if len(begins) != len(ends):
        errors.append(f"Unbalanced environments: {len(begins)} \\begin, {len(ends)} \\end")
    
    # Check for common unclosed commands
    unclosed = []
    for match in re.finditer(r'\\(section|subsection|subsubsection)\{', content):
        # Check if there's a closing brace
        pos = match.end()
        remaining = content[pos:]
        if not remaining.startswith('}'):
            # Find next brace
            next_brace = remaining.find('}')
            if next_brace == -1:
                unclosed.append(f"Unclosed {match.group(1)} at position {match.start()}")
    
    errors.extend(unclosed)
    return errors

def check_figure_paths(content):
    """Check if figure paths are valid"""
    errors = []
    figures_dir = PAPER_DIR / "figures"
    
    # Find all includegraphics
    figure_paths = re.findall(r'\\includegraphics.*?\{([^}]+)\}', content)
    
    for fig_path in figure_paths:
        # Remove path prefix if present
        fig_name = Path(fig_path).name
        fig_file = figures_dir / fig_name
        
        if not fig_file.exists():
            errors.append(f"Missing figure: {fig_path}")
    
    return errors

def check_table_paths(content):
    """Check if table paths are valid"""
    errors = []
    tables_dir = PAPER_DIR / "tables"
    
    # Find all \input{table/...}
    table_paths = re.findall(r'\\input\{tables/([^}]+)\}', content)
    
    for table_name in table_paths:
        table_file = tables_dir / table_name
        if not table_file.exists():
            errors.append(f"Missing table: {table_name}")
    
    return errors

def check_bibliography():
    """Check if bibliography file exists"""
    bib_file = PAPER_DIR / "references.bib"
    if not bib_file.exists():
        return ["Missing references.bib file"]
    return []

def main():
    print("=" * 60)
    print("LaTeX Syntax Validation")
    print("=" * 60)
    
    if not MANUSCRIPT.exists():
        print(f"❌ ERROR: Manuscript not found: {MANUSCRIPT}")
        return 1
    
    content = MANUSCRIPT.read_text(encoding='utf-8')
    
    all_errors = []
    
    # Run checks
    print("\n=== Checking Balanced Braces ===")
    brace_errors = check_balanced_braces(content)
    if brace_errors:
        all_errors.extend(brace_errors)
        for err in brace_errors:
            print(f"❌ {err}")
    else:
        print("✅ All braces balanced")
    
    print("\n=== Checking LaTeX Commands ===")
    cmd_errors = check_balanced_commands(content)
    if cmd_errors:
        all_errors.extend(cmd_errors)
        for err in cmd_errors:
            print(f"❌ {err}")
    else:
        print("✅ All commands balanced")
    
    print("\n=== Checking Figure Paths ===")
    fig_errors = check_figure_paths(content)
    if fig_errors:
        all_errors.extend(fig_errors)
        for err in fig_errors:
            print(f"❌ {err}")
    else:
        print("✅ All figure files exist")
    
    print("\n=== Checking Table Paths ===")
    table_errors = check_table_paths(content)
    if table_errors:
        all_errors.extend(table_errors)
        for err in table_errors:
            print(f"❌ {err}")
    else:
        print("✅ All table files exist")
    
    print("\n=== Checking Bibliography ===")
    bib_errors = check_bibliography()
    if bib_errors:
        all_errors.extend(bib_errors)
        for err in bib_errors:
            print(f"❌ {err}")
    else:
        print("✅ Bibliography file exists")
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    if all_errors:
        print(f"❌ Found {len(all_errors)} potential issues")
        print("\nPlease fix these before compilation.")
        return 1
    else:
        print("✅ No syntax errors detected")
        print("\nManuscript appears ready for LaTeX compilation.")
        print("Note: This is a basic check. Full compilation may reveal additional issues.")
        return 0

if __name__ == "__main__":
    exit(main())

