#!/usr/bin/env python3
"""
Convert markdown diagram files to HTML
Preserves ASCII art diagrams and formatting
"""

import os
import re
from pathlib import Path

# HTML template with styling
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
            line-height: 1.6;
        }}
        .container {{
            background-color: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-bottom: 2px solid #ecf0f1;
            padding-bottom: 5px;
        }}
        h3 {{
            color: #7f8c8d;
            margin-top: 20px;
        }}
        pre {{
            background-color: #f8f9fa;
            border: 1px solid #e1e4e8;
            border-radius: 6px;
            padding: 15px;
            overflow-x: auto;
            font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
            font-size: 13px;
            line-height: 1.4;
            white-space: pre;
        }}
        code {{
            background-color: #f1f3f5;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
            font-size: 0.9em;
        }}
        pre code {{
            background-color: transparent;
            padding: 0;
        }}
        ul, ol {{
            margin-left: 20px;
        }}
        li {{
            margin-bottom: 8px;
        }}
        p {{
            margin-bottom: 15px;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ecf0f1;
            color: #7f8c8d;
            font-size: 0.9em;
            text-align: center;
        }}
    </style>
</head>
<body>
    <div class="container">
        {content}
        <div class="footer">
            <p>Generated from: {source_file}</p>
        </div>
    </div>
</body>
</html>
"""

def markdown_to_html(md_content: str, source_file: str) -> str:
    """Convert markdown to HTML, preserving ASCII art"""
    
    # Extract title from first h1
    title_match = re.match(r'^#\s+(.+)$', md_content, re.MULTILINE)
    title = title_match.group(1) if title_match else "Diagram"
    
    html_parts = []
    lines = md_content.split('\n')
    i = 0
    
    while i < len(lines):
        line = lines[i]
        
        # H1
        if line.startswith('# '):
            html_parts.append(f'<h1>{line[2:]}</h1>')
            i += 1
        
        # H2
        elif line.startswith('## '):
            html_parts.append(f'<h2>{line[3:]}</h2>')
            i += 1
        
        # H3
        elif line.startswith('### '):
            html_parts.append(f'<h3>{line[4:]}</h3>')
            i += 1
        
        # Code block (preserves ASCII art)
        elif line.startswith('```'):
            code_lines = []
            i += 1
            # Optional language identifier
            lang = line[3:].strip() if len(line) > 3 else ''
            
            while i < len(lines) and not lines[i].startswith('```'):
                code_lines.append(lines[i])
                i += 1
            i += 1  # Skip closing ```
            
            code_content = '\n'.join(code_lines)
            html_parts.append(f'<pre><code>{code_content}</code></pre>')
        
        # Horizontal rule
        elif line.strip() == '---':
            html_parts.append('<hr>')
            i += 1
        
        # List items
        elif line.strip().startswith('- ') or line.strip().startswith('* '):
            html_parts.append('<ul>')
            while i < len(lines) and (lines[i].strip().startswith('- ') or 
                                      lines[i].strip().startswith('* ') or
                                      lines[i].strip().startswith('  ')):
                content = lines[i].strip()
                if content.startswith('- ') or content.startswith('* '):
                    # Remove bullet and clean
                    item_text = content[2:].strip()
                    # Handle inline code
                    item_text = re.sub(r'`([^`]+)`', r'<code>\1</code>', item_text)
                    # Handle bold
                    item_text = re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', item_text)
                    html_parts.append(f'<li>{item_text}</li>')
                i += 1
            html_parts.append('</ul>')
        
        # Numbered list
        elif re.match(r'^\d+\.\s+', line):
            html_parts.append('<ol>')
            while i < len(lines) and re.match(r'^\d+\.\s+', lines[i]):
                content = lines[i]
                item_text = re.sub(r'^\d+\.\s+', '', content)
                # Handle inline code
                item_text = re.sub(r'`([^`]+)`', r'<code>\1</code>', item_text)
                # Handle bold
                item_text = re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', item_text)
                html_parts.append(f'<li>{item_text}</li>')
                i += 1
            html_parts.append('</ol>')
        
        # Regular paragraph
        elif line.strip():
            para_lines = [line]
            i += 1
            # Collect consecutive non-empty lines
            while i < len(lines) and lines[i].strip() and not lines[i].startswith('#'):
                if lines[i].startswith('```') or lines[i].startswith('- ') or lines[i].startswith('* '):
                    break
                para_lines.append(lines[i])
                i += 1
            
            para_text = ' '.join(para_lines).strip()
            if para_text:
                # Handle inline code
                para_text = re.sub(r'`([^`]+)`', r'<code>\1</code>', para_text)
                # Handle bold
                para_text = re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', para_text)
                # Handle italic
                para_text = re.sub(r'\*([^*]+)\*', r'<em>\1</em>', para_text)
                html_parts.append(f'<p>{para_text}</p>')
        else:
            i += 1
    
    content = '\n'.join(html_parts)
    return HTML_TEMPLATE.format(title=title, content=content, source_file=source_file)


def convert_all_files():
    """Convert all .md files in current directory to .html"""
    current_dir = Path(__file__).parent
    
    md_files = list(current_dir.glob('*.md'))
    
    if not md_files:
        print("No .md files found in current directory")
        return
    
    print(f"Found {len(md_files)} markdown files to convert")
    
    for md_file in md_files:
        print(f"Converting: {md_file.name}")
        
        # Read markdown
        with open(md_file, 'r', encoding='utf-8') as f:
            md_content = f.read()
        
        # Convert to HTML
        html_content = markdown_to_html(md_content, md_file.name)
        
        # Write HTML
        html_file = md_file.with_suffix('.html')
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"  -> Created: {html_file.name}")
    
    print(f"\nConversion complete! {len(md_files)} HTML files created.")


if __name__ == '__main__':
    convert_all_files()

