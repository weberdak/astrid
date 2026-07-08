#!/usr/bin/env python3
"""
Generate a detailed per-module coverage markdown report from Cobertura XML.
"""

import xml.etree.ElementTree as ET
import sys
from pathlib import Path
from typing import List, Dict, Tuple


def parse_cobertura(xml_path: str) -> Dict:
    """Parse Cobertura XML and extract coverage data."""
    tree = ET.parse(xml_path)
    root = tree.getroot()
    
    # Get overall stats
    coverage = root.attrib
    overall_line_rate = float(coverage.get('line-rate', 0)) * 100
    
    # Extract per-package data
    packages = {}
    for package in root.findall('.//package'):
        pkg_name = package.attrib.get('name', 'unknown')
        pkg_line_rate = float(package.attrib.get('line-rate', 0)) * 100
        
        classes_data = []
        for cls in package.findall('.//class'):
            cls_name = cls.attrib.get('name', 'unknown')
            cls_filename = cls.attrib.get('filename', '')
            cls_line_rate = float(cls.attrib.get('line-rate', 0)) * 100
            
            # Extract line and method counts
            lines_covered = 0
            lines_total = 0
            for line in cls.findall('.//line'):
                lines_total += 1
                if line.attrib.get('hits', '0') != '0':
                    lines_covered += 1
            
            classes_data.append({
                'name': cls_name,
                'filename': cls_filename,
                'line_rate': cls_line_rate,
                'lines_covered': lines_covered,
                'lines_total': lines_total,
            })
        
        packages[pkg_name] = {
            'line_rate': pkg_line_rate,
            'classes': classes_data,
        }
    
    return {
        'overall_line_rate': overall_line_rate,
        'packages': packages,
    }


def extract_module_name(filename: str) -> str:
    """Extract module name from filename."""
    if not filename:
        return 'unknown'
    # Convert path like 'src/matrix_a.rs' to 'matrix_a'
    path = Path(filename)
    return path.stem


def generate_markdown(coverage_data: Dict) -> str:
    """Generate markdown report from coverage data."""
    lines = []
    
    lines.append("# Rust Code Coverage Report\n")
    
    # Overall summary
    lines.append("## Overall Coverage\n")
    lines.append(f"- **Line Coverage**: {coverage_data['overall_line_rate']:.2f}%\n")
    
    # Per-module breakdown
    lines.append("## Per-Module Coverage\n")
    lines.append("| Module | Line Coverage | Lines |")
    lines.append("|--------|---|---|")
    
    # Sort modules by line coverage (descending)
    modules_data: List[Tuple[str, Dict]] = []
    
    for pkg_name, pkg_data in coverage_data['packages'].items():
        for cls in pkg_data['classes']:
            module_name = extract_module_name(cls['filename'])
            # Avoid duplicates by using module name as key
            found = False
            for i, (name, data) in enumerate(modules_data):
                if name == module_name:
                    # Merge data
                    modules_data[i] = (
                        name,
                        {
                            'line_rate': max(data['line_rate'], cls['line_rate']),
                            'lines_covered': data['lines_covered'] + cls['lines_covered'],
                            'lines_total': data['lines_total'] + cls['lines_total'],
                        }
                    )
                    found = True
                    break
            
            if not found:
                modules_data.append((module_name, {
                    'line_rate': cls['line_rate'],
                    'lines_covered': cls['lines_covered'],
                    'lines_total': cls['lines_total'],
                }))
    
    # Sort by line coverage
    modules_data.sort(key=lambda x: x[1]['line_rate'], reverse=True)
    
    for module_name, data in modules_data:
        lines.append(
            f"| `{module_name}` | "
            f"{data['line_rate']:.1f}% | "
            f"{data['lines_covered']}/{data['lines_total']} |"
        )
    
    lines.append("")
    
    # Coverage legend
    lines.append("## Coverage Legend\n")
    lines.append("- 🟢 **90-100%**: Excellent coverage")
    lines.append("- 🟡 **70-89%**: Good coverage")
    lines.append("- 🔴 **Below 70%**: Needs improvement\n")
    
    return "\n".join(lines)


def main():
    if len(sys.argv) < 2:
        print("Usage: python coverage_report.py <cobertura.xml> [output.md]")
        sys.exit(1)
    
    xml_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else "coverage-report.md"
    
    try:
        coverage_data = parse_cobertura(xml_path)
        markdown = generate_markdown(coverage_data)
        
        with open(output_path, 'w') as f:
            f.write(markdown)
        
        print(f"✓ Coverage report generated: {output_path}")
        print(markdown)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
