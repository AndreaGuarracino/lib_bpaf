#!/usr/bin/env python3
"""
Normalize PAF files by truncating float values to 3 decimal places.
Equivalent to normalize_paf.pl but with explicit UTF-8 handling.
"""
import sys
import re

def normalize_line(line):
    """
    Normalize float fields in PAF format:
    1. Normal decimal: :f:0.993724 → :f:0.993
    2. Leading dot: :f:.0549 → :f:0.054
    3. Integer: :f:0 → :f:0.000
    """
    # Handle normal decimals: :f:0.993724 → :f:0.993
    def truncate_decimal(match):
        prefix = match.group(1)
        integer_part = match.group(2)
        decimal_part = match.group(3)
        # Pad with zeros if needed, then take first 3 digits
        padded = (decimal_part + '000')[:3]
        return f"{prefix}{integer_part}.{padded}"

    line = re.sub(r'(:f:)(\d+)\.(\d+)', truncate_decimal, line)

    # Handle leading dots: :f:.0549 → :f:0.054
    def truncate_leading_dot(match):
        prefix = match.group(1)
        decimal_part = match.group(2)
        padded = (decimal_part + '000')[:3]
        return f"{prefix}0.{padded}"

    line = re.sub(r'(:f:)\.(\d+)', truncate_leading_dot, line)

    # Handle integers: :f:0 → :f:0.000
    def add_decimal(match):
        prefix = match.group(1)
        integer_part = match.group(2)
        suffix = match.group(3)
        return f"{prefix}{integer_part}.000{suffix}"

    line = re.sub(r'(:f:)(\d+)(\s|\t|$)', add_decimal, line)

    return line

def main():
    """Process stdin line by line, normalizing float values."""
    # Use binary mode with explicit UTF-8 to avoid locale issues
    stdin = sys.stdin.buffer
    stdout = sys.stdout.buffer

    for line in stdin:
        # Decode from UTF-8, normalize, encode back
        try:
            text = line.decode('utf-8')
            normalized = normalize_line(text)
            stdout.write(normalized.encode('utf-8'))
        except UnicodeDecodeError:
            # If line isn't valid UTF-8, pass through unchanged
            stdout.write(line)

    stdout.flush()

if __name__ == '__main__':
    main()
