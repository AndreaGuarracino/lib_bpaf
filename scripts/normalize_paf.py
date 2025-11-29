#!/usr/bin/env python3
"""
Normalize PAF files by rounding float values to 3 decimal places.
Adds a small epsilon before rounding to handle f32 precision noise where
values like 0.939999997616 and 0.94000000 should both round to 0.940.
"""
import sys
import re

def normalize_line(line):
    """
    Normalize float fields in PAF format by rounding to 3 decimal places.
    Adds epsilon (1e-7) to handle f32 precision noise at boundaries.
    """
    # Handle all float formats: :f:0.993724, :f:.0549, :f:.1500, :f:0
    def round_float(match):
        prefix = match.group(1)
        value_str = match.group(2)
        try:
            value = float(value_str)
            # Add small epsilon to push f32 boundary cases up before rounding
            # f32 has ~7 significant digits, so 1e-7 is appropriate
            adjusted = value + 1e-7
            # Round to 3 decimal places
            rounded = round(adjusted, 3)
            return f"{prefix}{rounded:.3f}"
        except Exception:
            return match.group(0)

    # Match floats with optional leading zero: 0.123, .123, 123, 123.456
    line = re.sub(r'(:f:)(\.[0-9]+|[0-9]+\.?[0-9]*)', round_float, line)

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
