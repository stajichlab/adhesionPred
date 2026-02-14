#!/usr/bin/env python3
"""
Filter tab-delimited file based on percent identity threshold.
Only outputs specified columns (1st, 2nd, 3rd, and 11th - evalue).
"""

import argparse
import sys


def process_file(input_file, percent_id):
    try:
        with open(input_file) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                fields = line.split("\t")

                # Check if we have enough columns
                if len(fields) < 11:
                    print(
                        f"Warning: {input_file}:{line_num} has fewer than 11 columns, skipping",
                        file=sys.stderr,
                    )
                    continue

                try:
                    # 3rd column (index 2) is percent identity
                    percent_identity = float(fields[2])

                    # Filter based on percent identity threshold
                    if percent_identity > percent_id:
                        # Print columns 1, 2, 3, and 11 (indices 0, 1, 2, 10)
                        output_fields = [fields[0], fields[1], fields[2], fields[10]]
                        print("\t".join(output_fields))

                except ValueError:
                    print(
                        f"Warning: {input_file}:{line_num} has invalid percent identity value '{fields[2]}', skipping",
                        file=sys.stderr,
                    )
                    continue

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found", file=sys.stderr)
        return 1
    except OSError as e:
        print(f"Error reading file '{input_file}': {e}", file=sys.stderr)
        return 1

    return 0


def main():
    parser = argparse.ArgumentParser(description="Filter tab-delimited file by percent identity")
    parser.add_argument("input_files", nargs="+", help="Input tab-delimited file(s)")
    parser.add_argument(
        "--percent_id",
        type=float,
        default=90.0,
        help="Minimum percent identity threshold (default: 90.0)",
    )

    args = parser.parse_args()

    exit_code = 0
    for input_file in args.input_files:
        result = process_file(input_file, args.percent_id)
        if result != 0:
            exit_code = result

    if exit_code != 0:
        sys.exit(exit_code)


if __name__ == "__main__":
    main()
