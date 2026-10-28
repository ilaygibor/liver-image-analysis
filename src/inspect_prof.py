"""Validate and summarize the archived two-block permeability profile format."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path


def read_block(lines: list[str], marker: str) -> list[float]:
    start = lines.index(marker) + 1
    end = lines.index(")", start)
    return [float(value) for value in lines[start:end]]


def inspect(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    expected_header = "((surface_body nodes "
    if not lines or not lines[0].startswith(expected_header):
        raise ValueError("Not a legacy surface_body nodes profile")
    declared = int(lines[0][len(expected_header) : -1])
    xs = read_block(lines, "(x")
    values = read_block(lines, "(y")
    if len(xs) != declared or len(values) != declared:
        raise ValueError(
            f"Declared {declared} nodes, found {len(xs)} x values and {len(values)} profile values"
        )
    print(f"file: {path}")
    print(f"nodes: {declared}")
    print(f"x range: {min(xs):.6e} .. {max(xs):.6e}")
    print(f"profile values: {dict(sorted(Counter(values).items()))}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", type=Path)
    args = parser.parse_args()
    inspect(args.profile)


if __name__ == "__main__":
    main()

