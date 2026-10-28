"""Map an image-derived permeability field onto ANSYS mesh coordinates.

This reconstructs the KDTree workflow recovered from the historical
``Lobule Collagen Extractor`` notebook. The historical hard-coded node count
and coordinate-file parser are replaced with values parsed from the mesh file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

try:
    from src.collagen_prof import build_overlay, overlay_to_arrays, phase_image_paths
except ModuleNotFoundError:  # Support direct execution as ``python src/mesh_profile.py``.
    from collagen_prof import build_overlay, overlay_to_arrays, phase_image_paths


@dataclass(frozen=True)
class MeshNode:
    entity: str
    x: float
    y: float
    z: float
    node_id: int
    body: str


def read_ansys_nodes(path: Path) -> list[MeshNode]:
    nodes: list[MeshNode] = []
    declared_count: int | None = None
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.reader(handle, delimiter="\t")
        header = next(rows, None)
        if not header or [value.strip() for value in header[:6]] != [
            "Entity", "X(m)", "Y(m)", "Z(m)", "Node ID", "Body"
        ]:
            raise ValueError("Unexpected ANSYS node-table header")
        summary = next(rows, None)
        if summary:
            match = re.match(r"(\d+)\s+Nodes", summary[0].strip())
            if match:
                declared_count = int(match.group(1))
        for line_number, row in enumerate(rows, 3):
            if not row or not row[0].strip():
                continue
            try:
                nodes.append(
                    MeshNode(
                        entity=row[0].strip(),
                        x=float(row[1]),
                        y=float(row[2]),
                        z=float(row[3]),
                        node_id=int(row[4]),
                        body=row[5].strip(),
                    )
                )
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Invalid mesh row {line_number}: {row}") from exc

    if not nodes:
        raise ValueError("Mesh table contains no nodes")
    if declared_count is not None and declared_count != len(nodes):
        raise ValueError(f"Mesh declares {declared_count} nodes but contains {len(nodes)}")
    ids = [node.node_id for node in nodes]
    if len(set(ids)) != len(ids):
        raise ValueError("Mesh contains duplicate Node IDs")
    xy = [(node.x, node.y) for node in nodes]
    if len(set(xy)) != len(xy):
        raise ValueError("Mesh contains duplicate XY coordinates")
    return nodes


def normalize_mesh_to_image_coordinates(
    nodes: list[MeshNode], source_x: np.ndarray, source_y: np.ndarray
) -> np.ndarray:
    """Map mesh bounds onto the historical image-coordinate bounds.

    Image x=0 is the right-most column, so mesh X is reversed. Image Y is
    positive at the top, matching increasing mesh Y. Independent axis scaling
    is required because the historical pipeline forced every overlay to
    977x500 regardless of its original aspect ratio.
    """
    target_x = np.array([node.x for node in nodes], dtype=np.float64)
    target_y = np.array([node.y for node in nodes], dtype=np.float64)
    tx_min, tx_max = float(target_x.min()), float(target_x.max())
    ty_min, ty_max = float(target_y.min()), float(target_y.max())
    if tx_min == tx_max or ty_min == ty_max:
        raise ValueError("Mesh must span non-zero X and Y ranges")

    sx_min, sx_max = float(source_x.min()), float(source_x.max())
    sy_min, sy_max = float(source_y.min()), float(source_y.max())
    mapped_x = sx_max - (target_x - tx_min) / (tx_max - tx_min) * (sx_max - sx_min)
    mapped_y = sy_min + (target_y - ty_min) / (ty_max - ty_min) * (sy_max - sy_min)
    return np.column_stack((mapped_x, mapped_y))


def map_permeability(
    nodes: list[MeshNode],
    source_x: np.ndarray,
    source_y: np.ndarray,
    source_permeability: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not (len(source_x) == len(source_y) == len(source_permeability)):
        raise ValueError("Source coordinate/value lengths do not match")
    source_coordinates = np.column_stack((source_x, source_y))
    mapped_targets = normalize_mesh_to_image_coordinates(nodes, source_x, source_y)
    distances, indices = cKDTree(source_coordinates).query(mapped_targets)
    return source_permeability[indices], distances, indices


def write_mesh_profile(path: Path, nodes: list[MeshNode], permeability: np.ndarray) -> None:
    if len(nodes) != len(permeability):
        raise ValueError("Node and permeability counts do not match")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"((surface_body nodes {len(nodes)})\n")
        handle.write("(x\n")
        handle.writelines(f"{node.x:.7e}\n" for node in nodes)
        handle.write(")\n")
        handle.write("(y\n")
        handle.writelines(f"{node.y:.7e}\n" for node in nodes)
        handle.write(")\n")
        handle.write("(permeability\n")
        handle.writelines(f"{value:.7e}\n" for value in permeability)
        handle.write(")\n")
        handle.write("))\n")


def validate_mesh_profile(path: Path, expected_nodes: list[MeshNode]) -> dict[str, object]:
    lines = path.read_text(encoding="utf-8").splitlines()
    match = re.fullmatch(r"\(\(surface_body nodes (\d+)\)", lines[0])
    if not match:
        raise ValueError("Invalid mesh-profile header")
    declared = int(match.group(1))

    def block(marker: str) -> list[float]:
        start = lines.index(marker) + 1
        end = lines.index(")", start)
        return [float(value) for value in lines[start:end]]

    xs, ys, values = block("(x"), block("(y"), block("(permeability")
    if not (declared == len(expected_nodes) == len(xs) == len(ys) == len(values)):
        raise ValueError("Profile block lengths do not match the mesh node count")
    expected_x = np.array([node.x for node in expected_nodes])
    expected_y = np.array([node.y for node in expected_nodes])
    if not np.allclose(xs, expected_x, rtol=0, atol=5e-11):
        raise ValueError("Profile X coordinates do not match the mesh")
    if not np.allclose(ys, expected_y, rtol=0, atol=5e-11):
        raise ValueError("Profile Y coordinates do not match the mesh")
    if lines[-1] != "))":
        raise ValueError("Profile has an invalid closing delimiter")
    return {
        "nodes": declared,
        "x_range": (min(xs), max(xs)),
        "y_range": (min(ys), max(ys)),
        "permeability_counts": dict(sorted(Counter(values).items())),
    }


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate_mesh_profile(
    phase: int, stain_dir: Path, mesh_path: Path, output_path: Path
) -> Path:
    nodes = read_ansys_nodes(mesh_path)
    overlay, _ = build_overlay(phase, phase_image_paths(stain_dir, phase))
    source_x, source_y, source_permeability = overlay_to_arrays(overlay)
    mapped, distances, indices = map_permeability(
        nodes, source_x, source_y, source_permeability
    )
    write_mesh_profile(output_path, nodes, mapped)
    summary = validate_mesh_profile(output_path, nodes)
    print(f"Wrote {output_path}")
    print(f"SHA-256: {sha256(output_path)}")
    print(f"Validation: {summary}")
    print(f"Unique source pixels selected: {len(set(indices.tolist()))} / {len(nodes)}")
    print(
        "Nearest-neighbor distance (image-coordinate meters): "
        f"max={distances.max():.3e}, mean={distances.mean():.3e}"
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", type=int, choices=range(5))
    parser.add_argument("mesh", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--stain-dir", type=Path, default=Path("inputs/stains"))
    args = parser.parse_args()
    generate_mesh_profile(args.phase, args.stain_dir, args.mesh, args.output)


if __name__ == "__main__":
    main()
