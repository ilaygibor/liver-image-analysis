"""Generate legacy-compatible permeability profiles from collagen stain images.

The default behavior intentionally reproduces the archived July 2025 pipeline.
Known quirks are preserved until profiles from additional phases are available
for comparison.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from pathlib import Path

import cv2
import numpy as np


STAINS = ("COL-I", "COL-III", "COL-IV", "COL-VI")
FILE_SUFFIXES = {
    "COL-I": "COL-IB.png",
    "COL-III": "COL-IIIB.png",
    "COL-IV": "COL-IVB.png",
    "COL-VI": "COL-VIB.png",
}

ACCEPTED_CPA_VALUES = {
    "phase0": {"COL-I": 0.64, "COL-III": 0.64, "COL-IV": 17.01, "COL-VI": 14.76},
    "phase1": {"COL-I": 0.26, "COL-III": 1.05, "COL-IV": 15.50, "COL-VI": 13.06},
    "phase2": {"COL-I": 3.08, "COL-III": 3.80, "COL-IV": 18.75, "COL-VI": 18.44},
    "phase3": {"COL-I": 9.72, "COL-III": 6.33, "COL-IV": 24.21, "COL-VI": 25.90},
    "phase4": {"COL-I": 10.17, "COL-III": 11.88, "COL-IV": 27.00, "COL-VI": 31.41},
}

RGB_THRESHOLDS = {
    "COL-I": ((0, 255), (0, 200), (0, 150)),
    "COL-III": ((0, 255), (0, 255), (0, 145)),
    "COL-IV": ((100, 255), (0, 161), (0, 255)),
    "COL-VI": ((100, 255), (0, 176), (0, 255)),
}

# Phase 4 thresholds calibrated against the expert CPA targets stored in the
# historical notebooks. Only the active upper bound for each stain is adjusted.
PHASE_RGB_THRESHOLD_OVERRIDES = {
    "phase4": {
        "COL-I": ((0, 255), (0, 200), (0, 131)),
        "COL-III": ((0, 255), (0, 255), (0, 152)),
        "COL-IV": ((100, 255), (0, 146), (0, 255)),
        "COL-VI": ((100, 255), (0, 154), (0, 255)),
    }
}

EXCLUDE_RGBS = ((30, 9, 8), (83, 64, 50), (157, 132, 118), (28, 21, 20))

OVERLAY_COLORS = {
    "COL-I": np.array((255, 0, 0), dtype=np.uint8),
    "COL-III": np.array((0, 255, 0), dtype=np.uint8),
    "COL-IV": np.array((0, 0, 0), dtype=np.uint8),
    "COL-VI": np.array((0, 0, 0), dtype=np.uint8),
}

PERMEABILITY = {
    "COL-I": 2.5e14,
    "COL-III": 1.0e14,
    "healthy": 6.0e13,
    "background": 6.0e13,
}

WIDTH = 977
HEIGHT = 500
NODE_COUNT = WIDTH * HEIGHT


def phase_image_paths(input_dir: Path, phase: int) -> dict[str, Path]:
    paths = {
        stain: input_dir / f"phase{phase}_{suffix}"
        for stain, suffix in FILE_SUFFIXES.items()
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing stain image(s):\n" + "\n".join(missing))
    return paths


def load_image_and_tissue_mask(path: Path) -> tuple[np.ndarray, np.ndarray]:
    image_bgr = cv2.imread(str(path))
    if image_bgr is None:
        raise ValueError(f"Failed to load image: {path}")
    image_rgb = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2RGB)
    hsv = cv2.cvtColor(image_rgb, cv2.COLOR_RGB2HSV)
    blurred = cv2.GaussianBlur(hsv, (5, 5), 0)

    # The value < 300 condition is intentionally preserved from the archive.
    tissue = ((blurred[:, :, 1] > 3) & (blurred[:, :, 2] < 300)).astype(np.uint8) * 255
    tissue[np.all(image_rgb < 30, axis=2)] = 0
    return image_rgb, tissue


def collagen_mask(
    image_rgb: np.ndarray,
    stain: str,
    thresholds: dict[str, tuple[tuple[int, int], tuple[int, int], tuple[int, int]]] = RGB_THRESHOLDS,
) -> np.ndarray:
    r_range, g_range, b_range = thresholds[stain]
    selected = (
        (image_rgb[:, :, 0] > r_range[0])
        & (image_rgb[:, :, 0] < r_range[1])
        & (image_rgb[:, :, 1] > g_range[0])
        & (image_rgb[:, :, 1] < g_range[1])
        & (image_rgb[:, :, 2] > b_range[0])
        & (image_rgb[:, :, 2] < b_range[1])
    )
    for color in EXCLUDE_RGBS:
        selected &= ~np.all(image_rgb == color, axis=2)
    return selected.astype(np.uint8) * 255


def build_overlay(
    phase: int, paths: dict[str, Path]
) -> tuple[np.ndarray, dict[str, float]]:
    phase_key = f"phase{phase}"
    thresholds = PHASE_RGB_THRESHOLD_OVERRIDES.get(phase_key, RGB_THRESHOLDS)
    masks: list[np.ndarray] = []
    images: list[np.ndarray] = []
    cpas: dict[str, float] = {}

    for stain in STAINS:
        image, tissue = load_image_and_tissue_mask(paths[stain])
        mask = cv2.bitwise_and(collagen_mask(image, stain, thresholds), tissue)
        collagen_area = int(np.count_nonzero(mask))
        tissue_area = int(np.count_nonzero(tissue))
        cpas[stain] = 100 * collagen_area / tissue_area if tissue_area else 0.0
        masks.append(mask)
        images.append(image)

    height = min(mask.shape[0] for mask in masks)
    width = min(mask.shape[1] for mask in masks)
    masks = [mask[:height, :width] for mask in masks]
    images = [image[:height, :width] for image in images]

    stacked = np.stack(masks, axis=2)
    overlay = np.zeros((height, width, 3), dtype=np.uint8)
    active_y, active_x = np.where(np.any(stacked > 0, axis=2))

    # Keep the archived conflict-resolution order and arithmetic exactly.
    for y, x in zip(active_y, active_x):
        active = np.where(stacked[y, x] > 0)[0]
        if len(active) == 1:
            overlay[y, x] = OVERLAY_COLORS[STAINS[active[0]]]
            continue
        active_stains = [STAINS[index] for index in active]
        if "COL-I" in active_stains:
            overlay[y, x] = OVERLAY_COLORS["COL-I"]
        elif "COL-III" in active_stains:
            overlay[y, x] = OVERLAY_COLORS["COL-III"]
        else:
            red, green, blue = images[0][y, x].astype(np.int64)
            best_index = min(
                active,
                key=lambda index: _threshold_distance(
                    red, green, blue, thresholds[STAINS[index]]
                ),
            )
            overlay[y, x] = OVERLAY_COLORS[STAINS[best_index]]

    targets = ACCEPTED_CPA_VALUES.get(phase_key, {})
    for stain, cpa in cpas.items():
        target = targets.get(stain)
        detail = f" target={target:.2f}%" if target is not None else ""
        print(f"{stain}: CPA={cpa:.2f}%{detail}")
    return overlay, cpas


def _threshold_distance(
    red: np.int64,
    green: np.int64,
    blue: np.int64,
    thresholds: tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
) -> np.int64:
    (r_min, r_max), (g_min, g_max), (b_min, b_max) = thresholds
    return (
        abs(red - np.int64(r_min))
        + abs(red - np.int64(r_max))
        + abs(green - np.int64(g_min))
        + abs(green - np.int64(g_max))
        + abs(blue - np.int64(b_min))
        + abs(blue - np.int64(b_max))
    )


def overlay_to_arrays(overlay: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if overlay.shape[:2] != (HEIGHT, WIDTH):
        overlay = cv2.resize(overlay, (WIDTH, HEIGHT), interpolation=cv2.INTER_AREA)

    xs = np.tile(np.arange(WIDTH - 1, -1, -1, dtype=np.float64) * 1e-6, HEIGHT)
    ys = np.repeat((HEIGHT // 2 - np.arange(HEIGHT, dtype=np.float64)) * 1e-6, WIDTH)
    pixels = overlay.reshape(-1, 3)
    col_i = np.all(pixels == OVERLAY_COLORS["COL-I"], axis=1)
    col_iii = np.all(pixels == OVERLAY_COLORS["COL-III"], axis=1)
    permeability = np.full(NODE_COUNT, PERMEABILITY["background"], dtype=np.float64)
    permeability[col_iii] = PERMEABILITY["COL-III"]
    permeability[col_i] = PERMEABILITY["COL-I"]
    return xs, ys, permeability


def write_table(path: Path, xs: np.ndarray, ys: np.ndarray, permeability: np.ndarray) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("x\ty\tpermeability\n")
        for x, y, value in zip(xs, ys, permeability):
            handle.write(f"{x:.6e}\t{y:.6e}\t{value:.6e}\n")


def write_legacy_profile(path: Path, xs: np.ndarray, permeability: np.ndarray) -> None:
    """Write the exact two-block layout used by the successful Phase 2 file."""
    if len(xs) != NODE_COUNT or len(permeability) != NODE_COUNT:
        raise ValueError(f"Expected {NODE_COUNT} nodes")
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"((surface_body nodes {NODE_COUNT})\n")
        handle.write("(x\n")
        handle.writelines(f"{value:.6e}\n" for value in xs)
        handle.write(")\n")
        # Despite its name, this archived y block contains permeability values.
        handle.write("(y\n")
        handle.writelines(f"{value:.6e}\n" for value in permeability)
        handle.write(")\n")
        handle.write(")) \n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def generate_phase(phase: int, input_dir: Path, output_dir: Path) -> Path:
    paths = phase_image_paths(input_dir, phase)
    output_dir.mkdir(parents=True, exist_ok=True)
    overlay, _ = build_overlay(phase, paths)
    xs, ys, permeability = overlay_to_arrays(overlay)

    table_path = output_dir / f"overlay_table_phase{phase}.txt"
    profile_path = output_dir / f"my_permeability_profile_phase{phase}.prof"
    overlay_path = output_dir / f"overlay_phase{phase}.png"
    write_table(table_path, xs, ys, permeability)
    write_legacy_profile(profile_path, xs, permeability)
    cv2.imwrite(str(overlay_path), cv2.cvtColor(overlay, cv2.COLOR_RGB2BGR))

    counts = Counter(f"{value:.6e}" for value in permeability)
    print(f"Wrote {profile_path} ({NODE_COUNT} nodes)")
    print(f"SHA-256: {sha256(profile_path)}")
    print(f"Permeability counts: {dict(counts)}")
    return profile_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", type=int, choices=range(5), help="phase number, 0 through 4")
    parser.add_argument("--input-dir", type=Path, default=Path("inputs/stains"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    generate_phase(args.phase, args.input_dir, args.output_dir)


if __name__ == "__main__":
    main()
