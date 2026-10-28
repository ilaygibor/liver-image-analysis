from pathlib import Path

import numpy as np

from src.collagen_prof import HEIGHT, NODE_COUNT, WIDTH, write_legacy_profile


def test_grid_size() -> None:
    assert WIDTH * HEIGHT == NODE_COUNT == 488_500


def test_legacy_writer_layout(tmp_path: Path) -> None:
    xs = np.tile(np.arange(WIDTH - 1, -1, -1) * 1e-6, HEIGHT)
    permeability = np.full(NODE_COUNT, 6e13)
    output = tmp_path / "test.prof"
    write_legacy_profile(output, xs, permeability)
    lines = output.read_text().splitlines()
    assert lines[0] == "((surface_body nodes 488500)"
    assert lines[1] == "(x"
    assert lines[488_503] == "(y"
    assert lines[-1] == ")) "
    assert "(permeability" not in lines

