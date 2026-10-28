# Existing Phase 4 mesh profile

[Phase4_F4_permeability.prof](Phase4_F4_permeability.prof) is an existing generated output from the project archive. It contains mesh X and Y coordinates and a mapped coefficient for each node. The source stain images and mesh input are not included in this extract.

## Observed file contents

| Check | Value |
| --- | --- |
| Declared nodes | 8,427 |
| Values in each of the X, Y and permeability blocks | 8,427 |
| X range | -0.0096526 to 0.013502 m |
| Y range | -0.0082432 to 0.0073281 m |
| Nodes with coefficient `6.0e13` | 7,331 |
| Nodes with coefficient `1.0e14` | 582 |
| Nodes with coefficient `2.5e14` | 514 |

SHA-256:

```text
b8d4b8f2886919182feb051dfe308ad0c63d6af49c717c732ab0c1e5cb887fae
```

The preparation check parsed every line, verified all three block lengths and coefficient counts, and checked this checksum. The file contains only numeric blocks and structural delimiters.

The archived workflow normalized the mesh bounds onto the 977 × 500 image field, reversed X orientation, and used nearest-neighbor KDTree assignment. The input data and execution environment are needed to regenerate the file and assess that mapping. The checks above establish the file's structure and identity, not solver acceptance or scientific accuracy.

The mesh profile's `(y ...)` block contains spatial Y coordinates. This differs from the separate legacy grid format, where the same block name contains permeability values. See [input and format notes](../docs/INPUTS_AND_FORMATS.md).

No ANSYS import or simulation run was performed while preparing this extract.
