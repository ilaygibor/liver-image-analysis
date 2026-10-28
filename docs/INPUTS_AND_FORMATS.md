# Inputs and profile conventions

These notes describe the preserved implementation in [collagen_prof.py](../src/collagen_prof.py) and [mesh_profile.py](../src/mesh_profile.py). They are an input contract and code interpretation, not a claim of biological validation.

## Stain images

For each phase `N` from 0 to 4, supply these files in one directory:

| Stain | Required filename |
| --- | --- |
| COL-I | `phaseN_COL-IB.png` |
| COL-III | `phaseN_COL-IIIB.png` |
| COL-IV | `phaseN_COL-IVB.png` |
| COL-VI | `phaseN_COL-VIB.png` |

The images are loaded through OpenCV and converted from BGR to RGB. A 5 × 5 Gaussian blur is applied to an HSV representation for the tissue mask. The historical saturation `> 3`, value `< 300`, and dark-pixel exclusion rules are retained. On the 8-bit HSV arrays produced here, the value `< 300` condition does not exclude pixels.

Collagen masks use strict RGB bounds, then remove four fixed RGB colors. CPA is calculated separately for each original image as `100 × selected tissue pixels / tissue pixels`, before the masks are cropped to their common upper-left extent. A zero-area tissue mask produces CPA 0.

Images must already correspond spatially. The code performs no image registration. It crops differing dimensions to their common minimum width and height; that operation alone does not align tissue sections.

## Overlay and coefficients

The active thresholds are defined in `RGB_THRESHOLDS`, with fixed Phase 4 replacements in `PHASE_RGB_THRESHOLD_OVERRIDES`. `ACCEPTED_CPA_VALUES` stores historical comparison targets; printing a target does not mean the script is optimizing toward it.

Overlapping masks give COL-I precedence, followed by COL-III. COL-IV, COL-VI, and unselected pixels are black in the overlay. The final field is resized with `cv2.INTER_AREA` to width 977 and height 500. Only exact red and green matches after this resizing receive collagen-specific values; intermediate colors receive the background value.

| Exact overlay color | Stored coefficient |
| --- | ---: |
| Red, COL-I | `2.5e14` |
| Green, COL-III | `1.0e14` |
| All other colors | `6.0e13` |

The source calls these values permeability. The original solver configuration is needed to determine their intended physical meaning and units.

## ANSYS mesh input

The mesh reader requires a UTF-8, tab-separated file. The first six header fields, after trimming, must be:

```text
Entity    X(m)    Y(m)    Z(m)    Node ID    Body
```

The actual file must use tab separators, not runs of spaces. The second row is always consumed as a summary row, so include it even when no node count is declared. A summary beginning with `N Nodes` is checked against the parsed number of rows.

Each remaining nonempty row supplies the six fields above. X, Y, and Z must parse as floats; Node ID must parse as an integer. The parser rejects an empty mesh, duplicate node IDs, duplicate XY pairs, and count mismatches. Mapping requires nonzero X and Y spans.

Z is parsed but is not used in this two-dimensional mapping. X and Y are scaled independently to the image coordinate bounds, which changes aspect ratio when necessary. X is reversed. Each mapped point receives the nearest source pixel's value. The output preserves the original mesh row order and original mesh X/Y coordinates.

## Two output layouts

| Property | Legacy grid profile | Mesh profile |
| --- | --- | --- |
| Generator | `collagen_prof.py` | `mesh_profile.py` |
| Node count | 488,500 | Number of parsed mesh rows |
| `(x ...)` | Repeated descending grid X | Original mesh X |
| `(y ...)` | Stored permeability values | Original mesh Y |
| `(permeability ...)` | Absent | Mapped coefficients |
| Inspector | `inspect_prof.py` | `validate_mesh_profile` during generation |

The legacy grid X runs from `976e-6` to `0` in each row. Its spatial Y array, written to the separate tabular output, runs from `250e-6` to `-249e-6`. The unusual legacy block naming and closing delimiters are retained for compatibility with the archived workflow.

Do not use `inspect_prof.py` to assess the mesh profile: it interprets its `y` block as profile values. The mesh generator validates its own node counts and coordinate round trip using the supplied mesh.

## Verification scope

The selected code is preserved from a reviewed project snapshot. The public extract includes only the existing self-contained profile-format tests. The original dataset-dependent tests cannot be reproduced from this extract alone.

Preparation checks compiled all included Python source and tests without execution and inspected every line of the example numeric profile. OpenCV and pytest were unavailable in the preparation environment, so no image-processing run, pytest execution, original-data regression, or ANSYS import was performed.
