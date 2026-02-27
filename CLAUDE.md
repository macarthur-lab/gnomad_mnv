# gnomad_mnv Project Reference

## Project Overview

MNV (Multi-Nucleotide Variant) discovery and annotation for gnomAD. The `v2/` directory contains the original GRCh37 pipeline; `v4/` is the active GRCh38 rewrite using a scan-based approach.

## Code Style

### Formatting

Code is formatted with **black** (preview mode, line length 88) and **isort** (profile
`"black"`). Config is in `pyproject.toml`. Pre-commit hooks are installed — run
`pre-commit run --all-files` to check, or let the git hook run on commit.

```bash
# Manual formatting
black v4/
isort --profile black -o gnomad v4/
```

### Docstrings

Use **Sphinx-style** (`:param:`, `:return:`) docstrings following the gnomad_methods
convention:

- **Module docstrings**: Describe the module's purpose, key concepts, and usage examples
  (with ``Usage::`` code blocks).
- **Function/method docstrings**: Start with a concise one-line summary, followed by
  extended description if needed. Use `:param name: Description.` and
  `:return: Description.` for all parameters and return values.
- **Inline code references**: Use double backticks (````name````).
- **Notes/warnings**: Use ``.. note::`` or ``.. warning::`` directives indented under
  the docstring.
- **Constants**: Document with a docstring on the line after the assignment.

Example:

```python
def my_function(
    ht: hl.Table,
    max_distance: int = 2,
) -> hl.Table:
    """Short summary of what the function does.

    Extended description with more detail about behavior, edge cases, or
    design decisions.

    .. note::

        Any important caveats go here.

    :param ht: Description of the table parameter.
    :param max_distance: Description with default info.
    :return: Description of what is returned.
    """
```

### Type Annotations

- **All functions** must have type annotations on parameters and return values.
- Use `typing.List`, `typing.Optional`, etc. for generic types.
- For Hail expression parameters, use `hl.expr.StructExpression`,
  `hl.expr.BooleanExpression`, etc.
- For Hail table/matrix types, use `hl.Table`, `hl.MatrixTable`,
  `hl.vds.VariantDataset`.

## Known Gotchas

- **`entries_to_keep` field naming**: Use post-split global names (GT, PGT, PID), not
  L-prefixed. The `_split_and_filter_variant_data_for_loading` remaps GT->LGT, AD->LAD,
  PL->LPL, PGT->LPGT before split, then `sparse_split_multi` renames back.
- **`filters` row field**: After `get_gnomad_v4_vds(split=True, entries_to_keep=...)`,
  the returned variant_data MT may not have a `filters` row field. Don't depend on it in
  the discovery step; add filters during annotation from the release HT instead.
- **PGT is a Call after split**: After `sparse_split_multi`, PGT is a Hail
  **CallExpression** (not a string). Use `pgt.phased`, `pgt[0]`, `pgt[1]` — not
  `.split()` or `.contains()`. Pre-split, LPGT is also a Call type.
- **MT entry scans cause IR bug**: Using `hl.scan._prev_nonnull` or `hl.scan.fold` in
  MT entry context + ANY downstream operation (write, checkpoint, `entries().group_by()`,
  `annotate_rows` with `hl.agg`) causes `KeyError: 'agg_capability'` in Hail's
  `CSEAnalysisPass`. **Working pattern**: `mt._localize_entries("__entries", "__cols")`
  → Table → `hl.scan.array_agg(lambda entry: ..., ht.__entries)` → then unlocalize if
  needed. This is the pattern used by gnomAD's `compute_last_ref_block_end()`.
- **`hl.scan.fold` initial_value**: `hl.scan.fold(initial_value, seqop, combop)` —
  `initial_value` must be a plain Expression (e.g. `hl.empty_array(type)`), NOT a lambda.
  `seqop` is `lambda acc: ...`, `combop` is `lambda a, b: ...`.
- **`public_release("joint")` lacks `filters`**: The joint release table only has
  `joint_freq`, `joint_faf`, `joint_fafmax` — no `filters`, `freq`, or `vep` fields.
  Use `public_release("exomes")` for per-variant filters and exome-specific freq.

## gnomad_qc API Reference

### Key Imports

| Resource | Import |
|----------|--------|
| VDS loader | `from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds` |
| Split helper | `gnomad_qc.v4.resources.basics._split_and_filter_variant_data_for_loading` (private) |
| VEP annotations | `from gnomad_qc.v4.resources.annotations import get_vep` |
| Sample QC | `from gnomad_qc.v4.resources.sample_qc import get_sample_qc` |
| Sample metadata | `from gnomad_qc.v4.resources.meta import meta` |

### get_gnomad_v4_vds Key Parameters

```python
get_gnomad_v4_vds(
    split=True,                    # Split multiallelics
    filter_intervals=None,         # List[str | hl.tinterval]
    high_quality_only=False,       # Filter to HQ samples
    release_only=False,            # Filter to release samples
    entries_to_keep=None,          # List[str] of global entry names (GT, PGT, PID, etc.)
    chrom=None,                    # Filter to chromosome(s)
    autosomes_only=False,
    filter_variant_ht=None,        # Keyed HT to semi-join on
)
```

### entries_to_keep Field Naming

Pre-split VDS uses L-prefixed (local) names. `entries_to_keep` uses **global** (post-split) names. The `_split_and_filter_variant_data_for_loading` function remaps:

| Global (user provides) | Local (pre-split) |
|------------------------|-------------------|
| GT | LGT |
| AD | LAD |
| PL | LPL |
| PGT | LPGT |

Other fields (PID, DP, GQ, RGQ, END) are the same pre- and post-split.

### get_vep

```python
from gnomad_qc.v4.resources.annotations import get_vep

vep_ht = get_vep(data_type="exomes").ht()  # VersionedTableResource
# Has: vep.transcript_consequences, vep.most_severe_consequence, etc.
```

## gnomad_methods (gnomad package)

### Key Imports

| Resource | Import |
|----------|--------|
| TableResource | `from gnomad.resources.resource_utils import TableResource` |
| VersionedTableResource | `from gnomad.resources.resource_utils import VersionedTableResource` |
| Public release | `from gnomad.resources.grch38.gnomad import public_release` |

### public_release

```python
from gnomad.resources.grch38.gnomad import public_release

# Available data_types: "exomes", "genomes", "joint"
exomes_ht = public_release("exomes").ht()  # v4.1 by default
# Has: freq (array<struct{AC, AF, AN, ...}>), filters, vep, region_flags, etc.

# CAUTION: "joint" only has joint_freq, joint_faf, joint_fafmax — no filters/vep.
joint_ht = public_release("joint").ht()
```

### TableResource Usage

```python
from gnomad.resources.resource_utils import TableResource

resource = TableResource(path="gs://gnomad/v4.1/mnv/exomes/output.ht")
resource.path       # The GCS path string
resource.ht()       # Read and return the Hail Table
```

## Hail VDS Field Reference

After `hl.experimental.sparse_split_multi()`:
- `LGT` -> `GT` (genotype call)
- `LPGT` -> `PGT` (physical genotype, Call type — use `.phased`, `[0]`, `[1]`)
- `LAD` -> `AD` (allele depths)
- `LPL` -> `PL` (phred-scaled likelihoods)
- `LA` is consumed/dropped during split
- `PID`, `DP`, `GQ`, `RGQ`, `END` are unchanged

## v4 MNV Pipeline

### Files

| File | Purpose |
|------|---------|
| `v4/__init__.py` | Package init |
| `v4/resources.py` | Resource paths, `get_gnomad_v4_vds()`, TableResource definitions |
| `v4/discover_mnv.py` | MNV discovery (scan-based) and annotation |

### Pipeline Steps

`--discover` uses a single-pass architecture with `hl.scan.fold` + `hl.case()` branching:

1. **Localize + scan**: `unfilter_entries()` → `_localize_entries` → annotate entries
   with row alleles → `hl.scan.array_agg` + `hl.scan.fold` with `hl.case()` branching
   to track a window of recent non-ref `(locus, entry)` tuples per sample. Window uses
   missing state (not empty array) to distinguish "no data" from "no nearby variants".
   Filter to rows where any sample has a candidate MNV pair.
2. **Classify**: For each candidate, use LGT/LPGT directly (valid in local allele space)
   for het-het / hom-hom / het-hom classification. Extract biallelic alleles via
   `entry.LA[max(LGT[0], LGT[1])]` indexing into the stored `_alleles`. Filter to pairs
   where both carried alleles are SNPs.
3. **Aggregate**: Checkpoint → explode → `group_by(locus, alleles, prev_locus,
   prev_alleles)` → aggregate n_hethet, n_homhom, n_hethom, n_total.

No `sparse_split_multi` or unlocalize step is needed — biallelic alleles are extracted
inline via LA indexing, and LGT/LPGT classification is valid in local allele space.

`--annotate`: Read discovery HT → join AC/AF/filters from `public_release("exomes")` +
VEP from `get_vep("exomes")` → write annotated HT.

The fold-based scan tracks an array (window) of recent non-ref values per sample, pruned
by distance and contig each row. This catches all pairs within `--max-distance` bp,
including the (P, P+2) edge case when P+1 is also non-ref.

### Output Paths

- Discovery: `gs://gnomad/v4.1/mnv/exomes/gnomad.exomes.v4.1.mnv_discovery.ht`
- Annotated: `gs://gnomad/v4.1/mnv/exomes/gnomad.exomes.v4.1.mnv_annotated.ht`
- Test: `gs://gnomad-tmp/gnomad_v4.1_testing/mnv/exomes/...`

### Dataproc Submission

**Important**: hailctl repackages `--pyfiles` into a temp zip using `os.walk` with
`os.path.relpath(file, os.path.join(entry, '..'))` as arcnames. This nests files under
the directory name (e.g., `gnomad_mnv/v4/`), breaking `from v4.resources import ...`.

**Workaround**: Build a single `.zip` file (hailctl passes single zips through as-is,
per line 29 of `hailtop/hailctl/dataproc/submit.py`). The zip must have `v4/` and
`gnomad_qc/` at the top level:

```bash
# Build a single zip with correct top-level package structure
cd <gnomad_mnv_root> && \
  rm -f /tmp/pyfiles.zip && \
  zip -r /tmp/pyfiles.zip v4/ -x '*.pyc' '*__pycache__*' && \
  cd <gnomad_qc_root> && \
  zip -r /tmp/pyfiles.zip gnomad_qc/ -x '*.pyc' '*__pycache__*' '*.DS_Store'

# Submit to cluster (single zip = used directly, not repackaged)
hailctl dataproc submit <CLUSTER> v4/discover_mnv.py \
  --pyfiles /tmp/pyfiles.zip \
  -- --discover --test --overwrite
```

**Why a single zip**: hailctl's submit.py line 29 checks
`pyfiles.endswith('.zip') and ',' not in pyfiles`. If true, it uses the zip directly.
Otherwise it repackages via `os.walk`, which nests packages under their directory names.
