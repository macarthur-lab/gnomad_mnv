# gnomAD v4 MNV Pipeline — Documentation

This document describes the **v4 MNV discovery and annotation pipeline** in `v4/`.
It operates on the gnomAD v4 exomes VDS (GRCh38) using a scan-based approach that
replaces v2's `hl.window_by_locus()` (removed in modern Hail).

For the original v2 (GRCh37) pipeline and paper reference, see the `v2/` directory
and the [paper](https://www.nature.com/articles/s41467-019-12438-5).

---

## 1. Repository Structure (v4)

```
v4/
├── __init__.py          # Package init
├── resources.py         # Data loading, GCS output paths, TableResource definitions
└── discover_mnv.py      # MNV discovery (scan-based) and frequency/VEP annotation
```

### Dependencies

- **gnomad_qc** (`gnomad_qc.v4.resources.basics`): `get_gnomad_v4_vds` — loads the
  gnomAD v4 exomes VDS with standard sample filtering (hard filters, UKB exclusions).
- **gnomad_methods** (`gnomad.resources`): `TableResource`, `public_release("exomes")` —
  frequency and filter annotations from the gnomAD v4.1 exomes release.
- **gnomad_qc** (`gnomad_qc.v4.resources.annotations`): `get_vep("exomes")` — VEP
  transcript consequences for v4 exomes.

---

## 2. Pipeline Overview

The pipeline has two steps, selected via CLI flags:

| Step | Flag | What it does |
|------|------|--------------|
| **Discovery** | `--discover` | Scan the unsplit VDS for MNV pairs, classify, aggregate |
| **Annotation** | `--annotate` | Join AC/AF/filters and VEP onto discovery output |

Both can be run together: `--discover --annotate`.

### Usage

```bash
# Submit to Dataproc (see CLAUDE.md for zip packaging details)
hailctl dataproc submit <CLUSTER> v4/discover_mnv.py \
  --pyfiles /tmp/pyfiles.zip \
  -- --discover --annotate --test --overwrite
```

CLI arguments:

| Argument | Description |
|----------|-------------|
| `--discover` | Run MNV discovery |
| `--annotate` | Run frequency + VEP annotation |
| `--test` | Subset to PCSK9 region (chr1:55039447-55064852) |
| `--overwrite` | Overwrite existing output files |
| `--max-distance N` | Max bp distance between SNVs (default: 2) |
| `--high-quality-only` | Filter to high-quality samples only |
| `--release-only` | Filter to release samples only |

---

## 3. Discovery Step (`--discover`)

### Algorithm: Single-Pass Scan

The discovery uses a single-pass architecture operating on the **unsplit** sparse VDS
variant_data. It works in local allele space (LGT/LPGT) throughout, avoiding the need
for `sparse_split_multi`.

**Pipeline flow:**

1. **Load**: `get_gnomad_v4_vds(split=False)` — unsplit VDS with LGT, LPGT, PID, LA.
2. **Filter rows**: Keep rows with at least one SNP alt allele.
3. **Select entries**: Keep LGT, LPGT, PID, LA (pre-split names).
4. **Unfilter**: `mt.unfilter_entries()` — fills ref-block-covered positions with
   LGT=0/0 so the scan sees every row for every sample.
5. **Localize**: `mt._localize_entries()` — convert to Table (required to avoid Hail's
   `KeyError: 'agg_capability'` IR bug with entry-level scans).
6. **Scan** (`_scan_for_candidates`): Per-sample `hl.scan.fold` tracks a sliding window
   of recent non-ref entries. Each entry stores `(locus, entry_struct)`. The window is
   pruned by distance and contig at each row.
7. **Classify** (`_classify_mnv_pairs`): For each entry with a defined prev window,
   classify all pairs and extract biallelic alleles via `LA` indexing.
8. **Aggregate** (`_aggregate_mnv_pairs`): `group_by(locus, alleles, prev_locus,
   prev_alleles)` across all samples.

### Why local allele space works

- Het-ref / hom-var status is the same in local and global representations.
- Phase orientation (LPGT 0|1 vs 1|0) is relative to local ref/alt.
- Biallelic alleles are extracted via `LA[max(LGT[0], LGT[1])]` to get the global
  alt allele index, then used as the aggregation key.

### Classification rules

| Category | Criteria |
|----------|----------|
| **het-het** | Both het-ref, both LPGT phased, same PID, same phase orientation (`LPGT[0] == prev.LPGT[0]`) |
| **hom-hom** | Both hom-var |
| **het-hom** | One het-ref and one hom-var (either direction) |

### Key functions

| Function | Location | Purpose |
|----------|----------|---------|
| `_is_nonref(e)` | Expression helper | Check `is_defined(LGT) & is_non_ref()` |
| `_get_biallelic_snv(entry)` | Expression helper | Extract `{alleles, is_snp}` via LA indexing |
| `_classify_genotype_pair(cur, prev)` | Expression helper | Return `{is_hethet, is_homhom, is_hethom}` |
| `_scan_for_candidates(ht, max_distance)` | Pipeline step | `hl.scan.fold` sliding window |
| `_classify_mnv_pairs(ht)` | Pipeline step | Classify + checkpoint + explode |
| `_aggregate_mnv_pairs(ht)` | Pipeline step | group_by + aggregate counts |
| `discover_mnv(vds, ...)` | Top-level | Orchestrates scan → classify → aggregate |

---

## 4. Annotation Step (`--annotate`)

Reads the discovery HT and joins:

1. **Frequency**: AC, AF, filters from `public_release("joint")` for both SNVs.
2. **VEP**: Transcript consequences from `get_vep(data_type="exomes")` for both SNVs.

### Function

`annotate_mnv(mnv_ht)` — joins frequency and VEP, returns annotated Table.

---

## 5. Output Schema

### Discovery output

| Field | Type | Description |
|-------|------|-------------|
| `locus` | `locus<GRCh38>` | Position of SNP2 (downstream) |
| `alleles` | `array<str>` | [ref, alt] of SNP2 |
| `prev_locus` | `locus<GRCh38>` | Position of SNP1 (upstream) |
| `prev_alleles` | `array<str>` | [ref, alt] of SNP1 |
| `dist` | `int32` | Distance in bp (1 or 2) |
| `n_hethet` | `int64` | Phased het-het count across samples |
| `n_homhom` | `int64` | Hom-hom count |
| `n_hethom` | `int64` | Het-hom count (either direction) |
| `n_total` | `int64` | Sum of above |

### Annotated output (additional fields)

| Field | Type | Description |
|-------|------|-------------|
| `AC`, `AF` | `int32`, `float64` | SNP2 allele count and frequency |
| `filters` | `set<str>` | SNP2 variant filters |
| `prev_AC`, `prev_AF` | `int32`, `float64` | SNP1 allele count and frequency |
| `prev_filters` | `set<str>` | SNP1 variant filters |
| `vep` | `struct` | SNP2 VEP consequences |
| `prev_vep` | `struct` | SNP1 VEP consequences |

---

## 6. Resources (`v4/resources.py`)

### Constants

- `CURRENT_VERSION = "4.1"`
- `MNV_ENTRIES_TO_KEEP = ["GT", "PGT", "PID"]` — global (post-split) names; the
  pipeline remaps to pre-split names (LGT, LPGT, PID) internally.

### Functions

| Function | Returns | Description |
|----------|---------|-------------|
| `mnv_discovery(test=False)` | `TableResource` | Discovery output path |
| `mnv_annotated(test=False)` | `TableResource` | Annotated output path |
| `get_gnomad_v4_vds(...)` | `VariantDataset` | Loads unsplit VDS via `gnomad_qc` |

### Output paths

| Resource | Production path | Test path |
|----------|-----------------|-----------|
| Discovery | `gs://gnomad/v4.1/mnv/exomes/gnomad.exomes.v4.1.mnv_discovery.ht` | `gs://gnomad-tmp/gnomad_v4.1_testing/mnv/exomes/...` |
| Annotated | `gs://gnomad/v4.1/mnv/exomes/gnomad.exomes.v4.1.mnv_annotated.ht` | `gs://gnomad-tmp/gnomad_v4.1_testing/mnv/exomes/...` |

---

## 7. Technical Notes

### Why `_localize_entries` + `hl.scan.array_agg`?

Using `hl.scan._prev_nonnull` or `hl.scan.fold` directly in MT entry context causes
`KeyError: 'agg_capability'` in Hail's `CSEAnalysisPass`. The workaround — used by
gnomAD's `compute_last_ref_block_end()` — is to localize entries to a Table first:

```python
ht = mt._localize_entries("__entries", "__cols")
scan_result = hl.scan.array_agg(lambda e: hl.scan.fold(...), ht.__entries)
```

### Why `hl.scan.fold` instead of `_prev_nonnull`?

`_prev_nonnull` tracks only the single most recent non-null entry. With max_distance=2,
this misses the (P, P+2) pair when a sample has non-ref variants at P, P+1, and P+2
(the scan at P+2 sees P+1, not P). The fold-based scan maintains a window of all recent
non-ref entries within range, catching all valid pairs.

### Why `unfilter_entries()`?

The sparse VDS only stores entries at variant sites (non-ref calls) and ref block
boundaries. Without unfiltering, the scan might carry stale state across gaps where a
sample has ref-block coverage but no stored entry. `unfilter_entries()` fills in
LGT=0/0 at ref-block-covered positions, ensuring the scan sees every row and can
correctly prune the window by distance.

### PGT is a Call after split (but we don't split)

After `sparse_split_multi`, PGT is a Hail CallExpression. In our unsplit pipeline, LPGT
is also a Call type. Use `.phased`, `[0]`, `[1]` — not string methods.

---

## 8. Test Results (PCSK9 region)

Region: `chr1:55039447-55064852` (PCSK9 gene, ~25 kb)

| Metric | Single-pass (fold) | Previous 4-pass | `_prev_nonnull` |
|--------|-------------------|-----------------|-----------------|
| MNV pairs | 523 | 523 | 514 |
| Runtime | ~25 min | ~80 min | ~65 min |
| Edge case (P, P+2) | Caught | Caught | Missed (9 pairs) |

The 9 additional pairs (vs `_prev_nonnull`) are (P, P+2) edge cases where an
intervening variant at P+1 caused `_prev_nonnull` to return P+1 instead of P.

---

## 9. v2 → v4 Key Differences

| Aspect | v2 | v4 |
|--------|----|----|
| Reference | GRCh37 | GRCh38 |
| Data format | Dense MatrixTable | Sparse VDS (VariantDataset) |
| Windowing | `hl.window_by_locus()` | `hl.scan.fold` sliding window |
| Splitting | `hl.split_multi_hts` | None (works in local allele space) |
| Phasing check | `GT.phased` + `GT == prev_GT` | `LPGT.phased` + PID match + phase orientation |
| Resources | `gnomad_hail` | `gnomad_methods` + `gnomad_qc` |
| Hail version | 0.2.11 | 0.2.134+ |
| Processing | Per-chromosome | Whole-genome (or per-interval for testing) |

---

## 10. Future Work

- **Genome support**: Add `data_type="genomes"` support (different VDS, potentially
  different entry fields).
- **Functional annotation**: Port v2's combined-codon VEP annotation
  (`v2/code/annotate_vep_mnv.py`, `v2/util/mnv_functions.py`) to classify MNV pairs
  by functional impact (e.g., "Gained PTV", "Rescued PTV").
- **Release table formatting**: Add fields and formatting required for public release
  (refer to v2 release scripts).
- **Full-genome run**: Validate on full chromosomes, then genome-wide.
- **Sex chromosome handling**: Special handling for male hemizygosity on X/Y.
