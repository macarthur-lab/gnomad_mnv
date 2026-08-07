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
  frequency and filter annotations from the gnomAD v4.1 exomes release. Read once in
  `main` and passed to both steps, so within one invocation the hotfix AF and the
  reported AC/AF/filters come from the same release. Running `--discover` and
  `--annotate` as separate invocations reads it once each.
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
# Small test run over one gene (--intervals implies --test).
# See CLAUDE.md for zip packaging details.
hailctl dataproc submit <CLUSTER> v4/discover_mnv.py \
  --pyfiles /tmp/pyfiles.zip \
  -- --discover --annotate --intervals PCNT=chr21:46324141-46445769 --overwrite

# Full-cohort, genome-wide run (no interval flags). Note that --test alone does
# NOT subset intervals; it only redirects output to the test bucket.
hailctl dataproc submit <CLUSTER> v4/discover_mnv.py \
  --pyfiles /tmp/pyfiles.zip \
  -- --discover --annotate --overwrite
```

CLI arguments:

| Argument | Description |
|----------|-------------|
| `--discover` | Run MNV discovery |
| `--annotate` | Run frequency + VEP annotation |
| `--test` | Write output to the test bucket instead of production. **Does not subset intervals** — on its own this is a full-cohort, genome-wide run |
| `--intervals GENE=INTERVAL ...` | Subset discovery to one or more gene intervals. Implies `--test`. With no values, defaults to `PCNT=chr21:46324141-46445769`. Mutually exclusive with `--chr` |
| `--chr N` | Run on a single whole chromosome (e.g. `--chr 21`, `--chr X`); writes to the production bucket. Mutually exclusive with `--intervals` |
| `--ukb-only` | Use the UKB-only VDS subset; adds `.ukb_only` to the output path |
| `--output-suffix S` | Free-form suffix appended last to the output path, so runs don't overwrite each other |
| `--classify-n-partitions N` | Repartition before classification. Only applied with `--intervals` |
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

1. **Load**: `get_gnomad_v4_vds(split=False)` — unsplit VDS with LGT, LPGT, PID, LA,
   GQ, DP, LAD.
2. **Drop junk sites**: drop unsplit sites where every alt is `AS_lowqual` (unsplit info
   HT `get_info(split=False)`, `AS_lowqual` is `array<bool>` → all-lowqual = `hl.all`).
   Row-only filter; info HT re-read with the variant_data's partition intervals for a
   no-shuffle join.
3. **Filter rows**: Keep rows with at least one SNP alt allele.
4. **Select entries**: Keep LGT, LPGT, PID, GQ, DP, LAD, LA (pre-split names).
5. **Adjust genotypes** (on the MT, before localizing), in gnomAD QC's canonical order:
   (a) flag `_het_non_ref = LGT.is_het_non_ref()`; (b) hom-alt depletion hotfix —
   reclassify a high-AB (>0.9) het-ref call at a common variant
   (`public_release("exomes").freq[0].AF` > 0.01, joined per-alt) to hom-alt, skipping
   het_non_ref and `fixed_homalt_model` samples; (c) `adjusted_sex_ploidy_expr` (no-op on
   autosomes); (d) `adj = get_adj_expr(LGT, GQ, DP, LAD)`. Sex karyotype and the hom-alt
   model flag come from `meta()`.
6. **Localize + build `_alts`**: `mt._localize_entries()` — convert to Table (required to
   avoid Hail's `KeyError: 'agg_capability'` IR bug with entry-level scans). No
   `unfilter_entries`/densify (see below). Each non-ref entry gets `_alts`: one record per
   distinct carried **SNV** alt (`locus`, `alleles`, `is_hom`, `hap`), ploidy-safe.
   Indel alts are dropped here, so they never enter the scan window. Keep only
   LGT, PID, adj, `_alts` per entry.
7. **Scan** (`_scan_for_candidates`): Per-sample `hl.scan.fold` tracks a sliding window
   of recent non-ref entries. Each entry stores `(locus, entry_struct)`. The window is
   pruned by distance and contig at each row.
8. **Classify** (`_classify_mnv_pairs`): Filter to candidate carriers (entries with a
   `prev` window), `explode` to one row per carrier, then classify every
   `cur._alts × prev._alts` alt pair; a pair's `adj` = both entries' adj. Exploding
   before building pair arrays bounds per-task memory (no classify repartition needed).
9. **Aggregate** (`_aggregate_mnv_pairs`): `group_by(locus, alleles, prev_locus,
   prev_alleles)` across all samples, emitting raw and `_adj` counts.

### Why local allele space works

- Het / hom status is the same in local and global representations.
- The haplotype an alt sits on (`hap`) is derived from LPGT and is relative to local
  ref/alt.
- Biallelic alleles are extracted per carried alt via `LA[li]` (global alt index), then
  used as the aggregation key. het_non_ref (`1/2`) yields two records, one per alt.

### Classification rules

Applied per `cur._alts × prev._alts` alt pair (records from `_get_carried_alts`):

| Category | Criteria |
|----------|----------|
| **het-het** | Both het (not hom), both `hap` defined (phased), same PID, same `hap` (cis) |
| **hom-hom** | Both hom (incl. hotfix-converted high-AB het-refs, and hemizygous non-PAR chrX/Y calls) |
| **het-hom** | Exactly one hom (a hom occupies both haplotypes → always cis; see KNOWN_ISSUES.md #2 for the hemizygous case, where that does not hold) |

**Sex chromosomes.** `adjusted_sex_ploidy_expr` sets non-PAR XY het calls to missing
(only non-het calls become haploid), and `_is_nonref` then drops them, so on non-PAR
chrX/Y an XY sample contributes only hom-hom pairs. XX calls on chrY are set to missing
outright, so XX contributes nothing there at all. `--chr X`/`Y` output is therefore
incomplete for XY samples, not merely unvalidated — and on chrY it is barely usable:

| | chrX non-PAR (G6PD) | chrY non-PAR (chrY:2781480-14000000) |
|---|---|---|
| non-ref XY calls, pre-ploidy | 1,230,965 | 222,672,048 |
| dropped as hets | 17.3% | 94.7% |
| dropped, counting only adj-quality calls | 5.5% | 94.5% |
| surviving haploid calls that are hom-var | 1,017,790 / 1,017,790 | 11,719,615 / 11,719,615 |
| end-to-end XY output | 405 / 405 hom-hom | 708,336 / 708,336 hom-hom (7,235 pairs) |
| end-to-end XX output, same region | 3,072 het-het | 0 pairs |

The chrY window includes ampliconic/palindromic MSY with known mapping problems, which
inflates the het rate; treat 94.7% as specific to that window, not a clean chrY estimate.
Either way, chrY discards nearly all quality calls as artifact hets and should not be
released without a decision on hemizygous modelling.

Validation scope: ploidy handling and bucket assignment are measured on both chromosomes
as above. No MNV-level accuracy check (are the emitted pairs real, is cis correct) has
been run off the autosomes — the source paper used autosomes.

### Key functions

| Function | Location | Purpose |
|----------|----------|---------|
| `_is_nonref(e)` | Expression helper | Check `is_defined(LGT) & is_non_ref()` |
| `_get_carried_alts(entry)` | Expression helper | Array of `{locus, alleles, is_hom, hap}`, one per carried **SNV** alt (ploidy-safe) |
| `_classify_alt_pair(ca, pa, same_phase_set)` | Expression helper | Return `{is_hethet, is_homhom, is_hethom}` |
| `_scan_for_candidates(ht, max_distance)` | Pipeline step | `hl.scan.fold` sliding window |
| `_classify_mnv_pairs(ht)` | Pipeline step | Classify + checkpoint + explode |
| `_aggregate_mnv_pairs(ht)` | Pipeline step | group_by + aggregate counts |
| `discover_mnv(vds, ...)` | Top-level | Orchestrates scan → classify → aggregate |

---

## 4. Annotation Step (`--annotate`)

Reads the discovery HT and joins:

1. **Frequency**: AC, AF, filters from `public_release("exomes")` for both SNVs.
2. **VEP**: Transcript consequences from `get_vep(data_type="exomes")` for both SNVs.

### Function

`annotate_mnv(mnv_ht, release_ht)` — joins frequency and VEP using the given release
table, returns annotated Table.

---

## 5. Output Schema

### Discovery output

| Field | Type | Description |
|-------|------|-------------|
| `locus` | `locus<GRCh38>` | Position of SNP2 (downstream) |
| `alleles` | `array<str>` | [ref, alt] of SNP2 |
| `prev_locus` | `locus<GRCh38>` | Position of SNP1 (upstream) |
| `prev_alleles` | `array<str>` | [ref, alt] of SNP1 |
| `dist` | `int32` | Distance in bp (1 or 2 in the common case; see KNOWN_ISSUES.md #1 for rare min-rep shift exceptions, which `main()` writes to a sibling `*.dist_outliers.ht`) |
| `n_hethet` | `int64` | Phased het-het count across samples (raw) |
| `n_homhom` | `int64` | Hom-hom count (raw); includes hemizygous pairs |
| `n_hethom` | `int64` | Het-hom count, either direction (raw) |
| `n_total` | `int64` | Sum of above (raw) |
| `n_hethet_adj` | `int64` | het-het count where both genotypes pass adj |
| `n_homhom_adj` | `int64` | hom-hom count where both genotypes pass adj |
| `n_hethom_adj` | `int64` | het-hom count where both genotypes pass adj |
| `n_total_adj` | `int64` | Total count where both genotypes pass adj |

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
- `MNV_ENTRIES_TO_KEEP = ["GT", "PGT", "PID", "GQ", "DP", "AD"]` — global (post-split)
  names; the pipeline remaps to pre-split names (LGT, LPGT, PID, GQ, DP, LAD) internally.
  GQ/DP/AD feed the hom-alt hotfix (allele balance) and adj.

### Functions

| Function | Returns | Description |
|----------|---------|-------------|
| `mnv_discovery(test=False)` | `TableResource` | Discovery output path |
| `mnv_annotated(test=False)` | `TableResource` | Annotated output path |
| `get_gnomad_v4_vds(...)` | `VariantDataset` | Loads unsplit VDS via `gnomad_qc` |
| `get_gnomad_v4_ukb_vds(...)` | `VariantDataset` | Loads the UKB-only VDS subset (`--ukb-only`) |
| `UKB_VDS` | `VariantDatasetResource` | Path to the pre-built UKB-only VDS subset |

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

### Why no `unfilter_entries()` / densify?

The sparse VDS only stores entries at variant sites (non-ref calls). Densifying hom-ref
positions to `LGT=0/0` is **not needed**: the scan classifies a *missing* entry exactly
like `0/0` (both fail `_is_nonref`), and the window is pruned by comparing each stored
entry's locus to the **current row's** locus — not by intervening rows — so a stale
entry is dropped at the first row beyond `max_distance` regardless of what lies between.
Verified on synthetic data that dense (`0/0`), sparse+unfilter, and sparse-without-unfilter
produce byte-identical MNV output. Skipping the densify avoids materializing a `0/0` for
every sample at every site.

### PGT is a Call after split (but we don't split)

After `sparse_split_multi`, PGT is a Hail CallExpression. In our unsplit pipeline, LPGT
is also a Call type. Use `.phased`, `[0]`, `[1]` — not string methods.

---

## 8. Test Results (PCSK9 region)

> **Historical — pre-dates the current algorithm.** These counts were measured before
> the hom-alt depletion hotfix, `adj` annotation, sex-ploidy adjustment, min-rep-keyed
> aggregation, and het_non_ref alt-splitting were added. They compare scan
> implementations against each other, not the pipeline's present output, and have not
> been re-run. Retained for the scan-architecture comparison only.

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
- **Sex chromosome handling**: Decide whether to model hemizygosity, rather than
  inheriting the current behavior described in §3.
