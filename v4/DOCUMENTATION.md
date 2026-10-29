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
├── discover_mnv.py      # MNV discovery (scan-based) and frequency/VEP annotation
├── create_test_data.py  # Builds the small test dataset used by tests/ (see §8)
└── KNOWN_ISSUES.md      # Documented edge cases and shelved optimizations
tests/                   # pytest suite (see §8)
```

### Dependencies

- **gnomad_qc** (`gnomad_qc.v4.resources.basics`): `get_gnomad_v4_vds` — loads the
  gnomAD v4 exomes VDS with standard sample filtering (hard filters, UKB exclusions).
- **gnomad_methods** (`gnomad.resources`): `TableResource`, `public_release("exomes")` —
  frequency and filter annotations from the gnomAD v4.1 exomes release. Read once in
  `main` and passed to both steps so the hotfix AF and the reported AC/AF/filters are
  guaranteed to come from the same read. Running `--discover` and `--annotate` as
  separate invocations reads it independently each time, so consistency across the two
  runs isn't guaranteed if the underlying release changes between them.
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
  -- --discover --annotate --intervals SGCA=chr17:50164114-50176028 --overwrite

# Named test gene set: default, autosomes, sex_chr or all. default writes to the
# test bucket; the larger sets write to gs://gnomad/v4.1/mnv/exomes/test_genes/.
hailctl dataproc submit <CLUSTER> v4/discover_mnv.py \
  --pyfiles /tmp/pyfiles.zip \
  -- --discover --annotate --test-genes all --overwrite

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
| `--intervals GENE=INTERVAL ...` | Subset discovery to one or more gene intervals. Implies `--test`. With no values, uses the `default` test gene set: `SGCA=chr17:50164114-50176028` (padded by 100 bp). Mutually exclusive with `--chr` and `--test-genes` |
| `--test-genes {default,autosomes,sex_chr,all}` | Subset discovery to a named test gene set, each gene padded by 100 bp: `default` (SGCA, the smallest autosomal test gene), `autosomes` (AHNAK2, ANO5, CAPN3, DYSF, SGCA, the variant co-occurrence test genes), `sex_chr` (SHOX in PAR1, XG across the PAR1 boundary, EMD on non-PAR X, DDX3Y on Y), or `all`. `default` writes to the test bucket like `--intervals`; the larger sets write to a `test_genes/` folder under the production path, unless `--test` is given. Mutually exclusive with `--chr` and `--intervals` |
| `--chr N` | Run on a single whole chromosome (e.g. `--chr 21`, `--chr X`); writes to the production bucket. Mutually exclusive with `--intervals` and `--test-genes` |
| `--ukb-only` | Use the UKB-only VDS subset; adds `.ukb_only` to the output path |
| `--output-suffix S` | Free-form suffix appended last to the output path, so runs don't overwrite each other |
| `--classify-n-partitions N` | Repartition before classification. Requires `--intervals` or `--test-genes` (rejected otherwise) |
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
   model flag come from the sample metadata Table (`meta()`), which `main` reads and
   passes in, as it does the info Table. Steps 2–5 are `prepare_genotypes`, so tests can
   check the adjusted genotypes on their own.
6. **Localize + build `_alts`**: `mt._localize_entries()` — convert to Table (required to
   avoid Hail's `KeyError: 'agg_capability'` IR bug with entry-level scans). No
   `unfilter_entries`/densify (see below). Each non-ref entry gets `_alts`: one record per
   distinct carried **SNV** alt (`locus`, `alleles`, `is_hom`, `hap`), ploidy-safe.
   Indel alts are dropped here, so they never enter the scan window. Keep only
   PID, adj, `_alts` per entry (LGT is consumed building `_alts` and dropped).
7. **Scan** (`_scan_for_candidates`): Per-sample `hl.scan.fold` tracks a sliding window
   of recent SNV-carrying entries (gated by `_carries_snv`). Each entry stores
   `(locus, entry_struct)`; the window is pruned by distance and contig at each row.
   Rows are then narrowed to their candidate carriers (`_cands`, entries with a `prev`
   window) and checkpointed — the widened all-sample entries array is never
   materialized.
8. **Classify** (`_classify_mnv_pairs`): `explode` to one row per candidate carrier,
   then classify every `cur._alts × prev._alts` alt pair; a pair's `adj` = both
   entries' adj. Exploding before building pair arrays bounds per-task memory (no
   classify repartition needed).
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
outright. `--chr X`/`Y` output is therefore incomplete for XY samples, not merely
unvalidated. Measured: chrX non-PAR (G6PD) drops 17.3% of non-ref XY calls as hets
(5.5% of adj-quality calls); chrY non-PAR drops ~95% — chrY discards nearly all quality
calls as artifact hets and should not be released without a decision on hemizygous
modelling. Every surviving haploid call classified hom-hom on both chromosomes, and XX
on chrY produced zero pairs, as expected.

No MNV-level accuracy check (are the emitted pairs real, is cis correct) has been run
off the autosomes — the source paper used autosomes.

### Key functions

| Function | Location | Purpose |
|----------|----------|---------|
| `_is_nonref(e)` | Expression helper | Check `is_defined(LGT) & is_non_ref()` |
| `_carries_snv(e)` | Expression helper | Gate into the scan window: entry has ≥1 SNV alt |
| `_get_carried_alts(entry)` | Expression helper | Array of `{locus, alleles, is_hom, hap}`, one per carried **SNV** alt (ploidy-safe) |
| `_classify_alt_pair(ca, pa, same_phase_set)` | Expression helper | Return `{is_hethet, is_homhom, is_hethom}` |
| `_scan_for_candidates(ht, max_distance)` | Pipeline step | `hl.scan.fold` sliding window; narrows to candidate carriers and checkpoints |
| `_classify_mnv_pairs(ht)` | Pipeline step | Explode carriers + classify alt pairs |
| `_aggregate_mnv_pairs(ht)` | Pipeline step | group_by + aggregate counts |
| `prepare_genotypes(vds, release_ht, info_ht, meta_ht)` | Top-level | Site filters + genotype adjustments (steps 2–5) |
| `discover_mnv(vds, release_ht, info_ht, meta_ht, ...)` | Top-level | Orchestrates `prepare_genotypes` → scan → classify → aggregate |

---

## 4. Annotation Step (`--annotate`)

Reads the discovery HT and joins:

1. **Frequency**: AC, AF, filters from `public_release("exomes")` for both SNVs.
2. **VEP**: Transcript consequences from `get_vep(data_type="exomes")` for both SNVs.

### Function

`annotate_mnv(mnv_ht, release_ht, vep_ht)` — joins frequency and VEP from the given
tables, returns annotated Table. Both are read in `main` and passed in. Fields are
missing for SNVs absent from the release (missing `filters` ≠ PASS; PASS is an empty
set).

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
- `TEST_GENE_SETS` — named test gene sets for `--test-genes`: `default` (SGCA),
  `autosomes` (AHNAK2, ANO5, CAPN3, DYSF, SGCA), `sex_chr` (SHOX, XG, EMD, DDX3Y) and
  `all`. Built from `AUTOSOME_TEST_GENES` and `SEX_CHR_TEST_GENES`, which hold GRCh38
  `(contig, start, end)` gene spans.
- `TEST_GENE_PADDING = 100` — bp added to each side of a test gene's span.

### Functions

| Function | Returns | Description |
|----------|---------|-------------|
| `get_discovered_mnvs(test=False)` | `TableResource` | Discovery output path |
| `get_annotated_mnvs(test=False)` | `TableResource` | Annotated output path |
| `get_gnomad_v4_vds(...)` | `VariantDataset` | Loads unsplit VDS via `gnomad_qc` |
| `get_gnomad_v4_ukb_vds(...)` | `VariantDataset` | Loads the UKB-only VDS subset (`--ukb-only`) |
| `get_test_gene_intervals(gene_set, padding=TEST_GENE_PADDING)` | `List[Tuple[str, str]]` | Padded `(gene_name, interval)` pairs for a test gene set (`--test-genes`, and `--intervals` with no values) |
| `UKB_VDS` | `VariantDatasetResource` | Path to the pre-built UKB-only VDS subset |

### Output paths

| Resource | Production path | Test path |
|----------|-----------------|-----------|
| Discovery | `gs://gnomad/v4.1/mnv/exomes/gnomad.exomes.v4.1.mnv_discovery.ht` | `gs://gnomad-tmp/gnomad_v4.1_testing/mnv/exomes/...` |
| Annotated | `gs://gnomad/v4.1/mnv/exomes/gnomad.exomes.v4.1.mnv_annotated.ht` | `gs://gnomad-tmp/gnomad_v4.1_testing/mnv/exomes/...` |

Gene-interval runs add the gene names to the filename (e.g.
`...mnv_discovery.SHOX_XG_EMD_DDX3Y.ht`). `--intervals` and `--test-genes default`
write to the test path. `--test-genes autosomes`, `sex_chr` and `all` take long
enough that they write to `gs://gnomad/v4.1/mnv/exomes/test_genes/` instead, so the
output outlasts the test bucket's cleanup; add `--test` to send them to the test path.

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

## 8. Tests

### Local test suite

`pytest` (from the repo root, with `gnomad_qc` importable and a JDK for local Hail)
runs two kinds of test:

- `tests/test_genotypes.py` — hand-built inputs for each genotype edge case
  (hotfix, `1/2`, alt-first phased hets, AS_lowqual, sex ploidy, PAR/non-PAR, chrY)
  and for pair classification (cis, trans, unphased, hom-hom, het-hom, the 10 bp
  limit). No GCS access needed.
- `tests/test_test_data.py` (marker `gcs`) — runs the pipeline on the small test
  dataset and checks that:
  - every release SNV is counted, and any extra alt sits at an AS_lowqual site;
  - each SNV alt's adj AC, hom count and XY AC equal the v4.1.1 release's, and its raw
    AC and hom count do too (an `xfail` until gnomad_production#2300, because the raw
    counts are where the hotfix order shows);
  - pairs are within 10 bp, every pair class is present, and pair counts fit the alt
    counts;
  - the pairs match `tests/data/expected_mnv_pairs.tsv` (rewrite with
    `--update-expected`).

  These take about 3 minutes locally and need an 8 GB Spark driver heap (set
  `MNV_TEST_DRIVER_MEMORY` to change it).

Tests marked `gcs` are skipped when `CI` is set; `pytest -m "not gcs"` runs only the
fast ones. Tests marked `xfail(strict=True)` describe the v4.1-matching behaviour
decided in gnomad_production#2300 and start failing once it is implemented, so the
marker gets removed.

The test dataset is a release-only VDS cut down to `TEST_DATA_WINDOWS` (one partition,
reference blocks split at the window edges), plus slices of the release, info, VEP and
metadata Tables, under `gs://gnomad/v4.1/mnv/exomes/test_data/`. It is built by
`v4/create_test_data.py`: `--compute-site-features` scores every site in the test genes
for each test case, the windows are picked from those scores, `--build-vds` writes the
VDS (checkpointed at its own partitioning, then coalesced to one partition), and
`--build-tables` writes the Table slices.

### PCSK9 region (historical)

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
| Hail version | 0.2.11 | 0.2.128+ |
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
