# Known issues — v4 MNV discovery

## 1. min-rep locus shift vs. site-locus scan window

**Status:** known, un-fixed by design (a post-write guard *warns* but does not correct
it). Not observed in validation to date. Documented here so it can be addressed if/when
it matters. Owner decides whether to fix.

### One-paragraph summary

`discover_mnv` scans variants in **site-locus order** and decides which variants are
close enough to pair using the **site locus** (`_in_window`), but it *emits* each
variant's **min-repped** `(locus, alleles)` (`_get_carried_alts` calls `hl.min_rep` so
the output keys join the split/min-repped release + VEP HTs). `hl.min_rep` can shift a
locus **forward**, so the distance the scan uses (site locus) and the distance in the
output (min-repped locus) can disagree. At multiallelic sites where a SNP alt sits inside
a *left-padded* reference allele, this produces two edge-case defects:

1. **Missed pairs (false negatives):** a valid MNV whose min-repped distance ≤
   `max_distance` but whose *site* distance is greater is never formed by the scan.
2. **Unfiltered out-of-range pairs (false positives):** a pair whose min-repped distance
   falls outside `[1, max_distance]` is still written; the guard only logs a warning and
   `--annotate` consumes the rows as-is.

**Why it hasn't bitten:** across all validation runs (PCNT and AHNAK2, ukb-only and full
cohort, `max_distance` 2 and 10) every multiallelic SNP alt differed at the **first** base
of the ref, so `hl.min_rep` shift = 0 and site distance == min-repped distance. The issue
only arises when a SNP is left-padded — i.e. co-located at a multiallelic site with an
indel that anchors the site's ref start to the *left* of the SNP.

### Background: two coordinate systems

- **Scan** (`_scan_for_candidates`): walks the localized table one **site** row at a
  time in locus order, maintaining a per-sample window of recent non-ref entries.
  `_in_window` keeps an entry while
  `current_site_locus.position - stored_site_locus.position <= max_distance`.
- **Emit** (`_get_carried_alts`): for each carried alt it builds a biallelic `[ref, alt]`
  and calls `hl.min_rep(site_locus, [ref, alt])`, storing the **min-repped** locus and
  alleles. This is deliberate — without it, alts at multiallelic sites carry site context
  (e.g. `["TCCGGG","GCCGGG"]`) and the release/VEP/hotfix-AF joins keyed on
  `(locus, alleles)` miss. See commit `1453681`.
- `hl.min_rep` trims common prefix/suffix. Trimming a common **prefix** shifts the locus
  **forward** by the prefix length; suffix trimming leaves it. So
  `min_repped_locus = site_locus + shift`, with `shift >= 0`.

### The window math (answers "what should the window be?")

For a candidate pair (prev = earlier site, cur = later site):

```
min_rep_dist = site_dist + shift_cur - shift_prev          (shift >= 0)
```

The scan only forms the pair when `site_dist <= window`. To never miss a pair with
`min_rep_dist <= max_distance`:

```
site_dist = min_rep_dist - shift_cur + shift_prev  <=  max_distance + shift_prev
  =>  window must be >= max_distance + max_shift
```

where `max_shift = max(len(ref) - 1)` over SNP-candidate sites (a SNP differing at the
last base of an L-length padded ref shifts by `L - 1`).

**So for `max_distance = 2` the scan window is NOT 2 — it must be `2 + max_shift`.** In the
common case (SNP at the first base → `shift = 0`) that is still 2, which is why validation
matched the oracle exactly. With left-padded SNPs it is larger, and one co-located large
deletion (long ref) can make `max_shift` big — which is why a naive "use the global max"
is impractical (it would inflate the per-sample window for the whole scan).

A rarer related consequence: min-rep can **reorder** two variants (a later site
min-repping *before* an earlier one when `shift_prev - shift_cur > site_dist`); a
site-order scan cannot form that pair at all. `max_distance = 2` bounds the real-world
impact of both.

### Concrete example (constructed)

Multiallelic site at position `P` with alleles `["GCCGGG", "GTCGGG", <some indel>]`. The
SNP alt `GTCGGG` differs from ref at position 1, so
`hl.min_rep(P, ["GCCGGG","GTCGGG"]) = (P+1, ["C","T"])` (shift 1). A clean SNP at site
`P+3`:

- Site distance `= 3 > max_distance = 2` → the scan never forms the pair → **missed**
  (defect 1), and the guard cannot see it (it was never emitted).
- If a shift instead pushed an *emitted* pair's distance above `max_distance`, the guard
  would warn but not drop it (defect 2), and it would flow into `--annotate`.

### Current mitigation

- `_aggregate_mnv_pairs` computes `dist` from the min-repped loci.
- A post-write guard in `main()` reads the written discovery HT back and
  `logger.warning`s the count of pairs with `dist` outside `[1, max_distance]`
  (commit `46c2031`). **It only detects already-emitted (false-positive) pairs; it does
  not filter them, and it cannot see never-formed (false-negative) pairs.**

### How it was found

An independent `split_multi` "oracle" (classic biallelic per-sample pairing, entirely
partition- and scan-independent) matched the pipeline exactly on the validated regions on
min-repped keys. A subsequent two-model code review (Opus + Sonnet) then independently
flagged both symptoms; both models agreed on both when scored on a neutral validation
prompt.

### Options to fix (if you decide to)

**A — widen scan window + filter (small, bounded-correct):**
1. Compute `max_shift = max(hl.len(ref) - 1)` over SNP-candidate rows, **capped** at e.g.
   20 bp (cheap row-level aggregation).
2. Set the `_in_window` threshold to `max_distance + max_shift`.
3. In classification (`_classify_mnv_pairs`), keep only pairs with
   `1 <= (min-repped) dist <= max_distance`. This *replaces* the warn-only guard —
   invalid pairs are dropped, not just flagged (closes defect 2).
4. `logger.warning` if any site's ref exceeds the cap (those far-shifted SNPs may still be
   missed).
   - No-op on first-base data (window stays `max_distance`). Residual gap: SNPs shifted
     beyond the cap, and the reordering case above.

**B — scan on true (min-repped) positions (larger, fully correct):**
Explode each site into one min-repped biallelic variant per carried alt **before** the
scan, re-key/sort by the min-repped locus, and scan on that. Then the window is
`max_distance` exactly — no shift math, no reordering gap. Cost: a re-sort/shuffle and one
row per `(site, alt)` instead of per site (closer to, but not, a full `split_multi`).

Recommendation from the analysis: **A** with a cap of ~20 — a ~10-line change that is a
no-op on the data validated so far and closes defect 2 cleanly; **B** is correct by
construction but a meaningful rewrite of the scan input for a gap that is essentially
theoretical at `max_distance = 2`.

### How to re-verify a fix

Rebuild the `split_multi` oracle: load the VDS for a region, `hl.vds.split_multi` +
`hl.vds.to_dense_mt`, apply the same preprocessing as discovery (all-lowqual-site drop,
hom-alt hotfix; sex ploidy is a no-op on autosomes; adj not needed for raw counts), then
pair carried SNVs per sample via `group_by(s)` + explicit within-`max_distance` pairing on
the **split (min-repped)** loci. Compare the min-repped-keyed counts to discovery. To
actually exercise defect 1, pick/construct a region containing a left-padded multiallelic
SNP (a SNP co-located with an indel), not just first-base SNPs.

### Code pointers

- `_in_window` (in `_scan_for_candidates`) — the site-locus window; root of the mismatch.
- `_get_carried_alts` — the `hl.min_rep` call that introduces the forward shift.
- `_aggregate_mnv_pairs` — `dist` computed from min-repped loci.
- the guard block in `main()` (immediately after the discovery `.write(...)`) — warn-only.
