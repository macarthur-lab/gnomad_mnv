# gnomAD v2 MNV pipeline

Scripts, utilities, and tutorials for the gnomAD v2.1 MNV paper (GRCh37, Hail 0.2).

- **`code/`** — All analysis scripts (discovery, annotation, mechanism, phasing, release). See `code/readme.md` for categories and usage.
- **`util/`** — Shared functions (e.g. `mnv_functions.py`: MNV categories, heatmaps, null model).
- **`tutorials/`** — Six Jupyter notebooks (identify MNV, annotate, functional impact, global/per-region mechanisms, phase sensitivity). Use Hail 0.2 and gnomAD 2.1 paths; see root README for paper figure mapping.

**User-runnable (your own VCF):**

- `code/get_mnv.py <vcf_path>` — Discover MNVs from a VCF.
- `code/annotate_vep_mnv.py <mnv_hail_table> <distance>` — Annotate MNVs with VEP and category (distance 1 or 2).

Most other scripts depend on gnomAD internal data and paths. See the root [README.md](../README.md) and [DOCUMENTATION.md](../DOCUMENTATION.md) for the full pipeline and paper mapping.
