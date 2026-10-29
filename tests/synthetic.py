"""Build tiny hand-written inputs for the MNV pipeline functions.

Each test describes a few sites and samples in plain Python; :func:`make_inputs`
turns them into the variant data, release, info and metadata Tables that
:func:`v4.discover_mnv.prepare_genotypes` and :func:`v4.discover_mnv.discover_mnv`
take, so edge cases can be checked without the real gnomAD data.
"""

from types import SimpleNamespace
from typing import Dict, List, Optional

import hail as hl

ENTRY_TYPE = hl.tstruct(
    LGT=hl.tcall,
    LPGT=hl.tcall,
    PID=hl.tstr,
    LA=hl.tarray(hl.tint32),
    LAD=hl.tarray(hl.tint32),
    DP=hl.tint32,
    GQ=hl.tint32,
)
"""Entry schema of the unsplit variant data the pipeline reads."""


def entry(
    gt: List[int],
    ad: List[int],
    phased: bool = False,
    pid: Optional[str] = None,
    gq: int = 40,
) -> hl.Struct:
    """
    Build one variant-data entry in local allele space.

    ``LA`` covers every allele at the site, so local and global indices match.

    :param gt: Allele indices of the call, e.g. ``[0, 1]`` or ``[1]``.
    :param ad: Allele depths, one per allele at the site.
    :param phased: Whether the call is phased (``LPGT`` is set to the same call).
    :param pid: Phase set ID. Default is None.
    :param gq: Genotype quality. Default is 40, which passes adj.
    :return: Entry struct matching :data:`ENTRY_TYPE`.
    """
    call = hl.Call(gt, phased=phased)
    return hl.Struct(
        LGT=call,
        LPGT=call if phased else None,
        PID=pid,
        LA=list(range(len(ad))),
        LAD=ad,
        DP=sum(ad),
        GQ=gq,
    )


def make_inputs(
    sites: List[dict],
    samples: Dict[str, dict],
) -> SimpleNamespace:
    """
    Turn site and sample descriptions into the pipeline's input Tables.

    :param sites: One dict per site with ``locus`` (e.g. ``"chr1:1000"``),
        ``alleles``, ``af`` (release AF per alt), optional ``lowqual`` (AS_lowqual
        per alt, default all False) and ``entries`` (sample ID to :func:`entry`).
        Every alt must already be in min-rep form against the ref.
    :param samples: Sample ID to ``{"sex_karyotype": ..., "fixed_homalt_model":
        ...}``.
    :return: Namespace with ``vds`` (only ``variant_data`` is set), ``release_ht``,
        ``info_ht`` and ``meta_ht``.
    """
    rows = []
    for site in sites:
        contig, pos = site["locus"].split(":")
        locus = hl.Locus(contig, int(pos), reference_genome="GRCh38")
        for s in samples:
            rows.append(
                hl.Struct(
                    locus=locus,
                    alleles=site["alleles"],
                    s=s,
                    e=site["entries"].get(s),
                )
            )
    ht = hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            s=hl.tstr,
            e=ENTRY_TYPE,
        ),
    )
    mt = ht.to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])
    mt = mt.select_entries(**mt.e)

    release_rows = []
    info_rows = []
    for site in sites:
        contig, pos = site["locus"].split(":")
        locus = hl.Locus(contig, int(pos), reference_genome="GRCh38")
        n_alts = len(site["alleles"]) - 1
        for alt, af in zip(site["alleles"][1:], site["af"]):
            release_rows.append(
                hl.Struct(
                    locus=locus,
                    alleles=[site["alleles"][0], alt],
                    freq=[hl.Struct(AF=af)],
                )
            )
        info_rows.append(
            hl.Struct(
                locus=locus,
                alleles=site["alleles"],
                AS_lowqual=site.get("lowqual", [False] * n_alts),
            )
        )
    release_ht = hl.Table.parallelize(
        release_rows,
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            freq=hl.tarray(hl.tstruct(AF=hl.tfloat64)),
        ),
        key=["locus", "alleles"],
    )
    info_ht = hl.Table.parallelize(
        info_rows,
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            AS_lowqual=hl.tarray(hl.tbool),
        ),
        key=["locus", "alleles"],
    )
    meta_ht = hl.Table.parallelize(
        [
            hl.Struct(
                s=s,
                sex_imputation=hl.Struct(sex_karyotype=m["sex_karyotype"]),
                project_meta=hl.Struct(fixed_homalt_model=m["fixed_homalt_model"]),
            )
            for s, m in samples.items()
        ],
        hl.tstruct(
            s=hl.tstr,
            sex_imputation=hl.tstruct(sex_karyotype=hl.tstr),
            project_meta=hl.tstruct(fixed_homalt_model=hl.tbool),
        ),
        key="s",
    )
    return SimpleNamespace(
        vds=SimpleNamespace(variant_data=mt),
        release_ht=release_ht,
        info_ht=info_ht,
        meta_ht=meta_ht,
    )
