#!/usr/bin/env python3
"""GRCh37 REF/ALT contract for pre-GWAS-SSF GIGASTROKE summary statistics.

The source provides *effect_allele* and *other_allele*, not REF and ALT.
Do not interpret historical canonical variant_id = chr:pos:other:effect
as a reference-anchored variant identifier.

This module never mutates source or canonical files and never performs
liftover. FASTA base must be supplied by an independent GRCh37 lookup.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping


class VariantContractError(ValueError):
    """Allele or coordinate cannot be safely FASTA-harmonized."""


@dataclass(frozen=True)
class HarmonizedSNV:
    chrom: str
    pos: int
    ref: str
    alt: str
    variant_id: str
    effect_allele: str
    other_allele: str
    beta_alt: float
    eaf_alt: float | None
    effect_orientation: str
    historical_id: str | None
    historical_id_status: str


def normalize_record(row: Mapping[str, str], fasta_ref: str) -> HarmonizedSNV:
    chrom=str(row.get("chr",row.get("chromosome",""))).removeprefix("chr").strip()
    if not chrom or not (chrom.isdigit() or chrom in ("X","Y","MT")):
        raise VariantContractError("Chromosome missing or invalid")
    try:
        pos=int(row.get("pos",row.get("base_pair_location","")))
    except (ValueError,TypeError) as exc:
        raise VariantContractError("Position missing or invalid") from exc
    if pos<=0:raise VariantContractError("Position must be positive")
    ea=str(row.get("effect_allele","")).upper()
    oa=str(row.get("other_allele","")).upper()
    ref=str(fasta_ref).upper()
    if any(len(x)!=1 or x not in ("A","C","G","T") for x in (ea,oa,ref)):
        raise VariantContractError("This contract handles biallelic ACGT SNVs only")
    if ea==oa or ref not in (ea,oa):
        raise VariantContractError("Neither effect nor other allele matches FASTA REF")
    if row.get("build","GRCh37")!="GRCh37":
        raise VariantContractError("Expected GRCh37, no liftover performed")
    beta=float(row.get("beta","nan"))
    if not math.isfinite(beta):
        raise VariantContractError("Beta missing/nonfinite")
    eaf_value=row.get("eaf",row.get("effect_allele_frequency",""))
    if eaf_value is None or str(eaf_value).strip()=="":
        eaf=None
    else:
        eaf=float(eaf_value)
        if not math.isfinite(eaf) or not 0<=eaf<=1:
            raise VariantContractError("EAF out of range")
    if ea==ref:
        alt=oa
        beta_alt=-beta
        eaf_alt=None if eaf is None else 1-eaf
        orientation="EFFECT_REF_FLIPPED"
    else:
        alt=ea
        beta_alt=beta
        eaf_alt=eaf
        orientation="EFFECT_ALT"
    validated=f"{chrom}:{pos}:{ref}:{alt}"
    historical=row.get("variant_id")
    historical_status=("NOT_SUPPLIED" if not historical else
                       "FASTA_REF_ALT" if historical==validated else
                       "OTHER_EFFECT_ORDER" if historical==f"{chrom}:{pos}:{oa}:{ea}" else
                       "UNEXPECTED_ID")
    if historical_status=="UNEXPECTED_ID":
        raise VariantContractError("Historical ID disagrees with source allele pair/position")
    return HarmonizedSNV(
        chrom=chrom,pos=pos,ref=ref,alt=alt,variant_id=validated,
        effect_allele=ea,other_allele=oa,beta_alt=beta_alt,eaf_alt=eaf_alt,
        effect_orientation=orientation,historical_id=historical,
        historical_id_status=historical_status)

__all__=["HarmonizedSNV","VariantContractError","normalize_record"]
