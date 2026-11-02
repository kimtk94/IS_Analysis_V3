#!/usr/bin/env python3
"""Local-repository citation/source and final claim audit for GSE313736 thesis.

This is a static manuscript/source integrity audit; it does not verify that
webpages are still online, reproduce a scientific experiment, or do visual
inspection of assembled PDFs.
"""
import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from urllib.parse import unquote, urlsplit

DOCS = [
    "THESIS_MANUSCRIPT_METHODS_RESULTS_DISCUSSION_DRAFT_20261101.md",
    "THESIS_PUBLICATION_READINESS_EVIDENCE_INDEX_20261102.md",
    "THESIS_FIGURES_AND_MANUSCRIPT_RESULTS_VERIFIED_20261031.md",
    "THESIS_NOVELTY_CLAIM_OVERLAP_AND_ERRATUM_AUDIT_20261101.md",
    "FIGURE_VISUAL_REVIEW_AND_DURABLE_DRIVE_ARCHIVE_20261102.md",
    "GSE313736_LIBRARY_TREATMENT_ALIASING_AND_POOL_WEIGHT_DECOMPOSITION_20261102.md",
    "GSE301517_7DAY_STRAIN_LOO_AND_EMT_STRICT_REPLICATION_20261102.md",
]
REQUIRED_SOURCE_LINKS = [
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736",
    "https://pubmed.ncbi.nlm.nih.gov/41636427/",
    "https://pubmed.ncbi.nlm.nih.gov/42895266/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC9117214/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC10272683/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC13015409/",
]
# Forbid these claims as assertions; the manuscript may discuss them as
# specifically refuted potential conclusions.
FORBIDDEN_UNQUALIFIED = [
    r"we\s+(?:identified|discovered)\s+(?:a\s+)?novel\s+(?:causal|therapeutic)\s+gene",
    r"we\s+predict(?:ed)?\s+(?:individual\s+)?steroid.induced\s+glaucoma",
    r"all\s+124\s+(?:EMT\s+)?genes\s+(?:were\s+)?(?:reproduced|replicated)\s+in\s+all\s+3\s+strains",
]


def links(text):
    # Markdown normal inline links; code-fenced examples are absent from
    # tracked publication docs, and any validated reference must be literal.
    return re.findall(r"\]\(([^)\n]+)\)", text)


def check(root):
    basedir=root/"docs"/"ocular"
    selected={}
    total_links=[]
    for name in DOCS:
        p=basedir/name
        if not p.is_file():
            raise ValueError("MISSING_REQUIRED_PUBLICATION_DOC "+str(p))
        content=p.read_text(encoding="utf8")
        if len(content)<1500:
            raise ValueError("SUSPICIOUS_EMPTY_PUBLICATION_DOC_"+name)
        selected[name]=content
        total_links.extend((name,ref) for ref in links(content))
    missing=[]
    local=external=0
    for name, ref in total_links:
        # URL in GitHub markdown is a clickable reference; check only local
        # repo paths for actual existence. Live webpage verification is outside
        # this fast deterministic CI.
        if ref.startswith(("https://","http://")):
            external+=1
            if any(c.isspace() for c in ref) or not urlsplit(ref).netloc:
                raise ValueError("INVALID_EXTERNAL_LINK_IN_"+name+" "+ref)
        elif ref.startswith(("mailto:","data:","#")):
            continue
        else:
            local+=1
            raw=unquote(ref.split("#",1)[0].split("?",1)[0])
            target=(basedir/raw).resolve()
            if not target.is_relative_to(root.resolve()):
                raise ValueError("PATH_ESCAPE_IN_MARKDOWN_REFERENCE_"+name)
            if not target.is_file():
                missing.append(name+" -> "+ref)
    if missing:
        raise ValueError("BROKEN_RELATIVE_MARKDOWN_LINKS: "+"; ".join(missing[:25]))
    m=selected[DOCS[0]]
    sources=m.split("## 6. Sources",1)
    if len(sources)!=2:
        raise ValueError("BIBLIOGRAPHY_SECTION_MISSING")
    bibliography=sources[1].split("## 7. Explicit manuscript validation gates",1)[0]
    links_to_sources=[ref for ref in links(bibliography) if ref.startswith("https://")]
    for u in REQUIRED_SOURCE_LINKS:
        if u not in links_to_sources:
            raise ValueError("BIBLIOGRAPHY_MISSING_VERIFIED_SOURCE_"+u)
    # A successful experiment-free methodology thesis does not imply the
    # original source authors' 80 named genes are available.
    must_have=[
        ("GSE313736", "GSE313736"),
        ("GSE301517", "GSE301517"),
        ("57 donor", "57 donor"),
        ("55", "55"),
        ("124/194", "124/194"),
        ("90/194", "90/194"),
        ("rank 18", "rank 18"),
        ("NON-BLOCKING", "NON-BLOCKING"),
        ("original library limit", "library"),
        ("corrected supplemental source", "41435268815"),
        ("Mehrotra erratum", "PRJNA1243143"),
    ]
    for label, needle in must_have:
        if needle.lower() not in m.lower():
            raise ValueError("RELEASE_CLAIM_MISSING_"+label)
    for pattern in FORBIDDEN_UNQUALIFIED:
        if re.search(pattern,m,re.IGNORECASE):
            raise ValueError("UNQUALIFIED_CAUSAL_OR_CLINICAL_CLAIM "+pattern)
    if "do not require author correspondence" not in m.lower():
        raise ValueError("NONBLOCKING_AUTHOR_DEPENDENCY_NOT_STATED")
    if "arbitrary" not in m.lower() or "library" not in m.lower():
        raise ValueError("TREATMENT_LIBRARY_CONFOUNDING_NOT_QUALIFIED")
    figdoc=selected["FIGURE_VISUAL_REVIEW_AND_DURABLE_DRIVE_ARCHIVE_20261102.md"]
    if "original visual fail" not in figdoc.lower() or "visual go" not in figdoc.lower():
        raise ValueError("VISUAL_REVIEW_FAILURE_HISTORY_MISSING")
    if "1Zc_0p4l-HI7digyzQmaCIu2jWBKzT6tG" not in figdoc:
        raise ValueError("CORRECTED_DRIVE_ARCHIVE_LINK_MISSING")
    index=selected["THESIS_PUBLICATION_READINESS_EVIDENCE_INDEX_20261102.md"]
    if not all(x in index.lower() for x in ("journal submission hold", "visual inspection", "bibliography")):
        raise ValueError("SUBMISSION_SCOPE_OR_EDIT_REVIEW_MISSING")
    result={
        "status":"PASS_REPOSITORY_LOCAL_MANUSCRIPT_CITATION_AND_SCOPE_GATE",
        "manuscript_source_reference_urls":links_to_sources,
        "required_scientific_sources":REQUIRED_SOURCE_LINKS,
        "documents_scanned":len(selected),
        "relative_links_resolved":local,
        "external_link_syntax_checked":external,
        "all_relative_links_exist_in_checkout":True,
        "original_author80_nonblocking":True,
        "supplementary_figure_visual_revision_present":True,
        "abstract_and_figure_counts_claims_checked":True,
        "actual_live_external_webpage_status_checked":False,
        "visual_PDF_proof_checked_by_this_automated_audit":False,
        "submission_ready":False,
        "document_source_sha256":{name:hashlib.sha256(s.encode()).hexdigest() for name,s in selected.items()},
    }
    return result


def selftest():
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as t:
        d=Path(t)
        assert links("[a](./okay.md) [b](https://example.com)") == ["./okay.md","https://example.com"]
        assert unquote("file%20name.md") == "file name.md"
        assert not (d/"../etc").resolve().is_relative_to(d.resolve())
    print("OCULAR_MANUSCRIPT_REFERENCE_GATE_SELFTEST_PASS",flush=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--repo",type=Path,default=Path("."))
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        selftest();return
    r=check(a.repo.resolve())
    r["source_git_head_sha"]=os.environ.get("GITHUB_SHA","UNSET_LOCAL")
    if a.out:
        a.out.mkdir(parents=True,exist_ok=True)
        (a.out/"OCULAR_MANUSCRIPT_REFERENCE_RELEASE_GATE.json").write_text(
            json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("OCULAR_MANUSCRIPT_REFERENCE_GATE_PASS",json.dumps({
        "n_docs":r["documents_scanned"],"relative_links":r["relative_links_resolved"],
        "external_links":r["external_link_syntax_checked"],"submission_ready":False
    }),flush=True)


if __name__=="__main__":
    main()
