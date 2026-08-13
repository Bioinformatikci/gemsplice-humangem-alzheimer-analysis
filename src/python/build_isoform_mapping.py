#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build evidence-tagged isoform/reaction mapping candidates for the current
iMAT common Human-GEM model.

This script does not modify GEMsplice or MATLAB tool files. It only writes
mapping, fallback, and coverage reports under isoform_mapping_current_imat_common.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import re
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
WORK = Path(__file__).resolve().parents[1]
OUT = WORK / "outputs"
LOGS = WORK / "logs"
SNAP = WORK / "inputs_snapshot"

GENES_TSV = ROOT / "genes.tsv"
REACTIONS_TSV = ROOT / "reactions.tsv"
TX2GENE_TSV = ROOT / "imat_gemsplice_workspace" / "imat" / "prep_inputs" / "tx2gene.tsv"
KEPT_RXNS_TSV = ROOT / "imat_gemsplice_workspace" / "imat" / "outputs" / "imat_common_kept_reactions.tsv"
REACTION_SUMMARY_TSV = (
    ROOT
    / "imat_gemsplice_workspace"
    / "gemsplice_imat_alz_control_hgemrules_recon_reference_isoform"
    / "mapping_outputs"
    / "humanGEM_reaction_expression_source_summary.tsv"
)

OBJECTIVE_RXNS = {
    "MAR06921": "Mitochondrial Complex I / NADH:ubiquinone oxidoreductase",
    "MAR04388": "Lactate dehydrogenase",
    "MAR03802": "Glutamate dehydrogenase",
}

UNIPROT_FIELDS = "accession,id,reviewed,protein_name,gene_names,ec,cc_catalytic_activity"
USER_AGENT = "BurakKeskin-GEMspliceMapping/1.0 contact: local-analysis"


def log(message: str) -> None:
    stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    text = f"[{stamp}] {message}"
    print(text, flush=True)
    with (LOGS / "run_log.txt").open("a", encoding="utf-8") as fh:
        fh.write(text + "\n")


def strip_version(value: str) -> str:
    return re.sub(r"\..*$", "", (value or "").strip())


def split_list(value: str) -> list[str]:
    if not value:
        return []
    parts = re.split(r"[;,]\s*|\s+", value.strip())
    return [strip_version(p) for p in parts if strip_version(p)]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: normalize_cell(row.get(k, "")) for k in fieldnames})


def normalize_cell(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(x) for x in value if x is not None and str(x) != "")
    return str(value)


def request_text(url: str, timeout: int = 45) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def parse_tsv_text(text: str) -> list[dict[str, str]]:
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return []
    return list(csv.DictReader(lines, delimiter="\t"))


def chunked(values: list[str], size: int) -> list[list[str]]:
    return [values[i : i + size] for i in range(0, len(values), size)]


def join_expr(ids: list[str], fn: str = "max") -> str:
    ids = [x for x in ids if x]
    if not ids:
        return "1.0"
    expr = f"/{ids[0]}/"
    for item in ids[1:]:
        expr = f"{fn}({expr},/{item}/)"
    return expr


def parse_grrule_to_expression(rule: str, gene_expr: dict[str, str]) -> str:
    if not rule or not rule.strip():
        return ""
    tokens = re.findall(r"\(|\)|\band\b|\bor\b|[A-Za-z0-9_.:-]+", rule)
    pos = 0

    def parse_or() -> str:
        nonlocal pos
        node = parse_and()
        while pos < len(tokens) and tokens[pos].lower() == "or":
            pos += 1
            rhs = parse_and()
            node = f"max({node},{rhs})"
        return node

    def parse_and() -> str:
        nonlocal pos
        node = parse_factor()
        while pos < len(tokens) and tokens[pos].lower() == "and":
            pos += 1
            rhs = parse_factor()
            node = f"min({node},{rhs})"
        return node

    def parse_factor() -> str:
        nonlocal pos
        if pos >= len(tokens):
            return "1.0"
        tok = tokens[pos]
        if tok == "(":
            pos += 1
            node = parse_or()
            if pos < len(tokens) and tokens[pos] == ")":
                pos += 1
            return node
        pos += 1
        key = strip_version(tok)
        return gene_expr.get(key, f"/{key}/")

    try:
        return parse_or()
    except Exception:
        return ""


def extract_genes(rule: str) -> list[str]:
    if not rule:
        return []
    genes = re.findall(r"ENSG[0-9]+", rule)
    return sorted(set(strip_version(g) for g in genes))


def extract_ec_values(text: str) -> set[str]:
    if not text:
        return set()
    if isinstance(text, (list, tuple, set)):
        text = ";".join(str(x) for x in text)
    return set(re.findall(r"(?:EC[:=])?\b(\d+\.\d+\.\d+\.\d+|-)\b", text))


def extract_rhea_ids(text: str) -> list[str]:
    if not text:
        return []
    ids = []
    for raw in re.split(r"[;,]\s*|\s+", text.strip()):
        raw = raw.strip()
        if not raw:
            continue
        m = re.search(r"(\d+)", raw)
        if m:
            ids.append(m.group(1))
    return sorted(set(ids))


def build_tx2gene_maps(rows: list[dict[str, str]]) -> tuple[dict[str, str], dict[str, list[str]]]:
    tx_type: dict[str, str] = {}
    gene_to_pc_tx: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        tx = strip_version(row.get("transcript_id", ""))
        gene = strip_version(row.get("gene_id", ""))
        transcript_type = row.get("transcript_type", "") or row.get("gene_type", "")
        if not tx or not gene:
            continue
        tx_type[tx] = transcript_type
        if transcript_type == "protein_coding":
            gene_to_pc_tx[gene].append(tx)
    return tx_type, {g: sorted(set(v)) for g, v in gene_to_pc_tx.items()}


def build_gene_maps(genes_rows: list[dict[str, str]], tx_type: dict[str, str]) -> tuple[dict, list[dict[str, object]]]:
    gene_map: dict[str, dict[str, object]] = {}
    base_rows: list[dict[str, object]] = []
    for row in genes_rows:
        gene = strip_version(row.get("genes", ""))
        if not gene:
            continue
        ensts = sorted(set(split_list(row.get("geneENSTID", ""))))
        ensps = sorted(set(split_list(row.get("geneENSPID", ""))))
        uniprots = sorted(set(split_list(row.get("geneUniProtID", ""))))
        symbol = row.get("geneSymbols", "")
        entrez = row.get("geneEntrezID", "")
        pc_ensts = [tx for tx in ensts if tx_type.get(tx) == "protein_coding"]
        if not pc_ensts:
            pc_ensts = ensts

        gene_map[gene] = {
            "gene": gene,
            "ensts": ensts,
            "pc_ensts": pc_ensts,
            "ensps": ensps,
            "uniprots": uniprots,
            "symbol": symbol,
            "entrez": entrez,
        }
        for enst in ensts:
            base_rows.append(
                {
                    "ENSG": gene,
                    "ENST": enst,
                    "ENSP_candidates": ensps,
                    "UniProtID_candidates": uniprots,
                    "gene_symbol": symbol,
                    "Entrez": entrez,
                    "transcript_type": tx_type.get(enst, ""),
                    "pairing_note": "gene-level candidate; ENST/ENSP order not assumed",
                }
            )
    return gene_map, base_rows


def fetch_uniprot(uniprot_ids: list[str]) -> list[dict[str, str]]:
    cache = OUT / "UniProt_to_EC.tsv"
    if cache.exists():
        log(f"Using cached UniProt table: {cache}")
        return read_tsv(cache)

    rows: list[dict[str, str]] = []
    ids = sorted(set(x for x in uniprot_ids if x))
    log(f"Fetching UniProt annotations for {len(ids)} accessions")
    for idx, group in enumerate(chunked(ids, 45), 1):
        query = "(" + " OR ".join(f"accession:{x}" for x in group) + ")"
        params = {
            "query": query,
            "fields": UNIPROT_FIELDS,
            "format": "tsv",
            "size": str(len(group)),
        }
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)
        try:
            parsed = parse_tsv_text(request_text(url))
        except urllib.error.HTTPError as exc:
            log(f"UniProt chunk {idx} HTTP error {exc.code}; falling back to one-by-one")
            parsed = []
            for one in group:
                params["query"] = f"accession:{one}"
                params["size"] = "1"
                one_url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)
                try:
                    parsed.extend(parse_tsv_text(request_text(one_url)))
                    time.sleep(0.08)
                except Exception as inner:
                    log(f"UniProt failed for {one}: {inner}")
        except Exception as exc:
            log(f"UniProt chunk {idx} failed: {exc}")
            parsed = []

        for item in parsed:
            ec_text = item.get("EC number", "")
            cat = item.get("Catalytic activity", "")
            rows.append(
                {
                    "UniProtID": item.get("Entry", ""),
                    "EntryName": item.get("Entry Name", ""),
                    "reviewed_status": item.get("Reviewed", ""),
                    "protein_name": item.get("Protein names", ""),
                    "gene_names": item.get("Gene Names", ""),
                    "EC_number": ec_text,
                    "EC_from_catalytic_activity": sorted(extract_ec_values(cat)),
                    "catalytic_activity": cat,
                }
            )
        log(f"UniProt chunk {idx}/{len(chunked(ids, 45))}: {len(rows)} total rows")
        time.sleep(0.15)

    write_tsv(
        cache,
        rows,
        [
            "UniProtID",
            "EntryName",
            "reviewed_status",
            "protein_name",
            "gene_names",
            "EC_number",
            "EC_from_catalytic_activity",
            "catalytic_activity",
        ],
    )
    return rows


def fetch_rhea(rhea_ids: list[str]) -> list[dict[str, str]]:
    cache = OUT / "Rhea_reaction_to_EC.tsv"
    if cache.exists():
        log(f"Using cached Rhea table: {cache}")
        return read_tsv(cache)

    ids = sorted(set(extract_rhea_ids(";".join(rhea_ids))))
    log(f"Fetching Rhea annotations for {len(ids)} Rhea IDs")
    rows: list[dict[str, str]] = []
    for i, rhea_id in enumerate(ids, 1):
        query = f"RHEA:{rhea_id}"
        params = {
            "query": query,
            "columns": "rhea-id,equation,ec,reaction-xref(KEGG)",
            "format": "tsv",
            "limit": "10",
        }
        url = "https://www.rhea-db.org/rhea?" + urllib.parse.urlencode(params)
        try:
            parsed = parse_tsv_text(request_text(url))
        except Exception as exc:
            log(f"Rhea failed for {query}: {exc}")
            parsed = []
        if not parsed:
            rows.append(
                {
                    "query_rhea_id": f"RHEA:{rhea_id}",
                    "rhea_id": "",
                    "equation": "",
                    "EC_number": "",
                    "KEGG_reaction": "",
                    "fetch_status": "not_found_or_failed",
                }
            )
        for item in parsed:
            rows.append(
                {
                    "query_rhea_id": f"RHEA:{rhea_id}",
                    "rhea_id": item.get("Reaction identifier", ""),
                    "equation": item.get("Equation", ""),
                    "EC_number": item.get("EC number", ""),
                    "KEGG_reaction": item.get("Cross-reference (KEGG)", ""),
                    "fetch_status": "ok",
                }
            )
        if i % 50 == 0:
            log(f"Rhea {i}/{len(ids)}")
        time.sleep(0.08)

    write_tsv(
        cache,
        rows,
        ["query_rhea_id", "rhea_id", "equation", "EC_number", "KEGG_reaction", "fetch_status"],
    )
    return rows


def build_annotation_and_mapping() -> None:
    LOGS.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    SNAP.mkdir(parents=True, exist_ok=True)
    (LOGS / "run_log.txt").write_text("", encoding="utf-8")

    required = [GENES_TSV, REACTIONS_TSV, TX2GENE_TSV, KEPT_RXNS_TSV, REACTION_SUMMARY_TSV]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(f"Required input not found: {path}")
        shutil.copy2(path, SNAP / path.name)

    log("Reading input TSV files")
    genes_rows = read_tsv(GENES_TSV)
    reactions_rows = read_tsv(REACTIONS_TSV)
    tx2gene_rows = read_tsv(TX2GENE_TSV)
    kept_rows = read_tsv(KEPT_RXNS_TSV)
    summary_rows = read_tsv(REACTION_SUMMARY_TSV)

    tx_type, gene_to_pc_tx = build_tx2gene_maps(tx2gene_rows)
    gene_map, base_rows = build_gene_maps(genes_rows, tx_type)
    write_tsv(
        OUT / "ENST_ENSP_UniProt_base.tsv",
        base_rows,
        [
            "ENSG",
            "ENST",
            "ENSP_candidates",
            "UniProtID_candidates",
            "gene_symbol",
            "Entrez",
            "transcript_type",
            "pairing_note",
        ],
    )
    log(f"Built base ENST/ENSP/UniProt table: {len(base_rows)} transcript rows")

    reaction_ann = {r.get("rxns", ""): r for r in reactions_rows}
    summary_by_rxn = {r.get("human_rxn", ""): r for r in summary_rows}
    kept_rxns = [r.get("rxn_id", "") for r in kept_rows if r.get("rxn_id", "")]
    kept_set = set(kept_rxns)
    kept_summary = [summary_by_rxn[rxn] for rxn in kept_rxns if rxn in summary_by_rxn]
    log(f"Kept reactions: {len(kept_rxns)}; with grRules summary: {len(kept_summary)}")

    curated_rules = []
    all_genes_in_curated = set()
    rhea_ids = []
    for row in kept_rows:
        rxn = row.get("rxn_id", "")
        summary = summary_by_rxn.get(rxn, {})
        gr = summary.get("human_grRule", "")
        genes = extract_genes(gr)
        all_genes_in_curated.update(genes)
        ann = reaction_ann.get(rxn, {})
        rhea_ids.extend(split_list(ann.get("rxnRheaID", "")))
        rhea_ids.extend(split_list(ann.get("rxnRheaMasterID", "")))
        curated_rules.append(
            {
                "rxn": rxn,
                "rxnName": row.get("rxn_name", "") or summary.get("human_name", ""),
                "grRule": gr,
                "genes": genes,
                "subSystem": summary.get("human_subsystem", ""),
                "model_eccodes": summary.get("human_eccodes", ""),
            }
        )
    write_tsv(
        OUT / "curated_reaction_gene_rules.tsv",
        curated_rules,
        ["rxn", "rxnName", "grRule", "genes", "subSystem", "model_eccodes"],
    )

    uniprots = []
    for gene in all_genes_in_curated:
        uniprots.extend(gene_map.get(gene, {}).get("uniprots", []))
    uniprot_rows = fetch_uniprot(uniprots)
    uniprot_by_id = {r.get("UniProtID", ""): r for r in uniprot_rows}

    rhea_rows = fetch_rhea(rhea_ids)
    rhea_ec_by_query = defaultdict(set)
    for row in rhea_rows:
        q = row.get("query_rhea_id", "")
        for ec in extract_ec_values(row.get("EC_number", "")):
            rhea_ec_by_query[q].add(ec)

    curated_annotation = []
    reaction_expression_rows = []
    enst_to_reaction_rows = []
    fallback_rows = []
    manual_rows = []
    objective_rows = []
    summary_counter = Counter()

    for rule_row in curated_rules:
        rxn = rule_row["rxn"]
        ann = reaction_ann.get(rxn, {})
        genes = rule_row["genes"]
        reaction_rhea_ids = split_list(ann.get("rxnRheaID", "")) + split_list(ann.get("rxnRheaMasterID", ""))
        reaction_ec = set(extract_ec_values(rule_row.get("model_eccodes", "")))
        for rid in reaction_rhea_ids:
            rid_num = extract_rhea_ids(rid)
            for num in rid_num:
                reaction_ec.update(rhea_ec_by_query.get(f"RHEA:{num}", set()))

        gene_expr = {}
        gene_evidence = {}
        reaction_specific_hits = 0
        pc_fallback_hits = 0
        gene_level_hits = 0
        default_hits = 0
        all_uniprots = []
        all_uniprot_ec = set()

        for gene in genes:
            info = gene_map.get(gene, {})
            gene_uniprots = list(info.get("uniprots", []))
            all_uniprots.extend(gene_uniprots)
            gene_ec = set()
            reviewed_status = []
            for up in gene_uniprots:
                urow = uniprot_by_id.get(up, {})
                reviewed_status.append(urow.get("reviewed_status", ""))
                gene_ec.update(extract_ec_values(urow.get("EC_number", "")))
                gene_ec.update(extract_ec_values(urow.get("EC_from_catalytic_activity", "")))
            all_uniprot_ec.update(gene_ec)
            matching_ec = sorted(gene_ec & reaction_ec)

            ensts = list(info.get("pc_ensts", [])) or gene_to_pc_tx.get(gene, []) or list(info.get("ensts", []))
            ensts = sorted(set(ensts))
            if matching_ec and ensts:
                gene_expr[gene] = join_expr(ensts)
                gene_evidence[gene] = "reaction_specific_ec_match"
                reaction_specific_hits += 1
            elif ensts:
                gene_expr[gene] = join_expr(ensts)
                gene_evidence[gene] = "protein_coding_enst_fallback"
                pc_fallback_hits += 1
            elif gene:
                gene_expr[gene] = f"/{gene}/"
                gene_evidence[gene] = "gene_level_fallback"
                gene_level_hits += 1
            else:
                gene_evidence[gene] = "default_1"
                default_hits += 1

            for enst in ensts:
                enst_to_reaction_rows.append(
                    {
                        "rxn": rxn,
                        "rxnName": rule_row["rxnName"],
                        "ENSG": gene,
                        "ENST": enst,
                        "ENSP_candidates": info.get("ensps", []),
                        "UniProtID_candidates": gene_uniprots,
                        "gene_symbol": info.get("symbol", ""),
                        "gene_EC": sorted(gene_ec),
                        "reaction_EC": sorted(reaction_ec),
                        "matched_EC": matching_ec,
                        "evidence": gene_evidence.get(gene, ""),
                        "reviewed_status": sorted(set(x for x in reviewed_status if x)),
                    }
                )

        if not genes:
            expr = "1.0"
            mapping_source = "default_1"
            default_hits += 1
        else:
            expr = parse_grrule_to_expression(rule_row.get("grRule", ""), gene_expr)
            if not expr:
                expr = join_expr([tx for g in genes for tx in gene_map.get(g, {}).get("pc_ensts", [])])
            if reaction_specific_hits:
                mapping_source = "reaction_specific_ec_match"
            elif pc_fallback_hits:
                mapping_source = "protein_coding_enst_fallback"
            elif gene_level_hits:
                mapping_source = "gene_level_fallback"
            else:
                mapping_source = "default_1"

        summary_counter[mapping_source] += 1
        needs_review = mapping_source != "reaction_specific_ec_match"
        if needs_review:
            manual_rows.append(
                {
                    "rxn": rxn,
                    "rxnName": rule_row["rxnName"],
                    "grRule": rule_row.get("grRule", ""),
                    "genes": genes,
                    "reaction_EC": sorted(reaction_ec),
                    "gene_UniProtIDs": sorted(set(all_uniprots)),
                    "gene_UniProt_EC": sorted(all_uniprot_ec),
                    "current_mapping_source": mapping_source,
                    "reason": "no reaction-specific EC match; fallback used",
                }
            )

        curated_annotation.append(
            {
                "rxn": rxn,
                "rxnName": rule_row["rxnName"],
                "grRule": rule_row.get("grRule", ""),
                "genes": genes,
                "UniProtIDs": sorted(set(all_uniprots)),
                "rxnKEGGID": ann.get("rxnKEGGID", ""),
                "rxnRheaID": ann.get("rxnRheaID", ""),
                "rxnRheaMasterID": ann.get("rxnRheaMasterID", ""),
                "rxnMetaNetXID": ann.get("rxnMetaNetXID", ""),
                "rxnRecon3DID": ann.get("rxnRecon3DID", ""),
                "model_eccodes": rule_row.get("model_eccodes", ""),
                "reaction_candidate_EC": sorted(reaction_ec),
                "gene_UniProt_EC": sorted(all_uniprot_ec),
                "mapping_source": mapping_source,
            }
        )
        reaction_expression_rows.append(
            {
                "rxn": rxn,
                "rxnName": rule_row["rxnName"],
                "mapping_source": mapping_source,
                "reaction_expression_candidate": expr,
                "reaction_specific_gene_hits": reaction_specific_hits,
                "protein_coding_fallback_gene_hits": pc_fallback_hits,
                "gene_level_fallback_hits": gene_level_hits,
                "default_hits": default_hits,
                "needs_manual_review": "1" if needs_review else "0",
            }
        )
        fallback_rows.append(
            {
                "rxn": rxn,
                "rxnName": rule_row["rxnName"],
                "mapping_source": mapping_source,
                "reaction_specific_gene_hits": reaction_specific_hits,
                "protein_coding_fallback_gene_hits": pc_fallback_hits,
                "gene_level_fallback_hits": gene_level_hits,
                "default_hits": default_hits,
                "needs_manual_review": "1" if needs_review else "0",
            }
        )
        if rxn in OBJECTIVE_RXNS:
            objective_rows.append(
                {
                    "rxn": rxn,
                    "objective": OBJECTIVE_RXNS[rxn],
                    "rxnName": rule_row["rxnName"],
                    "grRule": rule_row.get("grRule", ""),
                    "reaction_EC": sorted(reaction_ec),
                    "gene_UniProt_EC": sorted(all_uniprot_ec),
                    "mapping_source": mapping_source,
                    "reaction_expression_candidate": expr,
                    "needs_manual_review": "1" if needs_review else "0",
                }
            )

    write_tsv(
        OUT / "curated_reaction_annotation.tsv",
        curated_annotation,
        [
            "rxn",
            "rxnName",
            "grRule",
            "genes",
            "UniProtIDs",
            "rxnKEGGID",
            "rxnRheaID",
            "rxnRheaMasterID",
            "rxnMetaNetXID",
            "rxnRecon3DID",
            "model_eccodes",
            "reaction_candidate_EC",
            "gene_UniProt_EC",
            "mapping_source",
        ],
    )
    write_tsv(
        OUT / "ENST_to_HumanGEM_reaction.tsv",
        enst_to_reaction_rows,
        [
            "rxn",
            "rxnName",
            "ENSG",
            "ENST",
            "ENSP_candidates",
            "UniProtID_candidates",
            "gene_symbol",
            "gene_EC",
            "reaction_EC",
            "matched_EC",
            "evidence",
            "reviewed_status",
        ],
    )
    write_tsv(
        OUT / "reaction_expression_candidates.tsv",
        reaction_expression_rows,
        [
            "rxn",
            "rxnName",
            "mapping_source",
            "reaction_expression_candidate",
            "reaction_specific_gene_hits",
            "protein_coding_fallback_gene_hits",
            "gene_level_fallback_hits",
            "default_hits",
            "needs_manual_review",
        ],
    )
    write_tsv(
        OUT / "fallback_report.tsv",
        fallback_rows,
        [
            "rxn",
            "rxnName",
            "mapping_source",
            "reaction_specific_gene_hits",
            "protein_coding_fallback_gene_hits",
            "gene_level_fallback_hits",
            "default_hits",
            "needs_manual_review",
        ],
    )
    write_tsv(
        OUT / "manual_review_needed.tsv",
        manual_rows,
        [
            "rxn",
            "rxnName",
            "grRule",
            "genes",
            "reaction_EC",
            "gene_UniProtIDs",
            "gene_UniProt_EC",
            "current_mapping_source",
            "reason",
        ],
    )
    write_tsv(
        OUT / "objective_reaction_mapping_check.tsv",
        objective_rows,
        [
            "rxn",
            "objective",
            "rxnName",
            "grRule",
            "reaction_EC",
            "gene_UniProt_EC",
            "mapping_source",
            "reaction_expression_candidate",
            "needs_manual_review",
        ],
    )

    summary_rows = [
        {"metric": "kept_reactions", "value": len(kept_rxns)},
        {"metric": "genes_in_curated_grRules", "value": len(all_genes_in_curated)},
        {"metric": "base_transcript_rows", "value": len(base_rows)},
        {"metric": "unique_uniprots_in_curated_genes", "value": len(set(uniprots))},
        {"metric": "uniprot_rows_returned", "value": len(uniprot_rows)},
        {"metric": "rhea_rows_returned", "value": len(rhea_rows)},
        {"metric": "reaction_specific_ec_match", "value": summary_counter["reaction_specific_ec_match"]},
        {"metric": "protein_coding_enst_fallback", "value": summary_counter["protein_coding_enst_fallback"]},
        {"metric": "gene_level_fallback", "value": summary_counter["gene_level_fallback"]},
        {"metric": "default_1", "value": summary_counter["default_1"]},
        {"metric": "manual_review_needed", "value": len(manual_rows)},
    ]
    write_tsv(OUT / "mapping_summary.tsv", summary_rows, ["metric", "value"])
    (OUT / "mapping_summary.json").write_text(json.dumps(summary_rows, indent=2), encoding="utf-8")

    log("Pipeline completed")
    log(f"Summary: {dict(summary_counter)}; manual review: {len(manual_rows)}")


def main() -> int:
    try:
        build_annotation_and_mapping()
        return 0
    except Exception as exc:
        log(f"ERROR: {type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
