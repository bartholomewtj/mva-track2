#!/usr/bin/env python3
"""Gene-symbol drug screen for Track 2. No genome. Public APIs only.

Queries Open Targets and ChEMBL for approved / late-phase drugs on
BUB1B and neighbouring spindle-assembly-checkpoint (SAC) genes, then
pulls openFDA label flags for unique approved names.

Usage (from repo root):
    python scripts/screen.py --out scripts/out
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_URL = "https://www.ebi.ac.uk/chembl/api/data"
FDA_URL = "https://api.fda.gov/drug/label.json"
UA = "mva-track2/0.1 (gene-symbol screen; no genome)"

# Primary disease gene, SAC neighbours, mitosis kinases often dumped as
# "checkpoint drugs", other MVA genes for the reuse example.
# Disease gene, SAC neighbours, other MVA genes, then nodes that show up
# as tempting traps (mitosis-killers, microtubules, HSP90, mTOR, SIRT2).
GENES = [
    "BUB1B",
    "BUB1",
    "BUB3",
    "MAD2L1",
    "MAD1L1",
    "CDC20",
    "CEP57",
    "TRIP13",
    "MAD2L1BP",
    "AURKB",
    "PLK1",
    "TTK",
    "TUBB",
    "HSP90AA1",
    "MTOR",
    "SIRT2",
]

# Harm-filter defaults. Inhibitors of these nodes weaken the checkpoint
# or mitotic fidelity — the opposite of what biallelic BUB1B LoF needs.
WORSEN_GENES = {"AURKB", "PLK1", "TTK", "BUB1", "BUB1B", "HSP90AA1"}
COMPENSATE_GENES = {"TUBB", "CDC20"}
TRAP_GENES = {"MTOR"}
WORSEN_MOA_BITS = (
    "inhibitor",
    "antagonist",
    "degrader",
    "suppressor",
    "blocker",
)


def http_json(url: str, payload: dict | None = None, timeout: int = 60) -> dict:
    data = None
    headers = {"User-Agent": UA, "Accept": "application/json"}
    if payload is not None:
        data = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        body = e.read().decode(errors="replace")
        try:
            parsed = json.loads(body)
        except json.JSONDecodeError:
            parsed = None
        if parsed is not None:
            return parsed
        raise RuntimeError(f"HTTP {e.code} {url}: {body[:800]}") from e


def ot_gql(query: str, variables: dict) -> dict:
    return http_json(OT_URL, {"query": query, "variables": variables})


def resolve_ensembl(symbol: str) -> str | None:
    q = """
    query Search($q: String!) {
      search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 10}) {
        hits { id name entity }
      }
    }
    """
    data = ot_gql(q, {"q": symbol})
    hits = (data.get("data") or {}).get("search", {}).get("hits") or []
    for hit in hits:
        if hit.get("entity") == "target" and (
            hit.get("name", "").upper() == symbol.upper()
            or hit.get("id", "").startswith("ENSG")
        ):
            # Prefer exact symbol match.
            if hit.get("name", "").upper() == symbol.upper():
                return hit["id"]
    for hit in hits:
        if hit.get("id", "").startswith("ENSG"):
            return hit["id"]
    return None


def ot_target(ensembl_id: str) -> dict:
    q = """
    query TargetDrugs($id: String!) {
      target(ensemblId: $id) {
        id
        approvedSymbol
        approvedName
        tractability {
          label
          modality
          value
        }
        chemicalProbes {
          id
          mechanismOfAction
          drugId
          isHighQuality
        }
        drugAndClinicalCandidates {
          count
          rows {
            id
            maxClinicalStage
            drug {
              id
              name
              drugType
              maximumClinicalStage
            }
            diseases {
              diseaseFromSource
              disease { name }
            }
          }
        }
      }
    }
    """
    data = ot_gql(q, {"id": ensembl_id})
    errs = data.get("errors")
    if errs:
        raise RuntimeError(f"Open Targets GraphQL: {errs}")
    return (data.get("data") or {}).get("target") or {}


def chembl_targets(symbol: str) -> list[dict]:
    url = (
        f"{CHEMBL_URL}/target.json?"
        + urllib.parse.urlencode({"target_synonym": symbol, "limit": 20})
    )
    data = http_json(url)
    return data.get("targets") or []


def chembl_mechanisms(target_chembl_id: str) -> list[dict]:
    url = (
        f"{CHEMBL_URL}/mechanism.json?"
        + urllib.parse.urlencode(
            {"target_chembl_id": target_chembl_id, "limit": 100}
        )
    )
    data = http_json(url)
    return data.get("mechanisms") or []


def fda_label(name: str) -> dict:
    # Brand or generic. Limit fields. One try generic, one try brand.
    params = {
        "search": f'openfda.generic_name:"{name}"',
        "limit": 1,
    }
    url = FDA_URL + "?" + urllib.parse.urlencode(params)
    try:
        data = http_json(url, timeout=30)
    except RuntimeError:
        params["search"] = f'openfda.brand_name:"{name}"'
        url = FDA_URL + "?" + urllib.parse.urlencode(params)
        try:
            data = http_json(url, timeout=30)
        except RuntimeError:
            return {"query": name, "found": False}
    results = data.get("results") or []
    if not results:
        return {"query": name, "found": False}
    row = results[0]
    openfda = row.get("openfda") or {}
    return {
        "query": name,
        "found": True,
        "brand": (openfda.get("brand_name") or [None])[0],
        "generic": (openfda.get("generic_name") or [None])[0],
        "route": openfda.get("route") or [],
        "boxed_warning": bool(row.get("boxed_warning")),
        "pediatric_use": bool(row.get("pediatric_use")),
        "pediatric_snippet": " ".join(row.get("pediatric_use") or [])[:400],
        "warnings_snippet": " ".join(row.get("warnings_and_cautions") or row.get("warnings") or [])[:400],
    }


def is_approved_stage(stage: str | None) -> bool:
    if not stage:
        return False
    s = str(stage).upper().replace(" ", "_")
    return s in {"APPROVAL", "APPROVED", "4", "PHASE_4"}


def classify(gene: str, moa: str, stage: str | None) -> str:
    moa_l = (moa or "").lower()
    gene_u = gene.upper()
    if gene_u in TRAP_GENES:
        return "trap"
    if gene_u in COMPENSATE_GENES:
        return "compensate"
    if gene_u in {"AURKB", "PLK1", "TTK", "HSP90AA1"}:
        return "worsen"
    if gene_u in WORSEN_GENES and any(bit in moa_l for bit in WORSEN_MOA_BITS):
        return "worsen"
    if gene_u == "SIRT2" and "activat" in moa_l:
        return "restore"
    if gene_u == "SIRT2":
        return "unrelated"
    if gene_u == "BUB1B" and not any(bit in moa_l for bit in WORSEN_MOA_BITS):
        return "restore" if moa_l else "unrelated"
    return "unrelated"


def unique_rows(rows: list[dict], keys: tuple[str, ...]) -> list[dict]:
    seen: set[tuple] = set()
    out = []
    for row in rows:
        sig = tuple(row.get(k) for k in keys)
        if sig in seen:
            continue
        seen.add(sig)
        out.append(row)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--skip-fda", action="store_true")
    args = p.parse_args()
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)

    ot_raw: dict[str, dict] = {}
    chembl_raw: dict[str, dict] = {}
    table: list[dict] = []
    gene_rows: list[dict] = []

    for symbol in GENES:
        print(f"Open Targets {symbol}...", file=sys.stderr)
        ens = resolve_ensembl(symbol)
        time.sleep(0.2)
        target = ot_target(ens) if ens else {}
        time.sleep(0.2)
        ot_raw[symbol] = {"ensembl_id": ens, "target": target}
        kd = target.get("drugAndClinicalCandidates") or {}
        for row in kd.get("rows") or []:
            drug = row.get("drug") or {}
            name = drug.get("name") or ""
            diseases = []
            for d in row.get("diseases") or []:
                dn = ((d.get("disease") or {}).get("name")) or d.get("diseaseFromSource")
                if dn:
                    diseases.append(dn)
            stage = row.get("maxClinicalStage") or drug.get("maximumClinicalStage") or ""
            approved = is_approved_stage(stage)
            table.append(
                {
                    "source": "open_targets",
                    "gene": target.get("approvedSymbol") or symbol,
                    "ensembl_id": ens,
                    "drug": name,
                    "chembl_id": drug.get("id") or "",
                    "moa": "",
                    "drug_type": drug.get("drugType") or "",
                    "phase": stage,
                    "max_phase": drug.get("maximumClinicalStage") or stage,
                    "approved": approved,
                    "withdrawn": False,
                    "black_box": False,
                    "disease": "; ".join(diseases[:5]),
                    "class": classify(symbol, "", stage),
                }
            )
        n_probes = len(target.get("chemicalProbes") or [])
        n_clinical = int((kd.get("count") or 0) or 0)
        print(f"ChEMBL {symbol}...", file=sys.stderr)
        try:
            targets = chembl_targets(symbol)
        except RuntimeError as e:
            print(f"  ChEMBL target search failed: {e}", file=sys.stderr)
            targets = []
        time.sleep(0.2)
        mechs: list[dict] = []
        for t in targets[:5]:
            tid = t.get("target_chembl_id")
            if not tid:
                continue
            try:
                mechs.extend(chembl_mechanisms(tid))
            except RuntimeError as e:
                print(f"  ChEMBL mechanisms {tid} failed: {e}", file=sys.stderr)
            time.sleep(0.2)
        chembl_raw[symbol] = {"targets": targets, "mechanisms": mechs}
        for m in mechs:
            moa = m.get("mechanism_of_action") or ""
            max_phase = m.get("max_phase_for_ind")
            approved = (max_phase == 4) or (str(max_phase) == "4.0")
            table.append(
                {
                    "source": "chembl",
                    "gene": symbol,
                    "ensembl_id": ens,
                    "drug": m.get("molecule_chembl_id") or "",
                    "chembl_id": m.get("molecule_chembl_id") or "",
                    "moa": moa,
                    "drug_type": m.get("action_type") or "",
                    "phase": max_phase,
                    "max_phase": max_phase,
                    "approved": approved,
                    "withdrawn": False,
                    "black_box": False,
                    "disease": "",
                    "class": classify(symbol, moa, str(max_phase) if approved else None),
                }
            )
        gene_rows.append(
            {
                "gene": target.get("approvedSymbol") or symbol,
                "ensembl_id": ens or "",
                "n_clinical": n_clinical,
                "n_probes": n_probes,
                "n_chembl_targets": len(targets),
                "n_chembl_mechanisms": len(mechs),
            }
        )

    table = unique_rows(
        table, ("source", "gene", "drug", "chembl_id", "moa", "disease")
    )

    approved_names = sorted(
        {
            r["drug"]
            for r in table
            if r["source"] == "open_targets" and r["approved"] and r["drug"]
        }
    )
    fda_raw: dict[str, dict] = {}
    if not args.skip_fda:
        for name in approved_names[:40]:
            print(f"openFDA {name}...", file=sys.stderr)
            try:
                fda_raw[name] = fda_label(name)
            except RuntimeError as e:
                fda_raw[name] = {"query": name, "found": False, "error": str(e)}
            time.sleep(0.25)

    (out / "open_targets.json").write_text(
        json.dumps(ot_raw, indent=2), encoding="utf-8"
    )
    (out / "chembl.json").write_text(
        json.dumps(chembl_raw, indent=2), encoding="utf-8"
    )
    (out / "openfda.json").write_text(
        json.dumps(fda_raw, indent=2), encoding="utf-8"
    )

    fields = [
        "source",
        "gene",
        "ensembl_id",
        "drug",
        "chembl_id",
        "moa",
        "drug_type",
        "phase",
        "max_phase",
        "approved",
        "withdrawn",
        "black_box",
        "disease",
        "class",
    ]
    with (out / "hits.tsv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(table)

    gene_fields = [
        "gene",
        "ensembl_id",
        "n_clinical",
        "n_probes",
        "n_chembl_targets",
        "n_chembl_mechanisms",
    ]
    with (out / "genes.tsv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=gene_fields, delimiter="\t")
        w.writeheader()
        w.writerows(gene_rows)

    summary = {
        "n_genes": len(GENES),
        "n_rows": len(table),
        "n_approved_names": len(approved_names),
        "by_class": {},
        "approved_names": approved_names,
    }
    for r in table:
        summary["by_class"][r["class"]] = summary["by_class"].get(r["class"], 0) + 1
    (out / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
