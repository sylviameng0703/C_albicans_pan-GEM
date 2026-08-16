#!/usr/bin/env bash
set -euo pipefail

PROJECT=${PROJECT:-$PWD}
RUN=${RUN:-$PROJECT/results/pangenome_relaxed_80}
OUT=${OUT:-$RUN/344_final_panGEM_detail_closed_reaction_projection_for_ssGEM}

MODEL=${MODEL:-$RUN/342_final_panGEM_detail_closed_release/Calb_panGEM_final_release.xml}
RELEASE_OVERVIEW=${RELEASE_OVERVIEW:-$RUN/342_final_panGEM_detail_closed_release/final_panGEM_detail_closed_release_overview.tsv}
GENE_CATALOG=${GENE_CATALOG:-$RUN/calb_pangem_v1_evidence_package/calb_pangem_v1_gene_catalog.tsv}
SAMPLES=${SAMPLES:-$RUN/00_lists/public_main_relaxed_80.samples}
ORTHOGROUPS=${ORTHOGROUPS:-$RUN/03_orthofinder/Results_pangenome_relaxed_80/Orthogroups/Orthogroups.tsv}
PYTHON_BIN=${PYTHON_BIN:-$(command -v python3)}
EXPECTED_RELEASE_GATE=${EXPECTED_RELEASE_GATE:-pass_final_detail_closed_panGEM_release}

mkdir -p "$OUT"

for f in "$MODEL" "$RELEASE_OVERVIEW" "$GENE_CATALOG" "$SAMPLES" "$ORTHOGROUPS"; do
  [ -s "$f" ] || { echo "[ERROR] Missing: $f" >&2; exit 1; }
done

"$PYTHON_BIN" - "$MODEL" "$RELEASE_OVERVIEW" "$GENE_CATALOG" "$SAMPLES" "$ORTHOGROUPS" "$OUT" "$EXPECTED_RELEASE_GATE" <<'PY'
import csv
import hashlib
import re
import sys
from pathlib import Path
from cobra.io import read_sbml_model

model_path, release_overview_path, gene_catalog_path, samples_path, orthogroups_path, outdir = map(Path, sys.argv[1:7])
expected_release_gate = sys.argv[7]
outdir.mkdir(parents=True, exist_ok=True)

def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def write_tsv(path, rows, header):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

release = {r["metric"]: r.get("value", "") for r in read_tsv(release_overview_path)}
if release.get("release_gate") != expected_release_gate:
    raise SystemExit(f"[ERROR] final pan-GEM release gate did not pass: {release.get('release_gate', '')}")
model_sha256 = sha256(model_path)
if model_sha256 != release.get("frozen_model_sha256"):
    raise SystemExit("[ERROR] pan-GEM SHA256 does not match the frozen release overview")

samples = [x.strip() for x in samples_path.read_text().splitlines() if x.strip()]
sample_set = set(samples)

catalog = read_tsv(gene_catalog_path)
catalog_by_gene = {r["model_gene_id"]: r for r in catalog if r.get("model_gene_id")}

orthogroups = read_tsv(orthogroups_path)
og_by_id = {}
for r in orthogroups:
    og = r.get("Orthogroup") or r.get("orthogroup") or r.get("orthogroup_id")
    if og:
        og_by_id[og] = r

model = read_sbml_model(str(model_path))

gene_presence = {}
gene_rows = []

for gene in sorted(g.id for g in model.genes):
    present = set()
    source = ""
    best_og = ""

    if gene in catalog_by_gene:
        row = catalog_by_gene[gene]
        source = "gene_catalog"
        best_og = row.get("best_orthogroup", "")
        support = row.get("support_samples", "")
        if support:
            present = {s for s in support.split(";") if s in sample_set}
        elif best_og and best_og in og_by_id:
            ogrow = og_by_id[best_og]
            present = {s for s in samples if str(ogrow.get(s, "")).strip()}
    elif re.fullmatch(r"OG\d+", gene):
        source = "orthogroup_proxy_gene"
        best_og = gene
        if gene in og_by_id:
            ogrow = og_by_id[gene]
            present = {s for s in samples if str(ogrow.get(s, "")).strip()}
    elif gene.startswith("CM_"):
        source = "retained_template_mitochondrial_gene"
        present = set(samples)

    gene_presence[gene] = present
    gene_rows.append({
        "model_gene_id": gene,
        "presence_source": source,
        "best_orthogroup": best_og,
        "present_sample_count": len(present),
        "present_sample_percent": round(100 * len(present) / len(samples), 4),
        "projection_status": "resolved" if present else "unresolved_or_absent",
    })

def eval_gpr(rule, sample):
    if not rule.strip():
        return True
    toks = re.findall(r"\(|\)|\band\b|\bor\b|[A-Za-z0-9_.-]+", rule)
    expr = []
    for t in toks:
        if t in {"(", ")", "and", "or"}:
            expr.append(t)
        else:
            expr.append("True" if sample in gene_presence.get(t, set()) else "False")
    return bool(eval(" ".join(expr), {"__builtins__": {}}, {}))

pa_rows = []
catalog_rows = []
presence_by_reaction = {}

for rxn in model.reactions:
    has_gpr = bool(rxn.gene_reaction_rule.strip())
    present_samples = []
    for sample in samples:
        present = eval_gpr(rxn.gene_reaction_rule, sample) if has_gpr else True
        if present:
            present_samples.append(sample)

    presence_by_reaction[rxn.id] = set(present_samples)
    present_count = len(present_samples)

    if not has_gpr:
        category = "no_gpr_scaffold"
        scope = "template_no_gpr_scaffold_not_genome_projected"
    elif present_count == len(samples):
        category = "strict_core_gpr"
        scope = "genome_projected_core"
    else:
        category = "shell_gpr"
        scope = "genome_projected_variable"

    catalog_rows.append({
        "reaction_id": rxn.id,
        "reaction_name": rxn.name,
        "gene_reaction_rule": rxn.gene_reaction_rule,
        "has_gpr": str(has_gpr),
        "present_sample_count": present_count,
        "present_sample_percent": round(100 * present_count / len(samples), 4),
        "pan_reaction_category": category,
        "reaction_pan_category": category,
        "reaction_category": category,
        "projection_category": category,
        "evidence_level": category,
        "reaction_pan_scope": scope,
    })

for sample in samples:
    row = {"sample": sample}
    for rxn in model.reactions:
        row[rxn.id] = 1 if sample in presence_by_reaction[rxn.id] else 0
    pa_rows.append(row)

write_tsv(outdir / "locked_sample_reaction_presence_absence.tsv", pa_rows, ["sample"] + [r.id for r in model.reactions])
write_tsv(outdir / "locked_reaction_projection_catalog.tsv", catalog_rows, [
    "reaction_id", "reaction_name", "gene_reaction_rule", "has_gpr",
    "present_sample_count", "present_sample_percent",
    "pan_reaction_category", "reaction_pan_category", "reaction_category", "projection_category",
    "evidence_level", "reaction_pan_scope",
])
write_tsv(outdir / "locked_gene_projection_used_for_reaction_PA.tsv", gene_rows, [
    "model_gene_id", "presence_source", "best_orthogroup",
    "present_sample_count", "present_sample_percent", "projection_status",
])

overview = [
    ("input_model", str(model_path)),
    ("input_model_sha256", model_sha256),
    ("input_release_overview", str(release_overview_path)),
    ("input_release_gate", release.get("release_gate", "")),
    ("samples", len(samples)),
    ("model_reactions", len(model.reactions)),
    ("model_metabolites", len(model.metabolites)),
    ("model_genes", len(model.genes)),
    ("gpr_reactions", sum(bool(r.gene_reaction_rule.strip()) for r in model.reactions)),
    ("no_gpr_scaffold_reactions", sum(not bool(r.gene_reaction_rule.strip()) for r in model.reactions)),
    ("strict_core_gpr_reactions", sum(r["reaction_pan_category"] == "strict_core_gpr" for r in catalog_rows)),
    ("variable_shell_gpr_reactions", sum(r["reaction_pan_category"] == "shell_gpr" for r in catalog_rows)),
    ("unresolved_or_absent_genes", sum(r["projection_status"] != "resolved" for r in gene_rows)),
    ("output_reaction_PA", str(outdir / "locked_sample_reaction_presence_absence.tsv")),
    ("output_reaction_catalog", str(outdir / "locked_reaction_projection_catalog.tsv")),
    ("projection_gate", "pass" if all(r["projection_status"] == "resolved" for r in gene_rows) else "review_unresolved_genes"),
    ("interpretation_boundary", "genome_projected_GPR_plus_no_GPR_scaffold_template_for_final_ssGEM_generation"),
]
with open(outdir / "locked_reaction_projection_overview.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["metric", "value"])
    w.writerows(overview)

print("[DONE] locked pan-GEM reaction projection for ssGEM")
print(f"  overview: {outdir / 'locked_reaction_projection_overview.tsv'}")
PY

echo
echo "[DONE] Inspect:"
echo "  column -t -s \$'\\t' $OUT/locked_reaction_projection_overview.tsv"
