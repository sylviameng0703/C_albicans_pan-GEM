#!/usr/bin/env python3
"""Generate final draft ssGEMs from the evidence-gated C. albicans pan-GEM.

Policy:
1. Start from the curated C. albicans pan-GEM template.
2. Remove sample-absent reactions in genome-projected removable classes.
3. If biomass growth is lost, restore the smallest greedy shell_gpr rescue set
   needed for nonzero growth. Restored reactions are explicitly labeled as
   growth rescue/gapfill, not as genome-supported strain presence.
4. Keep absent_gpr_projection removed by default; do not rescue it unless
   explicitly requested.

Strain-specific pruning is allowed, while nonzero growth and network consistency
are enforced and documented.
"""

from __future__ import annotations

import argparse
import hashlib
import warnings
from pathlib import Path

import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model


def clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "null"}:
        return ""
    return text


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input table: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).fillna("")


def read_model(path: Path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return read_sbml_model(str(path))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_set(value: str) -> set[str]:
    return {x.strip() for x in clean(value).replace(",", ";").split(";") if x.strip()}


def parse_samples(value: str, available: list[str]) -> list[str]:
    text = clean(value)
    if not text or text.lower() in {"all", "*"}:
        return available
    requested = [x.strip() for x in text.replace(",", ";").split(";") if x.strip()]
    missing = sorted(set(requested) - set(available))
    if missing:
        raise ValueError(f"Requested samples not found in reaction PA matrix: {missing}")
    return requested


def parse_float(value: object) -> float:
    try:
        return float(clean(value))
    except Exception:
        return 0.0


def is_absent(value: object) -> bool:
    return clean(value).lower() in {"", "0", "0.0", "false", "no", "absent", "na"}


def optimize(model, growth_eps: float) -> tuple[str, float, str]:
    try:
        sol = model.optimize()
        status = clean(sol.status)
        value = float(sol.objective_value) if status == "optimal" and sol.objective_value is not None else 0.0
    except Exception as exc:
        status = f"solver_error:{type(exc).__name__}"
        value = 0.0
    if status == "optimal" and value > growth_eps:
        call = "growth"
    elif status == "optimal":
        call = "zero_growth"
    else:
        call = "solver_or_infeasible"
    return status, value, call


def close_uptakes(model) -> None:
    for reaction in model.reactions:
        is_exchange = reaction.boundary or reaction.id.startswith(("Ex_", "EX_", "calb_EX_", "DM_", "SK_"))
        if is_exchange and reaction.lower_bound < 0:
            reaction.lower_bound = 0.0


def max_closed_uptake_atpm(model, atpm_id: str) -> tuple[str, float, str]:
    if atpm_id not in model.reactions:
        return "not_tested", 0.0, "atpm_reaction_missing"
    with model:
        close_uptakes(model)
        atpm = model.reactions.get_by_id(atpm_id)
        atpm.lower_bound = 0.0
        atpm.upper_bound = max(float(atpm.upper_bound), 1000.0)
        model.objective = atpm
        try:
            sol = model.optimize()
            status = clean(sol.status)
            value = float(sol.objective_value) if status == "optimal" and sol.objective_value is not None else 0.0
        except Exception as exc:
            return f"solver_error:{type(exc).__name__}", 0.0, "solver_error"
    call = "pass_no_closed_uptake_atpm_cycle" if status == "optimal" and value <= 1e-9 else "review_possible_cycle"
    return status, value, call


def reaction_prevalence(catalog_map: dict[str, dict[str, str]], rid: str) -> tuple[float, float]:
    row = catalog_map.get(rid, {})
    count = parse_float(row.get("present_sample_count", row.get("present_sample_count_numeric", "")))
    percent = parse_float(row.get("present_sample_percent", row.get("present_sample_percent_numeric", "")))
    return count, percent


def sorted_rescue_candidates(
    removed_ids: list[str],
    catalog_map: dict[str, dict[str, str]],
    rescue_categories: set[str],
) -> list[str]:
    candidates = []
    for rid in removed_ids:
        category = clean(catalog_map.get(rid, {}).get("pan_reaction_category", ""))
        if category not in rescue_categories:
            continue
        count, percent = reaction_prevalence(catalog_map, rid)
        candidates.append((rid, count, percent))
    candidates.sort(key=lambda x: (-x[1], -x[2], x[0]))
    return [x[0] for x in candidates]


def restore_reaction(model, template, rid: str) -> None:
    if rid not in model.reactions and rid in template.reactions:
        model.add_reactions([template.reactions.get_by_id(rid).copy()])


def remove_reaction(model, rid: str) -> None:
    if rid in model.reactions:
        model.remove_reactions([model.reactions.get_by_id(rid)], remove_orphans=False)


def cumulative_addback_until_growth(
    pruned,
    template,
    candidate_ids: list[str],
    growth_eps: float,
    max_addback: int,
) -> tuple[list[str], list[dict[str, object]], str, float, str]:
    model = pruned.copy()
    trace = []
    status, value, call = optimize(model, growth_eps)
    if call == "growth":
        return [], trace, status, value, call
    limit = min(max_addback, len(candidate_ids)) if max_addback > 0 else len(candidate_ids)
    restored = []
    for rank, rid in enumerate(candidate_ids[:limit], start=1):
        restore_reaction(model, template, rid)
        restored.append(rid)
        status, value, call = optimize(model, growth_eps)
        trace.append(
            {
                "addback_rank": rank,
                "reaction_id": rid,
                "objective_status": status,
                "objective_value": value,
                "growth_call": call,
            }
        )
        if call == "growth":
            return restored, trace, status, value, call
    return restored, trace, status, value, call


def greedy_minimize_rescue_set(
    pruned,
    template,
    rescue_ids: list[str],
    growth_eps: float,
) -> tuple[list[str], list[dict[str, object]], str, float, str]:
    model = pruned.copy()
    for rid in rescue_ids:
        restore_reaction(model, template, rid)
    status, value, call = optimize(model, growth_eps)
    if call != "growth":
        return rescue_ids, [], status, value, call

    kept = set(rescue_ids)
    trace = []
    # Reverse prevalence order first: the last addbacks are most likely removable.
    for rid in reversed(rescue_ids):
        if rid not in kept:
            continue
        test = model.copy()
        remove_reaction(test, rid)
        test_status, test_value, test_call = optimize(test, growth_eps)
        if test_call == "growth":
            model = test
            kept.remove(rid)
            decision = "remove_still_growth"
        else:
            decision = "keep_required_for_growth"
        trace.append(
            {
                "reaction_id": rid,
                "objective_status_without_reaction": test_status,
                "objective_value_without_reaction": test_value,
                "growth_call_without_reaction": test_call,
                "decision": decision,
            }
        )
    status, value, call = optimize(model, growth_eps)
    minimized = [rid for rid in rescue_ids if rid in kept]
    return minimized, trace, status, value, call


def build_sample_model(template, removed_ids: list[str], restored_ids: set[str]):
    model = template.copy()
    final_remove = [rid for rid in removed_ids if rid not in restored_ids and rid in model.reactions]
    model.remove_reactions([model.reactions.get_by_id(rid) for rid in final_remove], remove_orphans=True)
    return model


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--reaction-presence-absence", type=Path, required=True)
    parser.add_argument("--reaction-catalog", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--samples", default="all")
    parser.add_argument("--removable-categories", default="shell_gpr;absent_gpr_projection")
    parser.add_argument("--rescue-categories", default="shell_gpr")
    parser.add_argument("--max-addback", type=int, default=0, help="0 means all candidate reactions may be tested.")
    parser.add_argument("--growth-eps", type=float, default=1e-6)
    parser.add_argument("--atpm-id", default="calb_ATPM__c")
    parser.add_argument("--model-label", default="Calb")
    args = parser.parse_args()

    model_out = args.outdir / "models"
    qc_out = args.outdir / "qc"
    model_out.mkdir(parents=True, exist_ok=True)
    qc_out.mkdir(parents=True, exist_ok=True)

    template = read_model(args.model)
    reaction_pa = read_tsv(args.reaction_presence_absence)
    catalog = read_tsv(args.reaction_catalog)
    if "sample" not in reaction_pa.columns:
        raise ValueError("Reaction presence/absence matrix must contain a sample column.")
    if "reaction_id" not in catalog.columns:
        raise ValueError("Reaction catalog must contain reaction_id column.")

    available_samples = reaction_pa["sample"].map(clean).tolist()
    selected_samples = parse_samples(args.samples, available_samples)
    reaction_pa = reaction_pa.loc[reaction_pa["sample"].isin(selected_samples)].copy()

    catalog = catalog.copy()
    catalog["reaction_id"] = catalog["reaction_id"].map(clean)
    catalog_map = catalog.set_index("reaction_id").to_dict(orient="index")
    removable_categories = parse_set(args.removable_categories)
    rescue_categories = parse_set(args.rescue_categories)

    template_reactions = {rxn.id for rxn in template.reactions}
    objective_reactions = {rxn.id for rxn in template.reactions if rxn.objective_coefficient != 0}
    reaction_cols = [col for col in reaction_pa.columns if col != "sample"]
    usable_reaction_cols = [rid for rid in reaction_cols if rid in template_reactions]
    missing_pa_reactions = sorted(template_reactions - set(reaction_cols))
    missing_catalog_reactions = sorted(template_reactions - set(catalog_map))
    if missing_pa_reactions:
        raise ValueError(f"Template reactions missing from reaction PA matrix: {missing_pa_reactions[:20]}")
    if missing_catalog_reactions:
        raise ValueError(f"Template reactions missing from reaction catalog: {missing_catalog_reactions[:20]}")
    if args.atpm_id not in template.reactions:
        raise ValueError(f"ATPM reaction is missing from the template model: {args.atpm_id}")
    if not objective_reactions:
        raise ValueError("Template model has no objective reaction")

    qc_rows = []
    removed_rows = []
    rescue_rows = []
    addback_trace_rows = []
    minimize_trace_rows = []
    reaction_presence_rows = []

    for _, sample_row in reaction_pa.iterrows():
        sample = clean(sample_row["sample"])
        absent_removable = []
        for rid in usable_reaction_cols:
            if rid in objective_reactions or not is_absent(sample_row[rid]):
                continue
            category = clean(catalog_map.get(rid, {}).get("pan_reaction_category", ""))
            if category in removable_categories:
                absent_removable.append(rid)

        pruned = template.copy()
        pruned.remove_reactions(
            [pruned.reactions.get_by_id(rid) for rid in absent_removable if rid in pruned.reactions],
            remove_orphans=True,
        )
        initial_status, initial_value, initial_call = optimize(pruned, args.growth_eps)

        rescue_candidates = sorted_rescue_candidates(absent_removable, catalog_map, rescue_categories)
        prefix_rescue_ids, addback_trace, prefix_status, prefix_value, prefix_call = cumulative_addback_until_growth(
            pruned, template, rescue_candidates, args.growth_eps, args.max_addback
        )
        for row in addback_trace:
            addback_trace_rows.append({"sample": sample, **row, **catalog_map.get(clean(row["reaction_id"]), {})})

        minimized_ids = []
        minimize_trace = []
        minimized_status = prefix_status
        minimized_value = prefix_value
        minimized_call = prefix_call
        if prefix_call == "growth":
            minimized_ids, minimize_trace, minimized_status, minimized_value, minimized_call = greedy_minimize_rescue_set(
                pruned, template, prefix_rescue_ids, args.growth_eps
            )
        for row in minimize_trace:
            minimize_trace_rows.append({"sample": sample, **row, **catalog_map.get(clean(row["reaction_id"]), {})})

        restored_set = set(minimized_ids)
        final_model = build_sample_model(template, absent_removable, restored_set)
        final_model.id = f"{args.model_label}_{sample}_draft_ssGEM"
        final_model.name = f"Candida albicans {sample} draft ssGEM"
        final_status, final_value, final_call = optimize(final_model, args.growth_eps)
        egc_status, max_atpm, egc_call = max_closed_uptake_atpm(final_model, args.atpm_id)

        out_xml = model_out / f"{args.model_label}_{sample}_ssGEM.xml"
        write_sbml_model(final_model, str(out_xml))
        reread = read_model(out_xml)
        reread_status, reread_value, reread_call = optimize(reread, args.growth_eps)

        final_present = {rxn.id for rxn in reread.reactions}
        reaction_presence = {"sample": sample}
        for rid in usable_reaction_cols:
            reaction_presence[rid] = 1 if rid in final_present else 0
        reaction_presence_rows.append(reaction_presence)

        for rid in absent_removable:
            category = clean(catalog_map.get(rid, {}).get("pan_reaction_category", ""))
            removed_rows.append(
                {
                    "sample": sample,
                    "reaction_id": rid,
                    "initially_removed": True,
                    "restored_as_growth_rescue": rid in restored_set,
                    "final_model_presence": rid in final_present,
                    "rescue_interpretation": "growth_rescue_gapfill_not_genome_presence" if rid in restored_set else "genome_absent_removed",
                    **catalog_map.get(rid, {}),
                }
            )

        for rank, rid in enumerate(minimized_ids, start=1):
            rescue_rows.append(
                {
                    "sample": sample,
                    "rescue_rank": rank,
                    "reaction_id": rid,
                    "rescue_interpretation": "growth_rescue_gapfill_not_genome_presence",
                    **catalog_map.get(rid, {}),
                }
            )

        qc_rows.append(
            {
                "sample": sample,
                "model_file": str(out_xml),
                "base_reactions": len(template.reactions),
                "initial_removed_reactions": len(absent_removable),
                "initial_removed_shell_gpr": sum(clean(catalog_map.get(rid, {}).get("pan_reaction_category", "")) == "shell_gpr" for rid in absent_removable),
                "initial_removed_absent_gpr_projection": sum(clean(catalog_map.get(rid, {}).get("pan_reaction_category", "")) == "absent_gpr_projection" for rid in absent_removable),
                "initial_objective_status": initial_status,
                "initial_objective_value": initial_value,
                "initial_growth_call": initial_call,
                "prefix_rescue_reactions": len(prefix_rescue_ids),
                "minimized_rescue_reactions": len(minimized_ids),
                "final_removed_reactions": len(absent_removable) - len(restored_set),
                "final_reactions": len(reread.reactions),
                "final_metabolites": len(reread.metabolites),
                "final_genes": len(reread.genes),
                "objective_status": reread_status,
                "objective_value": reread_value,
                "growth_call": reread_call,
                "closed_uptake_atpm_status": egc_status,
                "closed_uptake_max_atpm": max_atpm,
                "egc_call": egc_call,
                "read_write_status": "pass",
                "builder_policy": "remove_absent_shell_and_absent_projection_then_minimal_shell_growth_rescue",
            }
        )

    qc = pd.DataFrame(qc_rows)
    removed = pd.DataFrame(removed_rows)
    rescued = pd.DataFrame(rescue_rows)
    addback_trace_df = pd.DataFrame(addback_trace_rows)
    minimize_trace_df = pd.DataFrame(minimize_trace_rows)
    reaction_presence_df = pd.DataFrame(reaction_presence_rows)

    objective_numeric = pd.to_numeric(qc["objective_value"], errors="coerce").fillna(0.0)
    summary = pd.DataFrame(
        [
            {"metric": "input_model", "value": str(args.model)},
            {"metric": "input_model_sha256", "value": sha256(args.model)},
            {"metric": "reaction_presence_absence", "value": str(args.reaction_presence_absence)},
            {"metric": "reaction_presence_absence_sha256", "value": sha256(args.reaction_presence_absence)},
            {"metric": "reaction_catalog", "value": str(args.reaction_catalog)},
            {"metric": "reaction_catalog_sha256", "value": sha256(args.reaction_catalog)},
            {"metric": "atpm_reaction", "value": args.atpm_id},
            {"metric": "samples_requested", "value": args.samples},
            {"metric": "samples_written", "value": len(qc)},
            {"metric": "initial_growth_samples_after_strict_pruning", "value": int((qc["initial_growth_call"] == "growth").sum())},
            {"metric": "final_growth_samples_after_rescue", "value": int((qc["growth_call"] == "growth").sum())},
            {"metric": "final_zero_or_failed_growth_samples", "value": int((qc["growth_call"] != "growth").sum())},
            {"metric": "egc_pass_samples", "value": int((qc["egc_call"] == "pass_no_closed_uptake_atpm_cycle").sum())},
            {"metric": "median_initial_removed_reactions", "value": float(qc["initial_removed_reactions"].median()) if len(qc) else ""},
            {"metric": "median_minimized_rescue_reactions", "value": float(qc["minimized_rescue_reactions"].median()) if len(qc) else ""},
            {"metric": "median_final_removed_reactions", "value": float(qc["final_removed_reactions"].median()) if len(qc) else ""},
            {"metric": "median_final_reactions", "value": float(qc["final_reactions"].median()) if len(qc) else ""},
            {"metric": "median_objective_value", "value": float(objective_numeric.median()) if len(qc) else ""},
            {"metric": "removable_categories", "value": ";".join(sorted(removable_categories))},
            {"metric": "rescue_categories", "value": ";".join(sorted(rescue_categories))},
            {"metric": "claim_boundary", "value": "final_ssGEMs_rescue_reactions_not_genome_presence"},
        ]
    )

    if not removed.empty:
        removed_summary = (
            removed.groupby(["pan_reaction_category", "restored_as_growth_rescue"], dropna=False)
            .size()
            .reset_index(name="sample_reaction_events")
            .sort_values("sample_reaction_events", ascending=False)
        )
    else:
        removed_summary = pd.DataFrame(columns=["pan_reaction_category", "restored_as_growth_rescue", "sample_reaction_events"])

    qc.to_csv(qc_out / "v18_afumigatus_style_ssgem_qc.tsv", sep="\t", index=False)
    removed.to_csv(qc_out / "v18_afumigatus_style_removed_and_rescued_reactions.tsv", sep="\t", index=False)
    rescued.to_csv(qc_out / "v18_afumigatus_style_growth_rescue_reactions.tsv", sep="\t", index=False)
    addback_trace_df.to_csv(qc_out / "v18_afumigatus_style_cumulative_addback_trace.tsv", sep="\t", index=False)
    minimize_trace_df.to_csv(qc_out / "v18_afumigatus_style_rescue_minimization_trace.tsv", sep="\t", index=False)
    removed_summary.to_csv(qc_out / "v18_afumigatus_style_removed_rescue_summary.tsv", sep="\t", index=False)
    reaction_presence_df.to_csv(qc_out / "v18_afumigatus_style_final_reaction_presence_absence.tsv", sep="\t", index=False)
    summary.to_csv(qc_out / "v18_afumigatus_style_generation_summary.tsv", sep="\t", index=False)

    print("[DONE] final ssGEM generation")
    print(f"  models:   {model_out}")
    print(f"  qc:       {qc_out / 'v18_afumigatus_style_ssgem_qc.tsv'}")
    print(f"  rescued:  {qc_out / 'v18_afumigatus_style_growth_rescue_reactions.tsv'}")
    print(f"  removed:  {qc_out / 'v18_afumigatus_style_removed_and_rescued_reactions.tsv'}")
    print(f"  matrix:   {qc_out / 'v18_afumigatus_style_final_reaction_presence_absence.tsv'}")
    print(f"  summary:  {qc_out / 'v18_afumigatus_style_generation_summary.tsv'}")
    print()
    print(summary.to_string(index=False))
    print()
    print(removed_summary.to_string(index=False))


if __name__ == "__main__":
    main()
