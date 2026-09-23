#!/usr/bin/env python3
"""Collect reviewer-simulation outputs into compact CSV summaries."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_key_value(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return {row["setting"]: row["value"] for row in csv.DictReader(handle)}


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def collect_rrt(root: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    scenario_rows: list[dict[str, object]] = []
    parameter_rows: list[dict[str, object]] = []
    for settings_path in sorted(root.rglob("simulation_settings.csv")):
        run_dir = settings_path.parent
        settings = read_key_value(settings_path)
        identity = {
            "run_id": settings.get("run_id", run_dir.name),
            "article_case": settings.get("article_case", ""),
            "case_model": settings.get("case_model", ""),
            "prevalence": settings.get("prevalence", ""),
            "n": settings.get("n", ""),
            "replications": settings.get("fixed_replications_no_replacement", ""),
            "p_true": settings.get("p_true", settings.get("p", "")),
            "c_true": settings.get("c_true", settings.get("c", "")),
            "p_fitted": settings.get("p_fitted", settings.get("p", "")),
            "c_fitted": settings.get("c_fitted", settings.get("c", "")),
            "n_starts": settings.get("number_of_EM_starts", ""),
            "result_directory": str(run_dir.resolve()),
        }
        row: dict[str, object] = dict(identity)
        for metric in read_rows(run_dir / "convergence_summary.csv"):
            name = metric.get("metric", "")
            row[f"conv_{name}_count"] = metric.get("count", "")
            row[f"conv_{name}_denominator"] = metric.get("denominator", "")
            row[f"conv_{name}_proportion"] = metric.get("proportion", "")
        for metric in read_rows(run_dir / "test_performance_summary.csv"):
            name = metric.get("metric", "")
            for field in ("value", "n_valid", "n_total", "MCSE", "undefined_rate"):
                row[f"test_{name}_{field}"] = metric.get(field, "")
        scenario_rows.append(row)

        source = run_dir / "parameter_summary_full_comparison.csv"
        if not source.exists():
            source = run_dir / "parameter_summary_full.csv"
        for parameter in read_rows(source):
            parameter_rows.append({**identity, **parameter})
    return scenario_rows, parameter_rows


def collect_npo(root: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    scenario_rows: list[dict[str, object]] = []
    parameter_rows: list[dict[str, object]] = []
    for settings_path in sorted(root.rglob("settings.csv")):
        run_dir = settings_path.parent
        settings = read_key_value(settings_path)
        identity = {
            "run_id": run_dir.name,
            "p": settings.get("p", ""),
            "c": settings.get("c", ""),
            "replications": settings.get("replications", ""),
            "n_starts": settings.get("n_starts", ""),
            "result_directory": str(run_dir.resolve()),
        }
        for scenario in read_rows(run_dir / "scenario_summary.csv"):
            scenario_rows.append({**identity, **scenario})
        for parameter in read_rows(run_dir / "parameter_summary.csv"):
            parameter_rows.append({**identity, **parameter})
    return scenario_rows, parameter_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("results_root", type=Path)
    parser.add_argument("summary_dir", type=Path)
    args = parser.parse_args()

    rrt_scenarios, rrt_parameters = collect_rrt(args.results_root / "rrt_sensitivity")
    npo_scenarios, npo_parameters = collect_npo(args.results_root / "npo_formal")
    write_rows(args.summary_dir / "rrt_scenario_summary.csv", rrt_scenarios)
    write_rows(args.summary_dir / "rrt_parameter_summary.csv", rrt_parameters)
    write_rows(args.summary_dir / "npo_scenario_summary.csv", npo_scenarios)
    write_rows(args.summary_dir / "npo_parameter_summary.csv", npo_parameters)
    print(
        f"Wrote {len(rrt_scenarios)} RRT scenario rows, "
        f"{len(rrt_parameters)} RRT parameter rows, "
        f"{len(npo_scenarios)} NPO scenario rows, and "
        f"{len(npo_parameters)} NPO parameter rows."
    )


if __name__ == "__main__":
    main()
