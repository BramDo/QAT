"""Aggregate the available 49Q Kingston ZNE/readout spotchecks.

This script intentionally labels the output as a diagnostic subset, not a final
tracker result. The current rows combine the earlier A/B mitigation records for
uniform indices 0 and 5 with the later index-1 ZNE+readout spotcheck.
"""

from __future__ import annotations

import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
AB_PATH = ROOT / "results" / "49q_small_budget_ab_mitigation_summary.json"
INDEX1_PATH = ROOT / "results" / "49q_zne_readout_spotcheck_index1_20260503.json"
OUTPUT_PATH = ROOT / "results" / "49q_zne_readout_subset3_summary.json"


def summarize(values: list[float]) -> dict[str, Any]:
    mean = statistics.fmean(values)
    standard_deviation = statistics.stdev(values) if len(values) > 1 else 0.0
    standard_error = standard_deviation / math.sqrt(len(values)) if len(values) > 1 else 0.0
    return {
        "n": len(values),
        "mean": mean,
        "standard_deviation": standard_deviation,
        "standard_error": standard_error,
        "normal_95_interval": [
            mean - 1.96 * standard_error,
            mean + 1.96 * standard_error,
        ],
    }


def load_rows() -> list[dict[str, Any]]:
    ab = json.loads(AB_PATH.read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    for row in ab["rows"]:
        zne = row["zne_delta_spotcheck_1024_shots"]
        rows.append(
            {
                "uniform_index": row["uniform_index"],
                "bitstring": row["bitstring"],
                "baseline_ratio": row["baseline"]["global_rescaled_ratio"],
                "raw_zne_ratio": zne["raw_zne_over_existing_delta0_ratio"],
                "full_assignment_readout_zne_ratio": zne["full_readout_zne_over_existing_delta0_ratio"],
                "local_tensor_readout_zne_ratio": zne["tensor_readout_zne_over_existing_delta0_ratio"],
                "zne_job_id": zne["hardware_job_id"],
                "readout_calibration_job_id": zne["calibration_job_id"],
                "source": str(AB_PATH.relative_to(ROOT)),
            }
        )

    index1 = json.loads(INDEX1_PATH.read_text(encoding="utf-8"))
    rows.append(
        {
            "uniform_index": index1["uniform_index"],
            "bitstring": index1["bitstring"],
            "baseline_ratio": index1["baseline"]["global_rescaled_ratio"],
            "raw_zne_ratio": index1["zne_result"]["zne_over_existing_delta0_ratio"],
            "full_assignment_readout_zne_ratio": index1["readout_mitigation"]["zne_result"]["full_assignment"][
                "zne_over_existing_delta0_ratio"
            ],
            "local_tensor_readout_zne_ratio": index1["readout_mitigation"]["zne_result"]["local_tensor"][
                "zne_over_existing_delta0_ratio"
            ],
            "zne_job_id": index1["hardware_job_id"],
            "readout_calibration_job_id": index1["readout_mitigation"]["calibration_job_id"],
            "source": str(INDEX1_PATH.relative_to(ROOT)),
        }
    )
    return sorted(rows, key=lambda row: int(row["uniform_index"]))


def main() -> None:
    rows = load_rows()
    payload = {
        "captured_at_utc": datetime.now(timezone.utc).isoformat(),
        "instance_id": "operator_loschmidt_echo_49x648",
        "backend": "ibm_kingston",
        "scope": {
            "uniform_indices": [row["uniform_index"] for row in rows],
            "claim_boundary": (
                "Diagnostic subset over available mitigated strings. It is not a final N_int=8 "
                "mitigated aggregate because the subset was assembled after inspecting earlier A/B results."
            ),
            "completion_path": (
                "For a clean complete small-budget result, run the same ZNE+readout pipeline on "
                "the remaining original uniform indices 2,3,4,6,7 and average all eight rows."
            ),
        },
        "mitigation": {
            "zne": "local-fold delta numerator, factors 1,3,5, 1024 shots per folded circuit",
            "readout": "3-bit observable readout calibration, full assignment and local tensor correction",
            "global_rescaling": "ZNE/readout-corrected delta numerator divided by existing matching delta0 denominator",
        },
        "rows": rows,
        "aggregate": {
            "baseline_ratio": summarize([row["baseline_ratio"] for row in rows]),
            "raw_zne_ratio": summarize([row["raw_zne_ratio"] for row in rows]),
            "full_assignment_readout_zne_ratio": summarize(
                [row["full_assignment_readout_zne_ratio"] for row in rows]
            ),
            "local_tensor_readout_zne_ratio": summarize(
                [row["local_tensor_readout_zne_ratio"] for row in rows]
            ),
        },
    }
    OUTPUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(f"rows={len(rows)}")
    print(f"indices={payload['scope']['uniform_indices']}")
    print(f"baseline_mean={payload['aggregate']['baseline_ratio']['mean']}")
    print(f"tensor_readout_zne_mean={payload['aggregate']['local_tensor_readout_zne_ratio']['mean']}")
    print(f"tensor_readout_zne_se={payload['aggregate']['local_tensor_readout_zne_ratio']['standard_error']}")
    print(f"output_json={OUTPUT_PATH}")


if __name__ == "__main__":
    main()
