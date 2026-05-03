"""49Q IBM Kingston runner matching the small-budget ~0.209 setup.

Default behavior is safe: it writes a plan JSON and does not submit hardware
jobs. Add ``--execute`` to submit IBM Runtime Sampler jobs.

This matches the original small-budget run:

- circuit: operator_loschmidt_echo_49x648, L=3, delta=0.15
- backend: ibm_kingston
- sample: N_int=8 predeclared uniform bitstrings, seed 20260430
- shots: 21 * 251 = 5271 per delta/delta0 circuit
- mitigation: global rescaling, DD XY4, active Pauli twirling
- no readout mitigation, no ZNE, no postselection

Important distinction: the reported N_int=8 value 0.2093896573 was produced
without readout mitigation or ZNE. Readout and ZNE were tested later as a small
two-string A/B diagnostic and did not produce a replacement N_int=8 aggregate.
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from qiskit import QuantumCircuit, qasm3, transpile


ROOT = Path(__file__).resolve().parents[1]
QASM_PATH = ROOT / "data" / "49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm"
SUMMARY_PATH = ROOT / "results" / "49q_small_budget_kingston_summary.json"
PLAN_OUTPUT_PATH = ROOT / "results" / "49q_kingston_small_budget_plan.json"
RUN_OUTPUT_PATH = ROOT / "results" / "49q_kingston_small_budget_run.json"

BACKEND = "ibm_kingston"
DELTA = 0.15
SHOTS = 21 * 251
TWIRLS = 21
SHOTS_PER_TWIRL = 251
SEED_TRANSPILER = 424242
OPTIMIZATION_LEVEL = 1

# Fixed calibration-selected compressed-49Q -> Kingston physical layout used by
# the original small-budget 0.209 run.
INITIAL_LAYOUT = (
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    16,
    17,
    18,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    36,
    37,
    38,
    41,
    42,
    43,
    44,
    45,
    46,
    47,
    48,
    49,
    50,
    51,
    56,
    57,
    58,
    63,
    64,
    65,
    66,
    67,
    68,
    69,
    70,
    71,
)

# For the compressed active-register circuit.
OBSERVABLE_ACTIVE = (33, 39, 45)  # tracker observable Z52 Z59 Z72
PERTURBATION_SUPPORT_ACTIVE = (
    4,
    5,
    6,
    7,
    8,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    26,
    27,
    28,
    29,
    30,
)


def load_full_qasm_circuit() -> QuantumCircuit:
    return qasm3.loads(QASM_PATH.read_text(encoding="utf-8"))


def active_qubit_indices(circuit: QuantumCircuit) -> list[int]:
    return sorted(
        {
            circuit.find_bit(qubit).index
            for instruction in circuit.data
            if instruction.operation.name != "barrier"
            for qubit in instruction.qubits
        }
    )


def compress_to_active_register(circuit: QuantumCircuit) -> tuple[QuantumCircuit, list[int]]:
    """Compress the public physical-register QASM to its 49 active qubits."""

    active = active_qubit_indices(circuit)
    mapping = {physical: logical for logical, physical in enumerate(active)}
    compressed = QuantumCircuit(len(active), name="operator_loschmidt_echo_49x648_active49")
    compressed.global_phase = circuit.global_phase
    for instruction in circuit.data:
        if instruction.operation.name == "barrier":
            continue
        qargs = [mapping[circuit.find_bit(qubit).index] for qubit in instruction.qubits]
        compressed.append(instruction.operation, qargs)
    return compressed, active


def zero_delta_rotations(circuit: QuantumCircuit) -> tuple[QuantumCircuit, int]:
    """Remove rz(2*delta) rotations on the perturbation support."""

    support = set(PERTURBATION_SUPPORT_ACTIVE)
    target_angle = 2.0 * DELTA
    zeroed = QuantumCircuit(circuit.num_qubits, name=f"{circuit.name}_delta0")
    zeroed.global_phase = circuit.global_phase
    removed = 0
    for instruction in circuit.data:
        operation = instruction.operation
        qargs = [circuit.find_bit(qubit).index for qubit in instruction.qubits]
        if (
            operation.name == "rz"
            and len(qargs) == 1
            and qargs[0] in support
            and len(operation.params) == 1
            and abs(float(operation.params[0]) - target_angle) <= 1e-9
        ):
            removed += 1
            continue
        zeroed.append(operation, qargs)
    return zeroed, removed


def build_measurement_circuit(circuit: QuantumCircuit, *, bitstring: str) -> QuantumCircuit:
    measured = QuantumCircuit(circuit.num_qubits, len(OBSERVABLE_ACTIVE))
    for qubit, bit in enumerate(bitstring):
        if bit == "1":
            measured.x(qubit)
    measured.compose(circuit, inplace=True)
    for clbit, qubit in enumerate(OBSERVABLE_ACTIVE):
        measured.measure(qubit, clbit)
    return measured


def input_parity(bitstring: str) -> int:
    ones = sum(1 for qubit in OBSERVABLE_ACTIVE if bitstring[qubit] == "1")
    return -1 if ones % 2 else 1


def counts_to_z_product_expectation(counts: Mapping[Any, Any]) -> float:
    total = sum(int(count) for count in counts.values())
    if total == 0:
        raise ValueError("Cannot compute expectation from empty counts.")
    weighted = 0
    for raw_key, raw_count in counts.items():
        key = str(raw_key).replace(" ", "").zfill(len(OBSERVABLE_ACTIVE))
        parity = -1 if key.count("1") % 2 else 1
        weighted += parity * int(raw_count)
    return weighted / total


def normalize_counts(counts: Mapping[Any, Any]) -> dict[str, int]:
    output: dict[str, int] = {}
    for key, value in counts.items():
        if isinstance(key, int):
            bitstring = format(key, f"0{len(OBSERVABLE_ACTIVE)}b")
        else:
            bitstring = str(key).replace(" ", "").zfill(len(OBSERVABLE_ACTIVE))
        output[bitstring] = int(round(float(value)))
    return output


def extract_counts_list_from_sampler_result(result: Any, *, n_items: int) -> list[dict[str, int]]:
    if hasattr(result, "quasi_dists"):
        quasi_dists = getattr(result, "quasi_dists")
        return [
            normalize_counts({key: float(probability) * SHOTS for key, probability in quasi_dist.items()})
            for quasi_dist in quasi_dists[:n_items]
        ]

    output: list[dict[str, int]] = []
    for item_index in range(n_items):
        item = result[item_index] if hasattr(result, "__getitem__") else None
        if item is None:
            raise RuntimeError(f"Could not index sampler result at item {item_index}.")
        counts = None
        data = getattr(item, "data", None)
        if data is not None:
            for register_name in ("c", "meas"):
                register = getattr(data, register_name, None)
                if register is not None and hasattr(register, "get_counts"):
                    counts = register.get_counts()
                    break
            if counts is None and hasattr(data, "get_counts"):
                counts = data.get_counts()
        if counts is None and hasattr(item, "get_counts"):
            counts = item.get_counts()
        if counts is None:
            raise RuntimeError(f"Could not extract counts for sampler item {item_index}.")
        output.append(normalize_counts(counts))
    return output


def sampler_options_payload() -> dict[str, Any]:
    return {
        "default_shots": SHOTS,
        "dynamical_decoupling": {
            "enable": True,
            "sequence_type": "XY4",
            "scheduling_method": "alap",
            "extra_slack_distribution": "middle",
            "skip_reset_qubits": True,
        },
        "twirling": {
            "enable_gates": True,
            "enable_measure": True,
            "strategy": "active",
            "num_randomizations": TWIRLS,
            "shots_per_randomization": SHOTS_PER_TWIRL,
        },
    }


def build_plan(*, sample_count: int) -> dict[str, Any]:
    full = load_full_qasm_circuit()
    active_circuit, active_physical = compress_to_active_register(full)
    delta0_circuit, removed = zero_delta_rotations(active_circuit)
    summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
    bitstrings = [record["bitstring"] for record in summary["records"][:sample_count]]

    return {
        "captured_at_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "plan",
        "execute_to_submit_hardware": False,
        "instance_id": "operator_loschmidt_echo_49x648",
        "qasm_path": str(QASM_PATH.relative_to(ROOT)),
        "backend": BACKEND,
        "shots": SHOTS,
        "twirling_num_randomizations": TWIRLS,
        "twirling_shots_per_randomization": SHOTS_PER_TWIRL,
        "optimization_level": OPTIMIZATION_LEVEL,
        "seed_transpiler": SEED_TRANSPILER,
        "initial_layout": list(INITIAL_LAYOUT),
        "sampler_options": sampler_options_payload(),
        "mitigation": {
            "global_rescaling": True,
            "dynamical_decoupling": "XY4",
            "pauli_twirling": "active",
            "readout_mitigation": False,
            "zne": False,
            "postselection": False,
        },
        "circuit": {
            "declared_qubits": full.num_qubits,
            "active_qubits": active_circuit.num_qubits,
            "active_physical_qubits": active_physical,
            "delta_depth": active_circuit.depth(),
            "delta0_depth": delta0_circuit.depth(),
            "removed_delta_rotations": removed,
            "observable_active": list(OBSERVABLE_ACTIVE),
            "perturbation_support_active": list(PERTURBATION_SUPPORT_ACTIVE),
        },
        "sampling": {
            "source": str(SUMMARY_PATH.relative_to(ROOT)),
            "strategy": "predeclared_uniform_random_basis_states_no_forced_zero",
            "random_seed": summary["sampling"]["random_seed"],
            "sample_count": sample_count,
            "bitstrings": bitstrings,
        },
        "baseline_reference": {
            "mean_hardware_global_rescaled_ratio": summary["aggregate"]["mean_hardware_global_rescaled_ratio"],
            "standard_error": summary["aggregate"]["hardware_standard_error"],
            "normal_95_interval": [
                summary["aggregate"]["hardware_lower95_normal_approx"],
                summary["aggregate"]["hardware_upper95_normal_approx"],
            ],
            "note": "This runner is configured to reproduce this small-budget setup, not the tracker N_int=500 protocol.",
            "readout_zne_note": "The 0.2093896573 aggregate did not include readout mitigation or ZNE; those were later two-string A/B diagnostics only.",
        },
    }


def build_runtime_service() -> Any:
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
    except ImportError as exc:
        raise RuntimeError("qiskit-ibm-runtime is required for --execute.") from exc

    token = os.getenv("QCAPI_TOKEN") or os.getenv("QISKIT_IBM_TOKEN") or os.getenv("IBM_QUANTUM_TOKEN")
    instance = os.getenv("QISKIT_IBM_INSTANCE")
    if token:
        kwargs: dict[str, str] = {"token": token}
        if instance:
            kwargs["instance"] = instance
        try:
            return QiskitRuntimeService(channel="ibm_quantum", **kwargs)
        except TypeError:
            return QiskitRuntimeService(**kwargs)
    return QiskitRuntimeService()


def run_sampler_job(backend: Any, circuits: list[QuantumCircuit], *, submit_only: bool) -> tuple[str, Any | None]:
    from qiskit_ibm_runtime import SamplerV2

    sampler = SamplerV2(mode=backend, options=sampler_options_payload())
    job = sampler.run(circuits, shots=SHOTS)
    if submit_only:
        return str(job.job_id()), None
    return str(job.job_id()), job.result()


def run_hardware(*, sample_count: int, submit_only: bool) -> dict[str, Any]:
    plan = build_plan(sample_count=sample_count)
    full = load_full_qasm_circuit()
    active_circuit, _ = compress_to_active_register(full)
    delta0_circuit, _ = zero_delta_rotations(active_circuit)

    service = build_runtime_service()
    backend = service.backend(BACKEND)
    records: list[dict[str, Any]] = []

    for index, bitstring in enumerate(plan["sampling"]["bitstrings"]):
        delta_meas = build_measurement_circuit(active_circuit, bitstring=bitstring)
        delta0_meas = build_measurement_circuit(delta0_circuit, bitstring=bitstring)
        transpiled = transpile(
            [delta_meas, delta0_meas],
            backend=backend,
            optimization_level=OPTIMIZATION_LEVEL,
            seed_transpiler=SEED_TRANSPILER,
            initial_layout=list(INITIAL_LAYOUT),
        )
        circuits = list(transpiled) if isinstance(transpiled, list) else [transpiled]
        job_id, result = run_sampler_job(backend, circuits, submit_only=submit_only)
        record: dict[str, Any] = {
            "uniform_index": index,
            "bitstring": bitstring,
            "hardware_job_id": job_id,
            "submit_only": submit_only,
            "transpiled": [
                {"depth": int(circuit.depth()), "size": int(circuit.size())}
                for circuit in circuits
            ],
        }
        if result is not None:
            counts_list = extract_counts_list_from_sampler_result(result, n_items=2)
            sigma = input_parity(bitstring)
            delta_expectation = counts_to_z_product_expectation(counts_list[0])
            delta0_expectation = counts_to_z_product_expectation(counts_list[1])
            delta_weighted = sigma * delta_expectation
            delta0_weighted = sigma * delta0_expectation
            record.update(
                {
                    "delta_counts": counts_list[0],
                    "delta0_counts": counts_list[1],
                    "input_parity": sigma,
                    "delta_weighted_term": delta_weighted,
                    "delta0_weighted_term": delta0_weighted,
                    "global_rescaled_ratio": delta_weighted / delta0_weighted,
                }
            )
        records.append(record)

    ratios = [record["global_rescaled_ratio"] for record in records if "global_rescaled_ratio" in record]
    aggregate = None
    if ratios:
        aggregate = {
            "completed_terms": len(ratios),
            "mean_global_rescaled_ratio": statistics.fmean(ratios),
            "standard_deviation": statistics.stdev(ratios) if len(ratios) > 1 else 0.0,
            "standard_error": (statistics.stdev(ratios) / (len(ratios) ** 0.5)) if len(ratios) > 1 else 0.0,
        }

    payload = dict(plan)
    payload.update(
        {
            "captured_at_utc": datetime.now(timezone.utc).isoformat(),
            "mode": "hardware-run",
            "execute_to_submit_hardware": True,
            "submit_only": submit_only,
            "records": records,
            "aggregate": aggregate,
        }
    )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-count", type=int, default=8, choices=range(1, 9))
    parser.add_argument("--execute", action="store_true", help="Submit IBM Runtime jobs.")
    parser.add_argument("--submit-only", action="store_true", help="Submit jobs and do not wait for results.")
    parser.add_argument("--output-json", default=None)
    args = parser.parse_args()

    output_path = Path(args.output_json) if args.output_json else (RUN_OUTPUT_PATH if args.execute else PLAN_OUTPUT_PATH)
    payload = run_hardware(sample_count=args.sample_count, submit_only=args.submit_only) if args.execute else build_plan(sample_count=args.sample_count)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    print(f"mode={payload['mode']}")
    print(f"backend={payload['backend']}")
    print(f"sample_count={payload['sampling']['sample_count']}")
    print(f"shots_per_circuit={payload['shots']}")
    print(f"mitigation={payload['mitigation']}")
    print(f"removed_delta_rotations={payload['circuit']['removed_delta_rotations']}")
    if payload.get("aggregate") is not None:
        print(f"mean_global_rescaled_ratio={payload['aggregate']['mean_global_rescaled_ratio']}")
    print(f"output_json={output_path}")
    if not args.execute:
        print("dry_run=true; add --execute to submit IBM Runtime jobs.")


if __name__ == "__main__":
    main()
