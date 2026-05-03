"""One-string 49Q ZNE spotcheck for the small-budget Kingston setup.

Default behavior is a dry-run plan. Add ``--execute`` to submit one IBM
Runtime job with folded delta circuits.

The ratio reported by this script is a diagnostic hybrid:

    ZNE(delta numerator) / existing baseline delta0 denominator

It is not a replacement for the full N_int=8 aggregate.
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from qiskit import QuantumCircuit, qasm3, transpile

from run_49q_kingston_small_budget import (
    BACKEND,
    DELTA,
    INITIAL_LAYOUT,
    OBSERVABLE_ACTIVE,
    OPTIMIZATION_LEVEL,
    QASM_PATH,
    ROOT,
    SEED_TRANSPILER,
    SUMMARY_PATH,
    active_qubit_indices,
    build_measurement_circuit,
    compress_to_active_register,
    counts_to_z_product_expectation,
    input_parity,
    load_full_qasm_circuit,
    normalize_counts,
    sampler_options_payload,
)


SPOTCHECK_OUTPUT_PATH = ROOT / "results" / "49q_zne_spotcheck_result.json"
PLAN_OUTPUT_PATH = ROOT / "results" / "49q_zne_spotcheck_plan.json"
ZNE_FACTORS = (1, 3, 5)
ZNE_SHOTS = 1024
READOUT_CAL_SHOTS = 2048


def fold_transpiled_two_qubit_local(circuit: QuantumCircuit, factor: int) -> QuantumCircuit:
    """Apply odd-factor local folding to transpiled two-qubit gates only."""

    if factor < 1 or factor % 2 == 0:
        raise ValueError(f"ZNE factor must be a positive odd integer, got {factor}.")
    if factor == 1:
        return circuit.copy()

    folded = QuantumCircuit(circuit.num_qubits, circuit.num_clbits, name=f"{circuit.name}_lfold{factor}")
    folded.global_phase = circuit.global_phase
    extra_pairs = (factor - 1) // 2
    for instruction in circuit.data:
        qargs = [circuit.find_bit(qubit).index for qubit in instruction.qubits]
        cargs = [circuit.find_bit(clbit).index for clbit in instruction.clbits]
        folded.append(instruction.operation, qargs, cargs)
        if len(instruction.qubits) == 2:
            inverse_operation = instruction.operation.inverse()
            for _ in range(extra_pairs):
                folded.append(inverse_operation, qargs, [])
                folded.append(instruction.operation, qargs, [])
    return folded


def extract_counts_list_from_sampler_result(result: Any, *, n_items: int) -> list[dict[str, int]]:
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


def extrapolate_zero_noise(*, factors: tuple[int, ...], values: list[float]) -> dict[str, Any]:
    coefficients = np.polyfit(np.asarray(factors, dtype=float), np.asarray(values, dtype=float), deg=1)
    prediction = np.polyval(coefficients, np.asarray(factors, dtype=float))
    return {
        "input_factors": [float(factor) for factor in factors],
        "input_values": values,
        "used_order": 1,
        "coefficients_high_to_low": [float(value) for value in coefficients.tolist()],
        "zero_noise_estimate": float(np.polyval(coefficients, 0.0)),
        "residuals": [float(value) for value in (np.asarray(values, dtype=float) - prediction).tolist()],
    }


def bitstring_to_index(bitstring: str, *, num_bits: int) -> int:
    normalized = str(bitstring).replace(" ", "").zfill(num_bits)
    if len(normalized) != num_bits or any(character not in "01" for character in normalized):
        raise ValueError(f"Expected a {num_bits}-bit binary string, got {bitstring!r}.")
    return int(normalized, 2)


def counts_probability_vector(counts: Mapping[Any, Any], *, num_bits: int) -> np.ndarray:
    normalized = normalize_counts(counts)
    total = float(sum(normalized.values()))
    if total <= 0.0:
        raise ValueError("Counts are empty; cannot form probability vector.")
    probabilities = np.zeros(2**num_bits, dtype=float)
    for bitstring, count in normalized.items():
        probabilities[bitstring_to_index(bitstring, num_bits=num_bits)] += float(count) / total
    return probabilities


def expectation_from_probability_vector(probabilities: np.ndarray, *, num_bits: int) -> float:
    expectation = 0.0
    for index, probability in enumerate(probabilities):
        bitstring = format(index, f"0{num_bits}b")
        parity = -1.0 if bitstring.count("1") % 2 else 1.0
        expectation += parity * float(probability)
    return float(expectation)


def vector_to_distribution(probabilities: np.ndarray, *, num_bits: int) -> dict[str, float]:
    return {
        format(index, f"0{num_bits}b"): float(probability)
        for index, probability in enumerate(probabilities)
    }


def tensor_assignment_matrix(local_assignments: Sequence[np.ndarray], *, num_bits: int) -> np.ndarray:
    if len(local_assignments) != num_bits:
        raise ValueError(f"Expected {num_bits} local assignment matrices.")
    tensor = np.zeros((2**num_bits, 2**num_bits), dtype=float)
    states = [format(index, f"0{num_bits}b") for index in range(2**num_bits)]
    for observed_index, observed_state in enumerate(states):
        for prepared_index, prepared_state in enumerate(states):
            probability = 1.0
            for display_pos, assignment in enumerate(local_assignments):
                observed_bit = int(observed_state[display_pos])
                prepared_bit = int(prepared_state[display_pos])
                probability *= float(assignment[observed_bit, prepared_bit])
            tensor[observed_index, prepared_index] = probability
    return tensor


def mitigate_counts_with_assignment(
    counts: Mapping[Any, Any],
    assignment_matrix: np.ndarray,
    *,
    num_bits: int,
) -> dict[str, Any]:
    observed = counts_probability_vector(counts, num_bits=num_bits)
    quasi = np.linalg.pinv(assignment_matrix) @ observed
    clipped = np.clip(quasi, 0.0, None)
    clipped_total = float(np.sum(clipped))
    if clipped_total > 0.0:
        clipped = clipped / clipped_total
    return {
        "expectation": expectation_from_probability_vector(quasi, num_bits=num_bits),
        "clipped_expectation": expectation_from_probability_vector(clipped, num_bits=num_bits),
        "quasi_distribution": vector_to_distribution(quasi, num_bits=num_bits),
        "clipped_distribution": vector_to_distribution(clipped, num_bits=num_bits),
    }


def build_readout_assignment_payload(
    calibration_counts: Sequence[Mapping[Any, Any]],
    *,
    prepared_states: Sequence[str],
    num_bits: int,
) -> dict[str, Any]:
    if len(calibration_counts) != len(prepared_states):
        raise ValueError("calibration_counts and prepared_states must have equal length.")
    if len(prepared_states) != 2**num_bits:
        raise ValueError(f"Expected {2**num_bits} prepared calibration states for {num_bits} bits.")

    full_assignment = np.zeros((2**num_bits, 2**num_bits), dtype=float)
    local_numerators = [np.zeros((2, 2), dtype=float) for _ in range(num_bits)]
    local_denominators = [np.zeros(2, dtype=float) for _ in range(num_bits)]

    for prepared_state, raw_counts in zip(prepared_states, calibration_counts, strict=True):
        prepared_index = bitstring_to_index(prepared_state, num_bits=num_bits)
        normalized_counts = normalize_counts(raw_counts)
        total = float(sum(normalized_counts.values()))
        if total <= 0.0:
            raise ValueError(f"Calibration state {prepared_state} has empty counts.")
        for observed_state, count in normalized_counts.items():
            observed_index = bitstring_to_index(observed_state, num_bits=num_bits)
            probability = float(count) / total
            full_assignment[observed_index, prepared_index] += probability

            for display_pos, (prepared_bit_text, observed_bit_text) in enumerate(
                zip(prepared_state, observed_state, strict=True)
            ):
                prepared_bit = int(prepared_bit_text)
                observed_bit = int(observed_bit_text)
                local_numerators[display_pos][observed_bit, prepared_bit] += float(count)
                local_denominators[display_pos][prepared_bit] += float(count)

    local_assignments: list[np.ndarray] = []
    for numerator, denominator in zip(local_numerators, local_denominators, strict=True):
        assignment = np.zeros((2, 2), dtype=float)
        for prepared_bit in (0, 1):
            if denominator[prepared_bit] <= 0.0:
                raise ValueError("Local readout calibration denominator is zero.")
            assignment[:, prepared_bit] = numerator[:, prepared_bit] / denominator[prepared_bit]
        local_assignments.append(assignment)

    tensor_assignment = tensor_assignment_matrix(local_assignments, num_bits=num_bits)
    return {
        "prepared_states": list(prepared_states),
        "full_assignment_matrix": full_assignment.tolist(),
        "local_assignment_matrices": [assignment.tolist() for assignment in local_assignments],
        "tensor_assignment_matrix": tensor_assignment.tolist(),
    }


def measured_physical_qubits_by_clbit(transpiled_circuit: QuantumCircuit, *, measured_bits: int) -> list[int]:
    measured: list[int | None] = [None] * measured_bits
    for instruction in transpiled_circuit.data:
        if instruction.operation.name != "measure":
            continue
        if not instruction.qubits or not instruction.clbits:
            continue
        clbit_index = transpiled_circuit.find_bit(instruction.clbits[0]).index
        if 0 <= clbit_index < measured_bits:
            measured[clbit_index] = transpiled_circuit.find_bit(instruction.qubits[0]).index
    if any(qubit is None for qubit in measured):
        raise RuntimeError(f"Could not extract {measured_bits} measured physical qubits from transpiled circuit.")
    return [int(qubit) for qubit in measured if qubit is not None]


def build_readout_calibration_circuits(
    transpiled_circuit: QuantumCircuit,
    *,
    measured_bits: int,
) -> tuple[list[QuantumCircuit], tuple[str, ...], list[int]]:
    measured_physical = measured_physical_qubits_by_clbit(transpiled_circuit, measured_bits=measured_bits)
    prepared_states = tuple(format(index, f"0{measured_bits}b") for index in range(2**measured_bits))
    circuits: list[QuantumCircuit] = []
    for prepared_state in prepared_states:
        calibration = QuantumCircuit(
            transpiled_circuit.num_qubits,
            measured_bits,
            name=f"readout_cal_{prepared_state}",
        )
        for clbit_index, physical_qubit in enumerate(measured_physical):
            display_position = measured_bits - 1 - clbit_index
            if prepared_state[display_position] == "1":
                calibration.x(physical_qubit)
            calibration.measure(physical_qubit, clbit_index)
        circuits.append(calibration)
    return circuits, prepared_states, measured_physical


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

    options = sampler_options_payload()
    options["default_shots"] = ZNE_SHOTS
    options["twirling"]["num_randomizations"] = "auto"
    options["twirling"]["shots_per_randomization"] = "auto"
    sampler = SamplerV2(mode=backend, options=options)
    job = sampler.run(circuits, shots=ZNE_SHOTS)
    if submit_only:
        return str(job.job_id()), None
    return str(job.job_id()), job.result()


def run_readout_calibration_job(
    backend: Any,
    calibration_circuits: list[QuantumCircuit],
    *,
    submit_only: bool,
) -> tuple[str, Any | None]:
    from qiskit_ibm_runtime import SamplerV2

    options = sampler_options_payload()
    options["default_shots"] = READOUT_CAL_SHOTS
    sampler = SamplerV2(mode=backend, options=options)
    job = sampler.run(calibration_circuits, shots=READOUT_CAL_SHOTS)
    if submit_only:
        return str(job.job_id()), None
    return str(job.job_id()), job.result()


def baseline_record(index: int) -> dict[str, Any]:
    summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
    for record in summary["records"]:
        if int(record["uniform_index"]) == index:
            return record
    raise ValueError(f"No baseline record for uniform index {index}.")


def build_plan(*, index: int) -> dict[str, Any]:
    full = load_full_qasm_circuit()
    active_circuit, active_physical = compress_to_active_register(full)
    record = baseline_record(index)
    return {
        "captured_at_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "plan",
        "execute_to_submit_hardware": False,
        "instance_id": "operator_loschmidt_echo_49x648",
        "qasm_path": str(QASM_PATH.relative_to(ROOT)),
        "backend": BACKEND,
        "uniform_index": index,
        "bitstring": record["bitstring"],
        "baseline": {
            "global_rescaled_ratio": record["hardware_global_rescaled_ratio"],
            "delta_weighted_term": record["hardware_delta_weighted_term"],
            "delta0_weighted_term": record["hardware_delta0_weighted_term"],
            "hardware_job_id": record["hardware_job_id"],
        },
        "zne": {
            "target": "delta numerator only",
            "factors": list(ZNE_FACTORS),
            "shots": ZNE_SHOTS,
            "ratio_formula": "zne_delta_weighted_term / existing_baseline_delta0_weighted_term",
            "note": "Diagnostic spotcheck only; not a replacement for N_int=8 aggregate.",
        },
        "mitigation": {
            "dynamical_decoupling": "XY4",
            "pauli_twirling": "active",
            "zne": "local-fold factors 1,3,5",
            "readout_mitigation": False,
            "postselection": False,
        },
        "circuit": {
            "declared_qubits": full.num_qubits,
            "active_qubits": active_circuit.num_qubits,
            "active_physical_qubits": active_physical,
            "observable_active": list(OBSERVABLE_ACTIVE),
        },
    }


def run_zne_spotcheck(*, index: int, submit_only: bool) -> dict[str, Any]:
    plan = build_plan(index=index)
    full = load_full_qasm_circuit()
    active_circuit, _ = compress_to_active_register(full)
    measurement = build_measurement_circuit(active_circuit, bitstring=plan["bitstring"])

    service = build_runtime_service()
    backend = service.backend(BACKEND)
    transpiled_base = transpile(
        measurement,
        backend=backend,
        optimization_level=OPTIMIZATION_LEVEL,
        seed_transpiler=SEED_TRANSPILER,
        initial_layout=list(INITIAL_LAYOUT),
    )
    folded = [fold_transpiled_two_qubit_local(transpiled_base, factor) for factor in ZNE_FACTORS]
    job_id, result = run_sampler_job(backend, folded, submit_only=submit_only)

    payload = dict(plan)
    payload.update(
        {
            "captured_at_utc": datetime.now(timezone.utc).isoformat(),
            "mode": "hardware-run",
            "execute_to_submit_hardware": True,
            "submit_only": submit_only,
            "hardware_job_id": job_id,
            "transpiled": {
                "base_depth": int(transpiled_base.depth()),
                "base_size": int(transpiled_base.size()),
                "folded": [
                    {"factor": factor, "depth": int(circuit.depth()), "size": int(circuit.size())}
                    for factor, circuit in zip(ZNE_FACTORS, folded, strict=True)
                ],
            },
        }
    )
    if result is not None:
        sigma = input_parity(plan["bitstring"])
        counts_list = extract_counts_list_from_sampler_result(result, n_items=len(ZNE_FACTORS))
        weighted_terms = [
            sigma * counts_to_z_product_expectation(counts)
            for counts in counts_list
        ]
        zne = extrapolate_zero_noise(factors=ZNE_FACTORS, values=weighted_terms)
        zne_ratio = zne["zero_noise_estimate"] / float(plan["baseline"]["delta0_weighted_term"])
        payload["records"] = [
            {"factor": factor, "counts": counts, "weighted_term": weighted}
            for factor, counts, weighted in zip(ZNE_FACTORS, counts_list, weighted_terms, strict=True)
        ]
        payload["zne_result"] = {
            "weighted_term": zne,
            "zne_over_existing_delta0_ratio": zne_ratio,
            "shift_vs_baseline_ratio": zne_ratio - float(plan["baseline"]["global_rescaled_ratio"]),
        }
    return payload


def transpiled_base_for_payload(payload: Mapping[str, Any], backend: Any) -> QuantumCircuit:
    full = load_full_qasm_circuit()
    active_circuit, _ = compress_to_active_register(full)
    measurement = build_measurement_circuit(active_circuit, bitstring=str(payload["bitstring"]))
    return transpile(
        measurement,
        backend=backend,
        optimization_level=OPTIMIZATION_LEVEL,
        seed_transpiler=SEED_TRANSPILER,
        initial_layout=list(INITIAL_LAYOUT),
    )


def apply_readout_mitigation(payload: dict[str, Any], *, submit_only: bool) -> dict[str, Any]:
    if "records" not in payload:
        raise ValueError("Readout mitigation needs ZNE records with counts; run without --submit-only first.")

    service = build_runtime_service()
    backend = service.backend(BACKEND)
    transpiled_base = transpiled_base_for_payload(payload, backend)
    measured_bits = len(OBSERVABLE_ACTIVE)
    calibration_circuits, prepared_states, measured_physical = build_readout_calibration_circuits(
        transpiled_base,
        measured_bits=measured_bits,
    )
    calibration_job_id, calibration_result = run_readout_calibration_job(
        backend,
        calibration_circuits,
        submit_only=submit_only,
    )

    readout_payload: dict[str, Any] = {
        "enabled": True,
        "cal_shots": READOUT_CAL_SHOTS,
        "calibration_job_id": calibration_job_id,
        "submit_only": submit_only,
        "prepared_states": list(prepared_states),
        "measured_physical_qubits_by_clbit": measured_physical,
        "note": "3-bit observable readout calibration applied to the existing ZNE folded delta counts.",
    }
    if calibration_result is None:
        payload["readout_mitigation"] = readout_payload
        return payload

    calibration_counts = extract_counts_list_from_sampler_result(
        calibration_result,
        n_items=len(calibration_circuits),
    )
    assignment_payload = build_readout_assignment_payload(
        calibration_counts,
        prepared_states=prepared_states,
        num_bits=measured_bits,
    )
    readout_payload.update(assignment_payload)
    readout_payload["calibration_counts"] = [
        {"prepared_state": prepared_state, "counts": counts}
        for prepared_state, counts in zip(prepared_states, calibration_counts, strict=True)
    ]

    sigma = input_parity(str(payload["bitstring"]))
    full_matrix = np.asarray(assignment_payload["full_assignment_matrix"], dtype=float)
    tensor_matrix = np.asarray(assignment_payload["tensor_assignment_matrix"], dtype=float)
    corrected_records: list[dict[str, Any]] = []
    full_weighted_terms: list[float] = []
    tensor_weighted_terms: list[float] = []
    for record in payload["records"]:
        counts = record["counts"]
        full_mitigated = mitigate_counts_with_assignment(counts, full_matrix, num_bits=measured_bits)
        tensor_mitigated = mitigate_counts_with_assignment(counts, tensor_matrix, num_bits=measured_bits)
        for mitigated in (full_mitigated, tensor_mitigated):
            mitigated["weighted_term"] = sigma * float(mitigated["expectation"])
            mitigated["clipped_weighted_term"] = sigma * float(mitigated["clipped_expectation"])
        full_weighted_terms.append(float(full_mitigated["weighted_term"]))
        tensor_weighted_terms.append(float(tensor_mitigated["weighted_term"]))
        corrected_records.append(
            {
                "factor": record["factor"],
                "full_assignment": full_mitigated,
                "local_tensor": tensor_mitigated,
            }
        )

    full_zne = extrapolate_zero_noise(factors=ZNE_FACTORS, values=full_weighted_terms)
    tensor_zne = extrapolate_zero_noise(factors=ZNE_FACTORS, values=tensor_weighted_terms)
    baseline_ratio = float(payload["baseline"]["global_rescaled_ratio"])
    denominator = float(payload["baseline"]["delta0_weighted_term"])
    readout_payload["records"] = corrected_records
    readout_payload["zne_result"] = {
        "full_assignment": {
            "weighted_term": full_zne,
            "zne_over_existing_delta0_ratio": float(full_zne["zero_noise_estimate"]) / denominator,
            "shift_vs_baseline_ratio": (float(full_zne["zero_noise_estimate"]) / denominator) - baseline_ratio,
        },
        "local_tensor": {
            "weighted_term": tensor_zne,
            "zne_over_existing_delta0_ratio": float(tensor_zne["zero_noise_estimate"]) / denominator,
            "shift_vs_baseline_ratio": (float(tensor_zne["zero_noise_estimate"]) / denominator) - baseline_ratio,
        },
    }
    payload["readout_mitigation"] = readout_payload
    payload["mitigation"]["readout_mitigation"] = "3-bit full assignment and local tensor calibration"
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=int, default=1, choices=range(8))
    parser.add_argument("--execute", action="store_true", help="Submit one IBM Runtime ZNE job.")
    parser.add_argument("--readout", action="store_true", help="Also submit a readout-calibration job and apply it.")
    parser.add_argument("--input-json", default=None, help="Apply --readout to an existing ZNE result JSON instead of rerunning ZNE.")
    parser.add_argument("--submit-only", action="store_true")
    parser.add_argument("--output-json", default=None)
    args = parser.parse_args()

    if args.input_json:
        payload = json.loads(Path(args.input_json).read_text(encoding="utf-8"))
        if args.execute and args.readout:
            payload = apply_readout_mitigation(payload, submit_only=args.submit_only)
        else:
            payload["mode"] = "plan"
            payload["readout_plan"] = {
                "execute_to_submit_hardware": bool(args.execute and args.readout),
                "cal_shots": READOUT_CAL_SHOTS,
                "note": "Add --execute --readout to submit only the readout calibration for this existing ZNE result.",
            }
    else:
        payload = run_zne_spotcheck(index=args.index, submit_only=args.submit_only) if args.execute else build_plan(index=args.index)
        if args.execute and args.readout and not args.submit_only:
            payload = apply_readout_mitigation(payload, submit_only=False)
    output_path = Path(args.output_json) if args.output_json else (SPOTCHECK_OUTPUT_PATH if args.execute else PLAN_OUTPUT_PATH)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    print(f"mode={payload['mode']}")
    print(f"backend={payload['backend']}")
    print(f"uniform_index={payload['uniform_index']}")
    print(f"baseline_ratio={payload['baseline']['global_rescaled_ratio']}")
    print(f"zne_factors={payload['zne']['factors']}")
    print(f"shots={payload['zne']['shots']}")
    if "zne_result" in payload:
        print(f"zne_ratio={payload['zne_result']['zne_over_existing_delta0_ratio']}")
        print(f"shift_vs_baseline={payload['zne_result']['shift_vs_baseline_ratio']}")
    if "readout_mitigation" in payload:
        readout = payload["readout_mitigation"]
        print(f"readout_calibration_job_id={readout['calibration_job_id']}")
        if "zne_result" in readout:
            print(f"readout_full_zne_ratio={readout['zne_result']['full_assignment']['zne_over_existing_delta0_ratio']}")
            print(f"readout_tensor_zne_ratio={readout['zne_result']['local_tensor']['zne_over_existing_delta0_ratio']}")
    print(f"output_json={output_path}")
    if not args.execute:
        print("dry_run=true; add --execute to submit IBM Runtime jobs.")


if __name__ == "__main__":
    main()
