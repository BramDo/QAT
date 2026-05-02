"""Matrix-product-state simulator runner for the six-qubit OLE surrogate.

This uses Qiskit Aer with ``method="matrix_product_state"``. It follows the
hardware-style measurement flow: prepare a basis state, run the delta and
delta=0 circuits, measure the Z observable support, and estimate the
global-rescale ratio from sampled counts.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Mapping

from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

from run_6q_statevector import (
    METADATA_PATH,
    OUTPUT_PATH,
    input_parity,
    load_circuit,
    zero_delta_rotations,
)


ROOT = Path(__file__).resolve().parents[1]
MPS_OUTPUT_PATH = ROOT / "results" / "6q_mps_simulator_result.json"


def build_measurement_circuit(
    circuit: QuantumCircuit,
    *,
    bitstring: str,
    observable_qubits: list[int],
) -> QuantumCircuit:
    """Prepare a bitstring, apply the OLE circuit, and measure O-support qubits."""

    measured = QuantumCircuit(circuit.num_qubits, len(observable_qubits))
    for qubit, bit in enumerate(bitstring):
        if bit == "1":
            measured.x(qubit)
    measured.compose(circuit, inplace=True)
    for clbit, qubit in enumerate(observable_qubits):
        measured.measure(qubit, clbit)
    return measured


def counts_to_z_product_expectation(counts: Mapping[str, int]) -> float:
    """Convert measured computational-basis counts to a Z-product expectation."""

    total = sum(int(count) for count in counts.values())
    if total == 0:
        raise ValueError("Cannot compute expectation from empty counts.")
    weighted = 0
    for raw_key, count in counts.items():
        key = raw_key.replace(" ", "")
        parity = -1 if key.count("1") % 2 else 1
        weighted += parity * int(count)
    return weighted / total


def run_mps_simulator(*, shots: int, seed: int) -> dict:
    metadata = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
    circuit = load_circuit()
    observable_qubits = [int(q) for q in metadata["local_observable_declared_zero_based"]]
    support_qubits = {0, 1, 2}
    delta = float(metadata["delta"])
    bitstring = "0" * circuit.num_qubits
    sigma = input_parity(bitstring, observable_qubits)

    delta0_circuit, removed = zero_delta_rotations(
        circuit,
        support_qubits=support_qubits,
        target_angle=2.0 * delta,
    )
    measurement_circuits = [
        build_measurement_circuit(circuit, bitstring=bitstring, observable_qubits=observable_qubits),
        build_measurement_circuit(delta0_circuit, bitstring=bitstring, observable_qubits=observable_qubits),
    ]
    simulator = AerSimulator(method="matrix_product_state", seed_simulator=seed)
    job = simulator.run(measurement_circuits, shots=shots)
    result = job.result()
    delta_counts = result.get_counts(0)
    delta0_counts = result.get_counts(1)

    delta_expectation = counts_to_z_product_expectation(delta_counts)
    delta0_expectation = counts_to_z_product_expectation(delta0_counts)
    delta_weighted = sigma * delta_expectation
    delta0_weighted = sigma * delta0_expectation
    ratio = delta_weighted / delta0_weighted

    return {
        "instance_id": metadata["instance_id"],
        "simulator": "qiskit-aer matrix_product_state",
        "shots": shots,
        "seed_simulator": seed,
        "bitstring": bitstring,
        "observable_qubits": observable_qubits,
        "perturbation_support_qubits": sorted(support_qubits),
        "removed_delta_rotations": removed,
        "source_gate_count": len(circuit.data),
        "depth": circuit.depth(),
        "delta_counts": delta_counts,
        "delta0_counts": delta0_counts,
        "delta_weighted_term": delta_weighted,
        "delta0_weighted_term": delta0_weighted,
        "global_rescale_ratio": ratio,
        "statevector_reference_output": str(OUTPUT_PATH.relative_to(ROOT)),
        "notes": [
            "This is a sampled MPS simulator run, so small shot noise is expected.",
            "It mirrors the hardware-style measurement path more closely than the exact statevector runner.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shots", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20260502)
    args = parser.parse_args()

    payload = run_mps_simulator(shots=args.shots, seed=args.seed)
    MPS_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    MPS_OUTPUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    print(f"instance: {payload['instance_id']}")
    print(f"simulator: {payload['simulator']}")
    print(f"shots: {payload['shots']}")
    print(f"removed delta rotations: {payload['removed_delta_rotations']}")
    print(f"delta weighted term: {payload['delta_weighted_term']:.12f}")
    print(f"delta0 weighted term: {payload['delta0_weighted_term']:.12f}")
    print(f"global-rescale ratio: {payload['global_rescale_ratio']:.12f}")
    print(f"wrote: {MPS_OUTPUT_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
