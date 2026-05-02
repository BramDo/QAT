"""Run the six-qubit OLE surrogate with a transparent statevector calculation."""

from __future__ import annotations

import json
from pathlib import Path

from qiskit import QuantumCircuit, qasm3
from qiskit.quantum_info import SparsePauliOp, Statevector


ROOT = Path(__file__).resolve().parents[1]
QASM_PATH = ROOT / "data" / "operator_loschmidt_echo_49x648_surrogate_6q.qasm"
METADATA_PATH = ROOT / "data" / "operator_loschmidt_echo_49x648_surrogate_6q.json"
OUTPUT_PATH = ROOT / "results" / "6q_statevector_result.json"


def load_circuit() -> QuantumCircuit:
    return qasm3.loads(QASM_PATH.read_text(encoding="utf-8"))


def zero_delta_rotations(
    circuit: QuantumCircuit,
    *,
    support_qubits: set[int],
    target_angle: float,
    tolerance: float = 1e-9,
) -> tuple[QuantumCircuit, int]:
    """Remove the inferred rz(2*delta) perturbation rotations."""

    zeroed = QuantumCircuit(circuit.num_qubits, name=f"{circuit.name}_delta0")
    zeroed.global_phase = circuit.global_phase
    removed = 0
    for instruction in circuit.data:
        operation = instruction.operation
        qargs = [circuit.find_bit(qubit).index for qubit in instruction.qubits]
        if (
            operation.name == "rz"
            and len(qargs) == 1
            and qargs[0] in support_qubits
            and len(operation.params) == 1
            and abs(float(operation.params[0]) - target_angle) <= tolerance
        ):
            removed += 1
            continue
        zeroed.append(operation, qargs)
    return zeroed, removed


def z_observable_expectation(circuit: QuantumCircuit, *, bitstring: str, observable_qubits: list[int]) -> float:
    """Prepare a bitstring, apply the circuit, and measure the Z-product observable."""

    prep = QuantumCircuit(circuit.num_qubits)
    for qubit, bit in enumerate(bitstring):
        if bit == "1":
            prep.x(qubit)
    state = Statevector.from_instruction(prep.compose(circuit))

    pauli = ["I"] * circuit.num_qubits
    for qubit in observable_qubits:
        # Qiskit Pauli strings are written most-significant qubit first.
        pauli[circuit.num_qubits - 1 - qubit] = "Z"
    return float(state.expectation_value(SparsePauliOp("".join(pauli))).real)


def input_parity(bitstring: str, observable_qubits: list[int]) -> int:
    ones = sum(1 for qubit in observable_qubits if bitstring[qubit] == "1")
    return -1 if ones % 2 else 1


def main() -> None:
    metadata = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
    circuit = load_circuit()
    observable_qubits = [int(q) for q in metadata["local_observable_declared_zero_based"]]
    support_qubits = {0, 1, 2}
    delta = float(metadata["delta"])
    bitstring = "0" * circuit.num_qubits

    delta0_circuit, removed = zero_delta_rotations(
        circuit,
        support_qubits=support_qubits,
        target_angle=2.0 * delta,
    )
    sigma = input_parity(bitstring, observable_qubits)
    delta_expectation = z_observable_expectation(
        circuit,
        bitstring=bitstring,
        observable_qubits=observable_qubits,
    )
    delta0_expectation = z_observable_expectation(
        delta0_circuit,
        bitstring=bitstring,
        observable_qubits=observable_qubits,
    )
    delta_weighted = sigma * delta_expectation
    delta0_weighted = sigma * delta0_expectation
    ratio = delta_weighted / delta0_weighted

    result = {
        "instance_id": metadata["instance_id"],
        "bitstring": bitstring,
        "observable_qubits": observable_qubits,
        "perturbation_support_qubits": sorted(support_qubits),
        "removed_delta_rotations": removed,
        "source_gate_count": len(circuit.data),
        "depth": circuit.depth(),
        "delta_weighted_term": delta_weighted,
        "delta0_weighted_term": delta0_weighted,
        "global_rescale_ratio": ratio,
        "notes": [
            "This is a transparent statevector calculation for the six-qubit surrogate.",
            "It is not the full tracker protocol and does not include hardware noise or sampling error.",
        ],
    }
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")

    print(f"instance: {result['instance_id']}")
    print(f"qubits: {circuit.num_qubits}")
    print(f"source gates: {len(circuit.data)}")
    print(f"depth: {circuit.depth()}")
    print(f"removed delta rotations: {removed}")
    print(f"delta weighted term: {delta_weighted:.12f}")
    print(f"delta0 weighted term: {delta0_weighted:.12f}")
    print(f"global-rescale ratio: {ratio:.12f}")
    print(f"wrote: {OUTPUT_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
