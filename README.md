# QAT: Quantum Advantage Tracker OLE Mini Repro

This repository is a small companion repo for a LinkedIn article about the
[Quantum Advantage Tracker](https://quantum-advantage-tracker.github.io/trackers/observable-estimations)
and the `operator_loschmidt_echo` circuit family.

It is intentionally not a full reproduction of the public tracker submissions.
The full IBM tracker protocol used hundreds of random initializations and
specific IBM processors such as `ibm_boston` and `ibm_pittsburgh`. This repo
keeps only the files needed for a reader to understand and run the small
six-qubit version.

## What Is Included

- `data/operator_loschmidt_echo_49x648_surrogate_6q.qasm`
  - a six-qubit surrogate induced from the public `49x648` OLE circuit.
- `data/operator_loschmidt_echo_49x648_surrogate_6q.json`
  - metadata for the six-qubit circuit.
- `data/49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm`
  - the public 49Q `operator_loschmidt_echo_49x648` OpenQASM 3.0 circuit.
- `data/operator_loschmidt_echo_49x648_49q.json`
  - metadata for the public 49Q circuit.
- `scripts/run_6q_statevector.py`
  - standalone Qiskit script that computes the all-zero six-qubit OLE signal.
- `results/reference_6q_and_scaling_results.json`
  - compact reference values from the local study.
- `docs/article.md`
  - longer article draft.
- `docs/linkedin_quantum_advantage_tracker_post_plain.txt`
  - LinkedIn-ready plain text.
- `docs/assets/*.png`
  - formula cards for LinkedIn.

## Quick Start

Use Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_6q_statevector.py
```

Expected local output is close to:

```text
delta weighted term: 0.969064735901
delta0 weighted term: 1.000000000000
global-rescale ratio: 0.969064735901
```

The local statevector result is close to the BP-TN BD64 reference
`0.9602546583`, but the methods are not identical. The statevector script is
included because it is transparent and easy to run.

## Six-Qubit Reference Results

| Method | Value |
|---|---:|
| Local BP-TN BD64 | `0.9602546583` |
| IBM Kingston raw | `0.8828125` |
| IBM Kingston global rescaling + tensor readout | `0.8389109357` |
| IBM Kingston ZNE + readout tensor | `0.9154327081` |
| Local statevector all-zero ratio from this repo | about `0.9690647359` |

The six-qubit circuit has:

- `6` qubits
- observable `Z4`
- `591` source gates
- circuit depth `149`
- `60` CZ gates

## Optional 49Q QASM

The repo also includes the public tracker `49x648` QASM:

```text
data/49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm
```

Important details:

- OpenQASM version: `3.0`
- declared qubits in the physical-register file: `156`
- active non-barrier qubits: `49`
- circuit depth: `149`
- CZ gates: `648`
- observable in the tracker instance: `Z52 Z59 Z72`
- parameters: `L=3`, `b=0.25`, `delta=0.15`

This file is useful for readers who want to inspect the full tracker circuit or
try their own transpilation/hardware workflow. It is not used by the quick local
statevector demo, because simulating the full 49Q circuit exactly is not
practical on a normal laptop.

## Full Tracker Context

The public tracker includes larger `operator_loschmidt_echo` instances:

| Instance | Public method | Public value | Resource |
|---|---:|---:|---|
| `49x648`, `L=3` | global rescaling | `0.824` | `ibm_boston` |
| `49x1296`, `L=6` | global rescaling | `0.649` | `ibm_pittsburgh` |
| `49x648`, `L=3` | BP-TN BD512 | `0.821658489` | Intel Xeon |
| `49x1296`, `L=6` | Single-path Monte Carlo | `0.619` | Apple M1 Max |

My small-budget 49-qubit Kingston run is included only as context:

```text
49Q Kingston small-budget global-rescaling estimate:
0.2093896573 +/- 0.0066239080 standard error
```

That result used `N_int=8` predeclared uniform basis states, not the public
tracker-scale `N_int=500` protocol. It should be described as a finite-budget
hardware diagnostic, not as a tracker reproduction.

## Sources

- Quantum Advantage Tracker:
  https://quantum-advantage-tracker.github.io/trackers/observable-estimations
- IBM blog:
  https://www.ibm.com/quantum/blog/quantum-advantage-tracker
- OLE circuit model:
  https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/tree/main/data/observable-estimations/circuit-models/operator_loschmidt_echo
- Tracker issue #38:
  https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/38
- Tracker issue #78:
  https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/78

## Attribution

The public Quantum Advantage Tracker repository is Apache-2.0 licensed. The
six-qubit QASM in this repo is a local surrogate derived from the public
tracker OLE circuit family, with metadata included in `data/`.
