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
- `scripts/qiskit_statevector_runner.py`
  - explicit Qiskit statevector runner entrypoint.
- `scripts/tn_mps_runner.py`
  - Qiskit Aer matrix-product-state simulator runner with sampled measurements.
- `scripts/run_49q_kingston_small_budget.py`
  - optional IBM Runtime runner configured like the small-budget `~0.209` 49Q
    Kingston result.
- `scripts/run_49q_zne_spotcheck.py`
  - optional one-string IBM Runtime ZNE diagnostic for the same Kingston layout.
- `scripts/aggregate_49q_zne_readout_subset.py`
  - aggregates the currently available 49Q ZNE/readout spotchecks.
- `results/reference_6q_and_scaling_results.json`
  - compact reference values from the local study.
- `results/49q_zne_spotcheck_index1_20260503.json`
  - one extra 49Q Kingston ZNE spotcheck on uniform string index `1`.
- `results/49q_zne_readout_spotcheck_index1_20260503.json`
  - the same spotcheck with a 3-bit readout calibration applied afterward.
- `results/49q_zne_readout_subset3_summary.json`
  - diagnostic average over the three available ZNE/readout strings.
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
python scripts/qiskit_statevector_runner.py
python scripts/tn_mps_runner.py --shots 20000
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

## Simulator Runners

There are two local simulator entrypoints:

```bash
python scripts/qiskit_statevector_runner.py
```

This computes the six-qubit all-zero signal exactly with Qiskit's statevector
tools. It should print:

```text
global-rescale ratio: 0.969064735901
```

The second runner uses Qiskit Aer with a matrix-product-state simulator:

```bash
python scripts/tn_mps_runner.py --shots 20000
```

This runner samples counts from the delta and delta=0 circuits and estimates the
same global-rescale ratio from measured Z-observable parity. Because it is a
shot-based simulator path, its value will fluctuate slightly around the exact
statevector value.

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

## Optional 49Q Hardware Runner

The repo includes a 49Q hardware runner configured to match the small-budget
Kingston run that produced the `0.2093896573 +/- 0.0066239080` result:

```bash
python scripts/run_49q_kingston_small_budget.py
```

By default this is a dry run. It parses the 49Q QASM, compresses the public
156-qubit physical-register circuit to its 49 active qubits, verifies the
`delta=0` transform, and writes a plan JSON. It does not submit IBM jobs unless
you explicitly add `--execute`.

Baseline configuration:

- backend: `ibm_kingston`
- circuit: `operator_loschmidt_echo_49x648`, `L=3`, `delta=0.15`
- sample: `N_int=8` predeclared uniform bitstrings, seed `20260430`
- shots: `21 * 251 = 5271` per delta/delta0 circuit
- transpiler: optimization level `1`, seed `424242`
- layout: fixed 49Q Kingston layout from the small-budget run
- mitigation: global rescaling, DD `XY4`, active Pauli twirling
- not used: readout mitigation, ZNE, postselection

Important: the `0.2093896573` N_int=8 aggregate was not obtained with readout
mitigation or ZNE. Those were tested later on two selected bitstrings as an A/B
diagnostic. Readout-only did not improve the two tested strings, and ZNE was
mixed: one string stayed flat while another improved. That A/B test is recorded
in `results/49q_small_budget_ab_mitigation_summary.json`, but it is not a new
full N_int=8 estimate.

To submit real IBM Runtime jobs:

```bash
python scripts/run_49q_kingston_small_budget.py --execute
```

That command spends quantum budget. For a cheaper smoke submission, reduce the
number of initializations:

```bash
python scripts/run_49q_kingston_small_budget.py --execute --sample-count 1
```

Use the same Python environment that has `qiskit-ibm-runtime` configured for
your IBM Quantum account.

## Optional 49Q ZNE Spotcheck

ZNE was the only mitigation variant that showed a hint of improvement in the
small A/B diagnostics. For that reason the repo includes a cheaper one-string
spotcheck runner:

```bash
python scripts/run_49q_zne_spotcheck.py
```

The default is again a dry run. The hardware command used for the recorded
spotcheck was:

```bash
python scripts/run_49q_zne_spotcheck.py --execute --index 1 --output-json results/49q_zne_spotcheck_index1_20260503.json
```

Readout calibration can be added without rerunning the folded delta circuits:

```bash
python scripts/run_49q_zne_spotcheck.py --execute --readout --input-json results/49q_zne_spotcheck_index1_20260503.json --output-json results/49q_zne_readout_spotcheck_index1_20260503.json
```

Recorded result on `ibm_kingston`, ZNE job `d7rmj5vljm6s73bah3ig`, readout
calibration job `d7rmsvkf3ras73b743b0`:

| Quantity | Value |
|---|---:|
| uniform string index | `1` |
| baseline global-rescale ratio | `0.2286591607` |
| raw ZNE ratio, factors `1,3,5` | `0.2834528131` |
| ZNE + full-assignment readout | `0.2969777829` |
| ZNE + local-tensor readout | `0.2964000464` |
| tensor-readout shift vs baseline | `+0.0677408857` |

This is a diagnostic result for one additional string. It supports the narrow
statement that ZNE, and on this string ZNE plus readout correction, can improve
some strings in this setup. It is not a new full 49Q aggregate.

The current mitigated subset can be aggregated with:

```bash
python scripts/aggregate_49q_zne_readout_subset.py
```

Current three-string diagnostic subset, indices `0,1,5`:

| Quantity | Mean | Standard error |
|---|---:|---:|
| baseline on same three strings | `0.2190697289` | `0.0149071606` |
| raw ZNE | `0.2660511423` | `0.0400938721` |
| ZNE + full-assignment readout | `0.2680347994` | `0.0396077227` |
| ZNE + local-tensor readout | `0.2704104168` | `0.0422483017` |

This subset was assembled after the earlier A/B diagnostics, so it should not be
presented as an unbiased final estimate. The clean small-budget completion path
is to run the same ZNE+readout pipeline on the remaining original uniform
indices `2,3,4,6,7`, then average all eight predeclared strings.

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
