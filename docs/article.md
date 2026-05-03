# Understanding the Quantum Advantage Tracker Through a Six-Qubit Loschmidt Echo

Draft status: working article draft.

## Why This Tracker Matters

The [Quantum Advantage Tracker](https://quantum-advantage-tracker.github.io/trackers/observable-estimations) is an open community effort to compare quantum and classical methods on concrete benchmark instances. IBM describes the motivation well in its [Quantum Advantage Tracker blog](https://www.ibm.com/quantum/blog/quantum-advantage-tracker): advantage should not be treated as a one-time announcement, but as a process where quantum and classical methods pressure-test each other in public.

The part I am studying is the observable-estimation tracker. These submissions report expectation values of observables, together with the method, runtime, and compute resources. My goal here is modest: explain the theory and circuit structure through a six-qubit version before discussing the full 49-qubit hardware runs.

## The Problem: Operator Loschmidt Echo

The circuit family is called `operator_loschmidt_echo`. The tracker's circuit README defines the Operator Loschmidt Echo as

```text
f_delta(O) = 2^-n Tr(U O U^\dagger V_delta^\dagger U O U^\dagger V_delta).
```

Here:

- `O` is a Pauli-Z observable, such as `Z52 Z59 Z72` in the public 49-qubit instance.
- `U` is a Trotterized time-evolution circuit.
- `V_delta = exp(-i delta G)` is a small perturbation.
- `f_delta(O)` measures how much the evolved observable changes when the system is perturbed.

For small `delta`, the OLE is related to an out-of-time-order correlator. In plain language: it probes how quickly local operator information spreads through the circuit.

The tracker also gives a measurement-friendly identity:

```text
f_delta(O) = 2^-n sum_z sigma_z <z| U^\dagger V_delta^\dagger U O U^\dagger V_delta |z>,
```

where `z` ranges over computational basis states and `sigma_z = <z|O|z>` is just the parity of the input bitstring on the observable qubits. Instead of summing all `2^n` basis states, the hardware protocol samples a number of random initializations `N_int`.

## Public Tracker Reference Points

The full tracker instances are large:

| Instance | Method | Value | Runtime | Resource |
|---|---:|---:|---:|---|
| `49x648`, `L=3` | IBM global rescaling | `0.824` | `4260s` | `ibm_boston` |
| `49x1296`, `L=6` | IBM global rescaling | `0.649` | `4120s` | `ibm_pittsburgh` |
| `49x648`, `L=3` | BP-TN BD512, classical | `0.821658489` | `149s` | Intel Xeon |
| `49x1296`, `L=6` | Single-path Monte Carlo | `0.619` | `1492s` | Apple M1 Max |

The detailed IBM issue for `49x648` says the hardware run used `N_int=500`, `21` twirls per initialization, and `251` shots per twirl. The later `49x1296` Pittsburgh issue uses the same `N_int=500`, `21 x 251` shot structure, but with `L=6`.

## Why Start With Six Qubits?

The full `49x648` circuit is already too large to understand visually. So I use a local six-qubit surrogate induced from the canonical `49x648` active-register circuit.

Local six-qubit surrogate:

| Property | Value |
|---|---:|
| Instance | `operator_loschmidt_echo_49x648_surrogate_6q` |
| Parent | `operator_loschmidt_echo_49x648` |
| Qubits | `6` |
| Observable | `Z4` |
| Kept source gates | `591` |
| Circuit depth | `149` |
| CZ gates | `60` |
| Surviving perturbation-support qubits | `3` |

This is not the public tracker instance. It is a teaching and debugging version. But it keeps the important shape: Trotterized dynamics, a perturbation region, and a Pauli-Z observable.

The first gates already show the structure. Each qubit gets one-qubit rotations such as `rx` and `rz`, then non-overlapping pairs receive `cz` gates. Repeated layers build up an interacting circuit. The perturbation is later represented by small `rz(2*delta)` rotations on selected support qubits in the active-register representation.

## What We Measured on the Six-Qubit Version

For the six-qubit circuit, the local and hardware results are close enough to be educational:

| Method | Value |
|---|---:|
| Local BP-TN BD64 | `0.9602546583` |
| IBM Kingston raw | `0.8828125` |
| IBM Kingston global rescaling + tensor readout | `0.8389109357` |
| IBM Kingston ZNE + readout tensor | `0.9154327081` |

This is the cleanest story for an article:

1. The local tensor-network result gives a high reference value.
2. Raw hardware is lower because the circuit is noisy.
3. Global rescaling compares the `delta=0.15` signal to a `delta=0` signal.
4. ZNE plus readout mitigation recovers part of the lost signal.

Even at six qubits, the circuit is not just a toy Bell-state experiment. It has 591 gates and 60 CZ gates, so the mitigation story is visible.

## Global Rescaling in One Sentence

Global rescaling estimates

```text
mitigated signal at delta = measured signal at delta / measured signal at delta=0.
```

The idea is that the `delta=0` circuit should retain the same broad hardware attenuation but without the perturbation. Dividing by it can remove part of the global signal loss. This is useful, but it is not magic: if the deep `delta=0.15` numerator loses too much information, global rescaling cannot fully recover it.

## From Six Qubits to Forty-Nine

The public `49x648` circuit has 49 active qubits and 648 CZ gates. My small-budget 49-qubit Kingston run used the same `49x648` family, but only `N_int=8` predeclared uniform basis states, not the tracker-scale `N_int=500`.

Small-budget 49Q result:

| Quantity | Value |
|---|---:|
| Backend | `ibm_kingston` |
| Circuit | `operator_loschmidt_echo_49x648` |
| Sample count | `N_int=8` |
| Shots per basis term | `5271` |
| Mitigation | global rescaling, DD `XY4`, active Pauli twirling |
| Mean hardware ratio | `0.2093896573` |
| Standard error | `0.0066239080` |
| 95% normal interval | `[0.1964, 0.2224]` |

This is not a reproduction of the tracker value `0.824`. It is better described as:

```text
49Q Kingston small-budget global-rescaling hardware estimate, L=3, N_int=8, no postselection.
```

For a hobby-shot-budget experiment, I still think it is a useful result. It is a real 49-qubit hardware measurement, with a fixed predeclared sample and no postselection. The result is low, but it teaches us where the difficulty is: the deep `delta=0.15` numerator is heavily attenuated, while the shallow `delta=0` denominator stays high.

## What The Mitigation A/B Test Suggests

I also ran a small A/B mitigation test on two of the existing 49Q uniform samples.

| Test | Uniform index 0 | Uniform index 5 |
|---|---:|---:|
| Baseline global rescaling | `0.189827` | `0.238723` |
| Readout tensor | `0.187700` | `0.238978` |
| ZNE + tensor-readout, delta-only spotcheck | `0.187787` | `0.327044` |

Readout mitigation alone did not help. ZNE was mixed: one sample stayed flat, one improved strongly. That is not enough to update the full `N_int=8` estimate, but it is enough to justify a small follow-up test.

I then ran one extra ZNE spotcheck on uniform index `1`, using the same Kingston
layout. The baseline ratio for that string was `0.228659`. Raw ZNE with folding
factors `1,3,5` moved it to `0.283453`. Adding a separate 3-bit readout
calibration moved it slightly further: `0.296978` with full-assignment
correction and `0.296400` with local-tensor correction. This is encouraging,
but it is still one string, not a new aggregate result.

Combining the three currently available mitigated strings, indices `0,1,5`,
gives a diagnostic mean of `0.270410 +/- 0.042248` standard error for
ZNE + local-tensor readout. The baseline mean on the same three strings is
`0.219070 +/- 0.014907`. This is a useful subset check, but not an unbiased
final estimate because the subset was assembled after looking at the mitigation
diagnostics. A clean small-budget completion would run exactly the same pipeline
on the remaining original indices `2,3,4,6,7` and average all eight predeclared
strings.

## Claim Boundary

The tracker results are serious benchmark submissions with large initialization counts and specific IBM processors. My results are a learning and reproduction effort from a limited-budget environment.

The six-qubit circuit is the right entry point for explanation. The 49-qubit result is the right entry point for humility: running the real family on real hardware is possible, but reproducing the public tracker value requires the full protocol, hardware layout, and budget.

## Sources

- [Quantum Advantage Tracker: observable estimations](https://quantum-advantage-tracker.github.io/trackers/observable-estimations)
- [Quantum Advantage Tracker GitHub repository](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io)
- [Operator Loschmidt Echo circuit README](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/tree/main/data/observable-estimations/circuit-models/operator_loschmidt_echo)
- [IBM blog: Quantum Advantage Tracker](https://www.ibm.com/quantum/blog/quantum-advantage-tracker)
- [Tracker issue #38: 49x648, ibm_boston](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/38)
- [Tracker issue #78: 49x1296, ibm_pittsburgh](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/78)
