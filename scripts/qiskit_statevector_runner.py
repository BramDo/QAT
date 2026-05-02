"""Qiskit statevector runner for the six-qubit OLE surrogate.

This is a thin, reader-friendly entrypoint around ``run_6q_statevector.py``.
It computes the all-zero delta and delta=0 signals exactly with Qiskit's
statevector tools.
"""

from __future__ import annotations

from run_6q_statevector import main


if __name__ == "__main__":
    main()
