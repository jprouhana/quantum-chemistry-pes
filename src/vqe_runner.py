"""
VQE execution for molecular ground state energy calculation.

Supports multiple ansatz types and tracks optimization convergence.
"""

import numpy as np
from qiskit.circuit.library import RealAmplitudes, EfficientSU2
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA
from qiskit.primitives import Estimator


def run_vqe(hamiltonian, ansatz_type='RealAmplitudes', reps=3,
            maxiter=200, seed=42):
    """
    Run VQE to find the ground state energy of a qubit Hamiltonian.

    Args:
        hamiltonian: SparsePauliOp qubit Hamiltonian
        ansatz_type: 'RealAmplitudes' or 'EfficientSU2'
        reps: number of ansatz repetitions
        maxiter: maximum COBYLA iterations
        seed: random seed for initial point

    Returns:
        dict with 'energy', 'optimal_params', 'energy_history',
            'num_params', 'ansatz_type'
    """
    num_qubits = hamiltonian.num_qubits

    # build ansatz
    if ansatz_type == 'RealAmplitudes':
        ansatz = RealAmplitudes(
            num_qubits=num_qubits,
            reps=reps,
            entanglement='linear'
        )
    elif ansatz_type == 'EfficientSU2':
        ansatz = EfficientSU2(
            num_qubits=num_qubits,
            reps=reps,
            entanglement='circular'
        )
    else:
        raise ValueError(f"Unknown ansatz type: {ansatz_type}")

    # track convergence
    energy_history = []

    def callback(eval_count, params, value, std):
        energy_history.append(value)

    optimizer = COBYLA(maxiter=maxiter)
    estimator = Estimator()

    vqe = VQE(
        estimator=estimator,
        ansatz=ansatz,
        optimizer=optimizer,
        callback=callback,
    )

    # random initial point
    rng = np.random.default_rng(seed)
    initial_point = rng.uniform(-np.pi, np.pi, ansatz.num_parameters)
    vqe.initial_point = initial_point

    result = vqe.compute_minimum_eigenvalue(hamiltonian)

    return {
        'energy': result.eigenvalue.real,
        'optimal_params': result.optimal_parameters,
        'energy_history': energy_history,
        'num_params': ansatz.num_parameters,
        'ansatz_type': ansatz_type,
    }
