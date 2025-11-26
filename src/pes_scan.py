"""
Potential energy surface scanning utilities.

Sweeps molecular geometries and computes VQE + exact energies
at each point to build a full PES with error analysis.
"""

import numpy as np

from .vqe_runner import run_vqe
from .exact_solver import exact_energy


def scan_pes(molecule_builder, distances, ansatz_type='RealAmplitudes',
             reps=3, maxiter=200, seed=42):
    """
    Scan the potential energy surface over a range of geometries.

    Computes both VQE and exact (FCI) energies at each geometry
    point. The molecule_builder function should accept a single
    geometric parameter (distance or angle) and return
    (qubit_op, problem).

    Args:
        molecule_builder: callable(distance) -> (qubit_op, problem)
        distances: array-like of geometric parameter values
        ansatz_type: 'RealAmplitudes' or 'EfficientSU2'
        reps: ansatz repetitions
        maxiter: maximum VQE optimizer iterations
        seed: random seed

    Returns:
        dict with 'distances', 'vqe_energies', 'exact_energies',
            'errors', 'nuclear_repulsions'
    """
    distances = np.array(distances)
    vqe_energies = []
    exact_energies = []
    errors = []
    nuclear_repulsions = []

    for i, d in enumerate(distances):
        print(f"  [{i + 1}/{len(distances)}] d = {d:.3f} ...", end=" ")

        qubit_op, problem = molecule_builder(d)
        nuc_repulsion = problem.nuclear_repulsion_energy

        # VQE energy (electronic only — add nuclear repulsion)
        vqe_result = run_vqe(
            qubit_op, ansatz_type=ansatz_type, reps=reps,
            maxiter=maxiter, seed=seed + i
        )
        vqe_e = vqe_result['energy'] + nuc_repulsion

        # exact energy
        exact_e = exact_energy(qubit_op) + nuc_repulsion

        error = abs(vqe_e - exact_e)

        vqe_energies.append(vqe_e)
        exact_energies.append(exact_e)
        errors.append(error)
        nuclear_repulsions.append(nuc_repulsion)

        print(f"VQE={vqe_e:.6f}, Exact={exact_e:.6f}, "
              f"Error={error * 1000:.2f} mHa")

    return {
        'distances': distances.tolist(),
        'vqe_energies': vqe_energies,
        'exact_energies': exact_energies,
        'errors': errors,
        'nuclear_repulsions': nuclear_repulsions,
    }
