"""
Exact classical solver for reference ground state energies.

Uses NumPyMinimumEigensolver for full diagonalization — equivalent
to full configuration interaction (FCI) in the given basis.
"""

from qiskit_algorithms import NumPyMinimumEigensolver


def exact_energy(hamiltonian):
    """
    Compute the exact ground state energy via full diagonalization.

    This is equivalent to FCI in the qubit-mapped Hilbert space.
    Used as the reference for evaluating VQE accuracy.

    Args:
        hamiltonian: SparsePauliOp qubit Hamiltonian

    Returns:
        float: exact ground state energy (in Hartree)
    """
    solver = NumPyMinimumEigensolver()
    result = solver.compute_minimum_eigenvalue(hamiltonian)
    return result.eigenvalue.real
