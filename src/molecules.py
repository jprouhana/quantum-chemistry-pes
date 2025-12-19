"""
Molecular Hamiltonian construction for potential energy surface studies.

Uses PySCF driver with Jordan-Wigner mapping to convert molecular
geometries into qubit Hamiltonians suitable for VQE.
"""

import numpy as np
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.units import DistanceUnit


def build_h2(distance=0.735):
    """
    Build the qubit Hamiltonian for H2 at a given bond distance.

    Uses STO-3G minimal basis set -> 4 spin-orbitals -> 4 qubits
    with Jordan-Wigner mapping.

    Args:
        distance: H-H bond distance in Angstroms

    Returns:
        qubit_op: SparsePauliOp qubit Hamiltonian
        problem: ElectronicStructureProblem (contains nuclear repulsion etc)
    """
    atom_string = f"H 0.0 0.0 0.0; H 0.0 0.0 {distance}"

    driver = PySCFDriver(
        atom=atom_string,
        basis='sto3g',
        charge=0,
        spin=0,
        unit=DistanceUnit.ANGSTROM
    )

    problem = driver.run()
    mapper = JordanWignerMapper()
    second_q_op = problem.hamiltonian.second_q_op()
    qubit_op = mapper.map(second_q_op)

    return qubit_op, problem


def build_lih(distance=1.546):
    """
    Build the qubit Hamiltonian for LiH at a given bond distance.

    Uses STO-3G basis set. LiH has 4 electrons in 6 spatial orbitals,
    giving 12 spin-orbitals -> 12 qubits with Jordan-Wigner.

    For tractability we use an active space of 2 electrons in 2 orbitals
    (4 qubits) when available, otherwise the full space.

    Args:
        distance: Li-H bond distance in Angstroms

    Returns:
        qubit_op: SparsePauliOp qubit Hamiltonian
        problem: ElectronicStructureProblem
    """
    atom_string = f"Li 0.0 0.0 0.0; H 0.0 0.0 {distance}"

    driver = PySCFDriver(
        atom=atom_string,
        basis='sto3g',
        charge=0,
        spin=0,
        unit=DistanceUnit.ANGSTROM
    )

    problem = driver.run()

    # apply active space reduction for tractability
    # keep 2 electrons in 2 spatial orbitals -> 4 qubits
    from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
    transformer = ActiveSpaceTransformer(
        num_electrons=2,
        num_spatial_orbitals=3
    )
    problem = transformer.transform(problem)

    mapper = JordanWignerMapper()
    second_q_op = problem.hamiltonian.second_q_op()
    qubit_op = mapper.map(second_q_op)

    return qubit_op, problem


def build_beh2(angle=180.0):
    """
    Build the qubit Hamiltonian for BeH2 at a given H-Be-H bond angle.

    The Be-H bond distance is fixed at 1.326 Angstroms (equilibrium).
    Scanning the angle from ~100 to 180 degrees traces the bending PES.

    Uses STO-3G basis with active space reduction for tractability.

    Args:
        angle: H-Be-H bond angle in degrees

    Returns:
        qubit_op: SparsePauliOp qubit Hamiltonian
        problem: ElectronicStructureProblem
    """
    # Be at origin, two H atoms at fixed distance, varying angle
    d = 1.326  # Be-H equilibrium distance in Angstroms
    half_angle = np.radians(angle / 2.0)

    # place H atoms symmetrically about the z-axis
    h1_y = d * np.sin(half_angle)
    h1_z = d * np.cos(half_angle)

    atom_string = (
        f"Be 0.0 0.0 0.0; "
        f"H 0.0 {h1_y:.6f} {h1_z:.6f}; "
        f"H 0.0 {-h1_y:.6f} {h1_z:.6f}"
    )

    driver = PySCFDriver(
        atom=atom_string,
        basis='sto3g',
        charge=0,
        spin=0,
        unit=DistanceUnit.ANGSTROM
    )

    problem = driver.run()

    # active space: 2 electrons in 2 spatial orbitals -> 4 qubits
    from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
    transformer = ActiveSpaceTransformer(
        num_electrons=2,
        num_spatial_orbitals=3
    )
    problem = transformer.transform(problem)

    mapper = JordanWignerMapper()
    second_q_op = problem.hamiltonian.second_q_op()
    qubit_op = mapper.map(second_q_op)

    return qubit_op, problem
# BeH2 bond angle scan from 90 to 180 degrees
