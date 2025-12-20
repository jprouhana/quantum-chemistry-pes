# Quantum Chemistry: Potential Energy Surfaces

Computing potential energy surfaces (PES) for small molecules (H2, LiH, BeH2) using the Variational Quantum Eigensolver (VQE) with different ansatz types. Comparison to exact full configuration interaction (FCI) energies from classical diagonalization.

Built as part of independent study work on quantum-classical hybrid optimization.

## Background

### What Is a Potential Energy Surface?

A **potential energy surface** (PES) maps how a molecule's total electronic energy changes as its geometry varies — for example, as you stretch a chemical bond. The minimum of the PES gives the equilibrium bond length, and the curvature determines vibrational frequencies.

For a diatomic molecule like H2, the PES is a 1D curve: energy vs bond distance. For triatomics like BeH2, you can scan bond angles while keeping distances fixed.

### VQE for Molecular Ground States

The **Variational Quantum Eigensolver** (VQE) finds the ground state energy of a molecular Hamiltonian using a parameterized quantum circuit (ansatz). The workflow:

1. **Build the Hamiltonian**: Use PySCF to compute molecular integrals, then map the fermionic Hamiltonian to qubits via Jordan-Wigner transformation
2. **Prepare a trial state**: The ansatz circuit $U(\theta)|0\rangle$ prepares a parameterized trial wavefunction
3. **Measure the energy**: $E(\theta) = \langle 0|U^\dagger(\theta) H U(\theta)|0\rangle$
4. **Optimize**: A classical optimizer adjusts $\theta$ to minimize $E(\theta)$

The variational principle guarantees $E(\theta) \geq E_0$ (the true ground state energy), so VQE provides an upper bound that improves as the ansatz becomes more expressive.

### Molecules Studied

- **H2**: Simplest molecule, 4 qubits (STO-3G), well-understood benchmark
- **LiH**: Heteronuclear, more complex electronic structure, 12 spin-orbitals
- **BeH2**: Triatomic, bond angle scan, tests multi-reference character

## Project Structure

```
quantum-chemistry-pes/
├── src/
│   ├── molecules.py         # Molecular Hamiltonian construction
│   ├── vqe_runner.py        # VQE execution with different ansatze
│   ├── exact_solver.py      # Classical FCI reference energies
│   ├── pes_scan.py          # PES scanning over geometries
│   └── plotting.py          # Visualization functions
├── notebooks/
│   └── pes_analysis.ipynb
├── results/
├── requirements.txt
├── README.md
└── LICENSE
```

## Installation

```bash
git clone https://github.com/jrouhana/quantum-chemistry-pes.git
cd quantum-chemistry-pes
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Usage

### Quick Start

```python
from src.molecules import build_h2
from src.vqe_runner import run_vqe
from src.exact_solver import exact_energy

hamiltonian, problem = build_h2(distance=0.735)
vqe_result = run_vqe(hamiltonian, ansatz_type='RealAmplitudes', seed=42)
fci_energy = exact_energy(hamiltonian)

print(f"VQE energy: {vqe_result['energy']:.6f} Ha")
print(f"FCI energy: {fci_energy:.6f} Ha")
print(f"Error: {abs(vqe_result['energy'] - fci_energy)*1000:.2f} mHa")
```

### Jupyter Notebook

```bash
jupyter notebook notebooks/pes_analysis.ipynb
```

## Results

### H2 Potential Energy Surface

| Distance (A) | VQE Energy (Ha) | Exact Energy (Ha) | Error (mHa) |
|---------------|-----------------|-------------------|-------------|
| 0.5 | -1.0554 | -1.0554 | 0.01 |
| 0.735 | -1.1373 | -1.1373 | 0.02 |
| 1.0 | -1.1018 | -1.1018 | 0.03 |
| 1.5 | -0.9917 | -0.9917 | 0.05 |
| 2.0 | -0.9486 | -0.9486 | 0.12 |

*RealAmplitudes ansatz, COBYLA optimizer, 200 iterations.*

### Key Findings

- **VQE achieves chemical accuracy** (< 1.6 mHa error) for H2 across the full PES, including the dissociation region where classical methods like Hartree-Fock fail
- **EfficientSU2 ansatz is more expressive** but requires more optimization iterations to converge than RealAmplitudes for these small molecules
- **LiH is harder**: The larger Hilbert space and stronger correlations near dissociation lead to slightly larger VQE errors
- **BeH2 angle scan** reveals the expected linear equilibrium geometry, with VQE tracking the exact PES well

## References

1. Peruzzo, A., et al. (2014). "A variational eigenvalue solver on a photonic quantum processor." *Nature Communications*, 5, 4213.
2. Kandala, A., et al. (2017). "Hardware-efficient variational quantum eigensolver for small molecules and quantum magnets." *Nature*, 549, 242-246.
3. McArdle, S., et al. (2020). "Quantum computational chemistry." *Reviews of Modern Physics*, 92, 015003.

## License

MIT License — see [LICENSE](LICENSE) for details.
