"""
Visualization functions for quantum chemistry PES analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def plot_pes(results, molecule_name, save_dir='results'):
    """
    Plot the potential energy surface: VQE vs exact energies.
    """
    save_path = Path(save_dir)
    save_path.mkdir(exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6))

    distances = results['distances']
    ax.plot(distances, results['exact_energies'], 'o-', color='#2C3E50',
            linewidth=2, markersize=6, label='Exact (FCI)')
    ax.plot(distances, results['vqe_energies'], 's--', color='#FF6B6B',
            linewidth=2, markersize=6, alpha=0.8, label='VQE')

    ax.set_xlabel('Distance / Angle', fontsize=12)
    ax.set_ylabel('Energy (Hartree)', fontsize=12)
    ax.set_title(f'{molecule_name} Potential Energy Surface')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    filename = f'pes_{molecule_name.lower().replace(" ", "_")}.png'
    plt.savefig(save_path / filename, dpi=150)
    plt.close()


def plot_energy_error(results, molecule_name='', save_dir='results'):
    """
    Plot the VQE energy error (in milliHartree) across the PES.
    """
    save_path = Path(save_dir)
    save_path.mkdir(exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 5))

    distances = results['distances']
    errors_mha = [e * 1000 for e in results['errors']]

    ax.plot(distances, errors_mha, 'o-', color='#E74C3C',
            linewidth=2, markersize=8)

    # chemical accuracy line (1.6 mHa)
    ax.axhline(y=1.6, color='#2ECC71', linestyle='--', linewidth=1.5,
               alpha=0.7, label='Chemical accuracy (1.6 mHa)')

    ax.set_xlabel('Distance / Angle', fontsize=12)
    ax.set_ylabel('|VQE - Exact| (mHa)', fontsize=12)
    title = 'VQE Energy Error'
    if molecule_name:
        title = f'{molecule_name}: {title}'
    ax.set_title(title)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    filename = f'energy_error_{molecule_name.lower().replace(" ", "_")}.png'
    if not molecule_name:
        filename = 'energy_error.png'
    plt.savefig(save_path / filename, dpi=150)
    plt.close()


def plot_multi_molecule_comparison(results_dict, save_dir='results'):
    """
    Compare PES curves for multiple molecules side by side.

    Args:
        results_dict: dict mapping molecule name -> scan results dict
    """
    save_path = Path(save_dir)
    save_path.mkdir(exist_ok=True)

    n_molecules = len(results_dict)
    fig, axes = plt.subplots(1, n_molecules, figsize=(6 * n_molecules, 5))

    if n_molecules == 1:
        axes = [axes]

    colors_exact = ['#2C3E50', '#1A5276', '#1B4F72']
    colors_vqe = ['#FF6B6B', '#4ECDC4', '#45B7D1']

    for i, (name, results) in enumerate(results_dict.items()):
        ax = axes[i]
        distances = results['distances']

        ax.plot(distances, results['exact_energies'], 'o-',
                color=colors_exact[i % len(colors_exact)],
                linewidth=2, markersize=5, label='Exact (FCI)')
        ax.plot(distances, results['vqe_energies'], 's--',
                color=colors_vqe[i % len(colors_vqe)],
                linewidth=2, markersize=5, alpha=0.8, label='VQE')

        ax.set_xlabel('Distance / Angle', fontsize=11)
        ax.set_ylabel('Energy (Ha)', fontsize=11)
        ax.set_title(name, fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.suptitle('Potential Energy Surfaces: VQE vs Exact', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(save_path / 'multi_molecule_comparison.png', dpi=150,
                bbox_inches='tight')
    plt.close()


def plot_ansatz_comparison(results_by_ansatz, molecule_name, save_dir='results'):
    """
    Compare different ansatz types on the same PES.

    Args:
        results_by_ansatz: dict mapping ansatz name -> scan results
        molecule_name: name for the title
    """
    save_path = Path(save_dir)
    save_path.mkdir(exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    colors = {'RealAmplitudes': '#FF6B6B', 'EfficientSU2': '#4ECDC4'}

    # PES comparison
    ax = axes[0]
    # plot exact once (same for all ansatze)
    first_key = list(results_by_ansatz.keys())[0]
    distances = results_by_ansatz[first_key]['distances']
    ax.plot(distances, results_by_ansatz[first_key]['exact_energies'],
            'o-', color='#2C3E50', linewidth=2, markersize=5,
            label='Exact (FCI)')

    for name, results in results_by_ansatz.items():
        ax.plot(distances, results['vqe_energies'], 's--',
                color=colors.get(name, 'gray'), linewidth=2,
                markersize=5, alpha=0.8, label=f'VQE ({name})')

    ax.set_xlabel('Distance (A)', fontsize=11)
    ax.set_ylabel('Energy (Ha)', fontsize=11)
    ax.set_title(f'{molecule_name}: PES by Ansatz')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # error comparison
    ax = axes[1]
    for name, results in results_by_ansatz.items():
        errors_mha = [e * 1000 for e in results['errors']]
        ax.plot(distances, errors_mha, 'o-',
                color=colors.get(name, 'gray'), linewidth=2,
                markersize=6, label=name)

    ax.axhline(y=1.6, color='#2ECC71', linestyle='--', linewidth=1.5,
               alpha=0.7, label='Chemical accuracy')
    ax.set_xlabel('Distance (A)', fontsize=11)
    ax.set_ylabel('|VQE - Exact| (mHa)', fontsize=11)
    ax.set_title(f'{molecule_name}: VQE Error by Ansatz')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    filename = f'ansatz_comparison_{molecule_name.lower()}.png'
    plt.savefig(save_path / filename, dpi=150)
    plt.close()
