import os
import numpy as np
import pickle
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Structure import Structure
import mdtraj as md


def load_alphafold_results(result_dir, model_num=1):
    """Load AlphaFold results from pickle file."""
    pkl_file = os.path.join(result_dir, f"result_model_{model_num}.pkl")
    with open(pkl_file, "rb") as f:
        results = pickle.load(f)
    return results


def compare_structures(pred_pdb, exp_pdb, output_prefix):
    """Compare predicted structure to experimental structure."""
    # Load structures with mdtraj for easier analysis
    pred = md.load(pred_pdb)
    exp = md.load(exp_pdb)

    # Calculate RMSD
    # First we need to align the structures
    pred.superpose(
        exp,
        atom_indices=pred.topology.select("name CA"),
        ref_atom_indices=exp.topology.select("name CA"),
    )
    rmsd = md.rmsd(
        pred,
        exp,
        atom_indices=pred.topology.select("name CA"),
        ref_atom_indices=exp.topology.select("name CA"),
    )[0]

    # Plot per-residue RMSD
    pred_xyz = pred.xyz[0, pred.topology.select("name CA")]
    exp_xyz = exp.xyz[0, exp.topology.select("name CA")]
    per_res_rmsd = np.sqrt(np.sum((pred_xyz - exp_xyz) ** 2, axis=1))

    plt.figure(figsize=(10, 6))
    plt.plot(per_res_rmsd)
    plt.xlabel("Residue")
    plt.ylabel("RMSD (nm)")
    plt.title(
        f"Per-residue RMSD: {os.path.basename(pred_pdb)} vs {os.path.basename(exp_pdb)}"
    )
    plt.savefig(f"{output_prefix}_per_res_rmsd.png")

    # Calculate lDDT (simplified version)
    # In a real implementation, you'd use a proper lDDT calculator

    return {"rmsd": rmsd, "per_res_rmsd": per_res_rmsd}


def plot_plddt_comparison(default_results, custom_results, output_prefix):
    """Plot pLDDT comparison between default and custom MSA runs."""
    default_plddt = default_results["plddt"]
    custom_plddt = custom_results["plddt"]

    plt.figure(figsize=(12, 6))
    plt.plot(default_plddt, label="Default MSA")
    plt.plot(custom_plddt, label="Custom MSA")
    plt.xlabel("Residue")
    plt.ylabel("pLDDT")
    plt.title("pLDDT Comparison: Default vs Custom MSA")
    plt.legend()
    plt.savefig(f"{output_prefix}_plddt_comparison.png")

    # Calculate improvement statistics
    improvement = custom_plddt - default_plddt
    avg_improvement = np.mean(improvement)
    pct_improved = np.sum(improvement > 0) / len(improvement) * 100

    return {
        "avg_improvement": avg_improvement,
        "pct_improved": pct_improved,
        "default_avg_plddt": np.mean(default_plddt),
        "custom_avg_plddt": np.mean(custom_plddt),
    }


def plot_pae_comparison(default_results, custom_results, output_prefix):
    """Plot PAE comparison between default and custom MSA runs."""
    if (
        "predicted_aligned_error" not in default_results
        or "predicted_aligned_error" not in custom_results
    ):
        return None

    default_pae = default_results["predicted_aligned_error"]
    custom_pae = custom_results["predicted_aligned_error"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    im1 = ax1.imshow(default_pae, cmap="Blues_r", vmin=0, vmax=30)
    ax1.set_title("Default MSA PAE")
    fig.colorbar(im1, ax=ax1)

    im2 = ax2.imshow(custom_pae, cmap="Blues_r", vmin=0, vmax=30)
    ax2.set_title("Custom MSA PAE")
    fig.colorbar(im2, ax=ax2)

    plt.savefig(f"{output_prefix}_pae_comparison.png")

    # Calculate improvement statistics
    improvement = default_pae - custom_pae
    avg_improvement = np.mean(improvement)

    return {
        "avg_pae_improvement": avg_improvement,
        "default_avg_pae": np.mean(default_pae),
        "custom_avg_pae": np.mean(custom_pae),
    }


def analyze_protein(protein_name, default_dir, custom_dir, exp_pdb, output_dir):
    """Analyze a single protein with all metrics."""
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = os.path.join(output_dir, protein_name)

    # Load results
    default_results = load_alphafold_results(default_dir)
    custom_results = load_alphafold_results(custom_dir)

    # Compare confidence metrics
    plddt_stats = plot_plddt_comparison(default_results, custom_results, output_prefix)
    pae_stats = plot_pae_comparison(default_results, custom_results, output_prefix)

    # Compare to experimental structure
    default_pdb = os.path.join(default_dir, "ranked_0.pdb")
    custom_pdb = os.path.join(custom_dir, "ranked_0.pdb")

    default_vs_exp = compare_structures(
        default_pdb, exp_pdb, f"{output_prefix}_default_vs_exp"
    )
    custom_vs_exp = compare_structures(
        custom_pdb, exp_pdb, f"{output_prefix}_custom_vs_exp"
    )

    # Create summary report
    with open(f"{output_prefix}_summary.txt", "w") as f:
        f.write(f"Analysis for {protein_name}\n")
        f.write("=" * 50 + "\n\n")

        f.write("Confidence Metrics:\n")
        f.write(f"Default MSA average pLDDT: {plddt_stats['default_avg_plddt']:.2f}\n")
        f.write(f"Custom MSA average pLDDT: {plddt_stats['custom_avg_plddt']:.2f}\n")
        f.write(f"Average pLDDT improvement: {plddt_stats['avg_improvement']:.2f}\n")
        f.write(
            f"Percentage of residues with improved pLDDT: {plddt_stats['pct_improved']:.2f}%\n\n"
        )

        if pae_stats:
            f.write(f"Default MSA average PAE: {pae_stats['default_avg_pae']:.2f}\n")
            f.write(f"Custom MSA average PAE: {pae_stats['custom_avg_pae']:.2f}\n")
            f.write(
                f"Average PAE improvement: {pae_stats['avg_pae_improvement']:.2f}\n\n"
            )

        f.write("Comparison to Experimental Structure:\n")
        f.write(f"Default MSA model RMSD to exp: {default_vs_exp['rmsd']:.2f} nm\n")
        f.write(f"Custom MSA model RMSD to exp: {custom_vs_exp['rmsd']:.2f} nm\n")
        f.write(
            f"RMSD improvement: {default_vs_exp['rmsd'] - custom_vs_exp['rmsd']:.2f} nm\n\n"
        )

        # Add TM-score and GDT if available (would require additional tools)

    return {
        "plddt_stats": plddt_stats,
        "pae_stats": pae_stats,
        "default_vs_exp": default_vs_exp,
        "custom_vs_exp": custom_vs_exp,
    }


# Example usage
proteins = [
    {"name": "ERb1", "uniprot": "A0A8M2B359", "pdb": "path/to/erb1_exp.pdb"},
    {"name": "TMTC4", "uniprot": "Q9VF81", "pdb": "path/to/tmtc4_exp.pdb"},
    {"name": "T1044", "uniprot": "T1044", "pdb": "path/to/t1044_exp.pdb"},
    {"name": "T1064", "uniprot": "T1064", "pdb": "path/to/t1064_exp.pdb"},
    {"name": "T1050", "uniprot": "T1050", "pdb": "path/to/t1050_exp.pdb"},
]

results = {}
for protein in proteins:
    default_dir = f"path/to/{protein['name']}_default/alphafold_output"
    custom_dir = f"path/to/{protein['name']}_custom/alphafold_output"
    results[protein["name"]] = analyze_protein(
        protein["name"], default_dir, custom_dir, protein["pdb"], "comparison_results"
    )

# Generate overall comparison
# [...code to create overall comparison table...]
