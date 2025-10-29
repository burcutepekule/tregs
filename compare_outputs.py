"""
Diagnostic script to compare R and Python simulation outputs.
Helps identify where differences occur.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_and_compare(r_file='simulation_results_R.csv', py_file='simulation_results.csv'):
    """Load R and Python outputs and compare them."""

    print("="*70)
    print("COMPARING R AND PYTHON SIMULATION OUTPUTS")
    print("="*70)

    try:
        df_r = pd.read_csv(r_file)
        print(f"\n✓ Loaded R results: {r_file} ({len(df_r)} rows)")
    except FileNotFoundError:
        print(f"\n✗ R results file not found: {r_file}")
        return

    try:
        df_py = pd.read_csv(py_file)
        print(f"✓ Loaded Python results: {py_file} ({len(df_py)} rows)")
    except FileNotFoundError:
        print(f"\n✗ Python results file not found: {py_file}")
        return

    # Check dimensions
    print(f"\nDimensions: R={df_r.shape}, Python={df_py.shape}")

    # Find common columns
    common_cols = set(df_r.columns) & set(df_py.columns)
    print(f"\nCommon columns: {len(common_cols)}")

    # Compare key metrics
    key_metrics = [
        'epithelial_healthy', 'pathogen', 'commensal',
        'phagocyte_M0', 'treg_resting', 'treg_active'
    ]

    print("\n" + "="*70)
    print("COMPARING KEY METRICS AT FINAL TIMEPOINT")
    print("="*70)

    final_r = df_r.iloc[-1]
    final_py = df_py.iloc[-1]

    for metric in key_metrics:
        if metric in common_cols:
            r_val = final_r[metric]
            py_val = final_py[metric]
            diff = abs(r_val - py_val)
            pct_diff = 100 * diff / (r_val + 1e-10)
            match = "✓" if diff < 1 else "✗"
            print(f"{match} {metric:25s} R={r_val:8.1f}  Py={py_val:8.1f}  Diff={diff:8.1f} ({pct_diff:6.2f}%)")

    # Find first timestep where differences appear
    print("\n" + "="*70)
    print("FINDING FIRST DIVERGENCE POINT")
    print("="*70)

    for t in range(min(len(df_r), len(df_py))):
        for metric in key_metrics[:3]:  # Check first 3 metrics
            if metric in common_cols:
                if abs(df_r.iloc[t][metric] - df_py.iloc[t][metric]) > 0.1:
                    print(f"\nFirst significant difference at t={t} in '{metric}':")
                    print(f"  R:      {df_r.iloc[t][metric]}")
                    print(f"  Python: {df_py.iloc[t][metric]}")
                    print(f"  Diff:   {abs(df_r.iloc[t][metric] - df_py.iloc[t][metric])}")

                    # Show context
                    if t > 0:
                        print(f"\nPrevious timestep (t={t-1}):")
                        print(f"  R:      {df_r.iloc[t-1][metric]}")
                        print(f"  Python: {df_py.iloc[t-1][metric]}")

                    return t  # Return first divergence point

    print("\n✓ No significant differences found!")
    return None

def plot_comparison(r_file='simulation_results_R.csv', py_file='simulation_results.csv'):
    """Create side-by-side comparison plots."""

    try:
        df_r = pd.read_csv(r_file)
        df_py = pd.read_csv(py_file)
    except FileNotFoundError as e:
        print(f"Cannot create plots: {e}")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    metrics = [
        ('pathogen', 'Pathogen Count', 'red'),
        ('commensal', 'Commensal Count', 'cyan'),
        ('epithelial_healthy', 'Healthy Epithelial Cells', 'blue'),
        ('treg_active', 'Active Tregs', 'purple')
    ]

    for ax, (metric, title, color) in zip(axes.flat, metrics):
        if metric in df_r.columns and metric in df_py.columns:
            ax.plot(df_r['t'], df_r[metric], label='R', color=color, linewidth=2, alpha=0.7)
            ax.plot(df_py['t'], df_py[metric], label='Python', color=color, linewidth=2, alpha=0.7, linestyle='--')
            ax.set_xlabel('Time')
            ax.set_ylabel('Count')
            ax.set_title(title)
            ax.legend()
            ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('comparison_R_vs_Python.png', dpi=150)
    print(f"\n✓ Comparison plot saved: comparison_R_vs_Python.png")

if __name__ == "__main__":
    # Compare outputs
    divergence_point = load_and_compare()

    # Create comparison plot
    plot_comparison()

    if divergence_point is not None:
        print(f"\n{'='*70}")
        print("DIAGNOSIS")
        print("="*70)
        print(f"Trajectories diverge at timestep t={divergence_point}")
        print("\nPossible causes:")
        print("  1. Random number consumption order differs between R and Python")
        print("  2. Numerical precision differences in ppf functions")
        print("  3. Different iteration order in loops")
        print("  4. Floating point arithmetic differences")
        print("\nRecommendation:")
        print("  - Check if statistical properties (mean, variance) match over many runs")
        print("  - Consider accepting statistical equivalence rather than exact match")
    else:
        print(f"\n{'='*70}")
        print("SUCCESS!")
        print("="*70)
        print("R and Python outputs match closely!")
