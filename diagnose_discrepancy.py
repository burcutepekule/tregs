"""
Diagnostic script to investigate discrepancy between expected and actual outcome proportions
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the balanced parameters file
balanced_params = pd.read_csv('balanced_lhs_parameters.csv')

print("\n=== BALANCED PARAMETERS INFO ===")
print(f"Total parameter sets: {len(balanced_params)}")
print(f"\nColumns in balanced_lhs_parameters.csv:")
print(balanced_params.columns.tolist())

# Check if target_category exists
if 'target_category' in balanced_params.columns:
    print("\n=== EXPECTED PROPORTIONS (from R script balancing) ===")
    expected = balanced_params['target_category'].value_counts()
    expected_pct = 100 * balanced_params['target_category'].value_counts(normalize=True)
    print(expected)
    print("\nPercentages:")
    print(expected_pct)
else:
    print("\nWARNING: 'target_category' not found in balanced_lhs_parameters.csv")
    print("This column should have been created by A006_decision_tree_a_bias_class.R")

# Check if we have the actual results
try:
    import pickle
    # Try to load the RDS file would require rpy2, so let's check if CSV results exist
    import glob

    result_files = glob.glob('mass_sim_results_presampled/simulation_results_param_set_*.csv')
    print(f"\n=== SIMULATION RESULTS ===")
    print(f"Found {len(result_files)} simulation result files")

    if len(result_files) > 0:
        print(f"\nProcessing sample of results to check format...")
        sample_df = pd.read_csv(result_files[0])
        print(f"Columns in simulation results:")
        print(sample_df.columns.tolist())
        print(f"\nFirst few rows:")
        print(sample_df.head())

except Exception as e:
    print(f"\nError loading simulation results: {e}")

# Additional diagnostics
print("\n=== PARAMETER BOUNDS ANALYSIS ===")

param_names = [
    "th_ROS_microbe", "th_ROS_epith_recover", "epith_recovery_chance",
    "rat_com_pat_threshold", "diffusion_speed_DAMPs", "diffusion_speed_SAMPs",
    "diffusion_speed_ROS", "add_ROS", "add_DAMPs", "add_SAMPs",
    "ros_decay", "DAMPs_decay", "SAMPs_decay",
    "activation_threshold_DAMPs", "activation_threshold_SAMPs",
    "activity_engulf_M0_baseline", "activity_engulf_M1_baseline",
    "activity_engulf_M2_baseline", "activity_ROS_M1_baseline",
    "rate_leak_commensal_injury", "rate_leak_pathogen_injury",
    "rate_leak_commensal_baseline", "active_age_limit",
    "treg_discrimination_efficiency"
]

# Check parameter ranges by target category if available
if 'target_category' in balanced_params.columns:
    print("\nParameter ranges by target category:")
    for category in balanced_params['target_category'].unique():
        cat_data = balanced_params[balanced_params['target_category'] == category]
        print(f"\n{category.upper()} (n={len(cat_data)}):")

        # Show a few key parameters
        key_params = ['th_ROS_microbe', 'diffusion_speed_DAMPs', 'treg_discrimination_efficiency']
        for param in key_params:
            if param in balanced_params.columns:
                print(f"  {param}: [{cat_data[param].min():.3f}, {cat_data[param].max():.3f}], "
                      f"mean={cat_data[param].mean():.3f}")

# Check for source_node if available
if 'source_node' in balanced_params.columns:
    print("\n=== NODE CONTRIBUTION ANALYSIS ===")
    node_summary = balanced_params.groupby('source_node').agg({
        'param_set_id': 'count',
        'target_category': lambda x: x.value_counts().to_dict() if len(x) > 0 else {}
    }).rename(columns={'param_set_id': 'n_samples'})

    print("\nSamples per node and their target categories:")
    for node_id in sorted(balanced_params['source_node'].unique()):
        node_data = balanced_params[balanced_params['source_node'] == node_id]
        cat_dist = node_data['target_category'].value_counts()
        print(f"\nNode {node_id} (n={len(node_data)}):")
        for cat, count in cat_dist.items():
            print(f"  {cat}: {count} ({100*count/len(node_data):.1f}%)")

# Create visualization if possible
if 'target_category' in balanced_params.columns:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Expected proportions
    ax1 = axes[0]
    expected_counts = balanced_params['target_category'].value_counts()
    colors = {'better': 'forestgreen', 'worse': 'firebrick', 'drift': 'gray'}
    bar_colors = [colors.get(cat, 'blue') for cat in expected_counts.index]
    ax1.bar(expected_counts.index, expected_counts.values, color=bar_colors)
    ax1.set_title('Expected Proportions\n(from balanced_lhs_parameters.csv)')
    ax1.set_ylabel('Count')

    # Add percentages on bars
    for i, (cat, count) in enumerate(expected_counts.items()):
        pct = 100 * count / expected_counts.sum()
        ax1.text(i, count, f'{pct:.1f}%', ha='center', va='bottom')

    # Source node distribution
    if 'source_node' in balanced_params.columns:
        ax2 = axes[1]
        node_data = balanced_params.groupby(['source_node', 'target_category']).size().unstack(fill_value=0)
        node_data.plot(kind='bar', stacked=True, ax=ax2, color=[colors.get(c, 'blue') for c in node_data.columns])
        ax2.set_title('Samples by Node and Category')
        ax2.set_ylabel('Count')
        ax2.set_xlabel('Node ID')
        ax2.legend(title='Target Category')
        plt.setp(ax2.xaxis.get_majorticklabels(), rotation=0)

    plt.tight_layout()
    plt.savefig('expected_proportions.png', dpi=150, bbox_inches='tight')
    print("\n=== Saved visualization to expected_proportions.png ===")
    plt.close()

print("\n=== DIAGNOSTIC COMPLETE ===")
print("\nTo fully diagnose the issue, you need to:")
print("1. Run A008_analyse_py_read_raw_presampled.R to process simulation results")
print("2. Run A009_analyse_py_selection_score_cts_presampled.R to get actual proportions")
print("3. Compare the actual proportions with the expected ones shown above")
