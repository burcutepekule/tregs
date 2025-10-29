"""
Analysis utilities for mass simulation results.

Functions to:
- Load and merge results with parameters
- Calculate summary statistics per parameter-scenario combo
- Extract key metrics (pathogen clearance time, final states, etc.)
- Export for downstream analysis in R or Python
"""

import pandas as pd
import numpy as np
from pathlib import Path


def load_mass_sim_results(results_dir='mass_sim_results'):
    """
    Load all mass simulation results and merge with parameters/scenarios.

    Returns:
        tuple: (results_df, params_df, scenarios_df)
    """
    results_dir = Path(results_dir)

    print(f"Loading results from {results_dir}...")

    results = pd.read_csv(results_dir / 'all_simulation_results.csv')
    params = pd.read_csv(results_dir / 'sampled_parameters.csv')
    scenarios = pd.read_csv(results_dir / 'scenarios.csv')

    print(f"  ✓ Loaded {len(results):,} result rows")
    print(f"  ✓ Loaded {len(params)} parameter sets")
    print(f"  ✓ Loaded {len(scenarios)} scenarios")

    return results, params, scenarios


def calculate_summary_statistics(results):
    """
    Calculate summary statistics per parameter-scenario-replicate combination.

    Returns:
        DataFrame with one row per simulation run
    """
    print("\nCalculating summary statistics...")

    summaries = []

    # Group by unique simulation
    grouped = results.groupby(['param_set_id', 'scenario_id', 'replicate_id'])

    for (param_id, scenario_id, rep_id), sim_data in grouped:
        sim_data = sim_data.sort_values('t')

        summary = {
            'param_set_id': param_id,
            'scenario_id': scenario_id,
            'replicate_id': rep_id,

            # Pathogen dynamics
            'pathogen_clearance_time': get_clearance_time(sim_data, 'pathogen'),
            'pathogen_peak': sim_data['pathogen'].max(),
            'pathogen_final': sim_data['pathogen'].iloc[-1],

            # Commensal dynamics
            'commensal_final': sim_data['commensal'].iloc[-1],
            'commensal_mean': sim_data['commensal'].mean(),
            'commensal_min': sim_data['commensal'].min(),

            # Epithelial health
            'epithelial_healthy_final': sim_data['epithelial_healthy'].iloc[-1],
            'epithelial_healthy_min': sim_data['epithelial_healthy'].min(),
            'max_injury_reached': get_max_injury(sim_data),

            # Phagocyte phenotypes
            'M0_final': sim_data['phagocyte_M0'].iloc[-1],
            'M1_peak': get_M1_total(sim_data).max(),
            'M1_final': get_M1_total(sim_data).iloc[-1],
            'M2_peak': get_M2_total(sim_data).max(),
            'M2_final': get_M2_total(sim_data).iloc[-1],

            # Treg activity
            'treg_active_final': sim_data['treg_active'].iloc[-1],
            'treg_active_mean': sim_data['treg_active'].mean(),
            'treg_active_peak': sim_data['treg_active'].max(),

            # Kill counts
            'pathogens_killed_by_ROS': sim_data['P_ROS'].iloc[-1],
            'pathogens_killed_by_Mac': sim_data['P_Mac'].iloc[-1],
            'commensals_killed_by_ROS': sim_data['C_ROS'].iloc[-1],
            'commensals_killed_by_Mac': sim_data['C_Mac'].iloc[-1],

            # Overall health score (0=bad, 1=good)
            'health_score': calculate_health_score(sim_data)
        }

        summaries.append(summary)

    summary_df = pd.DataFrame(summaries)
    print(f"  ✓ Created {len(summary_df)} summary rows")

    return summary_df


def get_clearance_time(sim_data, agent='pathogen'):
    """Get time when agent count drops to zero (or -1 if never cleared)."""
    zero_times = sim_data[sim_data[agent] == 0]['t']
    if len(zero_times) > 0:
        return zero_times.iloc[0]
    return -1


def get_M1_total(sim_data):
    """Sum all M1 columns."""
    m1_cols = [col for col in sim_data.columns if col.startswith('phagocyte_M1')]
    return sim_data[m1_cols].sum(axis=1)


def get_M2_total(sim_data):
    """Sum all M2 columns."""
    m2_cols = [col for col in sim_data.columns if col.startswith('phagocyte_M2')]
    return sim_data[m2_cols].sum(axis=1)


def get_max_injury(sim_data):
    """Get maximum injury level reached across all epithelial cells."""
    injury_cols = [col for col in sim_data.columns if col.startswith('epithelial_inj')]
    if len(injury_cols) > 0:
        # Check if any cell reached level 5
        if 'epithelial_inj_5' in sim_data.columns and sim_data['epithelial_inj_5'].max() > 0:
            return 5
        # Check level 4, 3, 2, 1
        for level in [4, 3, 2, 1]:
            col = f'epithelial_inj_{level}'
            if col in sim_data.columns and sim_data[col].max() > 0:
                return level
    return 0


def calculate_health_score(sim_data):
    """
    Calculate overall health score (0=bad, 1=good).
    Based on pathogen clearance, epithelial recovery, and commensal preservation.
    """
    # Pathogen cleared?
    pathogen_cleared = sim_data['pathogen'].iloc[-1] == 0

    # Epithelial health recovered?
    epithelial_recovered = sim_data['epithelial_healthy'].iloc[-1] > sim_data['epithelial_healthy'].iloc[0] * 0.8

    # Commensals preserved?
    commensals_preserved = sim_data['commensal'].iloc[-1] > 5

    # Combined score
    score = (0.5 * pathogen_cleared +
             0.3 * epithelial_recovered +
             0.2 * commensals_preserved)

    return score


def merge_with_parameters(summary_df, params_df, scenarios_df):
    """Merge summary statistics with parameters and scenarios."""
    print("\nMerging with parameters and scenarios...")

    # Merge with parameters
    merged = summary_df.merge(params_df, on='param_set_id', how='left')

    # Merge with scenarios
    merged = merged.merge(scenarios_df, on='scenario_id', how='left')

    print(f"  ✓ Merged dataset has {len(merged)} rows and {len(merged.columns)} columns")

    return merged


def aggregate_across_replicates(merged_df):
    """
    Aggregate metrics across replicates to get mean ± std per parameter-scenario combo.
    """
    print("\nAggregating across replicates...")

    # Columns to aggregate
    metric_cols = [
        'pathogen_clearance_time', 'pathogen_peak', 'pathogen_final',
        'commensal_final', 'commensal_mean',
        'epithelial_healthy_final', 'epithelial_healthy_min',
        'M1_peak', 'M2_peak',
        'treg_active_final', 'treg_active_mean',
        'health_score'
    ]

    # Group by param_set_id and scenario_id
    grouped = merged_df.groupby(['param_set_id', 'scenario_id'])

    agg_results = []
    for (param_id, scenario_id), group in grouped:
        row = {
            'param_set_id': param_id,
            'scenario_id': scenario_id,
            'n_replicates': len(group)
        }

        # Calculate mean and std for each metric
        for metric in metric_cols:
            if metric in group.columns:
                row[f'{metric}_mean'] = group[metric].mean()
                row[f'{metric}_std'] = group[metric].std()

        agg_results.append(row)

    agg_df = pd.DataFrame(agg_results)

    # Merge back with parameter values and scenario settings
    params_df = merged_df[['param_set_id'] + [col for col in merged_df.columns
                           if col in merged_df.columns and col not in metric_cols
                           and col not in ['replicate_id', 'scenario_id']]].drop_duplicates('param_set_id')

    scenarios_df = merged_df[['scenario_id', 'sterile', 'allow_tregs_to_do_their_job',
                              'allow_tregs_to_suppress_cognate', 'randomize_tregs']].drop_duplicates('scenario_id')

    agg_df = agg_df.merge(params_df, on='param_set_id', how='left')
    agg_df = agg_df.merge(scenarios_df, on='scenario_id', how='left')

    print(f"  ✓ Aggregated to {len(agg_df)} rows (one per param-scenario combo)")

    return agg_df


def export_for_analysis(merged_df, agg_df, output_dir='mass_sim_results'):
    """Export processed data for downstream analysis."""
    output_dir = Path(output_dir)

    print(f"\nExporting processed data to {output_dir}...")

    # Full detailed results (all replicates)
    merged_file = output_dir / 'analysis_full_with_params.csv'
    merged_df.to_csv(merged_file, index=False)
    print(f"  ✓ Saved full results: {merged_file}")

    # Aggregated results (mean ± std across replicates)
    agg_file = output_dir / 'analysis_aggregated.csv'
    agg_df.to_csv(agg_file, index=False)
    print(f"  ✓ Saved aggregated results: {agg_file}")

    print("\n✓ Analysis complete! Ready for downstream processing.")


def main():
    """Run full analysis pipeline."""
    import argparse

    parser = argparse.ArgumentParser(description='Analyze mass simulation results')
    parser.add_argument('--results_dir', type=str, default='mass_sim_results',
                       help='Directory with simulation results')
    args = parser.parse_args()

    # Load results
    results, params, scenarios = load_mass_sim_results(args.results_dir)

    # Calculate summaries
    summary_df = calculate_summary_statistics(results)

    # Merge with parameters
    merged_df = merge_with_parameters(summary_df, params, scenarios)

    # Aggregate across replicates
    agg_df = aggregate_across_replicates(merged_df)

    # Export
    export_for_analysis(merged_df, agg_df, args.results_dir)


if __name__ == "__main__":
    main()
