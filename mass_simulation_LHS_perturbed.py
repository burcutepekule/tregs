"""
Mass Parallel Simulation Framework - Pre-sampled Parameters

Uses pre-sampled parameters from balanced_lhs_parameters_PERTURBED.csv,
generates scenarios, runs 100 replicates per parameter-scenario combination
using multiprocessing.

Usage:
    python mass_simulation_LHS_perturbed.py --csv_file balanced_lhs_parameters_PERTURBED.csv --n_cores 10 --output_dir results/
"""

import numpy as np
import pandas as pd
import os
from pathlib import Path
from multiprocessing import Pool, cpu_count
from itertools import product
import argparse
import time
from main_simulation import initialize_parameters, run_simulation
import random


def load_parameters(csv_file):
    """
    Load pre-sampled parameter sets from CSV file.

    Args:
        csv_file: Path to CSV file with parameters

    Returns:
        DataFrame with parameter sets
    """
    # Read the CSV file
    params_df = pd.read_csv(csv_file)

    # Drop the first unnamed column if it exists (row numbers from R)
    if params_df.columns[0] == '' or 'Unnamed: 0' in params_df.columns[0]:
        params_df = params_df.drop(params_df.columns[0], axis=1)

    # Verify required columns exist
    required_cols = [
        'param_set_id',
        'th_ROS_microbe',
        'th_ROS_epith_recover',
        'epith_recovery_chance',
        'rat_com_pat_threshold',
        'diffusion_speed_DAMPs',
        'diffusion_speed_SAMPs',
        'diffusion_speed_ROS',
        'add_ROS',
        'add_DAMPs',
        'add_SAMPs',
        'ros_decay',
        'DAMPs_decay',
        'SAMPs_decay',
        'activation_threshold_DAMPs',
        'activation_threshold_SAMPs',
        'activity_engulf_M0_baseline',
        'activity_engulf_M1_baseline',
        'activity_engulf_M2_baseline',
        'activity_ROS_M1_baseline',
        'rate_leak_commensal_injury',
        'rate_leak_pathogen_injury',
        'rate_leak_commensal_baseline',
        'active_age_limit',
        'treg_discrimination_efficiency'
    ]

    missing_cols = [col for col in required_cols if col not in params_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in CSV: {missing_cols}")

    print(f"Loaded {len(params_df)} parameter sets from {csv_file}")

    return params_df


def generate_scenarios():
    """
    Generate all valid scenario combinations.
    Matches the R expand.grid + filtering logic.
    """
    scenarios = []

    for sterile in [0, 1]:
        for allow_tregs in [0, 1]:  # Note: R uses FALSE/TRUE, we use 0/1
            for allow_suppress in [0]:  # Only FALSE as per R code
                for randomize in [0, 1]:
                    # Filter: doesn't make sense to randomize if tregs don't work
                    if not (allow_tregs == 0 and randomize == 1):
                        scenarios.append({
                            'scenario_id': len(scenarios),
                            'sterile': sterile,
                            'allow_tregs_to_do_their_job': allow_tregs,
                            'allow_tregs_to_suppress_cognate': allow_suppress,
                            'randomize_tregs': randomize
                        })

    return pd.DataFrame(scenarios)


def update_params_with_sampled_and_scenario(base_params, sampled_row, scenario_row):
    """Update base parameters with sampled values and scenario settings."""
    params = base_params.copy()

    # Update with sampled parameters
    for key in sampled_row.index:
        if key != 'param_set_id' and key in params:
            params[key] = sampled_row[key]

    # Update with scenario settings
    params['sterile'] = int(scenario_row['sterile'])
    params['allow_tregs_to_do_their_job'] = int(scenario_row['allow_tregs_to_do_their_job'])
    params['allow_tregs_to_suppress_cognate'] = bool(scenario_row['allow_tregs_to_suppress_cognate'])
    params['randomize_tregs'] = int(scenario_row['randomize_tregs'])

    return params


def run_single_simulation(args):
    """
    Run a single simulation and return results.

    Args:
        args: tuple of (param_set_id, scenario_id, replicate_id, params, seed)

    Returns:
        DataFrame with results + identifiers
    """
    param_set_id, scenario_id, replicate_id, params, seed = args

    try:
        # Run simulation
        results = run_simulation(params, random_seed=seed)

        # Add identifiers
        results['param_set_id'] = param_set_id
        results['scenario_id'] = scenario_id
        results['replicate_id'] = replicate_id
        results['seed'] = seed

        return results

    except Exception as e:
        print(f"ERROR in param_set={param_set_id}, scenario={scenario_id}, replicate={replicate_id}: {e}")
        return None


def run_mass_simulations(csv_file, n_replicates=10, n_cores=10,
                         output_dir='mass_sim_results_perturbed', base_seed=42,
                         combine_at_end=False):
    """
    Main function to run mass parallel simulations with pre-sampled parameters.
    Saves results incrementally per parameter set for safety.

    Args:
        csv_file: Path to CSV file with pre-sampled parameters
        n_replicates: Number of replicate runs per parameter-scenario combo
        n_cores: Number of CPU cores to use
        output_dir: Directory to save results
        base_seed: Base random seed for reproducibility
        combine_at_end: If True, combine all param_set files at the end
    """

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load parameters from CSV
    params_df = load_parameters(csv_file)

    # Save a copy of the parameters to the output directory
    params_file = output_path / 'loaded_parameters.csv'
    params_df.to_csv(params_file, index=False)

    # Generate scenarios
    scenarios_df = generate_scenarios()
    scenarios_file = output_path / 'scenarios.csv'
    scenarios_df.to_csv(scenarios_file, index=False)

    # Calculate total simulations
    n_param_sets = len(params_df)
    total_sims = n_param_sets * len(scenarios_df) * n_replicates
    sims_per_param_set = len(scenarios_df) * n_replicates

    print(f"\n{'='*60}")
    print(f"Starting mass simulations with pre-sampled parameters")
    print(f"{'='*60}")
    print(f"Parameter sets: {n_param_sets}")
    print(f"Scenarios: {len(scenarios_df)}")
    print(f"Replicates per combo: {n_replicates}")
    print(f"Total simulations: {total_sims}")
    print(f"CPU cores: {n_cores}")
    print(f"Output directory: {output_dir}")
    print(f"{'='*60}\n")

    base_params = initialize_parameters()
    start_time = time.time()

    # Process one parameter set at a time
    completed_param_sets = 0

    for _, param_row in params_df.iterrows():
        param_set_id = int(param_row['param_set_id'])

        # Check if this param set already exists (for resume capability)
        param_set_file = output_path / f'simulation_results_param_set_{param_set_id}.csv'
        if param_set_file.exists():
            print(f"⏭️  Parameter set {param_set_id} already exists, skipping...")
            completed_param_sets += 1
            continue

        # Build tasks for this parameter set only
        tasks = []
        for _, scenario_row in scenarios_df.iterrows():
            scenario_id = int(scenario_row['scenario_id'])

            # Update parameters
            params = update_params_with_sampled_and_scenario(base_params, param_row, scenario_row)

            # Generate replicates with different seeds
            for replicate_id in range(n_replicates):
                seed = random.randint(0, 2**32 - 1)
                tasks.append((param_set_id, scenario_id, replicate_id, params, seed))

        # Run this parameter set in parallel
        param_set_start = time.time()

        with Pool(processes=n_cores) as pool:
            results_list = []
            for i, result in enumerate(pool.imap_unordered(run_single_simulation, tasks), 1):
                if result is not None:
                    results_list.append(result)

        # Save this parameter set immediately
        print(f"Running parameter set {param_set_id}...")
        param_set_results = pd.concat(results_list, ignore_index=True)
        param_set_results.to_csv(param_set_file, index=False)

        completed_param_sets += 1
        param_set_time = time.time() - param_set_start
        elapsed_time = time.time() - start_time
        avg_time_per_set = elapsed_time / completed_param_sets
        remaining_sets = n_param_sets - completed_param_sets
        estimated_remaining = avg_time_per_set * remaining_sets

        print(f"✓ Saved {param_set_file.name}")
        print(f"  Time for this set: {param_set_time:.1f}s")
        print(f"  Progress: {completed_param_sets}/{n_param_sets} sets ({100*completed_param_sets/n_param_sets:.1f}%)")
        print(f"  Estimated time remaining: {estimated_remaining/60:.1f} minutes\n")

    # Final summary
    total_time = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"All simulations completed!")
    print(f"Total time: {total_time/60:.1f} minutes")
    print(f"Average time per parameter set: {total_time/n_param_sets:.1f} seconds")
    print(f"{'='*60}\n")

    # Optionally combine all files
    if combine_at_end:
        print("Combining all results into single file...")
        all_files = sorted(output_path.glob('simulation_results_param_set_*.csv'))
        all_results = pd.concat([pd.read_csv(f) for f in all_files], ignore_index=True)
        combined_file = output_path / 'all_simulation_results.csv'
        all_results.to_csv(combined_file, index=False)
        print(f"✓ Combined results saved to {combined_file}")
    else:
        all_results = None

    return all_results


def main():
    """Command-line interface for mass simulations with pre-sampled parameters."""
    parser = argparse.ArgumentParser(
        description='Run mass parallel simulations with pre-sampled parameters'
    )
    parser.add_argument('--csv_file', type=str,
                       default='balanced_lhs_parameters_PERTURBED.csv',
                       help='CSV file with pre-sampled parameters (default: balanced_lhs_parameters_PERTURBED.csv)')
    parser.add_argument('--n_replicates', type=int, default=10,
                       help='Number of replicates per parameter-scenario combo (default: 10)')
    parser.add_argument('--n_cores', type=int, default=min(10, cpu_count()),
                       help='Number of CPU cores to use (default: 10 or max available)')
    parser.add_argument('--output_dir', type=str, default='mass_sim_results_perturbed',
                       help='Output directory (default: mass_sim_results_perturbed)')
    parser.add_argument('--base_seed', type=int, default=42,
                       help='Base random seed (default: 42)')
    parser.add_argument('--combine', action='store_true',
                       help='Combine all param_set files into single file at end (requires memory)')

    args = parser.parse_args()

    # Verify CSV file exists
    if not Path(args.csv_file).exists():
        print(f"ERROR: CSV file not found: {args.csv_file}")
        return

    # Set random seed
    random.seed(args.base_seed)
    np.random.seed(args.base_seed)

    # Run mass simulations
    run_mass_simulations(
        csv_file=args.csv_file,
        n_replicates=args.n_replicates,
        n_cores=args.n_cores,
        output_dir=args.output_dir,
        base_seed=args.base_seed,
        combine_at_end=args.combine
    )


if __name__ == "__main__":
    main()
