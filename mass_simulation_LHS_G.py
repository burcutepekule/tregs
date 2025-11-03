"""
Mass Parallel Simulation Framework - GOOGLE COLAB VERSION
Optimized for cloud execution with maximum CPU parallelization

Samples parameters using Latin Hypercube Sampling (LHS), generates scenarios,
runs multiple replicates per parameter-scenario combination using multiprocessing.

Usage:
    python mass_simulation_LHS_G.py --n_param_sets 100 --n_replicates 10 --output_dir results/
"""

import numpy as np
import pandas as pd
import os
from pathlib import Path
from multiprocessing import Pool, cpu_count
from itertools import product
import argparse
import time
from datetime import datetime
from main_simulation_G import initialize_parameters, run_simulation
import random
from scipy.stats import qmc


def sample_parameters(n_sets, random_seed=42):
    """
    Sample parameter sets using Latin Hypercube Sampling.
    Ensures good coverage of the parameter space.
    """
    np.random.seed(random_seed)

    # Define the parameter bounds
    param_bounds = [
        (0, 1),        # th_ROS_microbe
        (0, 1),        # th_ROS_epith_recover (will be adjusted)
        (0, 1),        # epith_recovery_chance
        (0.5, 1),      # rat_com_pat_threshold
        (0, 0.12),     # diffusion_speed_DAMPs
        (0, 0.12),     # diffusion_speed_SAMPs
        (0, 0.12),     # diffusion_speed_ROS
        (0, 1),        # add_ROS
        (0, 1),        # add_DAMPs
        (0, 1),        # add_SAMPs
        (0, 1),        # DAMPs_decay
        (0, 1),        # SAMPs_decay
        (0, 1),        # ros_decay
        (0, 1),        # activation_threshold_DAMPs
        (0, 1),        # activation_threshold_SAMPs
        (0, 0.5),      # activity_engulf_M0_baseline
        (0, 0.5),      # activity_engulf_M1_baseline
        (0, 0.5),      # activity_engulf_M2_baseline
        (0, 0.5),      # activity_ROS_M1_baseline
        (0.5, 1),      # rate_leak_commensal_injury
        (0.5, 1),      # rate_leak_pathogen_injury
        (0, 0.25),     # rate_leak_commensal_baseline
        (0, 1),        # treg_discrimination_efficiency
    ]

    # Create LHS sampler
    sampler = qmc.LatinHypercube(d=len(param_bounds), seed=random_seed)

    # Generate samples in [0,1]^d
    sample = sampler.random(n=n_sets)

    # Scale to actual parameter ranges
    lower_bounds = np.array([b[0] for b in param_bounds])
    upper_bounds = np.array([b[1] for b in param_bounds])
    scaled_sample = qmc.scale(sample, lower_bounds, upper_bounds)

    # Extract individual parameters
    idx = 0
    th_ROS_microbe = scaled_sample[:, idx]; idx += 1
    th_ROS_epith_recover = scaled_sample[:, idx]; idx += 1
    epith_recovery_chance = scaled_sample[:, idx]; idx += 1
    rat_com_pat_threshold = scaled_sample[:, idx]; idx += 1
    diffusion_speed_DAMPs = scaled_sample[:, idx]; idx += 1
    diffusion_speed_SAMPs = scaled_sample[:, idx]; idx += 1
    diffusion_speed_ROS = scaled_sample[:, idx]; idx += 1
    add_ROS = scaled_sample[:, idx]; idx += 1
    add_DAMPs = scaled_sample[:, idx]; idx += 1
    add_SAMPs = scaled_sample[:, idx]; idx += 1
    DAMPs_decay = scaled_sample[:, idx]; idx += 1
    SAMPs_decay = scaled_sample[:, idx]; idx += 1
    ros_decay = scaled_sample[:, idx]; idx += 1
    activation_threshold_DAMPs = scaled_sample[:, idx]; idx += 1
    activation_threshold_SAMPs = scaled_sample[:, idx]; idx += 1
    activity_engulf_M0_baseline = scaled_sample[:, idx]; idx += 1
    activity_engulf_M1_baseline = scaled_sample[:, idx]; idx += 1
    activity_engulf_M2_baseline = scaled_sample[:, idx]; idx += 1
    activity_ROS_M1_baseline = scaled_sample[:, idx]; idx += 1
    rate_leak_commensal_injury = scaled_sample[:, idx]; idx += 1
    rate_leak_pathogen_injury = scaled_sample[:, idx]; idx += 1
    rate_leak_commensal_baseline = scaled_sample[:, idx]; idx += 1
    treg_discrimination_efficiency = scaled_sample[:, idx]; idx += 1

    # Handle discrete parameter: active_age_limit
    age_sampler = qmc.LatinHypercube(d=1, seed=random_seed+1 if random_seed else None)
    active_age_sample = age_sampler.random(n=n_sets)[:, 0]
    active_age_limit = np.floor(3 + active_age_sample * 28).astype(int)

    # Create DataFrame
    params_df = pd.DataFrame({
        'param_set_id': range(n_sets),
        'th_ROS_microbe': th_ROS_microbe,
        'th_ROS_epith_recover': th_ROS_epith_recover,
        'epith_recovery_chance': epith_recovery_chance,
        'rat_com_pat_threshold': rat_com_pat_threshold,
        'diffusion_speed_DAMPs': diffusion_speed_DAMPs,
        'diffusion_speed_SAMPs': diffusion_speed_SAMPs,
        'diffusion_speed_ROS': diffusion_speed_ROS,
        'add_ROS': add_ROS,
        'add_DAMPs': add_DAMPs,
        'add_SAMPs': add_SAMPs,
        'ros_decay': ros_decay,
        'DAMPs_decay': DAMPs_decay,
        'SAMPs_decay': SAMPs_decay,
        'activation_threshold_DAMPs': activation_threshold_DAMPs,
        'activation_threshold_SAMPs': activation_threshold_SAMPs,
        'activity_engulf_M0_baseline': activity_engulf_M0_baseline,
        'activity_engulf_M1_baseline': activity_engulf_M1_baseline,
        'activity_engulf_M2_baseline': activity_engulf_M2_baseline,
        'activity_ROS_M1_baseline': activity_ROS_M1_baseline,
        'rate_leak_commensal_injury': rate_leak_commensal_injury,
        'rate_leak_pathogen_injury': rate_leak_pathogen_injury,
        'rate_leak_commensal_baseline': rate_leak_commensal_baseline,
        'active_age_limit': active_age_limit,
        'treg_discrimination_efficiency': treg_discrimination_efficiency
    })

    return params_df


def generate_scenarios():
    """
    Generate all valid scenario combinations.
    """
    scenarios = []

    for sterile in [0, 1]:
        for allow_tregs in [0, 1]:
            for allow_suppress in [0]:  # Only FALSE as per original
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


def run_mass_simulations(n_param_sets, n_replicates=10, n_cores=None,
                         output_dir='mass_sim_results', base_seed=42,
                         combine_at_end=False, verbose=True):
    """
    Main function to run mass parallel simulations optimized for Google Colab.

    Args:
        n_param_sets: Number of parameter sets to sample
        n_replicates: Number of replicate runs per parameter-scenario combo
        n_cores: Number of CPU cores to use (None = auto-detect all available)
        output_dir: Directory to save results
        base_seed: Base random seed for reproducibility
        combine_at_end: If True, combine all param_set files at the end
        verbose: Print detailed progress information
    """

    # Detect available CPUs
    if n_cores is None:
        n_cores = cpu_count()
        if verbose:
            print(f"🖥️  Auto-detected {n_cores} CPU cores available")
    else:
        if verbose:
            print(f"🖥️  Using {n_cores} CPU cores")

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Sample parameters
    if verbose:
        print(f"\n📊 Sampling {n_param_sets} parameter sets using Latin Hypercube Sampling...")
    params_df = sample_parameters(n_param_sets, random_seed=base_seed)
    params_file = output_path / 'sampled_parameters.csv'
    params_df.to_csv(params_file, index=False)
    if verbose:
        print(f"✓ Parameters saved to {params_file}")

    # Generate scenarios
    scenarios_df = generate_scenarios()
    scenarios_file = output_path / 'scenarios.csv'
    scenarios_df.to_csv(scenarios_file, index=False)
    if verbose:
        print(f"✓ {len(scenarios_df)} scenarios generated and saved to {scenarios_file}")

    # Calculate total simulations
    total_sims = n_param_sets * len(scenarios_df) * n_replicates
    sims_per_param_set = len(scenarios_df) * n_replicates

    if verbose:
        print(f"\n🎯 Total simulations to run: {total_sims:,}")
        print(f"   ({n_param_sets} param sets × {len(scenarios_df)} scenarios × {n_replicates} replicates)")
        print(f"   = {sims_per_param_set} simulations per parameter set\n")

    base_params = initialize_parameters()
    start_time = time.time()

    # Process one parameter set at a time
    completed_param_sets = 0
    total_simulations_completed = 0

    for param_idx, param_row in params_df.iterrows():
        param_set_id = int(param_row['param_set_id'])

        # Check if this param set already exists (for resume capability)
        param_set_file = output_path / f'simulation_results_param_set_{param_set_id}.csv'
        if param_set_file.exists():
            if verbose:
                print(f"⏭️  Parameter set {param_set_id}/{n_param_sets-1} already exists, skipping...")
            completed_param_sets += 1
            total_simulations_completed += sims_per_param_set
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

        if verbose:
            print(f"⚙️  Processing parameter set {param_set_id}/{n_param_sets-1} ({sims_per_param_set} simulations)...")

        with Pool(processes=n_cores) as pool:
            results_list = []
            for i, result in enumerate(pool.imap_unordered(run_single_simulation, tasks), 1):
                if result is not None:
                    results_list.append(result)

        # Save this parameter set immediately
        param_set_results = pd.concat(results_list, ignore_index=True)
        param_set_results.to_csv(param_set_file, index=False)

        completed_param_sets += 1
        total_simulations_completed += sims_per_param_set

        param_set_time = time.time() - param_set_start
        total_elapsed = time.time() - start_time

        if verbose:
            # Calculate progress statistics
            progress_pct = (completed_param_sets / n_param_sets) * 100
            avg_time_per_param_set = total_elapsed / completed_param_sets
            remaining_param_sets = n_param_sets - completed_param_sets
            estimated_remaining = avg_time_per_param_set * remaining_param_sets

            print(f"✓ Saved {param_set_file.name}")
            print(f"   Time: {param_set_time:.1f}s for this set | "
                  f"Total: {total_elapsed:.1f}s elapsed | "
                  f"ETA: {estimated_remaining:.1f}s remaining")
            print(f"   Progress: {completed_param_sets}/{n_param_sets} param sets "
                  f"({progress_pct:.1f}%) | {total_simulations_completed:,}/{total_sims:,} total sims\n")

    # Optionally combine all files
    if combine_at_end:
        if verbose:
            print("\n📦 Combining all results into single file...")
        all_files = sorted(output_path.glob('simulation_results_param_set_*.csv'))
        all_results = pd.concat([pd.read_csv(f) for f in all_files], ignore_index=True)
        combined_file = output_path / 'all_simulation_results.csv'
        all_results.to_csv(combined_file, index=False)
        if verbose:
            print(f"✓ Combined results saved to {combined_file}")
    else:
        all_results = None

    total_time = time.time() - start_time

    if verbose:
        print("\n" + "="*60)
        print("✅ MASS SIMULATION COMPLETED!")
        print("="*60)
        print(f"Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
        print(f"Total simulations: {total_sims:,}")
        print(f"Average time per simulation: {total_time/total_sims:.2f}s")
        print(f"Results saved in: {output_path.absolute()}")
        print("="*60)

    return all_results


def main():
    """Command-line interface for mass simulations."""
    parser = argparse.ArgumentParser(
        description='Run mass parallel simulations - Google Colab Optimized',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--n_param_sets', type=int, default=10,
                       help='Number of parameter sets to sample')
    parser.add_argument('--n_replicates', type=int, default=10,
                       help='Number of replicates per parameter-scenario combo')
    parser.add_argument('--n_cores', type=int, default=None,
                       help='Number of CPU cores to use (default: auto-detect all)')
    parser.add_argument('--output_dir', type=str, default='mass_sim_results',
                       help='Output directory')
    parser.add_argument('--base_seed', type=int, default=42,
                       help='Base random seed for reproducibility')
    parser.add_argument('--combine', action='store_true',
                       help='Combine all param_set files into single file at end')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress progress output')

    args = parser.parse_args()

    # Display header
    if not args.quiet:
        print("\n" + "="*60)
        print("🚀 MASS SIMULATION - GOOGLE COLAB VERSION")
        print("="*60)
        print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("="*60 + "\n")

    # Run mass simulations
    run_mass_simulations(
        n_param_sets=args.n_param_sets,
        n_replicates=args.n_replicates,
        n_cores=args.n_cores,
        output_dir=args.output_dir,
        base_seed=args.base_seed,
        combine_at_end=args.combine,
        verbose=not args.quiet
    )


if __name__ == "__main__":
    main()
