"""
Mass Parallel Simulation Framework

Samples parameters uniformly, generates scenarios, runs 100 replicates
per parameter-scenario combination using multiprocessing.

Usage:
    python mass_simulation.py --n_param_sets 100 --n_cores 10 --output_dir results/
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

def sample_parameters(n_sets, random_seed=42):
    """
    Sample parameter sets uniformly.
    Matches the R sampling scheme.
    """
    np.random.seed(random_seed)

    # Sample all parameters
    th_ROS_microbe = np.random.uniform(0, 0.5, n_sets)
    th_ROS_epith_recover = np.random.uniform(th_ROS_microbe, 1.0)  # Must be > th_ROS_microbe
    epith_recovery_chance = np.random.uniform(0, 1, n_sets)
    rat_com_pat_threshold = np.random.uniform(0.5, 1, n_sets)  # At least above half

    # Diffusion speeds (MUST be < 0.125 for stability!)
    diffusion_speed_DAMPs = np.random.uniform(0, 0.12, n_sets)
    diffusion_speed_SAMPs = np.random.uniform(0, 0.12, n_sets)
    diffusion_speed_ROS = np.random.uniform(0, 0.12, n_sets)

    # Signal production rates
    add_ROS = np.random.uniform(0, 1, n_sets)
    add_DAMPs = np.random.uniform(0, 1, n_sets)
    add_SAMPs = np.random.uniform(0, 1, n_sets)

    # Decay rates
    DAMPs_decay = np.random.uniform(0, 1, n_sets)
    SAMPs_decay = np.random.uniform(0, 1, n_sets)
    ros_decay = np.random.uniform(0, 1, n_sets)

    # Activation thresholds
    activation_threshold_DAMPs = np.random.uniform(0, 1, n_sets)
    activation_threshold_SAMPs = np.random.uniform(0, 1, n_sets)

    # Engulfing and ROS rates - baseline is max 50% of capacity
    activity_engulf_M0_baseline = np.random.uniform(0, 0.5, n_sets)
    activity_engulf_M1_baseline = np.random.uniform(0, 0.5, n_sets)
    activity_engulf_M2_baseline = np.random.uniform(0, 0.5, n_sets)
    activity_ROS_M1_baseline = np.random.uniform(0, 0.5, n_sets)

    # Leakage rates
    rate_leak_commensal_injury = np.random.uniform(0.5, 1, n_sets)
    rate_leak_pathogen_injury = np.random.uniform(0.5, 1, n_sets)
    rate_leak_commensal_baseline = np.random.uniform(0.0, 0.25, n_sets)

    # Active age limit (discrete values 3-30)
    active_age_limit = np.random.randint(3, 31, n_sets)

    # Treg discrimination efficiency
    treg_discrimination_efficiency = np.random.uniform(0, 1, n_sets)

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


def run_mass_simulations(n_param_sets, n_replicates=100, n_cores=10,
                         output_dir='mass_sim_results', base_seed=42):
    """
    Main function to run mass parallel simulations.

    Args:
        n_param_sets: Number of parameter sets to sample
        n_replicates: Number of replicate runs per parameter-scenario combo
        n_cores: Number of CPU cores to use
        output_dir: Directory to save results
        base_seed: Base random seed for reproducibility
    """

    print("="*70)
    print("MASS SIMULATION FRAMEWORK")
    print("="*70)

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Sample parameters
    print(f"\n1. Sampling {n_param_sets} parameter sets...")
    params_df = sample_parameters(n_param_sets, random_seed=base_seed)
    params_file = output_path / 'sampled_parameters.csv'
    params_df.to_csv(params_file, index=False)
    print(f"   ✓ Saved to {params_file}")

    # Generate scenarios
    print(f"\n2. Generating scenarios...")
    scenarios_df = generate_scenarios()
    scenarios_file = output_path / 'scenarios.csv'
    scenarios_df.to_csv(scenarios_file, index=False)
    print(f"   ✓ Generated {len(scenarios_df)} scenarios")
    print(f"   ✓ Saved to {scenarios_file}")

    # Calculate total simulations
    total_sims = n_param_sets * len(scenarios_df) * n_replicates
    print(f"\n3. Setting up {total_sims:,} simulations:")
    print(f"   - {n_param_sets} parameter sets")
    print(f"   - {len(scenarios_df)} scenarios")
    print(f"   - {n_replicates} replicates each")
    print(f"   - Using {n_cores} CPU cores")

    # Build task list
    print(f"\n4. Building task list...")
    base_params = initialize_parameters()
    tasks = []

    for _, param_row in params_df.iterrows():
        param_set_id = int(param_row['param_set_id'])

        for _, scenario_row in scenarios_df.iterrows():
            scenario_id = int(scenario_row['scenario_id'])

            # Update parameters
            params = update_params_with_sampled_and_scenario(base_params, param_row, scenario_row)

            # Generate replicates with different seeds
            for replicate_id in range(n_replicates):
                # Create unique seed: base_seed + param_id*1000000 + scenario_id*1000 + replicate_id
                seed = base_seed + param_set_id * 1000000 + scenario_id * 1000 + replicate_id

                tasks.append((param_set_id, scenario_id, replicate_id, params, seed))

    print(f"   ✓ Created {len(tasks):,} tasks")

    # Run simulations in parallel
    print(f"\n5. Running simulations on {n_cores} cores...")
    print(f"   (This will take a while...)")

    start_time = time.time()

    with Pool(processes=n_cores) as pool:
        # Use imap for progress tracking
        results_list = []
        for i, result in enumerate(pool.imap_unordered(run_single_simulation, tasks), 1):
            if result is not None:
                results_list.append(result)

            # Progress update every 10 simulations
            if i % 10 == 0 or i == len(tasks):
                elapsed = time.time() - start_time
                rate = i / elapsed
                remaining = (len(tasks) - i) / rate if rate > 0 else 0
                print(f"   Progress: {i}/{len(tasks)} ({100*i/len(tasks):.1f}%) | "
                      f"Elapsed: {elapsed/60:.1f}min | "
                      f"ETA: {remaining/60:.1f}min | "
                      f"Rate: {rate:.1f} sim/s")

    elapsed_time = time.time() - start_time

    print(f"\n6. Combining results...")
    all_results = pd.concat(results_list, ignore_index=True)

    # Save combined results
    results_file = output_path / 'all_simulation_results.csv'
    all_results.to_csv(results_file, index=False)
    print(f"   ✓ Saved {len(all_results)} rows to {results_file}")

    # Summary
    print(f"\n{'='*70}")
    print("SIMULATION COMPLETE")
    print("="*70)
    print(f"Total time: {elapsed_time/60:.1f} minutes ({elapsed_time/3600:.2f} hours)")
    print(f"Average time per simulation: {elapsed_time/len(tasks):.2f} seconds")
    print(f"Throughput: {len(tasks)/elapsed_time:.2f} simulations/second")
    print(f"\nResults saved to: {output_dir}/")
    print(f"  - sampled_parameters.csv")
    print(f"  - scenarios.csv")
    print(f"  - all_simulation_results.csv")

    return all_results


def main():
    """Command-line interface for mass simulations."""
    parser = argparse.ArgumentParser(description='Run mass parallel simulations')
    parser.add_argument('--n_param_sets', type=int, default=10,
                       help='Number of parameter sets to sample (default: 10)')
    parser.add_argument('--n_replicates', type=int, default=100,
                       help='Number of replicates per parameter-scenario combo (default: 100)')
    parser.add_argument('--n_cores', type=int, default=min(10, cpu_count()),
                       help='Number of CPU cores to use (default: 10 or max available)')
    parser.add_argument('--output_dir', type=str, default='mass_sim_results',
                       help='Output directory (default: mass_sim_results)')
    parser.add_argument('--base_seed', type=int, default=42,
                       help='Base random seed (default: 42)')

    args = parser.parse_args()

    # Run mass simulations
    run_mass_simulations(
        n_param_sets=args.n_param_sets,
        n_replicates=args.n_replicates,
        n_cores=args.n_cores,
        output_dir=args.output_dir,
        base_seed=args.base_seed
    )


if __name__ == "__main__":
    main()
