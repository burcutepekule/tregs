"""
Quick test of mass simulation framework.
Runs a small batch to verify everything works before full run.
"""

from mass_simulation import run_mass_simulations

def main():
    # Test with small numbers
    print("Running test with 2 parameter sets, 2 replicates, 2 cores...")
    print("This should complete in ~1-2 minutes\n")

    results = run_mass_simulations(
        n_param_sets=2,       # Just 2 parameter sets
        n_replicates=2,       # Just 2 replicates
        n_cores=2,            # Just 2 cores
        output_dir='test_mass_sim_results',
        base_seed=42
    )

    print("\n" + "="*70)
    print("TEST RESULTS")
    print("="*70)
    print(f"Total rows: {len(results)}")
    print(f"Columns: {list(results.columns)}")
    print(f"\nUnique param_set_ids: {results['param_set_id'].nunique()}")
    print(f"Unique scenario_ids: {results['scenario_id'].nunique()}")
    print(f"Unique replicate_ids: {results['replicate_id'].nunique()}")
    print(f"\nFirst few rows:")
    print(results[['param_set_id', 'scenario_id', 'replicate_id', 't', 'pathogen', 'commensal']].head(10))

    print("\n✓ Test successful! Ready for full run.")

if __name__ == '__main__':
    main()
