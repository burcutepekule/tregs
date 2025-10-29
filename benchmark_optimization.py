"""
Benchmark script to compare v0 (baseline) vs optimized version.
Tests both performance and output equivalence.
"""

import time
import numpy as np
import pandas as pd
from main_simulation_v0 import initialize_parameters, run_simulation as run_v0
from main_simulation import run_simulation as run_optimized

def benchmark_simulations(n_runs=3, seed=42):
    """Run both versions multiple times and compare."""

    print("="*70)
    print("BENCHMARK: BASELINE (v0) vs OPTIMIZED")
    print("="*70)

    params = initialize_parameters()

    # Warm-up run (JIT compilation for numba)
    print("\nWarm-up run...")
    np.random.seed(seed)
    _ = run_v0(params, random_seed=seed)

    # Baseline version (v0)
    print(f"\n{'='*70}")
    print("BASELINE VERSION (v0)")
    print("="*70)

    v0_times = []
    for i in range(n_runs):
        np.random.seed(seed + i)
        start = time.time()
        df_v0 = run_v0(params, random_seed=seed + i)
        elapsed = time.time() - start
        v0_times.append(elapsed)
        print(f"Run {i+1}: {elapsed:.2f}s")

    v0_mean = np.mean(v0_times)
    v0_std = np.std(v0_times)

    # Optimized version
    print(f"\n{'='*70}")
    print("OPTIMIZED VERSION")
    print("="*70)

    opt_times = []
    for i in range(n_runs):
        np.random.seed(seed + i)
        start = time.time()
        df_opt = run_optimized(params, random_seed=seed + i)
        elapsed = time.time() - start
        opt_times.append(elapsed)
        print(f"Run {i+1}: {elapsed:.2f}s")

    opt_mean = np.mean(opt_times)
    opt_std = np.std(opt_times)

    # Results
    print(f"\n{'='*70}")
    print("RESULTS")
    print("="*70)
    print(f"\nBaseline (v0):  {v0_mean:.2f}s ± {v0_std:.2f}s")
    print(f"Optimized:      {opt_mean:.2f}s ± {opt_std:.2f}s")
    print(f"\nSpeedup:        {v0_mean/opt_mean:.2f}x")
    print(f"Time saved:     {v0_mean - opt_mean:.2f}s ({100*(v0_mean-opt_mean)/v0_mean:.1f}%)")

    # Compare outputs (using same seed)
    print(f"\n{'='*70}")
    print("OUTPUT COMPARISON (same seed)")
    print("="*70)

    np.random.seed(seed)
    df_v0 = run_v0(params, random_seed=seed)

    np.random.seed(seed)
    df_opt = run_optimized(params, random_seed=seed)

    # Compare key metrics at final timestep
    metrics = ['epithelial_healthy', 'pathogen', 'commensal',
               'phagocyte_M0', 'treg_resting', 'treg_active']

    all_match = True
    for metric in metrics:
        if metric in df_v0.columns and metric in df_opt.columns:
            v0_val = df_v0[metric].iloc[-1]
            opt_val = df_opt[metric].iloc[-1]
            diff = abs(v0_val - opt_val)
            match = "✓" if diff < 1 else "✗"
            if diff >= 1:
                all_match = False
            print(f"{match} {metric:20s} v0={v0_val:6.0f}  opt={opt_val:6.0f}  diff={diff:6.1f}")

    if all_match:
        print(f"\n✓ All metrics match! Optimization preserves logic.")
    else:
        print(f"\n✗ WARNING: Some metrics differ. Check for bugs in optimization.")

    return {
        'v0_mean': v0_mean,
        'opt_mean': opt_mean,
        'speedup': v0_mean / opt_mean,
        'outputs_match': all_match
    }

if __name__ == "__main__":
    results = benchmark_simulations(n_runs=3, seed=42)

    print(f"\n{'='*70}")
    print(f"Speedup achieved: {results['speedup']:.2f}x")
    print(f"Outputs match: {'Yes' if results['outputs_match'] else 'No'}")
    print("="*70)
