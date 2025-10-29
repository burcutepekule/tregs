"""
Simple script to run the ABM simulation.
This demonstrates how to use the simulation code and customize parameters.
"""

from main_simulation import initialize_parameters, run_simulation, plot_results
import time

def main():
    """Run the simulation with custom parameters if desired."""
    
    print("=" * 60)
    print("AGENT-BASED MODEL SIMULATION")
    print("Python implementation optimized with Numba")
    print("=" * 60)
    
    # Initialize default parameters
    params = initialize_parameters()
    
    # You can modify parameters here if needed
    # Examples:
    # params['sterile'] = 1  # Run sterile version (no pathogens)
    # params['t_max'] = 100  # Shorter simulation for testing
    # params['grid_size'] = 15  # Smaller grid for faster testing
    
    print(f"\nSimulation parameters:")
    print(f"  Grid size: {params['grid_size']}x{params['grid_size']}")
    print(f"  Time steps: {params['t_max']}")
    print(f"  Phagocytes: {params['n_phagocytes']}")
    print(f"  Tregs: {params['n_tregs']}")
    print(f"  Sterile: {'Yes' if params['sterile'] else 'No'}")
    print("-" * 60)
    
    # Start timer
    start_time = time.time()
    
    # Run the simulation with random stream file
    print("\nStarting simulation...")
    results = run_simulation(params, random_stream_file='random_numbers_seed_42.txt')
    
    # Calculate runtime
    runtime = time.time() - start_time
    
    print(f"\n✓ Simulation completed in {runtime:.2f} seconds")
    print(f"  Average time per step: {runtime/params['t_max']*1000:.2f} ms")
    
    # Save results
    output_file = 'simulation_results.csv'
    results.to_csv(output_file, index=False)
    print(f"\n✓ Results saved to: {output_file}")
    
    # Generate plots
    print("\n✓ Generating visualization plots...")
    plot_results(results)
    
    # Print final state summary
    print("\n" + "=" * 60)
    print("SIMULATION SUMMARY")
    print("=" * 60)
    
    final_state = results.iloc[-1]
    
    print(f"\nFinal Epithelial State:")
    print(f"  Healthy cells: {final_state['epithelial_healthy']:.0f}")
    for i in range(1, 6):
        if f'epithelial_inj_{i}' in final_state:
            print(f"  Injury level {i}: {final_state[f'epithelial_inj_{i}']:.0f}")
    
    print(f"\nFinal Microbe Counts:")
    print(f"  Commensals: {final_state['commensal']:.0f}")
    print(f"  Pathogens: {final_state['pathogen']:.0f}")
    
    print(f"\nFinal Phagocyte State:")
    print(f"  M0 (resting): {final_state['phagocyte_M0']:.0f}")
    
    # Sum M1 and M2 across all activation levels
    m1_total = sum(final_state[f'phagocyte_M1_L_{i}'] for i in range(params['cc_phagocyte'] + 1))
    m2_total = sum(final_state[f'phagocyte_M2_L_{i}'] for i in range(params['cc_phagocyte'] + 1))
    print(f"  M1 (pro-inflammatory): {m1_total:.0f}")
    print(f"  M2 (anti-inflammatory): {m2_total:.0f}")
    
    print(f"\nFinal Treg State:")
    print(f"  Resting: {final_state['treg_resting']:.0f}")
    print(f"  Active: {final_state['treg_active']:.0f}")
    
    print(f"\nMicrobe Deaths (cumulative):")
    print(f"  Commensals killed by ROS: {final_state['C_ROS']:.0f}")
    print(f"  Commensals killed by Macrophages: {final_state['C_Mac']:.0f}")
    print(f"  Pathogens killed by ROS: {final_state['P_ROS']:.0f}")
    print(f"  Pathogens killed by Macrophages: {final_state['P_Mac']:.0f}")
    
    print("\n" + "=" * 60)
    print("Simulation complete! Check the output files for detailed results.")
    print("=" * 60)


if __name__ == "__main__":
    main()