"""
Main Agent-Based Model Simulation
Optimized Python implementation with Numba for maximum performance.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from time import time
from simulation_utils import *
import numba

# Set random seed for reproducibility at module level
np.random.seed(42)

def initialize_parameters():
    """Initialize all simulation parameters."""
    params = {}
    
    # Grid parameters
    params['grid_size'] = 25
    params['t_max'] = 500
    
    # Agent counts
    params['n_phagocytes'] = round(params['grid_size'] * params['grid_size'] * 0.35)
    params['n_tregs'] = round(params['grid_size'] * params['grid_size'] * 0.35)
    
    # Initial microbe counts
    params['n_commensals_lp'] = 20
    
    # Injury parameters
    params['injury_percentage'] = 60
    params['max_level_injury'] = 5
    
    # Diffusion parameters (must be < 0.125 for stability)
    params['diffusion_speed_DAMPs'] = 0.07209796
    params['diffusion_speed_SAMPs'] = 0.08999438
    params['diffusion_speed_ROS'] = 0.05232508
    
    # Decay rates
    params['DAMPs_decay'] = 0.6555613
    params['SAMPs_decay'] = 0.2551983
    params['ros_decay'] = 0.02983629
    
    # Signal production rates
    params['add_DAMPs'] = 0.4806108
    params['add_SAMPs'] = 0.8993819
    params['add_ROS'] = 0.8972294
    
    # Activation thresholds
    params['activation_threshold_DAMPs'] = 0.03546454
    params['activation_threshold_SAMPs'] = 0.7854582
    
    # Killing thresholds
    params['th_ROS_microbe'] = 0.2691068
    params['th_ROS_epith_recover'] = 0.8875888
    params['epith_recovery_chance'] = 0.473025
    
    # Engulfment activities (baseline rates)
    params['activity_engulf_M0_baseline'] = 0.07882888
    params['activity_engulf_M1_baseline'] = 0.1559077
    params['activity_engulf_M2_baseline'] = 0.2466745
    params['activity_engulf_max'] = 0.99
    
    # ROS production activities
    params['activity_ROS_M0_baseline'] = 0.00
    params['activity_ROS_M1_baseline'] = 0.4768265
    params['activity_ROS_M2_baseline'] = 0.00
    params['activity_ROS_max'] = 0.99
    
    # Treg parameters
    params['treg_vicinity_effect'] = 1  # if 0, that means has to be at the very same pixel
    params['treg_discrimination_efficiency'] = 0.03508874
    params['rat_com_pat_threshold'] = 0.7703095
    params['allow_tregs_to_do_their_job'] = 1
    params['allow_tregs_to_suppress_cognate'] = False
    params['randomize_tregs'] = 1
    
    # Phagocyte parameters
    params['cc_phagocyte'] = 5
    params['active_age_limit'] = 15
    params['digestion_time'] = 1
    
    # Infection parameters
    params['rate_leak_pathogen_injury'] = 0.7240775
    params['rate_leak_commensal_injury'] = 0.9778595
    params['rate_leak_commensal_baseline'] = 0.2162822
    params['sterile'] = 0
    
    # Logistic function parameters
    params['k_in'] = 0.044
    params['x0_in'] = 50
    
    # Action radii
    params['act_radius_ROS'] = 1
    params['act_radius_DAMPs'] = 1
    params['act_radius_SAMPs'] = 1
    
    # Max cell values
    params['max_cell_value_ROS'] = 1
    params['max_cell_value_DAMPs'] = 1
    params['max_cell_value_SAMPs'] = 1
    
    return params


def run_simulation(params, random_seed=42):
    """Run the main ABM simulation."""
    
    # Set random seed for reproducibility
    np.random.seed(random_seed)
    
    # Extract frequently used parameters
    grid_size = params['grid_size']
    t_max = params['t_max']
    n_phagocytes = params['n_phagocytes']
    n_tregs = params['n_tregs']
    
    # Calculate derived parameters
    injury_site = get_middle_percent(np.arange(grid_size), params['injury_percentage'])
    precision = 10 * np.exp(5 * params['treg_discrimination_efficiency'])
    
    if params['sterile'] == 1:
        params['rate_leak_pathogen_injury'] = 0.0
    
    n_pathogens_lp = round(params['rate_leak_pathogen_injury'] * len(injury_site))
    
    # Initialize signal fields
    DAMPs = np.zeros((grid_size, grid_size), dtype=np.float32)
    SAMPs = np.zeros((grid_size, grid_size), dtype=np.float32)
    ROS   = np.zeros((grid_size, grid_size), dtype=np.float32)
    
    # Calculate activity steps
    activity_engulf_M1_step = (params['activity_engulf_max'] - params['activity_engulf_M1_baseline']) / params['cc_phagocyte']
    activity_engulf_M2_step = (params['activity_engulf_max'] - params['activity_engulf_M2_baseline']) / params['cc_phagocyte']
    activity_ROS_M1_step = (params['activity_ROS_max'] - params['activity_ROS_M1_baseline']) / params['cc_phagocyte']
    
    # Initialize epithelium
    epithelium_injury = np.zeros(grid_size, dtype=np.int32)
    epithelium_injury[injury_site] = 1
    
    # Initialize phagocytes
    phagocyte_x = np.random.randint(0, grid_size, n_phagocytes)
    phagocyte_y = np.random.randint(1, grid_size, n_phagocytes)  # Not on epithelium row
    phagocyte_pathogens_engulfed = np.zeros(n_phagocytes, dtype=np.int32)
    phagocyte_commensals_engulfed = np.zeros(n_phagocytes, dtype=np.int32)
    phagocyte_num_times_activated = np.zeros(n_phagocytes, dtype=np.int32)
    phagocyte_phenotype = np.zeros(n_phagocytes, dtype=np.int32)  # 0=M0, 1=M1, 2=M2
    phagocyte_active_age = np.zeros(n_phagocytes, dtype=np.int32)
    phagocyte_digestion_counter = np.zeros(n_phagocytes, dtype=np.int32)
    phagocyte_activity_engulf = np.full(n_phagocytes, params['activity_engulf_M0_baseline'], dtype=np.float32)
    phagocyte_activity_ROS = np.full(n_phagocytes, params['activity_ROS_M0_baseline'], dtype=np.float32)
    phagocyte_bacteria_registry = np.zeros((n_phagocytes, params['cc_phagocyte']), dtype=np.int32)
    
    # Initialize Tregs
    treg_x = np.random.randint(0, grid_size, n_tregs)
    treg_y = np.random.randint(1, grid_size, n_tregs)  # Not on epithelium row
    
    treg_phenotype = np.zeros(n_tregs, dtype=np.int32)  # 0=resting, 1=active
    treg_active_age = np.zeros(n_tregs, dtype=np.int32)
    treg_activity_SAMPs_binary = np.zeros(n_tregs, dtype=np.int32)
    
    # Initialize bacteria coordinates
    # Pathogens start at y=1 (one row above epithelium which is at y=0)
    pathogen_x = np.random.choice(injury_site, n_pathogens_lp, replace=True) if n_pathogens_lp > 0 else np.array([])
    pathogen_y = np.ones(n_pathogens_lp, dtype=np.int32)  # y=1, not y=0!
    pathogen_coords = np.column_stack([pathogen_x, pathogen_y]) if n_pathogens_lp > 0 else np.empty((0, 2), dtype=np.int32)
    
    # Commensals start at random positions throughout the grid
    commensal_x = np.random.randint(0, grid_size, params['n_commensals_lp'])
    commensal_y = np.random.randint(0, grid_size, params['n_commensals_lp'])  # Random y position
    commensal_coords = np.column_stack([commensal_x, commensal_y])
    
    # Initialize tracking variables
    pathogens_killed_by_ROS = 0
    pathogens_killed_by_Mac = 0
    commensals_killed_by_ROS = 0
    commensals_killed_by_Mac = 0
    
    # Initialize data storage matrices
    epithelium_longitudinal = np.zeros((t_max, 6), dtype=np.int32)
    macrophages_longitudinal = np.zeros((t_max, 1 + 2 * (params['cc_phagocyte'] + 1)), dtype=np.int32)
    microbes_longitudinal = np.zeros((t_max, 2), dtype=np.int32)
    tregs_longitudinal = np.zeros((t_max, 2), dtype=np.int32)
    microbes_cumdeath_longitudinal = np.zeros((t_max, 8), dtype=np.int32)
    
    print(f"Starting simulation with {n_phagocytes} phagocytes and {n_tregs} Tregs...")
    print(f"Initial bacteria: {params['n_commensals_lp']} commensals, {n_pathogens_lp} pathogens")
    
    # Main simulation loop
    for t in range(t_max):
        if t % 100 == 0:
            print(f"Time step {t}/{t_max}")
        
        # Update SAMPs based on activated tregs
        for i in range(len(treg_x)):
            if treg_activity_SAMPs_binary[i] == 1:
             SAMPs[treg_y[i], treg_x[i]] += params['add_SAMPs']
             
        # Update ROS based on M1 phagocytes  
        for i in range(len(phagocyte_x)):
            if phagocyte_phenotype[i] == 1:  # M1
              ROS[phagocyte_y[i], phagocyte_x[i]] += phagocyte_activity_ROS[i]*params['add_ROS']
              
        #   # Move pathogens and commensals randomly

        if len(pathogen_coords) > 0:
            # When at y=0 (touching epithelium), can only move to LP (dy=1), not lumen
            dy = np.where(pathogen_coords[:, 1] == 0,
                          np.ones(len(pathogen_coords), dtype=np.int32),  # Must move up
                          np.random.choice([-1, 0, 1], size=len(pathogen_coords), replace=True))  # Random movement
            dx = iszero_coordinates(pathogen_coords[:, 0])
            pathogen_coords[:, 0] = np.clip(pathogen_coords[:, 0] + dx, 0, grid_size - 1)
            pathogen_coords[:, 1] = np.clip(pathogen_coords[:, 1] + dy, 0, grid_size - 1)
            
        if len(commensal_coords) > 0:
            # When at y=0 (touching epithelium), can only move uLP (dy=1), not lumen
            dy = np.where(commensal_coords[:, 1] == 0,
                          np.ones(len(commensal_coords), dtype=np.int32),  # Must move up
                          np.random.choice([-1, 0, 1], size=len(commensal_coords), replace=True))  # Random movement
            dx = iszero_coordinates(commensal_coords[:, 0])
            commensal_coords[:, 0] = np.clip(commensal_coords[:, 0] + dx, 0, grid_size - 1)
            commensal_coords[:, 1] = np.clip(commensal_coords[:, 1] + dy, 0, grid_size - 1)
        
        # Add DAMPs from bacteria at epithelium - this is calculated separately in R but summed in the end
        bacteria_at_epithelium = np.zeros(grid_size)
        if len(pathogen_coords) > 0:
            path_at_epi = pathogen_coords[pathogen_coords[:, 1] == 0]
            if len(path_at_epi) > 0:
                unique, counts = np.unique(path_at_epi[:, 0], return_counts=True)
                bacteria_at_epithelium[unique] += counts
        if len(commensal_coords) > 0:
            comm_at_epi = commensal_coords[commensal_coords[:, 1] == 0]
            if len(comm_at_epi) > 0:
                unique, counts = np.unique(comm_at_epi[:, 0], return_counts=True)
                bacteria_at_epithelium[unique] += counts
        
        # Add DAMPs from bacteria counts at epithelium
        for i in range(grid_size):
            DAMPs[0, i] += epithelium_injury[i] * params['add_DAMPs']
            if bacteria_at_epithelium[i] > 0:
                DAMPs[0, i] += logistic_scaled_0_to_5_quantized(bacteria_at_epithelium[i], params['k_in'], params['x0_in']) * params['add_DAMPs']

        # Diffuse signals BEFORE decay (like R code)
        DAMPs = diffuse_matrix(DAMPs, params['diffusion_speed_DAMPs'], params['max_cell_value_DAMPs'])
        SAMPs = diffuse_matrix(SAMPs, params['diffusion_speed_SAMPs'], params['max_cell_value_SAMPs'])
        ROS = diffuse_matrix(ROS, params['diffusion_speed_ROS'], params['max_cell_value_ROS'])
        
        # Decay AFTER diffusion (like R code)
        DAMPs = DAMPs * (1 - params['DAMPs_decay'])
        SAMPs = SAMPs * (1 - params['SAMPs_decay'])
        ROS = ROS * (1 - params['ros_decay'])
        
        # Create density matrix for Treg movement
        if params['randomize_tregs'] == 0:
            density_matrix_tregs = DAMPs  # Tregs follow DAMPs gradient
        else:
            density_matrix_tregs = np.zeros_like(DAMPs)  # All zeros = random movement
            
        density_matrix_phagocytes = DAMPs
        
        # Check if density matrix has variation
        if not np.all(density_matrix_tregs == density_matrix_tregs[0, 0]):
            # Move based on gradient (chemotaxis)
            for i in range(n_tregs):
                x, y = treg_x[i], treg_y[i]
                
                # Get 3x3 neighborhood
                x_start = max(0, x - 1)
                x_end = min(grid_size, x + 2)
                y_start = max(0, y - 1)  
                y_end = min(grid_size, y + 2)
                
                # Get density values in neighborhood
                neighbors_x = []
                neighbors_y = []
                neighbor_density = []
                
                for ny in range(y_start, y_end):
                    for nx in range(x_start, x_end):
                        neighbors_x.append(nx)
                        neighbors_y.append(ny)
                        neighbor_density.append(density_matrix_tregs[ny, nx])
                
                # Normalize to get probabilities
                neighbor_density = np.array(neighbor_density)
                if neighbor_density.sum() > 0:
                    probs = neighbor_density / neighbor_density.sum()
                else:
                    probs = np.ones(len(neighbor_density)) / len(neighbor_density)
                
                # Sample weighted by density concentration
                chosen_idx = np.random.choice(len(neighbors_x), p=probs)
                treg_x[i] = neighbors_x[chosen_idx]
                treg_y[i] = neighbors_y[chosen_idx]
        else:
            # Random movement if no gradient
            dy = np.where(treg_y == 0,
                 np.ones(len(treg_y), dtype=np.int32),  # Must move up
                 np.random.choice([-1, 0, 1], size=len(treg_y), replace=True))
                        
            dx = iszero_coordinates(dy)
            treg_x = np.clip(treg_x + dx, 0, grid_size - 1)
            treg_y = np.clip(treg_y + dy, 0, grid_size - 1)
        
        # Phagocyte movement (chemotaxis towards DAMPs if present)
        # Check if DAMPs field has variation
        if not np.all(density_matrix_phagocytes == density_matrix_phagocytes[0, 0]):
            # Move based on DAMP gradient (chemotaxis)
            for i in range(n_phagocytes):
                x, y = phagocyte_x[i], phagocyte_y[i]
                
                # Get 3x3 neighborhood
                x_start = max(0, x - 1)
                x_end = min(grid_size, x + 2)
                y_start = max(0, y - 1) 
                y_end = min(grid_size, y + 2)
                
                # Get DAMP values in neighborhood
                neighbors_x = []
                neighbors_y = []
                neighbor_damps = []
                
                for ny in range(y_start, y_end):
                    for nx in range(x_start, x_end):
                        neighbors_x.append(nx)
                        neighbors_y.append(ny)
                        neighbor_damps.append(DAMPs[ny, nx])
                
                # Normalize to get probabilities
                neighbor_damps = np.array(neighbor_damps)
                if neighbor_damps.sum() > 0:
                    probs = neighbor_damps / neighbor_damps.sum()
                else:
                    probs = np.ones(len(neighbor_damps)) / len(neighbor_damps)
                
                # Sample weighted by DAMP concentration
                chosen_idx = np.random.choice(len(neighbors_x), p=probs)
                phagocyte_x[i] = neighbors_x[chosen_idx]
                phagocyte_y[i] = neighbors_y[chosen_idx]
        else:
            # Random movement if no gradient
            dy = np.where(phagocyte_y == 0,
                 np.ones(len(phagocyte_y), dtype=np.int32),  # Must move up
                 np.random.choice([-1, 0, 1], size=len(phagocyte_y), replace=True))
                        
            dx = iszero_coordinates(dy)
            phagocyte_x = np.clip(phagocyte_x + dx, 0, grid_size - 1)
            phagocyte_y = np.clip(phagocyte_y + dy, 0, grid_size - 1) 

        # Leak bacteria based on injury and baseline rates
        # Calculate number of new bacteria to add (like R code)
        injured_sites = np.where(epithelium_injury > 0)[0]
        
        # Pathogens leak based on mean injury level
        if len(injured_sites) > 0:
            mean_injury = np.mean(epithelium_injury)
            n_pathogens_new = int(np.round(mean_injury * params['rate_leak_pathogen_injury'] * len(injured_sites)))
            
            if n_pathogens_new > 0:
                # Sample locations weighted by injury level
                probs = epithelium_injury / epithelium_injury.sum() if epithelium_injury.sum() > 0 else None
                new_x = np.random.choice(grid_size, n_pathogens_new, replace=True, p=probs)
                new_pathogens = np.column_stack([new_x, np.ones(n_pathogens_new, dtype=np.int32)])
                pathogen_coords = np.vstack([pathogen_coords, new_pathogens]) if len(pathogen_coords) > 0 else new_pathogens
        
        # Commensals leak from both injury and baseline
        mean_injury = np.mean(epithelium_injury)
        n_commensals_injury = int(np.round(mean_injury * params['rate_leak_commensal_injury'] * len(injured_sites)))
        n_commensals_baseline = int(np.round(params['rate_leak_commensal_baseline'] * grid_size))
        total_new_commensals = n_commensals_injury + n_commensals_baseline
        
        if total_new_commensals > 0:
            new_commensals_x = []
            
            # Baseline commensals - random locations
            if n_commensals_baseline > 0:
                baseline_x = np.random.randint(0, grid_size, n_commensals_baseline)
                new_commensals_x.extend(baseline_x)
            
            # Injury-site commensals - weighted by injury
            if n_commensals_injury > 0:
                probs = epithelium_injury / epithelium_injury.sum() if epithelium_injury.sum() > 0 else None
                injury_x = np.random.choice(grid_size, n_commensals_injury, replace=True, p=probs)
                new_commensals_x.extend(injury_x)
            
            if len(new_commensals_x) > 0:
                new_commensals = np.column_stack([
                    np.array(new_commensals_x),
                    np.ones(len(new_commensals_x), dtype=np.int32)
                ])
                commensal_coords = np.vstack([commensal_coords, new_commensals]) if len(commensal_coords) > 0 else new_commensals
        
        # Process phagocytes
        # First handle bacteria registry shifting (digestion) every digestion_time steps
        if t % params['digestion_time'] == 0:
            # Shift registry to the right (oldest bacteria are digested)
            phagocyte_bacteria_registry = np.column_stack([
                np.zeros((n_phagocytes, 1), dtype=np.int32),
                phagocyte_bacteria_registry[:, :-1]
            ])
        
        # Process M0 phagocytes first (check for activation)
        M0_indices = np.where(phagocyte_phenotype == 0)[0]
        for i in M0_indices:
            px, py = phagocyte_x[i], phagocyte_y[i]
            
            # Check for activation
            avg_DAMP = get_8n_avg_signal_fast(px, py, params['act_radius_DAMPs'], DAMPs, grid_size)
            avg_SAMP = get_8n_avg_signal_fast(px, py, params['act_radius_SAMPs'], SAMPs, grid_size)
            bacteria_count = np.sum(phagocyte_bacteria_registry[i, :] > 0)
            
            # Competition between signals determines activation
            if avg_DAMP >= params['activation_threshold_DAMPs'] and avg_DAMP > avg_SAMP:
                # Activate to M1
                phagocyte_phenotype[i] = 1
                phagocyte_active_age[i] = 1
                phagocyte_activity_ROS[i] = params['activity_ROS_M1_baseline'] + activity_ROS_M1_step * bacteria_count
                phagocyte_activity_engulf[i] = params['activity_engulf_M1_baseline'] + activity_engulf_M1_step * bacteria_count
            elif avg_SAMP >= params['activation_threshold_SAMPs'] and avg_SAMP > avg_DAMP:
                # Activate to M2
                phagocyte_phenotype[i] = 2
                phagocyte_active_age[i] = 1
                phagocyte_activity_ROS[i] = params['activity_ROS_M2_baseline']
                phagocyte_activity_engulf[i] = params['activity_engulf_M2_baseline'] + activity_engulf_M2_step * bacteria_count

        # Check aged phagocytes for state transitions
        # Process M1 and M2 phagocytes
        active_indices = np.where(phagocyte_phenotype != 0)[0]
        # Age all active phagocytes
        phagocyte_active_age[active_indices] += 1
        
        # Process active phagocytes for engulfment and check aged ones for state transitions
        for i in active_indices:
            px, py = phagocyte_x[i], phagocyte_y[i]
            
            # Check aged phagocytes for state transitions
            if phagocyte_active_age[i] >= params['active_age_limit']:
                avg_DAMP = get_8n_avg_signal_fast(px, py, params['act_radius_DAMPs'], DAMPs, grid_size)
                avg_SAMP = get_8n_avg_signal_fast(px, py, params['act_radius_SAMPs'], SAMPs, grid_size)
                bacteria_count = np.sum(phagocyte_bacteria_registry[i, :] > 0)
                
                if avg_DAMP >= params['activation_threshold_DAMPs'] and avg_DAMP > avg_SAMP:
                    # Transition/stay M1
                    phagocyte_phenotype[i] = 1
                    phagocyte_active_age[i] = 1
                    phagocyte_activity_ROS[i] = params['activity_ROS_M1_baseline'] + activity_ROS_M1_step * bacteria_count
                    phagocyte_activity_engulf[i] = params['activity_engulf_M1_baseline'] + activity_engulf_M1_step * bacteria_count
                elif avg_SAMP >= params['activation_threshold_SAMPs'] and avg_SAMP > avg_DAMP:
                    # Transition/stay M2
                    phagocyte_phenotype[i] = 2
                    phagocyte_active_age[i] = 1
                    phagocyte_activity_ROS[i] = params['activity_ROS_M2_baseline']
                    phagocyte_activity_engulf[i] = params['activity_engulf_M2_baseline'] + activity_engulf_M2_step * bacteria_count
                elif avg_SAMP < params['activation_threshold_SAMPs'] and avg_DAMP < params['activation_threshold_DAMPs']:
                    # Return to M0
                    phagocyte_phenotype[i] = 0  # Back to M0
                    phagocyte_active_age[i] = 0
                    phagocyte_activity_ROS[i] = params['activity_ROS_M0_baseline']
                    phagocyte_activity_engulf[i] = params['activity_engulf_M0_baseline']
                    # Note: Keep bacteria memory! - this is probably a bug though

        # Age Tregs
        active_tregs = (treg_phenotype == 1)
        old_tregs    = (treg_active_age >= params['active_age_limit'])
        young_tregs  = (treg_active_age < params['active_age_limit'])
        
        if len(active_tregs)>0:
            if len(young_tregs)>0:
                treg_active_age[young_tregs] += 1
            if len(old_tregs)>0:
                treg_phenotype[old_tregs] = 0
                treg_active_age[old_tregs] = 0
                treg_activity_SAMPs_binary[old_tregs] = 0
        

        # Process active phagocytes for engulfment and check aged ones for state transitions
        all_indices = np.where(np.isin(phagocyte_phenotype, [0, 1, 2]))[0]
        
        for i in all_indices:
            px, py = phagocyte_x[i], phagocyte_y[i]
            
            pathogen_indices = np.where((pathogen_coords[:, 0] == px) & (pathogen_coords[:, 1] == py))[0]
            if len(pathogen_indices) > 0:
                engulf_success = np.random.rand(len(pathogen_indices)) < phagocyte_activity_engulf[i]
                indices_to_engulf = np.array(pathogen_indices)[engulf_success]
    
                if len(indices_to_engulf) > 0:
                    # Update counters
                    phagocyte_pathogens_engulfed[i] += len(indices_to_engulf)
                    pathogens_killed_by_Mac += len(indices_to_engulf)
                    # Remove all engulfed pathogens at once
                    pathogen_coords = np.delete(pathogen_coords, indices_to_engulf, axis=0)

                    # Update registry
                    bacteria_count = np.sum(phagocyte_bacteria_registry[i, :] > 0)
                    n_insert = min(len(indices_to_engulf), params['cc_phagocyte'] - bacteria_count)
                    phagocyte_bacteria_registry[i, bacteria_count:bacteria_count+n_insert] = 2

            commensal_indices = np.where((commensal_coords[:, 0] == px) & (commensal_coords[:, 1] == py))[0]
            if len(commensal_indices) > 0:
                engulf_success = np.random.rand(len(commensal_indices)) < phagocyte_activity_engulf[i]
                indices_to_engulf = np.array(commensal_indices)[engulf_success]
                
                if len(indices_to_engulf) > 0:
                    # Update counters
                    phagocyte_commensals_engulfed[i] += len(indices_to_engulf)
                    commensals_killed_by_Mac += len(indices_to_engulf)
                    # Remove all engulfed pathogens at once
                    commensal_coords = np.delete(commensal_coords, indices_to_engulf, axis=0)

                    # Update registry
                    bacteria_count = np.sum(phagocyte_bacteria_registry[i, :] > 0)
                    n_insert = min(len(indices_to_engulf), params['cc_phagocyte'] - bacteria_count)
                    phagocyte_bacteria_registry[i, bacteria_count:bacteria_count+n_insert] = 2
        
        # Treg activation
        if params['allow_tregs_to_do_their_job'] == 1:
            for i in range(n_phagocytes):
                if phagocyte_phenotype[i] == 1:  # M1
                    px, py = phagocyte_x[i], phagocyte_y[i]
                    
                    # Find nearby Tregs using Manhattan distance
                    treg_distances_x = np.abs(treg_x - px)
                    treg_distances_y = np.abs(treg_y - py)
                    nearby_tregs = np.where(
                        (treg_distances_x <= params['treg_vicinity_effect']) & 
                        (treg_distances_y <= params['treg_vicinity_effect']) 
                    )[0]
                    
                    if len(nearby_tregs) > 0:
                        # Use the engulfed COUNTERS, not the registry!
                        num_pat_antigens = phagocyte_pathogens_engulfed[i]
                        num_com_antigens = phagocyte_commensals_engulfed[i]
                        
                        if (num_pat_antigens + num_com_antigens) > 0:
                            rat_com_pat_real = num_com_antigens / (num_com_antigens + num_pat_antigens)
                            
                            # Use GLOBAL treg_discrimination_efficiency for ALL Tregs
                            # This is the same for every Treg, not individualized
                            alpha = (1 - params['treg_discrimination_efficiency']) + params['treg_discrimination_efficiency'] * (rat_com_pat_real * precision)
                            beta = (1 - params['treg_discrimination_efficiency']) + params['treg_discrimination_efficiency'] * ((1 - rat_com_pat_real) * precision)
                            rat_com_pat = sample_beta_numba(alpha, beta)
                            
                            if rat_com_pat > params['rat_com_pat_threshold']:
                                # Activate ALL nearby Tregs at once (like in R code)
                                treg_phenotype[nearby_tregs] = 1
                                treg_activity_SAMPs_binary[nearby_tregs] = 1
                                treg_active_age[nearby_tregs] = 1
                                
                                if params['allow_tregs_to_suppress_cognate']:
                                    # Suppress the M1 phagocyte
                                    phagocyte_phenotype[i] = 2
                                    phagocyte_active_age[i] = 1
                                    bacteria_count = np.sum(phagocyte_bacteria_registry[i, :] > 0)
                                    phagocyte_activity_ROS[i] = params['activity_ROS_M2_baseline']
                                    phagocyte_activity_engulf[i] = params['activity_engulf_M2_baseline'] + activity_engulf_M2_step * bacteria_count
        
        # Kill bacteria with ROS
        if len(pathogen_coords) > 0:
            pathogen_ROS = np.array([get_8n_avg_signal_fast(p[0], p[1], params['act_radius_ROS'], ROS, grid_size) 
                                    for p in pathogen_coords])
            survivors = pathogen_ROS <= params['th_ROS_microbe']
            killed = np.sum(~survivors)
            pathogens_killed_by_ROS += killed
            pathogen_coords = pathogen_coords[survivors]
        
        if len(commensal_coords) > 0:
            commensal_ROS = np.array([get_8n_avg_signal_fast(c[0], c[1], params['act_radius_ROS'], ROS, grid_size) 
                                     for c in commensal_coords])
            survivors = commensal_ROS <= params['th_ROS_microbe']
            killed = np.sum(~survivors)
            commensals_killed_by_ROS += killed
            commensal_coords = commensal_coords[survivors]
        
        # Injure epithelium
        pathogen_epithelium_counts = np.zeros(grid_size, dtype=np.int32)
        if len(pathogen_coords) > 0:
            # Bacteria at y=0 are touching the epithelium
            epithelium_pathogens = pathogen_coords[pathogen_coords[:, 1] == 0]
            if len(epithelium_pathogens) > 0:
                unique, counts = np.unique(epithelium_pathogens[:, 0], return_counts=True)
                pathogen_epithelium_counts[unique] = counts
        
        for i in range(grid_size):
            # Get ROS in vicinity from row y=1 (above epithelium)
            x_start = max(0, i - params['act_radius_ROS'])
            x_end = min(grid_size, i + params['act_radius_ROS'] + 1)
            mean_ros = np.mean(ROS[0, x_start:x_end])  
            
            # Increase injury from pathogens
            count_pathogens = pathogen_epithelium_counts[i]
            epithelium_injury[i] += int(logistic_scaled_0_to_5_quantized(count_pathogens, params['k_in'], params['x0_in']))
            
            # Increase injury from ROS
            if mean_ros > params['th_ROS_epith_recover']:
                epithelium_injury[i] += 1
            
            # Cap at maximum
            epithelium_injury[i] = min(epithelium_injury[i], params['max_level_injury'])
            
            # Stochastic recovery
            if epithelium_injury[i] > 0 and np.random.random() < params['epith_recovery_chance']:
                epithelium_injury[i] = max(0, epithelium_injury[i] - 1)
        
        # Save longitudinal data
        epithelium_counts = np.bincount(epithelium_injury, minlength=6)[:6]
        epithelium_longitudinal[t, :] = epithelium_counts
        
        # Count phagocyte phenotypes
        phagocyte_counts = [np.sum(phagocyte_phenotype == 0)]  # M0
        
        # M1 by activation level
        m1_mask = phagocyte_phenotype == 1
        if np.any(m1_mask):
            m1_ages = phagocyte_active_age[m1_mask]
            m1_counts = np.bincount(m1_ages, minlength=params['cc_phagocyte']+1)[:params['cc_phagocyte']+1]
        else:
            m1_counts = np.zeros(params['cc_phagocyte']+1, dtype=np.int32)
        phagocyte_counts.extend(m1_counts)
        
        # M2 by activation level
        m2_mask = phagocyte_phenotype == 2
        if np.any(m2_mask):
            m2_ages = phagocyte_active_age[m2_mask]
            m2_counts = np.bincount(m2_ages, minlength=params['cc_phagocyte']+1)[:params['cc_phagocyte']+1]
        else:
            m2_counts = np.zeros(params['cc_phagocyte']+1, dtype=np.int32)
        phagocyte_counts.extend(m2_counts)
        
        macrophages_longitudinal[t, :] = phagocyte_counts
        
        # Microbe counts
        microbes_longitudinal[t, :] = [len(commensal_coords), len(pathogen_coords)]
        
        # Treg counts
        tregs_longitudinal[t, :] = [np.sum(treg_phenotype == 0), np.sum(treg_phenotype == 1)]
        
        # Death counts (cumulative for this implementation)
        microbes_cumdeath_longitudinal[t, :] = [
            commensals_killed_by_ROS, commensals_killed_by_Mac, 0, 0,  # Commensals (M0, M1, M2 not tracked separately)
            pathogens_killed_by_ROS, pathogens_killed_by_Mac, 0, 0     # Pathogens (M0, M1, M2 not tracked separately)
        ]
    
    print("Simulation complete!")
    
    # Create results DataFrame
    results = create_results_dataframe(
        epithelium_longitudinal, macrophages_longitudinal, microbes_longitudinal,
        tregs_longitudinal, microbes_cumdeath_longitudinal, params
    )
    
    return results


def create_results_dataframe(epithelium_longitudinal, macrophages_longitudinal, 
                            microbes_longitudinal, tregs_longitudinal, 
                            microbes_cumdeath_longitudinal, params):
    """Create a pandas DataFrame with all simulation results."""
    
    # Combine all data
    data = np.column_stack([
        epithelium_longitudinal,
        macrophages_longitudinal,
        microbes_longitudinal,
        tregs_longitudinal,
        microbes_cumdeath_longitudinal
    ])
    
    # Create column names
    column_names = [
        'epithelial_healthy', 'epithelial_inj_1', 'epithelial_inj_2', 
        'epithelial_inj_3', 'epithelial_inj_4', 'epithelial_inj_5',
        'phagocyte_M0'
    ]
    
    # M1 levels
    for i in range(params['cc_phagocyte'] + 1):
        column_names.append(f'phagocyte_M1_L_{i}')
    
    # M2 levels
    for i in range(params['cc_phagocyte'] + 1):
        column_names.append(f'phagocyte_M2_L_{i}')
    
    column_names.extend([
        'commensal', 'pathogen',
        'treg_resting', 'treg_active',
        'C_ROS', 'C_Mac', 'C_M1', 'C_M2',
        'P_ROS', 'P_Mac', 'P_M1', 'P_M2'
    ])
    
    # Create DataFrame
    df = pd.DataFrame(data, columns=column_names)
    
    # Add time and metadata
    df['t'] = np.arange(params['t_max'])
    df['sterile'] = params['sterile']
    df['allow_tregs_to_do_their_job'] = params['allow_tregs_to_do_their_job']
    df['allow_tregs_to_suppress_cognate'] = params['allow_tregs_to_suppress_cognate']
    df['randomize_tregs'] = params['randomize_tregs']
    
    # Reorder columns
    first_cols = ['t', 'sterile', 'allow_tregs_to_do_their_job', 
                  'allow_tregs_to_suppress_cognate', 'randomize_tregs']
    other_cols = [col for col in df.columns if col not in first_cols]
    df = df[first_cols + other_cols]
    
    return df


def plot_results(df):
    """Create visualization plots for the simulation results."""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Epithelial dynamics
    ax = axes[0, 0]
    epithelial_cols = [col for col in df.columns if col.startswith('epithelial_')]
    for col in epithelial_cols:
        ax.plot(df['t'], df[col], label=col.replace('epithelial_', ''), linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Count')
    ax.set_title('Epithelial Cell Dynamics')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Microbe dynamics
    ax = axes[0, 1]
    ax.plot(df['t'], df['commensal'], label='Commensals', color='turquoise', linewidth=2)
    ax.plot(df['t'], df['pathogen'], label='Pathogens', color='firebrick', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Count')
    ax.set_title('Microbe Dynamics')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Phagocyte dynamics (simplified)
    ax = axes[1, 0]
    # Sum M1 and M2 across all levels
    m1_cols = [col for col in df.columns if 'phagocyte_M1_' in col]
    m2_cols = [col for col in df.columns if 'phagocyte_M2_' in col]
    
    ax.plot(df['t'], df['phagocyte_M0'], label='M0', color='grey', linewidth=2)
    ax.plot(df['t'], df[m1_cols].sum(axis=1), label='M1 (total)', color='magenta', linewidth=2)
    ax.plot(df['t'], df[m2_cols].sum(axis=1), label='M2 (total)', color='green', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Count')
    ax.set_title('Phagocyte Phenotype Dynamics')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Treg dynamics
    ax = axes[1, 1]
    ax.plot(df['t'], df['treg_resting'], label='Resting', color='plum', linewidth=2)
    ax.plot(df['t'], df['treg_active'], label='Active', color='purple', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Count')
    ax.set_title('Treg Dynamics')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('simulation_results.png', dpi=150, bbox_inches='tight')
    #plt.show()
    
    print("Plots saved to simulation_results.png")


if __name__ == "__main__":
    # Start timer
    start_time = time()
    
    # Initialize parameters
    params = initialize_parameters()
    
    # Run simulation
    results_df = run_simulation(params)
    
    # Calculate runtime
    runtime = time() - start_time
    print(f"\nSimulation runtime: {runtime:.2f} seconds")
    
    # Save results
    results_df.to_csv('simulation_results.csv', index=False)
    print("Results saved to simulation_results.csv")
    
    # Create plots
    plot_results(results_df)
    
    # Display summary statistics
    print("\n=== Summary Statistics ===")
    print(f"Final epithelial health: {results_df['epithelial_healthy'].iloc[-1]} healthy cells")
    print(f"Final pathogen count: {results_df['pathogen'].iloc[-1]}")
    print(f"Final commensal count: {results_df['commensal'].iloc[-1]}")
    print(f"Total pathogens killed by ROS: {results_df['P_ROS'].iloc[-1]}")
    print(f"Total pathogens killed by Macrophages: {results_df['P_Mac'].iloc[-1]}")