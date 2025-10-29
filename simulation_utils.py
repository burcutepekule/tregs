"""
Utility functions for the Agent-Based Model simulation.
All computationally intensive functions are JIT-compiled with Numba for maximum speed.
"""

import numpy as np
from numba import jit, njit, prange
from scipy.special import expit as sigmoid

# Use njit (no-python mode) for maximum speed
@njit
def get_middle_percent(seq_vector, percent):
    """Get middle percentage of a sequence."""
    n_total = len(seq_vector)
    n_select = int(np.ceil(n_total * percent / 100))
    
    mid = n_total // 2
    half_window = n_select // 2
    
    start_idx = max(0, mid - half_window)# No +1 as in R (0-based indexing)
    end_idx = min(n_total, start_idx + n_select)# No -1 as in R (Python excludes the end)

    
    return seq_vector[start_idx:end_idx] 

#@njit
def iszero_coordinates(x):
    """Generate coordinate transformations with special handling for zeros."""
    n = len(x)
    
    # Initialize with default sampling: -1, 0, 1 (equal probability)
    y = np.random.choice([-1, 0, 1], size=n, replace=True)
    
    # Replace where x == 0 with sample from -1 or 1
    zero_idx = np.where(x == 0)[0]
    if len(zero_idx) > 0:
        y[zero_idx] = np.random.choice([-1, 1], size=len(zero_idx), replace=True)
    
    return y


def logistic_scaled_0_to_5_quantized(x, k, x0):
    """Scaled and quantized logistic function."""
    # Using scipy's sigmoid for vectorized operations
    return np.round(5 * sigmoid(k * (x - x0)))


@njit(parallel=True)
def diffuse_matrix(mat, D, max_cell_value):
    """
    Fast 8-neighbor diffusion using Moore neighborhood.
    Uses parallel processing for maximum speed.
    """
    nr, nc = mat.shape
    
    # Create padded matrix
    padded = np.zeros((nr + 2, nc + 2))
    padded[1:nr+1, 1:nc+1] = mat
    
    # Compute 8-neighbor Laplacian (Moore neighborhood)
    laplacian = np.zeros((nr, nc))
    
    # Parallel computation of Laplacian
    for i in prange(nr):
        for j in range(nc):
            laplacian[i, j] = (
                padded[i, j] +      # top-left
                padded[i, j+1] +    # top
                padded[i, j+2] +    # top-right
                padded[i+1, j] +    # left
                padded[i+1, j+2] +  # right
                padded[i+2, j] +    # bottom-left
                padded[i+2, j+1] +  # bottom
                padded[i+2, j+2] -  # bottom-right
                8 * mat[i, j]       # center subtraction
            )
    
    # Update matrix with diffusion
    mat_new = mat + D * laplacian
    
    # Apply maximum constraint
    mat_new = np.minimum(max_cell_value, mat_new)
    
    return mat_new


@njit
def get_8n_avg_signal_fast(x, y, act_radius_signal, signal_matrix, grid_size):
    """Get average signal in 8-neighborhood."""
    x_start = max(0, x - act_radius_signal)
    x_end = min(grid_size, x + act_radius_signal + 1)
    y_start = max(0, y - act_radius_signal)
    y_end = min(grid_size, y + act_radius_signal + 1)
    
    # Extract submatrix and compute mean
    total = 0.0
    count = 0
    for i in range(y_start, y_end):
        for j in range(x_start, x_end):
            total += signal_matrix[i, j]
            count += 1
    
    return total / count if count > 0 else 0.0

#@njit
def sample_beta_numba(alpha, beta):
    """
    Fast beta distribution sampling using Numba.
    Uses the relationship between Beta and Gamma distributions.
    """
    # Sample from Gamma distributions
    x = np.random.gamma(alpha, 1.0)
    y = np.random.gamma(beta, 1.0)
    
    # Beta is the ratio
    return x / (x + y)