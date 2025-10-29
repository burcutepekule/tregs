"""
Utility functions for the Agent-Based Model simulation.
All computationally intensive functions are JIT-compiled with Numba for maximum speed.
"""

import numpy as np
from numba import jit, njit, prange
from scipy.special import expit as sigmoid
from scipy import stats

# ============================================================================
# RANDOM NUMBER STREAM MANAGEMENT (matching R's approach)
# ============================================================================

class RandomStream:
    """
    Manages a pre-generated stream of uniform random numbers.
    This ensures exact reproducibility with R's custom RNG.
    """
    def __init__(self, stream_file=None):
        if stream_file is not None:
            print(f"Loading random stream from {stream_file}...")
            self.stream = np.loadtxt(stream_file, dtype=np.float64)
            print(f"Loaded {len(self.stream)} random numbers")
        else:
            self.stream = None
        self.index = 0

    def runif(self, n=1):
        """Sample n uniform random numbers from the stream."""
        if self.stream is None:
            # Fallback to numpy if no stream loaded
            return np.random.random(n)

        if self.index + n > len(self.stream):
            raise RuntimeError(f"Random stream exhausted! Requested {n} values at index {self.index}, but stream has {len(self.stream)} values")

        vals = self.stream[self.index:self.index + n]
        self.index += n

        return vals if n > 1 else vals[0]

    def sample(self, x, size, replace=False, prob=None):
        """
        Sample from array x, matching R's sample behavior.
        Uses inverse CDF method for exact reproducibility.
        """
        n = len(x)
        if not replace and size > n:
            raise ValueError("Cannot take a sample larger than the population when replace=False")

        # Normalize probabilities
        if prob is None:
            probs = np.ones(n) / n
        else:
            probs = np.array(prob) / np.sum(prob)

        # Get uniform random values
        u = self.runif(size) if size > 1 else np.array([self.runif(1)])

        # Compute cumulative probabilities
        cumprob = np.cumsum(probs)

        # Use searchsorted to find intervals (equivalent to R's findInterval + 1)
        indices = np.searchsorted(cumprob, u, side='right')

        return x[indices]

    def rbeta(self, n, shape1, shape2):
        """
        Sample from Beta distribution using quantile function.
        Matches R's qbeta approach for exact reproducibility.
        """
        u = self.runif(n) if n > 1 else np.array([self.runif(1)])
        return stats.beta.ppf(u, shape1, shape2)

    def rgamma(self, n, shape, rate=1.0):
        """
        Sample from Gamma distribution using quantile function.
        Matches R's qgamma approach for exact reproducibility.
        Note: R's rgamma uses 'rate' parameter, while scipy uses 'scale' = 1/rate
        """
        u = self.runif(n) if n > 1 else np.array([self.runif(1)])
        scale = 1.0 / rate
        return stats.gamma.ppf(u, shape, scale=scale)

    def choice(self, a, size=None, replace=True, p=None):
        """
        Equivalent to np.random.choice but using our stream.
        """
        if isinstance(a, int):
            x = np.arange(a)
        else:
            x = np.array(a)

        if size is None:
            size = 1
            return self.sample(x, size, replace=replace, prob=p)[0]
        else:
            return self.sample(x, size, replace=replace, prob=p)

    def randint(self, low, high, size=None):
        """
        Equivalent to np.random.randint but using our stream.
        """
        if size is None:
            u = self.runif(1)
            return int(low + u * (high - low))
        else:
            u = self.runif(size)
            return np.floor(low + u * (high - low)).astype(np.int32)

# Global random stream instance
_rng_stream = None

def initialize_random_stream(stream_file=None):
    """Initialize the global random stream."""
    global _rng_stream
    _rng_stream = RandomStream(stream_file)
    return _rng_stream

def get_rng():
    """Get the global random stream instance."""
    global _rng_stream
    if _rng_stream is None:
        # Initialize with default (no file = use numpy)
        _rng_stream = RandomStream(None)
    return _rng_stream

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

def iszero_coordinates(x):
    """
    Generate coordinate transformations with special handling for zeros.
    Uses the global random stream for reproducibility.
    """
    rng = get_rng()
    n = len(x)

    # Initialize with default sampling: -1, 0, 1 (equal probability)
    y = rng.choice([-1, 0, 1], size=n, replace=True)

    # Replace where x == 0 with sample from -1 or 1
    zero_idx = np.where(x == 0)[0]
    if len(zero_idx) > 0:
        y[zero_idx] = rng.choice([-1, 1], size=len(zero_idx), replace=True)

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

def sample_beta_from_gamma(alpha, beta):
    """
    Sample from Beta distribution using Gamma distributions.
    Matches R's sample_rbeta implementation exactly.
    Uses quantile-based gamma sampling for reproducibility.
    """
    rng = get_rng()

    # Sample from Gamma distributions using quantile functions
    x = rng.rgamma(1, shape=alpha, rate=1.0)
    y = rng.rgamma(1, shape=beta, rate=1.0)

    # Beta is the ratio - extract scalar from array
    result = x / (x + y)
    return float(result[0]) if hasattr(result, '__len__') else float(result)