import numpy as np

def _approximate_k(total_RAM_available_gb, num_samples):
    if num_samples == 0:
        raise ValueError("num_samples must be greater than 0.")
    
    if total_RAM_available_gb < 0.21:
        raise ValueError("total_RAM_available_gb must be greater than 0.21.")
    
    k = 10**7 * (total_RAM_available_gb - 0.21) / (2.35 * num_samples)
    return np.min([num_samples - 1, int(np.floor(k))])
    
    