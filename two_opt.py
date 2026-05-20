import numpy as np


def _path_length(path):
    """Compute total path length (closed loop)."""
    if path.shape[0] <= 1:
        return 0.0
    seg_lengths = np.linalg.norm(np.diff(path, axis=0), axis=1)
    closing = np.linalg.norm(path[-1, :] - path[0, :])
    return np.sum(seg_lengths) + closing


def _two_opt_swap(path, i, j):
    """Perform a 2-opt swap: reverse segment between indices i and j."""
    new_path = np.vstack([
        path[:i + 1, :],
        path[j:i:-1, :],   # reverse from j down to i+1
        path[j + 1:, :],
    ])
    return new_path


def two_opt(path, n):
    """
    2-opt optimization for a path.

    Parameters
    ----------
    path : np.ndarray (N, D)
        Input waypoints.
    n : int
        Maximum search window size.

    Returns
    -------
    optimized_path : np.ndarray (N, D)
        Optimized path.
    """
    num_points = path.shape[0]
    best_path = path.copy()
    improved = True

    while improved:
        improved = False
        for i in range(num_points - 2):
            j_max = min(i + n, num_points - 1)
            for j in range(i + 2, j_max + 1):
                new_path = _two_opt_swap(best_path, i, j)
                if _path_length(new_path) < _path_length(best_path):
                    best_path = new_path
                    improved = True

    return best_path
