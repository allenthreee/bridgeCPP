import time
import numpy as np
from scipy.spatial.distance import pdist, squareform


def ratio_compute(centroid_uv, centroid_3d):
    """
    Compute the UV-to-3D unit ratio from centroid distances.

    Uses scipy.spatial.distance.pdist for efficient O(N²) pairwise
    computation without broadcasting large intermediate matrices.

    Parameters
    ----------
    centroid_uv : np.ndarray (N, 2)
        UV centroids.
    centroid_3d : np.ndarray (N, 3)
        3D centroids.

    Returns
    -------
    uv_to_3d_ratio : float
        Mean ratio of UV distance to 3D distance.
    """
    t_start = time.time()

    if centroid_uv.shape[0] != centroid_3d.shape[0]:
        raise ValueError("centroid_uv and centroid_3d must have the same number of rows")

    # Compute condensed pairwise distances (avoids N×N intermediate arrays)
    c2c_dist_uv = squareform(pdist(centroid_uv))   # (N, N)
    c2c_dist_3d = squareform(pdist(centroid_3d))   # (N, N)

    # Avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = c2c_dist_uv / c2c_dist_3d
        ratio[c2c_dist_3d == 0] = np.nan

    uv_to_3d_ratio = np.nanmean(ratio)

    t_elapsed = time.time() - t_start
    print(f"ratio_compute elapsed: {t_elapsed:.2f}s")

    return uv_to_3d_ratio
