import numpy as np


def hermite_plot(points):
    """
    Generate Hermite spline interpolation for 3D points.

    Parameters
    ----------
    points : np.ndarray (N, 3)
        Input waypoints.

    Returns
    -------
    interp_points : np.ndarray (M, 3)
        Interpolated points along the Hermite spline.
    """
    n = points.shape[0]

    # Calculate tangent vectors using finite differencing
    tangents = np.zeros((n, 3))
    for i in range(1, n - 1):
        diff = points[i + 1, :] - points[i - 1, :]
        tangents[i, :] = diff / np.linalg.norm(diff)

    tangents[0, :] = (points[1, :] - points[0, :]) / (np.linalg.norm(points[1, :] - points[0, :]) + 1e-12)
    tangents[n - 1, :] = (points[n - 1, :] - points[n - 2, :]) / (np.linalg.norm(points[n - 1, :] - points[n - 2, :]) + 1e-12)

    # Parameter t from 0 to 1
    t = np.linspace(0, 1, n)

    # Interpolated points — pre-allocate based on actual point count
    pts_per_seg = 100
    interp_points = np.zeros(((n - 1) * pts_per_seg, 3))

    for i in range(n - 1):
        ti = np.linspace(t[i], t[i + 1], pts_per_seg)

        # Hermite blending functions
        h00 = 2 * ti**3 - 3 * ti**2 + 1
        h01 = ti**3 - 2 * ti**2 + ti
        h10 = -2 * ti**3 + 3 * ti**2
        h11 = ti**3 - ti**2

        start_idx = i * pts_per_seg
        end_idx = (i + 1) * pts_per_seg
        interp_points[start_idx:end_idx, :] = (
            np.outer(h00, points[i, :]) +
            np.outer(h01, tangents[i, :]) +
            np.outer(h10, points[i + 1, :]) +
            np.outer(h11, tangents[i + 1, :])
        )

    return interp_points
