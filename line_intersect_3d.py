import numpy as np


def line_intersect_3d(pa, pb):
    """
    Find intersection point of lines in 3D space, in the least squares sense.

    Parameters
    ----------
    pa : np.ndarray (N, 3)
        Starting points of N lines.
    pb : np.ndarray (N, 3)
        End points of N lines.

    Returns
    -------
    p_intersect : np.ndarray (3,)
        Best intersection point of the N lines, in least squares sense.
    distances : np.ndarray (N,)
        Distances from intersection point to the input lines.
    """
    si = pb - pa  # N lines as vectors
    norms = np.linalg.norm(si, axis=1, keepdims=True)
    norms[norms == 0] = 1.0  # avoid division by zero
    ni = si / norms
    nx, ny, nz = ni[:, 0], ni[:, 1], ni[:, 2]

    sxx = np.sum(nx**2 - 1)
    syy = np.sum(ny**2 - 1)
    szz = np.sum(nz**2 - 1)
    sxy = np.sum(nx * ny)
    sxz = np.sum(nx * nz)
    syz = np.sum(ny * nz)

    s = np.array([[sxx, sxy, sxz],
                  [sxy, syy, syz],
                  [sxz, syz, szz]])

    cx = np.sum(pa[:, 0] * (nx**2 - 1) + pa[:, 1] * (nx * ny) + pa[:, 2] * (nx * nz))
    cy = np.sum(pa[:, 0] * (nx * ny) + pa[:, 1] * (ny**2 - 1) + pa[:, 2] * (ny * nz))
    cz = np.sum(pa[:, 0] * (nx * nz) + pa[:, 1] * (ny * nz) + pa[:, 2] * (nz**2 - 1))
    c = np.array([cx, cy, cz])

    # Use least-squares in case S is singular (parallel lines)
    try:
        p_intersect = np.linalg.solve(s, c)
    except np.linalg.LinAlgError:
        p_intersect = np.linalg.lstsq(s, c, rcond=None)[0]

    n = pa.shape[0]
    distances = np.zeros(n)
    for i in range(n):
        distances[i] = np.linalg.norm(pa[i, :] - p_intersect)

    return p_intersect, distances
