import numpy as np


def detect_large_angle_of_incidence(surface_normal, viewing_direction, threshold):
    """
    Check if the angle of incidence is larger than a threshold.

    Parameters
    ----------
    surface_normal : np.ndarray (3,)
        Surface normal vector.
    viewing_direction : np.ndarray (3,)
        Viewing direction vector.
    threshold : float
        Angle threshold in degrees.

    Returns
    -------
    is_large : bool
        True if angle of incidence > threshold, or if dot product >= 0.
    """
    surface_normal = surface_normal / (np.linalg.norm(surface_normal) + 1e-12)
    viewing_direction = viewing_direction / (np.linalg.norm(viewing_direction) + 1e-12)

    if surface_normal.shape != viewing_direction.shape:
        raise ValueError("surface_normal and viewing_direction must have the same shape")

    dot_product = np.dot(surface_normal, viewing_direction)

    if dot_product >= 0:
        return True

    angle_of_incidence = np.rad2deg(np.arccos(-dot_product))
    return angle_of_incidence > threshold
