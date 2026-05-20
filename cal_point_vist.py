import numpy as np


def cal_point_vist(viewpoint_location, viewing_angle, point_location, fov_degrees):
    """
    Check if a point is within the view cone.

    Parameters
    ----------
    viewpoint_location : np.ndarray (3,)
        Position of the viewpoint.
    viewing_angle : np.ndarray (3,)
        Viewing direction vector.
    point_location : np.ndarray (3,)
        Position of the point to check.
    fov_degrees : float
        Field of view in degrees.

    Returns
    -------
    is_within_view_cone : bool
    """
    vector_to_point = point_location - viewpoint_location
    vector_to_point = vector_to_point / (np.linalg.norm(vector_to_point) + 1e-12)
    viewing_angle = viewing_angle / (np.linalg.norm(viewing_angle) + 1e-12)

    dot_product = np.dot(vector_to_point, viewing_angle)
    cos_half_fov = np.cos(np.deg2rad(fov_degrees / 2))

    return dot_product >= cos_half_fov
