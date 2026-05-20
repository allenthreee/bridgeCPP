import numpy as np


def get_grid_id(point, grid_size, x_lower, x_upper, y_lower, y_upper,
                z_lower, z_upper, dimension, find_neighbor):
    """
    Compute grid ID for a 3D or 2D point, and optionally its neighbour list.

    Parameters
    ----------
    point : np.ndarray (3,) or (2,)
        Point coordinates.
    grid_size : float
        Grid cell size.
    x_lower, x_upper : float
        X-axis bounds.
    y_lower, y_upper : float
        Y-axis bounds.
    z_lower, z_upper : float
        Z-axis bounds (ignored for 2D).
    dimension : int
        2 or 3.
    find_neighbor : bool
        Whether to compute neighbour list.

    Returns
    -------
    grid_id : int
        1-based grid index.
    neighbour_list : list
        List of neighbouring grid IDs (1-based). Only returned if find_neighbor=True.
    """
    neighbour_list = []

    if dimension == 3:
        point = point.copy()
        if point[2] < 0:
            point[2] = 0

        x_size = (x_upper - x_lower) / grid_size
        y_size = (y_upper - y_lower) / grid_size
        z_size = (z_upper - z_lower) / grid_size
        no_of_grid = int(x_size * y_size * z_size)

        grid_x = int(np.floor((point[0] - x_lower) / grid_size))
        grid_y = int(np.floor((point[1] - y_lower) / grid_size))
        grid_z = int(np.floor((point[2] - z_lower) / grid_size))

        grid_id = grid_x + grid_y * int(x_size) + grid_z * int(x_size) * int(y_size) + 1

        if grid_id < 0:
            raise ValueError("grid_id < 0")

        if find_neighbor:
            xs = int(x_size)
            ys = int(y_size)
            candidates = [
                grid_id + 1, grid_id - 1,
                grid_id + xs, grid_id + xs + 1, grid_id + xs - 1,
                grid_id - xs, grid_id - xs + 1, grid_id - xs - 1,
                grid_id + xs * ys,
                grid_id + xs * ys + 1,
                grid_id + xs * ys - 1,
                grid_id + xs * ys + xs,
                grid_id + xs * ys + xs + 1,
                grid_id + xs * ys + xs - 1,
                grid_id + xs * ys - xs,
                grid_id + xs * ys - xs + 1,
                grid_id + xs * ys - xs - 1,
                grid_id - xs * ys,
                grid_id - xs * ys + 1,
                grid_id - xs * ys - 1,
                grid_id - xs * ys + xs,
                grid_id - xs * ys + xs + 1,
                grid_id - xs * ys + xs - 1,
                grid_id - xs * ys - xs,
                grid_id - xs * ys - xs + 1,
                grid_id - xs * ys - xs - 1,
            ]
            neighbour_list = [n for n in candidates if 1 <= n <= no_of_grid]

    else:  # dimension == 2
        x_size = (x_upper - x_lower) / grid_size
        y_size = (y_upper - y_lower) / grid_size
        xs = int(x_size)

        grid_x = int(np.floor((point[0] - x_lower) / grid_size))
        grid_y = int(np.floor((point[1] - y_lower) / grid_size))
        grid_id = grid_x + grid_y * xs + 1

        if find_neighbor:
            # Left/Right neighbours
            if grid_x == 0:
                neighbour_list.append(grid_id + 1)
            elif grid_x == xs - 1:
                neighbour_list.append(grid_id - 1)
            else:
                neighbour_list.extend([grid_id - 1, grid_id + 1])

            if grid_y == 0:
                temp_id = grid_x + (grid_y + 1) * xs + 1
                neighbour_list.append(temp_id)
                if grid_x == 0:
                    neighbour_list.append(temp_id + 1)
                elif grid_x == xs - 1:
                    neighbour_list.append(temp_id - 1)
                else:
                    neighbour_list.extend([temp_id - 1, temp_id + 1])
            elif grid_y == int(y_size) - 1:
                temp_id = grid_x + (grid_y - 1) * xs + 1
                neighbour_list.append(temp_id)
                if grid_x == 0:
                    neighbour_list.append(temp_id + 1)
                elif grid_x == xs - 1:
                    neighbour_list.append(temp_id - 1)
                else:
                    neighbour_list.extend([temp_id - 1, temp_id + 1])
            else:
                temp_id1 = grid_x + (grid_y - 1) * xs + 1  # UP
                temp_id2 = grid_x + (grid_y + 1) * xs + 1  # DOWN
                neighbour_list.extend([temp_id1, temp_id2])
                if grid_x == 0:
                    neighbour_list.extend([temp_id1 + 1, temp_id2 + 1])
                elif grid_x == xs - 1:
                    neighbour_list.extend([temp_id1 - 1, temp_id2 - 1])
                else:
                    neighbour_list.extend([
                        temp_id1 - 1, temp_id1 + 1,
                        temp_id2 - 1, temp_id2 + 1,
                    ])

    return grid_id, neighbour_list
