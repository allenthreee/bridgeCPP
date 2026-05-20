import numpy as np
from get_grid_id import get_grid_id
from cal_point_vist import cal_point_vist
from detect_large_angle import detect_large_angle_of_incidence


def calculate_vist_3d(traj, viewing_dist, viewing_angle, fov, centroid_3d,
                       normal, threshold, all_grid, grid_info):
    """
    Compute 3D visibility along a trajectory.

    Parameters
    ----------
    traj : np.ndarray (T, 3)
        Trajectory waypoints.
    viewing_dist : float
        Maximum viewing distance.
    viewing_angle : np.ndarray (T-1, 3)
        Viewing direction for each segment.
    fov : float
        Field of view in degrees.
    centroid_3d : np.ndarray (N, 3)
        Face centroids.
    normal : np.ndarray (N, 3)
        Face normals.
    threshold : float
        Angle-of-incidence threshold (degrees).
    all_grid : list of dict
        Grid data (each dict has key 'point' with shape Mx4: x,y,z,face_id).
    grid_info : np.ndarray (7,)
        [Q1_x, Q1_y, Q1_z, Q3_x, Q3_y, Q3_z, grid_size].

    Returns
    -------
    visibility_3d : list of dict
        List with key 'Visibility' containing face IDs visible per segment.
    """
    q1_bound = grid_info[:3]
    q3_bound = grid_info[3:6]
    grid_size = grid_info[6]

    visibility_3d = []

    for i in range(traj.shape[0] - 1):
        point1 = traj[i, :]
        point2 = traj[i + 1, :]
        segment = np.column_stack([
            np.linspace(point1[0], point2[0], 5),
            np.linspace(point1[1], point2[1], 5),
            np.linspace(point1[2], point2[2], 5),
        ])

        vist = []
        angle = viewing_angle[i, :]

        for j in range(segment.shape[0]):
            temp_point = segment[j, :]
            grid_id, neighbour = get_grid_id(
                temp_point, grid_size,
                q3_bound[0], q1_bound[0],
                q3_bound[1], q1_bound[1],
                q3_bound[2], q1_bound[2],
                3, True,
            )

            part1 = all_grid[grid_id - 1]['point']  # 0-based index
            point_id_p1 = []
            if part1.size > 0:
                point_id_p1 = part1[:, -1].astype(int)
                part1 = part1[:, :3]
            else:
                part1 = np.empty((0, 3))

            part2_list = []
            for k in range(len(neighbour)):
                p = all_grid[neighbour[k] - 1]['point']
                if p.size > 0:
                    part2_list.append(p)
            if part2_list:
                part2 = np.vstack(part2_list)
                point_id_p2 = part2[:, -1].astype(int)
                part2 = part2[:, :3]
                point_ids = np.concatenate([point_id_p1, point_id_p2])
            else:
                part2 = np.empty((0, 3))
                point_ids = point_id_p1

            neighbour_centroids = np.vstack([part1, part2]) if part1.size > 0 or part2.size > 0 else np.empty((0, 3))

            for k in range(neighbour_centroids.shape[0]):
                temp_c = neighbour_centroids[k, :]
                temp_dist = np.linalg.norm(temp_point - temp_c)

                # Find matching centroid index
                matches = np.where(np.all(np.abs(centroid_3d - temp_c) < 1e-10, axis=1))[0]
                if len(matches) == 0:
                    continue
                cid = matches[0]
                temp_normal = normal[cid, :]

                is_large_angle = detect_large_angle_of_incidence(temp_normal, angle, threshold)

                if temp_dist > viewing_dist or is_large_angle:
                    continue

                is_visible = cal_point_vist(temp_point, angle, temp_c, fov)
                if is_visible:
                    vist.append(cid)

        vist = list(set(vist))
        visibility_3d.append({'Visibility': vist})

    return visibility_3d
