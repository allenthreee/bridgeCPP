import numpy as np
from line_intersect_3d import line_intersect_3d


def untrim_offset_3d(og_segment, offset_distance, normal, ccentroid):
    """
    Offset trajectory segments along surface normals, connecting
    adjacent segments with arcs when they don't intersect.

    Parameters
    ----------
    og_segment : list of dict
        Each dict has keys 'segment' (Mx3 np.ndarray) and 'FaceID' (int).
    offset_distance : float
        Offset distance.
    normal : np.ndarray (K, 3)
        Surface normals, indexed by FaceID.
    ccentroid : np.ndarray (K, 3)
        Face centroids, indexed by FaceID.

    Returns
    -------
    untrim_segments : list of dict
        Offset segments with same structure, plus arc connections.
    """
    n_segments = len(og_segment)
    untrim_segments = []

    for i in range(n_segments):
        temp_face_id = og_segment[i]['FaceID']
        temp_traj = og_segment[i]['segment']

        # Offset direction = normal of this face
        offset_dir = normal[temp_face_id, :]
        xy_angle = np.arctan2(offset_dir[1], offset_dir[0])
        z_angle = np.arctan2(offset_dir[2],
                             np.sqrt(offset_dir[0]**2 + offset_dir[1]**2))

        # Build offset trajectory
        temp_offset_traj = np.zeros_like(temp_traj)
        for ii in range(temp_traj.shape[0]):
            tx = temp_traj[ii, 0] + (offset_distance * np.cos(z_angle)) * np.cos(xy_angle)
            ty = temp_traj[ii, 1] + (offset_distance * np.cos(z_angle)) * np.sin(xy_angle)
            tz = temp_traj[ii, 2] + offset_distance * np.sin(z_angle)
            temp_offset_traj[ii, :] = [tx, ty, tz]

        untrim_segments.append({
            'segment': temp_offset_traj,
            'FaceID': temp_face_id,
        })

        if i == 0:
            continue

        # --- Connect adjacent segments (i-1 and i) ---
        segment1 = untrim_segments[i - 1]['segment']
        segment2 = untrim_segments[i]['segment']

        # Find intersections between line pairs
        intersection = []
        for j in range(segment1.shape[0] - 1):
            p1_start = segment1[j, :]
            p1_end = segment1[j + 1, :]
            dist1 = np.linalg.norm(p1_end - p1_start)

            for jj in range(segment2.shape[0] - 1):
                p2_start = segment2[jj, :]
                p2_end = segment2[jj + 1, :]
                dist2 = np.linalg.norm(p2_end - p2_start)

                start_points = np.vstack([p1_start, p2_start])
                end_points = np.vstack([p1_end, p2_end])
                point, dist = line_intersect_3d(start_points, end_points)

                # Valid intersection: both distances <= respective segment lengths
                if not (np.any(dist > dist1) or np.any(dist > dist2)) and np.any(~np.isnan(dist)):
                    intersection.append(point)

        # If no intersection, connect with arc
        if len(intersection) == 0:
            face_id_prev = untrim_segments[i - 1]['FaceID']
            face_id_curr = untrim_segments[i]['FaceID']

            # Mean normal of the two adjacent faces
            mean_normal = np.array([
                np.mean([normal[face_id_curr, 0], normal[face_id_prev, 0]]),
                np.mean([normal[face_id_curr, 1], normal[face_id_prev, 1]]),
                np.mean([normal[face_id_curr, 2], normal[face_id_prev, 2]]),
            ])

            # Shared endpoint (last point of previous original segment)
            seg_end = og_segment[i - 1]['segment'][-1, :]
            seg_end_x, seg_end_y, seg_end_z = seg_end

            # point1: offset along direction1 (normal of previous face)
            dir1 = normal[face_id_prev, :]
            a_xy = np.arctan2(dir1[1], dir1[0])
            a_z = np.arctan2(dir1[2], np.sqrt(dir1[0]**2 + dir1[1]**2))
            p1 = np.array([
                seg_end_x + (offset_distance * np.cos(a_z)) * np.cos(a_xy),
                seg_end_y + (offset_distance * np.cos(a_z)) * np.sin(a_xy),
                seg_end_z + offset_distance * np.sin(a_z),
            ])

            # point2: offset along mean_normal
            a_xy = np.arctan2(mean_normal[1], mean_normal[0])
            a_z = np.arctan2(mean_normal[2], np.sqrt(mean_normal[0]**2 + mean_normal[1]**2))
            p2 = np.array([
                seg_end_x + (offset_distance * np.cos(a_z)) * np.cos(a_xy),
                seg_end_y + (offset_distance * np.cos(a_z)) * np.sin(a_xy),
                seg_end_z + offset_distance * np.sin(a_z),
            ])

            # point3: offset along direction3 (normal of current face)
            dir3 = normal[face_id_curr, :]
            a_xy = np.arctan2(dir3[1], dir3[0])
            a_z = np.arctan2(dir3[2], np.sqrt(dir3[0]**2 + dir3[1]**2))
            p3 = np.array([
                seg_end_x + (offset_distance * np.cos(a_z)) * np.cos(a_xy),
                seg_end_y + (offset_distance * np.cos(a_z)) * np.sin(a_xy),
                seg_end_z + offset_distance * np.sin(a_z),
            ])

            # Build arc
            arc = p1.reshape(1, 3)  # start with point1

            if not np.allclose(p1, p2):
                temp_seg = np.column_stack([
                    np.linspace(p1[0], p2[0], 20),
                    np.linspace(p1[1], p2[1], 20),
                    np.linspace(p1[2], p2[2], 20),
                ])
                arc = np.vstack([arc[:-1, :], temp_seg])

            if not np.allclose(p1, p3) and not np.allclose(p2, p3):
                temp_seg = np.column_stack([
                    np.linspace(p2[0], p3[0], 20),
                    np.linspace(p2[1], p3[1], 20),
                    np.linspace(p2[2], p3[2], 20),
                ])
                arc = np.vstack([arc[:-1, :], temp_seg])

            # Prepend arc to current segment
            untrim_segments[i]['segment'] = np.vstack([arc, untrim_segments[i]['segment']])

    return untrim_segments
