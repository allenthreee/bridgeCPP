import numpy as np
from get_grid_id import get_grid_id


def _in_polyhedron(connectivity, vertex, point):
    """
    Check if a point is inside a 3D polyhedron (mesh).

    Uses the ray-casting method: count intersections of a ray from the
    point with the mesh faces. Odd count = inside.

    Parameters
    ----------
    connectivity : np.ndarray (F, 3)
        Face vertex indices.
    vertex : np.ndarray (V, 3)
        Vertex coordinates.
    point : np.ndarray (3,)
        Point to test.

    Returns
    -------
    inside : bool
    """
    # Simple implementation using ray casting along +X
    ray_dir = np.array([1.0, 0.0, 0.0])
    count = 0

    for face in connectivity:
        v0 = vertex[face[0], :]
        v1 = vertex[face[1], :]
        v2 = vertex[face[2], :]

        # Möller–Trumbore ray-triangle intersection
        edge1 = v1 - v0
        edge2 = v2 - v0
        h = np.cross(ray_dir, edge2)
        a = np.dot(edge1, h)

        if abs(a) < 1e-10:
            continue

        f = 1.0 / a
        s = point - v0
        u = f * np.dot(s, h)

        if u < 0.0 or u > 1.0:
            continue

        q = np.cross(s, edge1)
        v = f * np.dot(ray_dir, q)

        if v < 0.0 or u + v > 1.0:
            continue

        t = f * np.dot(edge2, q)
        if t > 1e-10:
            count += 1

    return count % 2 == 1


def _calculate_new_coordinate(old_point, vector, point_list, offset_distance, ground_lvl):
    """
    Calculate new coordinate by stepping along a direction while
    satisfying collision constraints. Uses binary-search over step size
    (fast, no scipy.optimize dependency).

    Parameters
    ----------
    old_point : np.ndarray (3,)
        Current position.
    vector : np.ndarray (3,)
        Direction vector.
    point_list : np.ndarray (M, 3)
        Collision check points.
    offset_distance : float
        Minimum allowed distance.
    ground_lvl : float
        Ground level (minimum Z).

    Returns
    -------
    new_coordinate : np.ndarray (3,)
    step_size : float
    """
    vector = vector / (np.linalg.norm(vector) + 1e-12)

    def _is_valid(step):
        displacement = vector * step
        new_pt = old_point + displacement
        if new_pt[2] < ground_lvl:
            return False
        if point_list.size > 0:
            dists = np.linalg.norm(point_list - new_pt, axis=1)
            if np.min(dists) < offset_distance:
                return False
        return True

    # Binary search for max valid step size
    lo, hi = 0.0, offset_distance * 4.0  # search up to 4x offset
    for _ in range(12):
        mid = (lo + hi) / 2.0
        if _is_valid(mid):
            lo = mid
        else:
            hi = mid

    step_size = lo
    if step_size < offset_distance * 0.1:
        step_size = offset_distance * 0.1  # minimum step

    displacement = vector * step_size
    new_coordinate = old_point + displacement
    return new_coordinate, step_size


def _distance_changed(old_point, vector, step_size):
    """Distance from old_point when moving along vector by step_size."""
    vector = vector / (np.linalg.norm(vector) + 1e-12)
    displacement = vector * step_size
    new_point = old_point + displacement
    return np.linalg.norm(new_point - old_point)


def pointwise_replan(old_point, vector, offset_distance, ground_lvl, point_list,
                     colli_type, vertex, connectivity, all_grid, grid_info,
                     ccentroid, normal):
    """
    Replan a waypoint to avoid collisions.

    Parameters
    ----------
    old_point : np.ndarray (3,)
        Starting waypoint.
    vector : np.ndarray (K, 3)
        Direction vectors: vector[0] = primary escape direction.
    offset_distance : float
        Safety distance.
    ground_lvl : float
        Ground clearance level.
    point_list : np.ndarray (M, 3)
        Collision check points.
    colli_type : int
        1 = ground collision, 2 = mesh collision.
    vertex : np.ndarray (V, 3)
        Mesh vertices.
    connectivity : np.ndarray (F, 3)
        Mesh face indices.
    all_grid : list of dict
        Grid spatial index.
    grid_info : np.ndarray (7,)
        Grid bounds and size.
    ccentroid : np.ndarray (N, 3)
        Face centroids.
    normal : np.ndarray (N, 3)
        Face normals.

    Returns
    -------
    update_point : np.ndarray (3,)
        Replanned waypoint.
    """
    q1_bound = grid_info[:3]
    q3_bound = grid_info[3:6]
    grid_size = grid_info[6]

    prim_dir = vector[0, :]
    second_dir = np.sum(vector[1:, :], axis=0)
    if np.linalg.norm(second_dir) > 1e-12:
        second_dir = second_dir / np.linalg.norm(second_dir)

    if colli_type == 1:
        # Ground collision
        closest_point = np.array([old_point[0], old_point[1], 0.0])
        update_point, _ = _calculate_new_coordinate(
            old_point, prim_dir, closest_point.reshape(1, 3), 0, ground_lvl,
        )
        no_of_iter = 1

        while True:
            print(f"Num of Iteration in Replanning = {no_of_iter}")

            if update_point[2] < ground_lvl:
                pass  # MATLAB: GIGA = 1 marker, continue anyway

            grid_id, neighbour = get_grid_id(
                update_point, grid_size,
                q3_bound[0], q1_bound[0],
                q3_bound[1], q1_bound[1],
                q3_bound[2], q1_bound[2],
                3, True,
            )

            part1 = all_grid[grid_id - 1]['point']
            point_id_p1 = []
            if part1.size > 0:
                point_id_p1 = part1[:, -1].astype(int)
                part1_pts = part1[:, :3]
            else:
                part1_pts = np.empty((0, 3))

            part2_list = []
            for k in range(len(neighbour)):
                p = all_grid[neighbour[k] - 1]['point']
                if p.size > 0:
                    part2_list.append(p)
            if part2_list:
                part2 = np.vstack(part2_list)
                point_id_p2 = part2[:, -1].astype(int)
                part2_pts = part2[:, :3]
                point_ids = np.concatenate([point_id_p1, point_id_p2])
            else:
                part2_pts = np.empty((0, 3))
                point_ids = point_id_p1

            neighbour_centroids = np.vstack([part1_pts, part2_pts]) if part1_pts.size > 0 or part2_pts.size > 0 else np.empty((0, 3))

            collided_point = []
            dist_list = []
            for j in range(neighbour_centroids.shape[0]):
                temp_c = neighbour_centroids[j, :]
                temp_dist = np.linalg.norm(update_point - temp_c)
                if np.round(temp_dist - offset_distance, 1) < 0:
                    collided_point.append(temp_c)
                    dist_list.append(temp_dist)

            if len(dist_list) > 0:
                dist_list = np.array(dist_list)
                closest_idx = np.argmin(dist_list)
                closest = dist_list[closest_idx]
                need_replan = True
            else:
                closest = None
                need_replan = False

            if no_of_iter > 20:
                need_replan = False
                if len(collided_point) > 0:
                    pass  # ERROR condition in original

            if not need_replan:
                if update_point[2] < ground_lvl:
                    pass  # MATLAB: GIGA = 1 marker, continue anyway
                break

            # Compute weighted escape direction from collided points
            collided_point = np.array(collided_point)
            dist_list = np.array(dist_list)
            normal_list = []
            weight = []
            weighted_normal = []

            elite = max(int(np.floor(collided_point.shape[0] * 0.1)), 1)

            for jj in range(collided_point.shape[0]):
                temp_normal = (update_point - collided_point[jj, :])
                temp_normal = temp_normal / (np.linalg.norm(temp_normal) + 1e-12)
                normal_list.append(temp_normal)
                temp_ratio = dist_list[jj] / np.sum(dist_list)
                temp_weight = 1.0 / temp_ratio if temp_ratio > 0 else 1.0
                weight.append(temp_weight)
                weighted_normal.append(temp_weight * temp_normal)

            normal_list = np.array(normal_list)
            weight = np.array(weight)
            weighted_normal = np.array(weighted_normal)

            sorted_idx = np.argsort(weight)[::-1]
            sorted_normal = normal_list[sorted_idx, :]
            sorted_weight = weight[sorted_idx]

            elite_weighted_normal = np.zeros(3)
            for jj in range(elite):
                elite_weighted_normal += sorted_weight[jj] * sorted_normal[jj, :]
            elite_norm = np.linalg.norm(elite_weighted_normal)
            if elite_norm > 1e-12:
                elite_weighted_normal = elite_weighted_normal / elite_norm

            escape_dir = elite_weighted_normal
            if escape_dir[2] < 0:
                escape_dir[2] = 0.0

            # Check if inside polyhedron (debug)
            # if _in_polyhedron(connectivity, vertex, update_point):
            #     pass

            prev_replan_point = update_point.copy()
            update_point, _ = _calculate_new_coordinate(
                update_point, escape_dir,
                neighbour_centroids,
                offset_distance, ground_lvl,
            )
            no_of_iter += 1

        if update_point[2] < ground_lvl:
            pass  # GIGA = 1 in original

    elif colli_type == 2:
        # Mesh collision
        update_point, _ = _calculate_new_coordinate(
            old_point, prim_dir, point_list, offset_distance, ground_lvl,
        )
        no_of_iter = 1

        while True:
            print(f"No of Iteration in Replanning = {no_of_iter}")

            grid_id, neighbour = get_grid_id(
                update_point, grid_size,
                q3_bound[0], q1_bound[0],
                q3_bound[1], q1_bound[1],
                q3_bound[2], q1_bound[2],
                3, True,
            )

            if grid_id < 0 or np.isnan(grid_id) or grid_id == 0:
                raise ValueError("Invalid grid_id")

            part1 = all_grid[grid_id - 1]['point']
            point_id_p1 = []
            if part1.size > 0:
                point_id_p1 = part1[:, -1].astype(int)
                part1_pts = part1[:, :3]
            else:
                part1_pts = np.empty((0, 3))

            part2_list = []
            for k in range(len(neighbour)):
                p = all_grid[neighbour[k] - 1]['point']
                if p.size > 0:
                    part2_list.append(p)
            if part2_list:
                part2 = np.vstack(part2_list)
                point_id_p2 = part2[:, -1].astype(int)
                part2_pts = part2[:, :3]
                point_ids = np.concatenate([point_id_p1, point_id_p2])
            else:
                part2_pts = np.empty((0, 3))
                point_ids = point_id_p1

            neighbour_centroids = np.vstack([part1_pts, part2_pts]) if part1_pts.size > 0 or part2_pts.size > 0 else np.empty((0, 3))

            collided_point = []
            dist_list = []
            for j in range(neighbour_centroids.shape[0]):
                temp_c = neighbour_centroids[j, :]
                temp_dist = np.linalg.norm(update_point - temp_c)
                if np.round(temp_dist - offset_distance, 1) < 0:
                    collided_point.append(temp_c)
                    dist_list.append(temp_dist)

            if len(dist_list) > 0:
                dist_list = np.array(dist_list)
                closest_idx = np.argmin(dist_list)
                closest = dist_list[closest_idx]
                need_replan = True
            else:
                closest = None
                need_replan = False

            if no_of_iter > 30:
                need_replan = False
                if len(collided_point) > 0:
                    pass

            if not need_replan:
                if update_point[2] < ground_lvl:
                    pass  # MATLAB: GIGA = 1 marker, continue anyway
                break

            # Initial escape direction from closest collided point
            if no_of_iter == 1:
                collided_arr = np.array(collided_point)
                escape_dir = (update_point - collided_arr[closest_idx, :])
                escape_dir = escape_dir / (np.linalg.norm(escape_dir) + 1e-12)

            if no_of_iter > 1:
                collided_point = np.array(collided_point)
                dist_list = np.array(dist_list)
                normal_list = []
                weighted_normal = []
                elite_weighted_normal = []
                weight = []

                multiplier = (no_of_iter // 5) + 1
                elite = max(int(np.floor(collided_point.shape[0] * multiplier * 0.1)), 1)

                for jj in range(collided_point.shape[0]):
                    temp_normal = (update_point - collided_point[jj, :])
                    temp_normal = temp_normal / (np.linalg.norm(temp_normal) + 1e-12)
                    normal_list.append(temp_normal)
                    temp_ratio = dist_list[jj] / np.sum(dist_list)
                    temp_weight = 1.0 / temp_ratio if temp_ratio > 0 else 1.0
                    weight.append(temp_weight)
                    weighted_normal.append(temp_weight * temp_normal)

                normal_list = np.array(normal_list)
                weight = np.array(weight)
                weighted_normal = np.array(weighted_normal)

                if weight.size == 0:
                    raise ValueError("Empty weight list")

                sorted_idx = np.argsort(weight)[::-1]
                sorted_normal = normal_list[sorted_idx, :]
                sorted_weight = weight[sorted_idx]

                elite_weighted_normal_accum = np.zeros(3)
                for jj in range(min(elite, sorted_weight.size)):
                    elite_weighted_normal_accum += sorted_weight[jj] * sorted_normal[jj, :]
                elite_norm = np.linalg.norm(elite_weighted_normal_accum)
                if elite_norm > 1e-12:
                    elite_weighted_normal_accum = elite_weighted_normal_accum / elite_norm

                escape_dir = elite_weighted_normal_accum
                if escape_dir.shape[0] < 3:
                    raise ValueError("Escape direction invalid size")
                if update_point[2] < ground_lvl:
                    escape_dir[2] = abs(escape_dir[2])

            # if _in_polyhedron(connectivity, vertex, update_point):
            #     pass

            prev_replan_point = update_point.copy()
            update_point, _ = _calculate_new_coordinate(
                update_point, escape_dir,
                neighbour_centroids,
                offset_distance, ground_lvl,
            )
            no_of_iter += 1

    return update_point
