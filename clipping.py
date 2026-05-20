import numpy as np
from get_grid_id import get_grid_id
from pointwise_replan import pointwise_replan


def clipping(untrim_segment, offset_distance, ccentroid, unclip_traj,
             connectivity, vertex, normal, all_grid, grid_info):
    """
    Detect collisions, replan waypoints, and return clipped segments.

    Parameters
    ----------
    untrim_segment : list of dict
        Each dict has keys:
        - 'segment': np.ndarray (M, 3)
        - 'FaceID': int
        Output is mutated in place with additional keys:
        - 'collisionID': list
        - 'GroundClear': list
        - 'closest': list of dict with keys 'id' and 'dist'
        - 'collisionFace': list of dict with key 'id'
    offset_distance : float
    ccentroid : np.ndarray (N, 3)
    unclip_traj : np.ndarray (T, 3)
    connectivity : np.ndarray (F, 3)
    vertex : np.ndarray (V, 3)
    normal : np.ndarray (N, 3)
    all_grid : list of dict
    grid_info : np.ndarray (7,)

    Returns
    -------
    untrim_segment : list of dict (mutated)
    replan_segment : list of dict
    replan_segment_end : list of dict
    """
    q1_bound = grid_info[:3]
    q3_bound = grid_info[3:6]
    grid_size = grid_info[6]
    ground_lvl = 1.0

    closest_pool_size = max(int(np.round(ccentroid.shape[0] / 10)), 1)

    # ---- Phase 1: Collision detection ----
    for i in range(len(untrim_segment)):
        temp_segment = untrim_segment[i]['segment']
        untrim_segment[i].setdefault('collisionID', [])
        untrim_segment[i].setdefault('GroundClear', [])
        untrim_segment[i]['closest'] = []

        for ii in range(temp_segment.shape[0]):
            temp_point = temp_segment[ii, :]

            # Initialize closest storage
            untrim_segment[i]['closest'].append({
                'id': np.zeros(closest_pool_size, dtype=int),
                'dist': np.full(closest_pool_size, np.inf),
            })

            if temp_point[2] < ground_lvl:
                untrim_segment[i]['GroundClear'].append(ii)

            # Grid-based collision check
            grid_id, neighbour = get_grid_id(
                temp_point, grid_size,
                q3_bound[0], q1_bound[0],
                q3_bound[1], q1_bound[1],
                q3_bound[2], q1_bound[2],
                3, True,
            )

            part1 = all_grid[grid_id - 1]['point']
            point_id_p1 = np.array([], dtype=int)
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

            for j in range(neighbour_centroids.shape[0]):
                temp_c = neighbour_centroids[j, :]
                temp_dist = np.linalg.norm(temp_point - temp_c)
                centroid_id = point_ids[j] if j < len(point_ids) else 0

                # Track closest centroids
                closest_rec = untrim_segment[i]['closest'][ii]
                if np.any(np.isinf(closest_rec['dist'])):
                    inf_idx = np.where(np.isinf(closest_rec['dist']))[0][0]
                    closest_rec['dist'][inf_idx] = temp_dist
                    closest_rec['id'][inf_idx] = centroid_id
                else:
                    mask = temp_dist < closest_rec['dist']
                    if np.any(mask):
                        idx = np.where(mask)[0][0]
                        closest_rec['dist'][idx] = temp_dist
                        closest_rec['id'][idx] = centroid_id

                # Collision check
                if temp_dist < offset_distance:
                    if f'collisionFace_{ii}' not in untrim_segment[i]:
                        untrim_segment[i][f'collisionFace_{ii}'] = {'id': []}
                    untrim_segment[i][f'collisionFace_{ii}']['id'].append(centroid_id)
                    if ii not in untrim_segment[i]['collisionID']:
                        untrim_segment[i]['collisionID'].append(ii)

        # Clean up closest records
        for ii in range(temp_segment.shape[0]):
            rec = untrim_segment[i]['closest'][ii]
            keep = rec['id'] != 0
            rec['id'] = rec['id'][keep]
            rec['dist'] = rec['dist'][keep]

        untrim_segment[i]['GroundClear'] = list(set(untrim_segment[i].get('GroundClear', [])))
        untrim_segment[i]['collisionID'] = list(set(untrim_segment[i].get('collisionID', [])))
        untrim_segment[i]['collisionID'].sort()
        untrim_segment[i]['GroundClear'].sort()

    # ---- Phase 2: Replanning ----
    replan_seg_id = 0
    seg_cont = False
    is_start = False
    prev_seg = 0
    insertion_pt = 0
    insertion_end = 0
    temp_replan_seg = []
    temp_replan_end = []
    ground_clear_idlist = []
    collision_idlist = []
    segment_idlist = []
    replan_segment = []
    replan_segment_end = []

    for i in range(len(untrim_segment)):
        print(f'Current Segment = {i + 1}')

        has_ground = len(untrim_segment[i].get('GroundClear', [])) > 0
        has_collision = len(untrim_segment[i].get('collisionID', [])) > 0

        if not has_ground and not has_collision:
            continue

        # ---- Ground clearance replanning ----
        if has_ground:
            # Determine start/end/continuation flags
            if i != 0:
                if not seg_cont:
                    is_start = True
                else:
                    is_start = False
            else:
                is_start = True

            gc_list = untrim_segment[i]['GroundClear']
            if gc_list[-1] == untrim_segment[i]['segment'].shape[0] - 1:  # collision till end
                if i != len(untrim_segment) - 1:
                    next_gc = untrim_segment[i + 1].get('GroundClear', [])
                    if len(next_gc) == 0:
                        end_here = True
                    elif next_gc[0] == 0:
                        seg_cont = True
                        end_here = False
                    else:
                        seg_cont = False
                        end_here = True
                else:
                    seg_cont = False
                    end_here = True
            else:
                seg_cont = False
                end_here = True

            for ii_idx, ii in enumerate(gc_list):
                point_id = ii
                if ii_idx < len(gc_list) - 1:
                    if gc_list[ii_idx + 1] != point_id + 1:
                        is_start = True
                        seg_cont = False
                        end_here = True

                old_point = untrim_segment[i]['segment'][point_id, :]

                # Build vector_list from closest faces
                closest_rec = untrim_segment[i]['closest'][point_id]
                vector_list = np.zeros((len(closest_rec['id']), 3))
                collision_id = []
                for j in range(len(closest_rec['id'])):
                    face_id = closest_rec['id'][j]
                    temp_dist = np.linalg.norm(old_point - ccentroid[face_id, :])
                    if temp_dist < offset_distance:
                        collision_id.append(j)
                    temp_vector = normal[face_id, :]
                    temp_vector = temp_vector * (1.0 / max(closest_rec['dist'][j], 1e-12))
                    vector_list[j, :] = temp_vector

                if len(collision_id) > 0:
                    vector_list = vector_list[collision_id, :]

                if vector_list.shape[0] < 2:
                    test_sum = vector_list / (np.linalg.norm(vector_list) + 1e-12)
                else:
                    s = np.sum(vector_list, axis=0)
                    test_sum = s / (np.linalg.norm(s) + 1e-12)

                test_sum[2] = abs(test_sum[2])
                sum_vector = test_sum / (np.linalg.norm(test_sum) + 1e-12)

                # Gather closest points for collision check
                cid_list = closest_rec['id']
                closest_centroid_list = ccentroid[cid_list, :]
                closest_vertex_list = connectivity[cid_list, :].reshape(-1)
                closest_vertex_list = vertex[closest_vertex_list, :]
                closest_pointlist = np.vstack([closest_centroid_list, closest_vertex_list])

                giga_vector = np.vstack([np.array([0, 0, 1.0]), vector_list])

                replanned = pointwise_replan(
                    old_point, giga_vector, offset_distance, ground_lvl,
                    closest_pointlist, 1, vertex, connectivity,
                    all_grid, grid_info, ccentroid, normal,
                )

                temp_replan_seg.append(replanned)
                ground_clear_idlist.append(ii_idx)
                segment_idlist.append(point_id)

            temp_replan_end.append(replanned)

            if is_start:
                prev_seg = i
                insertion_pt = segment_idlist[0] - 1
                if end_here:
                    insertion_end = segment_idlist[0]

            if end_here:
                next_seg = i
                replan_segment.append({
                    'segment': np.array(temp_replan_seg),
                    'Head': prev_seg,
                    'Tail': next_seg,
                })
                replan_segment_end.append({
                    'segment': np.array(temp_replan_end),
                    'Head': prev_seg,
                    'Tail': next_seg,
                    'insertion_pt': insertion_pt,
                    'insertion_end': insertion_end,
                })
                temp_replan_seg = []
                temp_replan_end = []
                replan_seg_id += 1

            # Remove replanned points from original segment
            untrim_segment[i]['GroundClear'] = [
                gc for idx, gc in enumerate(gc_list)
                if idx not in ground_clear_idlist
            ]
            untrim_segment[i]['segment'] = np.delete(
                untrim_segment[i]['segment'], segment_idlist, axis=0,
            )
            # Adjust remaining index references after deletion
            for d in sorted(segment_idlist, reverse=True):
                untrim_segment[i]['collisionID'] = [
                    c - 1 if c > d else c
                    for c in untrim_segment[i].get('collisionID', [])
                ]
            ground_clear_idlist = []
            segment_idlist = []

        # ---- Mesh collision replanning ----
        if has_collision:
            coll_list = untrim_segment[i]['collisionID']

            # Handle local segments (discontinuous collision groups)
            if len(coll_list) > 1:
                diffs = np.diff(coll_list)
                breakpoints = np.where(diffs > 1)[0]
                no_local_seg = len(breakpoints) + 1

                local_seg = []
                if no_local_seg == 1:
                    local_seg.append({'segment': coll_list})
                else:
                    local_seg_id = 0
                    temp_pointer = 0
                    for ii in range(no_local_seg - 1):
                        bp = breakpoints[ii]
                        local_seg.append({'segment': coll_list[temp_pointer:bp + 1]})
                        temp_pointer = bp + 1
                        if ii == no_local_seg - 2:
                            local_seg.append({'segment': coll_list[temp_pointer:]})
            else:
                local_seg = [{'segment': coll_list}]

            for seg_idx, ls in enumerate(local_seg):
                # Determine start/end/cont for this local segment
                if len(local_seg) > 1:
                    if seg_idx == 0 and seg_idx < len(local_seg) - 1:
                        seg_cont_l = False
                        end_here_l = True
                    elif seg_idx == len(local_seg) - 1:
                        if coll_list[-1] == untrim_segment[i]['segment'].shape[0] - 1:
                            if i != len(untrim_segment) - 1:
                                next_coll = untrim_segment[i + 1].get('collisionID', [])
                                if len(next_coll) == 0:
                                    is_start = True
                                    end_here_l = True
                                    seg_cont_l = False
                                elif len(next_coll) > 0 and next_coll[0] == 0:
                                    is_start = True
                                    seg_cont_l = True
                                    end_here_l = False
                                else:
                                    is_start = True
                                    seg_cont_l = False
                                    end_here_l = True
                            else:
                                is_start = True
                                seg_cont_l = False
                                end_here_l = True
                        else:
                            is_start = True
                            seg_cont_l = False
                            end_here_l = True
                    elif seg_idx == 0:
                        seg_cont_l = False
                        end_here_l = True
                    else:
                        seg_cont_l = False
                        end_here_l = True
                else:
                    seg_cont_l = False
                    end_here_l = True
                    is_start = True

                for ii_idx, ii in enumerate(ls['segment']):
                    point_id = ii

                    # Build vector_list from closest 10%
                    closest_rec = untrim_segment[i]['closest'][point_id]
                    vector_list = np.zeros((len(closest_rec['id']), 3))
                    for j in range(len(closest_rec['id'])):
                        face_id = closest_rec['id'][j]
                        temp_vector = normal[face_id, :]
                        temp_vector = temp_vector * (1.0 / max(closest_rec['dist'][j], 1e-12))
                        vector_list[j, :] = temp_vector

                    # Find colliding faces and escaping normal
                    coll_face_key = f'collisionFace_{point_id}'
                    coll_face_ids = untrim_segment[i].get(coll_face_key, {}).get('id', [])

                    escaping_normal = []
                    temp_vector_list = []
                    closest_collide_pt = None
                    closest_collide_dist = np.inf
                    old_point = untrim_segment[i]['segment'][point_id, :]

                    for j, face_id in enumerate(coll_face_ids):
                        temp_c = ccentroid[face_id, :]
                        temp_n = (old_point - temp_c)
                        temp_n = temp_n / (np.linalg.norm(temp_n) + 1e-12)
                        escaping_normal.append(temp_n)
                        temp_vector_list.append(normal[face_id, :])
                        temp_dist = np.linalg.norm(old_point - temp_c)
                        if closest_collide_pt is None or temp_dist < closest_collide_dist:
                            closest_collide_pt = face_id
                            closest_collide_dist = temp_dist

                    if closest_collide_pt is not None:
                        closest_collide_normal = normal[closest_collide_pt, :]
                        # Filter out faces on opposite side (> 135 deg)
                        opposite_side_id = []
                        for j in range(len(temp_vector_list)):
                            dot_p = np.dot(closest_collide_normal, temp_vector_list[j])
                            dot_p = np.rad2deg(np.arccos(
                                dot_p / (np.linalg.norm(closest_collide_normal) *
                                         np.linalg.norm(temp_vector_list[j]) + 1e-12)
                            ))
                            if dot_p > 135:
                                opposite_side_id.append(j)
                        for idx in sorted(opposite_side_id, reverse=True):
                            temp_vector_list.pop(idx)
                            if idx < len(escaping_normal):
                                escaping_normal.pop(idx)

                    if len(temp_vector_list) > 0:
                        temp_vector = np.sum(temp_vector_list, axis=0)
                        temp_vector = temp_vector / (np.linalg.norm(temp_vector) + 1e-12)
                    else:
                        temp_vector = np.zeros(3)

                    if len(escaping_normal) > 0:
                        escaping_normal_arr = np.sum(escaping_normal, axis=0)
                        escaping_normal_arr = escaping_normal_arr / (np.linalg.norm(escaping_normal_arr) + 1e-12)
                    else:
                        escaping_normal_arr = np.zeros(3)

                    # Gather closest points
                    cid_list = closest_rec['id']
                    closest_centroid_list = ccentroid[cid_list, :]
                    closest_vertex_list = connectivity[cid_list, :].reshape(-1)
                    closest_vertex_list = vertex[closest_vertex_list, :]
                    closest_pointlist = np.vstack([closest_centroid_list, closest_vertex_list])

                    giga_vector = np.vstack([
                        escaping_normal_arr.reshape(1, 3),
                        vector_list,
                    ])

                    replanned = pointwise_replan(
                        old_point, giga_vector, offset_distance, ground_lvl,
                        closest_pointlist, 2, vertex, connectivity,
                        all_grid, grid_info, ccentroid, normal,
                    )

                    temp_replan_seg.append(replanned)
                    collision_idlist.append(ii_idx)
                    segment_idlist.append(point_id)

                temp_replan_end.append(replanned)

                if is_start:
                    prev_seg = i
                    insertion_pt = segment_idlist[0] - 1
                    if end_here_l:
                        insertion_end = segment_idlist[0]

                if end_here_l:
                    next_seg = i
                    replan_segment.append({
                        'segment': np.array(temp_replan_seg),
                        'Head': prev_seg,
                        'Tail': next_seg,
                    })
                    replan_segment_end.append({
                        'segment': np.array(temp_replan_end),
                        'Head': prev_seg,
                        'Tail': next_seg,
                        'insertion_pt': insertion_pt if is_start else 0,
                        'insertion_end': insertion_end if is_start else segment_idlist[0],
                    })
                    temp_replan_seg = []
                    temp_replan_end = []
                    replan_seg_id += 1

                # Remove replanned points from original
                for idx in sorted(collision_idlist, reverse=True):
                    if idx < len(untrim_segment[i]['collisionID']):
                        del untrim_segment[i]['collisionID'][idx]
                for idx in sorted(segment_idlist, reverse=True):
                    untrim_segment[i]['segment'] = np.delete(
                        untrim_segment[i]['segment'], idx, axis=0,
                    )
                # Adjust remaining indices in coll_list and local_seg after deletion
                n_deleted = len(segment_idlist)
                if n_deleted > 0 and seg_idx < len(local_seg) - 1:
                    for kk in range(seg_idx + 1, len(local_seg)):
                        local_seg[kk]['segment'] = [s - n_deleted for s in local_seg[kk]['segment']]
                    coll_list = [c - n_deleted for c in coll_list]
                    untrim_segment[i]['collisionID'] = coll_list
                collision_idlist = []
                segment_idlist = []

    return untrim_segment, replan_segment, replan_segment_end
