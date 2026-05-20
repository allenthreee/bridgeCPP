"""
main.py — Canning parameterization-based inspection path planning.

Converted from main.m. Reads a mesh model (pillar default), generates
a 2D sweep path in UV space, maps it to 3D, offsets, clips, replans,
and outputs waypoints with camera commands.
"""

import os
import sys
import time

# ── Fix mpl_toolkits: system has old 3.5.1 nspkg.pth that injects broken paths ──
# Clear any stale mpl_toolkits cached by the system .pth file
for _k in list(sys.modules):
    if 'mpl_toolkits' in _k:
        del sys.modules[_k]

import numpy as np
from scipy.io import loadmat
from scipy.stats import norm as stats_norm
import matplotlib.pyplot as plt

# Local imports
from ratio_compute import ratio_compute
from get_grid_id import get_grid_id
from line_intersect_3d import line_intersect_3d
from untrim_offset_3d import untrim_offset_3d
from clipping import clipping
from pointwise_replan import pointwise_replan
from calculate_vist_3d import calculate_vist_3d
from cal_point_vist import cal_point_vist
from detect_large_angle import detect_large_angle_of_incidence
from two_opt import two_opt
from hermite_plot import hermite_plot


# ────────────────── Mesh Triangulation ──────────────────
class MeshTri:
    """Lightweight triangulation holding connectivity and vertex coordinates."""
    def __init__(self, simplices, points):
        self.simplices = simplices
        self.points = points


# ────────────────── Barycentric helpers ──────────────────
def _barycentric_to_cartesian(tri, face_id, bary):
    """
    Convert barycentric coordinates to Cartesian on a triangulation.

    Parameters
    ----------
    tri : Delaunay or similar (with `simplices` and `points`)
    face_id : int
        Face index.
    bary : np.ndarray (3,)
        Barycentric coordinates (sum to 1).

    Returns
    -------
    point : np.ndarray (3,)
    """
    simplex = tri.simplices[face_id]
    pts = tri.points[simplex]
    return bary @ pts


def _cartesian_to_barycentric(tri, face_id, point):
    """
    Convert a Cartesian point to barycentric coordinates on a face.

    Parameters
    ----------
    tri : Delaunay or similar
    face_id : int
    point : np.ndarray (2,)

    Returns
    -------
    bary : np.ndarray (3,)
    """
    simplex = tri.simplices[face_id]
    pts = tri.points[simplex]  # 3x2
    # Solve: [pts.T; 1 1 1] * [b0 b1 b2]^T = [point; 1]
    A = np.vstack([pts.T, np.ones(3)])
    b = np.append(point, 1.0)
    return np.linalg.solve(A, b)


# ────────────────── Main ──────────────────
def main():
    t_start = time.time()
    model = 'pillar'

    # Data paths
    base_dir = os.path.dirname(os.path.abspath(__file__))

    if model == 'pillar':
        mapping_file = os.path.join(base_dir, 'pillar', 'pillar_maping.txt')
        faces_file = os.path.join(base_dir, 'pillar', 'pillar_faces.txt')
        normal_file = os.path.join(base_dir, 'pillar', 'pillar_normal.txt')
        mat_path = os.path.join(base_dir, 'pillar', 'pillar_path.mat')
    else:
        raise ValueError(f"Unknown model: {model}")

    # ── Load and parse data ──
    with open(mapping_file, 'r') as f:
        lines = [line.strip().replace('(', '').replace(')', '').replace(',', '')
                 .replace('<', '').replace('>', '') for line in f]
    file_data = [line.split() for line in lines if line]
    # file_data format: Col 1=VertexID, 2=Vector, 3=Vx, 4=Vy, 5=Vz, 6=Vector, 7=UVx, 8=UVy

    with open(faces_file, 'r') as f:
        lines2 = [line.strip().replace(':', '') for line in f]
    file2_data = [line.split() for line in lines2 if line]
    # file2 format: Col 1=Face, 2=FaceID, 3=V1, 4=V2, 5=V3

    with open(normal_file, 'r') as f:
        lines3 = [line.strip().replace('(', '').replace(')', '').replace(',', '')
                  .replace('<', '').replace('>', '') for line in f]
    file3_data = [line.split() for line in lines3 if line]

    # Parse arrays
    Id = np.array([float(row[0]) for row in file_data])
    vertex_X = np.array([float(row[2]) for row in file_data])
    vertex_Y = np.array([float(row[3]) for row in file_data])
    vertex_Z = np.array([float(row[4]) for row in file_data])
    if model == 'pillar':
        Vertex = np.column_stack([vertex_X / 7.29, vertex_Y / 7.29, vertex_Z / 7.29])
    else:
        Vertex = np.column_stack([vertex_X, vertex_Y, vertex_Z])
    uv_x = np.array([float(row[6]) for row in file_data])
    uv_y = np.array([float(row[7]) for row in file_data])
    UV = np.column_stack([uv_x, uv_y])

    FaceID = np.array([int(float(row[1])) for row in file2_data])
    Faces = np.array([[int(float(row[2])), int(float(row[3])), int(float(row[4]))]
                      for row in file2_data])

    Normal = np.column_stack([
        np.array([float(row[1]) for row in file3_data]),
        np.array([float(row[2]) for row in file3_data]),
        np.array([float(row[3]) for row in file3_data]),
    ])

    # ── Build connectivity ──
    max_vertex_id = int(Faces.max())
    count_list = np.zeros(max_vertex_id + 1, dtype=int)
    connectivity = np.zeros_like(Faces, dtype=int)

    for i in range(Faces.shape[0]):
        for ii in range(Faces.shape[1]):
            times = count_list[Faces[i, ii]]
            # Find vertex_id matching this Face vertex
            vertex_ids = np.where(Id == Faces[i, ii])[0]
            vertex_ids_sorted = vertex_ids[times]
            count_list[Faces[i, ii]] += 1
            connectivity[i, ii] = vertex_ids_sorted

    # ── UV triangulation ──
    # Plot UV map (uses connectivity directly, not Delaunay)
    fig1 = plt.figure(figsize=(12, 5))
    ax_uv = fig1.add_subplot(1, 2, 1)
    ax_3d = fig1.add_subplot(1, 2, 2, projection='3d')
    ax_uv.triplot(UV[:, 0], UV[:, 1], connectivity, color='k')
    ax_uv.set_aspect('equal')

    # 3D mesh
    ax_3d_trisurf = ax_3d.plot_trisurf(
        Vertex[:, 0], Vertex[:, 1], Vertex[:, 2],
        triangles=connectivity, color=[0.7, 0.7, 0.7],
    )
    ax_3d.set_aspect('equal')
    ax_3d.grid(True)
    ax_3d.set_xlabel('X (m)')
    ax_3d.set_ylabel('Y (m)')
    ax_3d.set_zlabel('Z (m)')
    plt.draw()

    # ── Load sweep path ──
    mat_data = loadmat(mat_path)
    # Try common key names
    if 'pillar_path' in mat_data:
        test_traj = mat_data['pillar_path']
    elif 'path' in mat_data:
        test_traj = mat_data['path']
    else:
        # Use first non-meta array
        for k, v in mat_data.items():
            if not k.startswith('__') and isinstance(v, np.ndarray) and v.ndim == 2:
                test_traj = v
                break
        else:
            raise ValueError("Could not find trajectory data in .mat file")

    # ── Discretization ──
    offset_distance = 1.5
    grid_size = offset_distance

    Q1_bound = grid_size - np.remainder(Vertex.max(axis=0), grid_size) + Vertex.max(axis=0) + 3 * grid_size
    Q3_bound = Vertex.min(axis=0) - (grid_size + np.remainder(Vertex.min(axis=0), grid_size)) - 3 * grid_size
    Q3_bound[2] = 0.0

    no_of_grid = (Q1_bound - Q3_bound) / grid_size
    print('Number of grid:', no_of_grid)

    grid_size_2d = 0.05
    u_upper, u_lower = 1.0, 0.0
    v_upper, v_lower = 1.0, 0.0
    no_grid_2d = int(((u_upper - u_lower) / grid_size_2d) ** 2)
    ALL_2Dgrid = [{'point': np.empty((0, 2))} for _ in range(no_grid_2d)]

    Grid_INFO = np.array([
        Q1_bound[0], Q1_bound[1], Q1_bound[2],
        Q3_bound[0], Q3_bound[1], Q3_bound[2],
        grid_size,
    ])

    no_grid = int(np.prod((Q1_bound - Q3_bound) / grid_size))
    ALL_grid = [{'point': np.empty((0, 4))} for _ in range(no_grid)]  # x,y,z,face_id

    FOV = 90.0
    viewingDist = offset_distance / np.cos(np.deg2rad(FOV / 2))
    Vist_radius3D = offset_distance * np.tan(np.deg2rad(FOV / 2)) * 0.5
    viewingAngle_threshold = 60.0

    # ── Compute centroids and populate grids ──
    CCentroid = np.zeros((connectivity.shape[0], 2))
    CCentroid_3D = np.zeros((connectivity.shape[0], 3))

    for i in range(connectivity.shape[0]):
        uv_c = UV[connectivity[i, :], :]
        dd_c = Vertex[connectivity[i, :], :]
        CCentroid[i, :] = np.mean(uv_c, axis=0)
        CCentroid_3D[i, :] = np.mean(dd_c, axis=0)

        grid_id_2d, _ = get_grid_id(CCentroid[i, :], grid_size_2d,
                                     u_lower, u_upper, v_lower, v_upper,
                                     0, 0, 2, False)
        if ALL_2Dgrid[grid_id_2d - 1]['point'].size == 0:
            ALL_2Dgrid[grid_id_2d - 1]['point'] = CCentroid[i, :].reshape(1, 2)
        else:
            ALL_2Dgrid[grid_id_2d - 1]['point'] = np.vstack([
                ALL_2Dgrid[grid_id_2d - 1]['point'], CCentroid[i, :],
            ])

        grid_id_3d, _ = get_grid_id(CCentroid_3D[i, :], grid_size,
                                     Q3_bound[0], Q1_bound[0],
                                     Q3_bound[1], Q1_bound[1],
                                     Q3_bound[2], Q1_bound[2],
                                     3, False)
        entry = np.array([CCentroid_3D[i, 0], CCentroid_3D[i, 1], CCentroid_3D[i, 2], i])
        if ALL_grid[grid_id_3d - 1]['point'].size == 0:
            ALL_grid[grid_id_3d - 1]['point'] = entry.reshape(1, 4)
        else:
            ALL_grid[grid_id_3d - 1]['point'] = np.vstack([
                ALL_grid[grid_id_3d - 1]['point'], entry,
            ])

    # ── UV-to-3D ratio ──
    uv_to_3d_ratio = ratio_compute(CCentroid, CCentroid_3D)
    print(f"UV to 3D unit ratio = {uv_to_3d_ratio}")

    uv_radius = Vist_radius3D * uv_to_3d_ratio

    # ── UV-to-3D mapping & visibility ──
    traj_3D = np.zeros((test_traj.shape[0], 3))
    traj_FaceID = np.zeros(test_traj.shape[0], dtype=int)
    Visibility_UV = []

    for i in range(test_traj.shape[0]):
        uv_traj = test_traj[i, :]
        grid_id_2d, neighbour_2d = get_grid_id(uv_traj, grid_size_2d,
                                                u_lower, u_upper, v_lower, v_upper,
                                                0, 0, 2, True)

        list_of_2d = ALL_2Dgrid[grid_id_2d - 1]['point']
        for nid in neighbour_2d:
            p = ALL_2Dgrid[nid - 1]['point']
            if p.size > 0:
                if list_of_2d.size == 0:
                    list_of_2d = p
                else:
                    list_of_2d = np.vstack([list_of_2d, p])

        # Remove zero rows
        if list_of_2d.size > 0:
            list_of_2d = list_of_2d[~np.all(list_of_2d == 0, axis=1)]

        visible = []
        closest = None
        best_dist = np.inf

        for ii in range(list_of_2d.shape[0]):
            temp_dist = np.linalg.norm(uv_traj - list_of_2d[ii, :])
            if temp_dist < best_dist:
                best_dist = temp_dist
                closest = list_of_2d[ii, :]
            if temp_dist < uv_radius:
                matches = np.where(np.all(np.abs(CCentroid - list_of_2d[ii, :]) < 1e-10, axis=1))[0]
                if len(matches) > 0:
                    visible.append(matches[0])

        if closest is not None:
            face_matches = np.where(np.all(np.abs(CCentroid - closest) < 1e-10, axis=1))[0]
            if len(face_matches) > 0:
                traj_FaceID[i] = face_matches[0]

        Visibility_UV.append({'Visibility': visible})

        # Barycentric mapping UV -> 3D
        if traj_FaceID[i] < connectivity.shape[0]:
            uv_tri = MeshTri(connectivity, UV)
            model_tri = MeshTri(connectivity, Vertex)
            bary_uv = _cartesian_to_barycentric(uv_tri, traj_FaceID[i], uv_traj)
            mesh_point = _barycentric_to_cartesian(model_tri, traj_FaceID[i], bary_uv)
            traj_3D[i, :] = mesh_point

    # Visible faces in UV
    UV_Vist = np.unique(np.concatenate([v['Visibility'] for v in Visibility_UV if v['Visibility']]))

    # Plot UV visibility (green edges for visible faces)
    if len(UV_Vist) > 0:
        ax_uv.triplot(UV[:, 0], UV[:, 1], connectivity[UV_Vist, :],
                      color='g', linewidth=0.5)
    ax_uv.plot(test_traj[:, 0], test_traj[:, 1], 'r', linewidth=2.5)
    ax_3d.plot(traj_3D[:, 0], traj_3D[:, 1], traj_3D[:, 2], 'r')
    plt.draw()

    # ── New figure for 3D mesh ──
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.plot_trisurf(Vertex[:, 0], Vertex[:, 1], Vertex[:, 2],
                     triangles=connectivity, color=[0.7, 0.7, 0.7])
    ax2.view_init(elev=35, azim=25)
    ax2.set_aspect('equal')
    ax2.grid(True)

    # ── Build OG segments ──
    OG_segment = []
    prev_end = 0
    for i in range(traj_FaceID.shape[0] - 1):
        current = traj_FaceID[i]
        next_f = traj_FaceID[i + 1]
        if current != next_f:
            OG_segment.append({
                'segment': traj_3D[prev_end:i + 1, :],
                'FaceID': int(current),
            })
            prev_end = i + 1
        elif i == traj_FaceID.shape[0] - 2:
            OG_segment.append({
                'segment': traj_3D[prev_end:, :],
                'FaceID': int(current),
            })

    # ── Offset ──
    untrim_segment = untrim_offset_3d(OG_segment, offset_distance, Normal, CCentroid_3D)
    unclip_traj = np.vstack([s['segment'] for s in untrim_segment])
    ax2.plot(unclip_traj[:, 0], unclip_traj[:, 1], unclip_traj[:, 2],
             'r', linewidth=2.5)
    plt.draw()

    # ── Clipping & Replanning ──
    for s in untrim_segment:
        s.setdefault('collisionID', [])
        s.setdefault('GroundClear', [])

    untrim_segment, replan_segment, replan_segment_end = clipping(
        untrim_segment, offset_distance, CCentroid_3D, unclip_traj,
        connectivity, Vertex, Normal, ALL_grid, Grid_INFO,
    )

    # Plot non-colliding segments green
    for s in untrim_segment:
        if len(s.get('collisionID', [])) == 0 and len(s.get('GroundClear', [])) == 0:
            ax2.plot(s['segment'][:, 0], s['segment'][:, 1], s['segment'][:, 2], 'g')

    # ── Reassemble trajectory ──
    replan_segment_HEADs = np.array([r['Head'] for r in replan_segment_end])
    replan_segment_TAILs = np.array([r['Tail'] for r in replan_segment_end])

    # Plot replanned segments
    fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(14, 6),
                                       subplot_kw={'projection': '3d'})
    for ax_s in [ax3a, ax3b]:
        ax_s.plot_trisurf(Vertex[:, 0], Vertex[:, 1], Vertex[:, 2],
                          triangles=connectivity, color=[0.7, 0.7, 0.7])
        ax_s.view_init(elev=35, azim=25)
        ax_s.set_aspect('equal')
        ax_s.grid(True)

    ax3a.plot(traj_3D[:, 0], traj_3D[:, 1], traj_3D[:, 2], 'b')
    for s in untrim_segment:
        ax3a.plot(s['segment'][:, 0], s['segment'][:, 1], s['segment'][:, 2], 'r')
    for r in replan_segment_end:
        ax3b.plot(r['segment'][:, 0], r['segment'][:, 1], r['segment'][:, 2], 'b')

    # Reassemble replan trajectory
    replan_traj = []
    i = 0
    while i < len(untrim_segment):
        if i in replan_segment_HEADs or i in replan_segment_TAILs:
            head_match = np.where(replan_segment_HEADs == i)[0]
            if len(head_match) > 0:
                for hidx in head_match:
                    r = replan_segment_end[hidx]
                    if r['Head'] == r['Tail']:
                        start_pt = r['insertion_pt']
                        end_pt = r['insertion_end']
                        seg_i = untrim_segment[i]['segment']
                        if start_pt == 0:
                            replan_traj.append(r['segment'])
                            if end_pt < seg_i.shape[0]:
                                replan_traj.append(seg_i[end_pt:, :])
                        else:
                            if start_pt > 0:
                                replan_traj.append(seg_i[:start_pt, :])
                            replan_traj.append(r['segment'])
                            if end_pt < seg_i.shape[0]:
                                replan_traj.append(seg_i[end_pt:, :])
                    else:
                        start_pt = r['insertion_pt']
                        seg_i = untrim_segment[i]['segment']
                        if start_pt == 0:
                            replan_traj.append(r['segment'])
                        else:
                            replan_traj.append(seg_i[:start_pt, :])
                            replan_traj.append(r['segment'])
            else:
                # Tail match
                tail_match = np.where(replan_segment_TAILs == i)[0]
                if len(tail_match) > 0:
                    r = replan_segment_end[tail_match[0]]
                    if r['Head'] != r['Tail']:
                        replan_traj.append(r['segment'])
                        replan_traj.append(untrim_segment[i]['segment'])
        else:
            replan_traj.append(untrim_segment[i]['segment'])
        i += 1

    if replan_traj:
        replan_traj = np.vstack([arr.reshape(-1, 3) for arr in replan_traj if arr.size > 0])
    else:
        replan_traj = np.empty((0, 3))

    ax3b.plot(replan_traj[:, 0], replan_traj[:, 1], replan_traj[:, 2], 'r')
    plt.draw()
    print("Replan Finished")

    # ── Down-sampling ──
    start_point = replan_traj[0, :].copy()
    end_point = replan_traj[-1, :].copy()
    down_sampled_traj = replan_traj.copy()
    down_sample_goal = 3 * (replan_traj.max() * 0.01)

    while True:
        distances = np.linalg.norm(np.diff(down_sampled_traj, axis=0), axis=1)
        threshold = np.mean(distances)
        i = 0
        down_sampled_ID = []
        temp_dist_acc = 0.0
        while i < len(distances):
            temp_dist_acc += distances[i]
            if temp_dist_acc > threshold:
                down_sampled_ID.append(i)
                temp_dist_acc = 0.0
            i += 1
        if len(down_sampled_ID) > 0:
            down_sampled_traj = down_sampled_traj[down_sampled_ID, :]
        if threshold >= down_sample_goal:
            break

    down_sampled_traj = np.vstack([start_point, down_sampled_traj, end_point])
    print("Starting Two-opt")
    down_sampled_traj = two_opt(down_sampled_traj, 20)

    # Plot down-sampled trajectory
    fig4 = plt.figure()
    ax4 = fig4.add_subplot(111, projection='3d')
    ax4.plot_trisurf(Vertex[:, 0], Vertex[:, 1], Vertex[:, 2],
                     triangles=connectivity, color=[0.7, 0.7, 0.7])
    ax4.view_init(elev=35, azim=25)
    ax4.plot(down_sampled_traj[:, 0], down_sampled_traj[:, 1],
             down_sampled_traj[:, 2], 'r')
    ax4.plot(Q1_bound[0], Q1_bound[1], Q1_bound[2], 'm*', markersize=15)
    ax4.plot(Q3_bound[0], Q3_bound[1], Q3_bound[2], 'm*', markersize=15)
    ax4.set_aspect('equal')
    ax4.grid(True)

    # ── Compute viewing distance and spray direction ──
    dist2surf = np.zeros(down_sampled_traj.shape[0])
    target_point = np.zeros_like(down_sampled_traj)
    spraying_dir = np.zeros_like(down_sampled_traj)

    for i in range(down_sampled_traj.shape[0]):
        temp_point = down_sampled_traj[i, :]
        grid_id, neighbour = get_grid_id(temp_point, grid_size,
                                          Q3_bound[0], Q1_bound[0],
                                          Q3_bound[1], Q1_bound[1],
                                          Q3_bound[2], Q1_bound[2],
                                          3, True)

        part1 = ALL_grid[grid_id - 1]['point']
        if part1.size > 0:
            part1_pts = part1[:, :3]
        else:
            part1_pts = np.empty((0, 3))

        part2_list = []
        for k in range(len(neighbour)):
            p = ALL_grid[neighbour[k] - 1]['point']
            if p.size > 0:
                part2_list.append(p)
        if part2_list:
            part2 = np.vstack(part2_list)
            part2_pts = part2[:, :3]
        else:
            part2_pts = np.empty((0, 3))

        neighbour_c = np.vstack([part1_pts, part2_pts]) if part1_pts.size > 0 or part2_pts.size > 0 else np.empty((0, 3))

        closest_point = None
        current_shortest = np.inf
        for j in range(neighbour_c.shape[0]):
            temp_c = neighbour_c[j, :]
            temp_dist = np.linalg.norm(temp_point - temp_c)
            if temp_dist < current_shortest:
                current_shortest = temp_dist
                closest_point = temp_c

        if closest_point is not None:
            vector_to_surface = temp_point - closest_point
            matches = np.where(np.all(np.abs(CCentroid_3D - closest_point) < 1e-10, axis=1))[0]
            if len(matches) > 0:
                temp_normal = Normal[matches[0], :]
                perpendicular_dist = abs(np.dot(vector_to_surface, temp_normal))
                dist2surf[i] = perpendicular_dist
                target_point[i, :] = closest_point
                spraying_dir[i, :] = -temp_normal

    # Zero out upward spray directions
    spraying_dir[spraying_dir[:, 2] > 0, :] = 0.0

    # Remove collided points
    remove_id = []
    for i in range(down_sampled_traj.shape[0]):
        temp_point = down_sampled_traj[i, :]
        grid_id, neighbour = get_grid_id(temp_point, grid_size,
                                          Q3_bound[0], Q1_bound[0],
                                          Q3_bound[1], Q1_bound[1],
                                          Q3_bound[2], Q1_bound[2],
                                          3, True)
        if grid_id < 0 or np.isnan(grid_id) or grid_id == 0:
            continue

        part1 = ALL_grid[grid_id - 1]['point']
        if part1.size > 0:
            part1_pts = part1[:, :3]
        else:
            part1_pts = np.empty((0, 3))

        part2_list = []
        for k in range(len(neighbour)):
            p = ALL_grid[neighbour[k] - 1]['point']
            if p.size > 0:
                part2_list.append(p)
        if part2_list:
            part2 = np.vstack(part2_list)
            part2_pts = part2[:, :3]
        else:
            part2_pts = np.empty((0, 3))

        neighbour_c = np.vstack([part1_pts, part2_pts]) if part1_pts.size > 0 or part2_pts.size > 0 else np.empty((0, 3))

        collided = False
        for j in range(neighbour_c.shape[0]):
            temp_c = neighbour_c[j, :]
            temp_dist = np.linalg.norm(temp_point - temp_c)
            if np.round(temp_dist - offset_distance, 1) < 0:
                collided = True
                break
        if collided:
            remove_id.append(i)

    if remove_id:
        down_sampled_traj = np.delete(down_sampled_traj, remove_id, axis=0)
        spraying_dir = np.delete(spraying_dir, remove_id, axis=0)

    # ── Statistics ──
    print(f'Mean viewing distance = {np.mean(dist2surf):.4f}')
    print(f'Median viewing distance = {np.median(dist2surf):.4f}')
    print(f'Max viewing distance = {np.max(dist2surf):.4f}')
    print(f'Min viewing distance = {np.min(dist2surf):.4f}')
    print(f'Standard Deviation = {np.std(dist2surf):.4f}')

    # Distribution fit
    mu, sigma = stats_norm.fit(dist2surf)
    print(f'Normal fit: mu={mu:.4f}, sigma={sigma:.4f}')

    # Histogram with fit
    fig5, (ax5a, ax5b) = plt.subplots(1, 2, figsize=(12, 5))
    ax5a.hist(dist2surf, bins=30, density=True, alpha=0.6, edgecolor='k')
    x_fit = np.linspace(dist2surf.min(), dist2surf.max(), 100)
    ax5a.plot(x_fit, stats_norm.pdf(x_fit, mu, sigma), 'r-', linewidth=2)
    ax5a.set_title('Distribution of viewing distance across the inspection path')
    ax5a.set_xlabel('Viewing Distance (m)')
    ax5a.set_ylabel('Density')

    ax5b.hist(dist2surf, bins=np.arange(dist2surf.min(), dist2surf.max() + 0.025, 0.025),
              density=True, alpha=0.6, edgecolor='k')
    ax5b.plot(x_fit, stats_norm.pdf(x_fit, mu, sigma), 'r-', linewidth=2)
    ax5b.set_title('Distribution (bin width=0.025)')
    ax5b.set_xlabel('Viewing Distance (m)')
    ax5b.set_ylabel('Density')

    # UV visibility plot
    fig6, ax6 = plt.subplots(figsize=(8, 6))
    ax6.triplot(UV[:, 0], UV[:, 1], connectivity, color='k')
    if len(UV_Vist) > 0:
        ax6.triplot(UV[:, 0], UV[:, 1], connectivity[UV_Vist, :],
                    color='g', linewidth=0.5)
    ax6.plot(test_traj[:, 0], test_traj[:, 1], 'r', linewidth=2)
    ax6.legend(['', 'Expected coverage', 'Planned 2D sweep path'])
    ax6.set_title('UV map')
    ax6.set_aspect('equal')

    # ── 3D Visibility ──
    Visibility_3D = calculate_vist_3d(
        down_sampled_traj, viewingDist, spraying_dir, FOV,
        CCentroid_3D, Normal, viewingAngle_threshold, ALL_grid, Grid_INFO,
    )
    Vist3D = np.unique(np.concatenate([v['Visibility'] for v in Visibility_3D if v['Visibility']]))

    # Final plot with visibility
    fig7 = plt.figure()
    ax7 = fig7.add_subplot(111, projection='3d')
    ax7.plot_trisurf(Vertex[:, 0], Vertex[:, 1], Vertex[:, 2],
                     triangles=connectivity, color=[0.9, 0.9, 0.9],
                     edgecolor=[0.4, 0.4, 0.4])
    ax7.view_init(elev=35, azim=25)

    hermite_curve = hermite_plot(down_sampled_traj)
    ax7.plot(hermite_curve[:, 0], hermite_curve[:, 1], hermite_curve[:, 2],
             color=[0, 0.447, 0.741], linewidth=2)

    ax7.quiver(down_sampled_traj[:, 0], down_sampled_traj[:, 1], down_sampled_traj[:, 2],
               spraying_dir[:, 0], spraying_dir[:, 1], spraying_dir[:, 2],
               length=0.5, color=[0.635, 0.078, 0.184])

    if len(Vist3D) > 0:
        ax7.plot_trisurf(Vertex[:, 0], Vertex[:, 1], Vertex[:, 2],
                         triangles=connectivity[Vist3D, :],
                         color='g', edgecolor=[0, 0, 0])

    ax7.set_title('The planned 3D inspection path')
    ax7.legend(['', 'Waypoints', 'Inspection Path', 'Viewing Direction', 'Covered surface'])
    ax7.set_aspect('equal')

    # Coverage percentage
    if len(UV_Vist) > 0 and len(Vist3D) > 0:
        coverage = np.sum(np.isin(UV_Vist, Vist3D)) / len(UV_Vist)
        print(f"Coverage Percentage = {coverage:.4f}")

    # ── Output ──
    output_path = np.column_stack([
        -down_sampled_traj[:, 0],
        -down_sampled_traj[:, 1],
        down_sampled_traj[:, 2],
    ])
    output_viewing_dir = np.column_stack([
        -spraying_dir[:, 1],
        -spraying_dir[:, 0],
        spraying_dir[:, 2],
    ])

    output_cameraCMD = np.zeros((output_viewing_dir.shape[0], 3))
    for i in range(output_viewing_dir.shape[0]):
        temp_vector = output_viewing_dir[i, :]
        temp_vector = temp_vector / (np.linalg.norm(temp_vector) + 1e-12)
        yaw = np.arctan2(temp_vector[1], temp_vector[0])
        if yaw < 0:
            yaw += 2 * np.pi
        pitch = np.arcsin(temp_vector[2])
        yaw = np.rad2deg(yaw)
        if yaw > 180:
            yaw -= 360
        pitch = np.rad2deg(pitch)

        if (np.isnan(yaw) and np.isnan(pitch) and i != 0) or (yaw == 0 and pitch == 0 and i != 0):
            output_cameraCMD[i, :] = output_cameraCMD[i - 1, :]
        else:
            output_cameraCMD[i, :] = [pitch, 0, yaw]

    np.savetxt(os.path.join(base_dir, 'pillar_path.txt'), output_path, delimiter='\t')
    np.savetxt(os.path.join(base_dir, 'pillar_CMD.txt'), output_cameraCMD, delimiter='\t')

    t_elapsed = time.time() - t_start
    print(f"Total elapsed time: {t_elapsed:.2f}s")

    plt.show()


if __name__ == '__main__':
    main()
