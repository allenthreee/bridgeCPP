clc
close all 
clear
tic
model = 'pillar';

if strcmp(model, 'dice')
    fileID = fopen('3D 2D maping.txt');
    fileID2 = fopen("faces.txt");
    % load transform_sweep_path.mat
    load test_traj_simple.mat
    test_traj = test_traj_simple;n 
    fileID3 = fopen("surface_normal.txt");
elseif strcmp(model, 'portal')
    fileID = fopen('portal_mapping.txt');
    fileID2 = fopen("portal_faces.txt");
    load portal_path.mat
    fileID3 = fopen("portal_surface_normal.txt");
elseif strcmp(model, 'rounded_portal')
    fileID = fopen('rounded_portal.txt');
    fileID2 = fopen("rounded_portal_faces.txt");
    % load rounded_portal_path.mat
    load portal_uni_photo_path.mat
    fileID3 = fopen("rounded_portal_normal.txt");
elseif strcmp(model, 'HOA')
    fileID = fopen("HOA_mapping.txt");
    fileID2 = fopen("HOA_faces.txt");
    fileID3 = fopen("HOA_surface_normal.txt");
    load HOA_sweep_two.mat
elseif strcmp(model,'HOA_complete')
    fileID = fopen("HOA_complete_mapping.txt");
    fileID2 = fopen("HOA_complete_faces.txt");
    fileID3 = fopen("HOA_complete_normal.txt");
    load HOA_complete_sweep.mat
elseif strcmp(model, 'portal_reduced')
    fileID = fopen("round_portal_reduced_maping.txt");
    fileID2 = fopen("round_portal_reduced_faces.txt");
    fileID3 = fopen("round_portal_reduced_surface_normal.txt");
    load portal_reduced_sweep.mat
elseif strcmp(model, 'A380')
    fileID = fopen("A380_maping.txt");
    fileID2 = fopen("A380_faces.txt");
    fileID3 = fopen("A380_surface_normal.txt");
    % load A380_fuseNwing_detailed.mat
    load A380_photo1.mat
elseif strcmp(model, 'Jet')
    fileID = fopen("LearJet45_maping.txt");
    fileID2 = fopen("LearJet45_faces.txt");
    fileID3 = fopen("LearJet45_normal.txt");
    load jet_path.mat
    % load jet_photo.mat
elseif strcmp(model, 'pillar')
    fileID = fopen("./pillar/pillar_maping.txt");
    fileID2 = fopen("./pillar/pillar_faces.txt");
    fileID3 = fopen("./pillar/pillar_normal.txt");
    load ./pillar/pillar_path.mat
    % load jet_photo.mat
end

file = textscan(fileID, '%s %s %s %s %s %s %s %s');
file2 = textscan(fileID2, '%s %s %s %s %s');
file3 = textscan(fileID3, '%s %s %s %s');
for i = 1:size(file,2)
    file{1,i} = erase(file{1,i},'(');
    file{1,i} = erase(file{1,i},')');
    file{1,i} = erase(file{1,i},',');
    file{1,i} = erase(file{1,i},'<');
    file{1,i} = erase(file{1,i},'>');
end
% File 1 format: 
% Column 1: Vertex ID; 
% Column 2: 'Vector'
% Column 3: Vertex X;
% Column 4: Vertex Y;
% Column 5: Vertex Z;
% Column 6: 'Vector'
% Column 7: UV X;
% Column 8: UV Y;

for i = 1:size(file2,2)
    file2{1,i} = erase(file2{1,i},':');
end
%File 2 format:
% Colum 1: Face;
% Column 2: Face ID;
% Column 3: Vertex 1;
% Column 4: Vertex 2:
% Column 5: Vertex 3;

for i = 1:size(file3,2)
    file3{1,i} = erase(file3{1,i},'(');
    file3{1,i} = erase(file3{1,i},')');
    file3{1,i} = erase(file3{1,i},',');
    file3{1,i} = erase(file3{1,i},'<');
    file3{1,i} = erase(file3{1,i},'>');
end
% File 3 format:
% Column 1: 'Vector';
% Column 2: normal X;
% Column 3: normal Y;
% Column 4: normal Z;

Id = str2double(file{1,1});
vertex_X = str2double(file{1,3});
vertex_Y = str2double(file{1,4});
vertex_Z = str2double(file{1,5});
if strcmp(model, 'pillar')
    Vertex = [vertex_X/7.29 vertex_Y/7.29 vertex_Z/7.29];
else
    Vertex = [vertex_X vertex_Y vertex_Z];
end
uv_x = str2double(file{1,7});
uv_y = str2double(file{1,8});
UV = [uv_x uv_y];

FaceID = str2double(file2{1,2});
Faces = [str2double(file2{1,3}) str2double(file2{1,4}) str2double(file2{1,5})];

Normal = [str2double(file3{1,2}) str2double(file3{1,3}) str2double(file3{1,4})];

count_list = zeros(max(max(Faces))+1,1);
connectivity = zeros(size(Faces));
for i = 1:length(Faces)
    for ii = 1:size(Faces,2)
        times = count_list(Faces(i,ii)+1); % Find out how many times the vertex is selected
        vertex_id = find(Id == Faces(i,ii)); % Find out all the element with the corresponding vertex ID
        vertex_id = vertex_id(times+1);
        count_list(Faces(i,ii)+1) = count_list(Faces(i,ii)+1) + 1;
        connectivity(i,ii) = vertex_id;
    end
end

if strcmp(model, 'HOA') 
    Vertex = Vertex/60;
end

UV_map = triangulation(connectivity, UV);
boundaryPoints = freeBoundary(UV_map);
figure
subplot(1,2,1)
triplot(UV_map,'k')
axis equal
hold on
subplot(1,2,2)
model_3D = trisurf(connectivity, Vertex(:,1), Vertex(:,2), Vertex(:,3),'facecolor',[.7 .7 .7]);
model_triangulation = triangulation(connectivity, Vertex(:,1), Vertex(:,2), Vertex(:,3));
hold on
grid on
axis equal
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
drawnow

% Discretization
offset_distance = 1.5;
% if strcmp(model, 'pillar') 
%     grid_size = 0.5;
% else
    grid_size = offset_distance;
% end
Q1_bound = grid_size-rem(max(Vertex),grid_size)+max(Vertex)+3*grid_size;
% Q1_bound = (max(Vertex)+offset_distance)*1.25;
Q3_bound = min(Vertex)-(grid_size+rem(min(Vertex),grid_size))-3*grid_size;
Q3_bound(3) = 0;
plot3(Q1_bound(1),Q1_bound(2),Q1_bound(3),'m*','MarkerSize',30)
plot3(Q3_bound(1),Q3_bound(2),Q3_bound(3),'m*','MarkerSize',30)
no_of_grid = (Q1_bound-Q3_bound)/grid_size;
disp('Number of grid:')
disp((Q1_bound-Q3_bound)/grid_size)

grid_size_2D = 0.05;
u_upper = 1;
u_lower = 0;
v_upper = u_upper;
v_lower = u_lower;
no_grid_2D = ((u_upper-u_lower)/grid_size_2D)*((v_upper-v_lower)/grid_size_2D);
ALL_2Dgrid(no_grid_2D).point = [];
Grid_INFO = [Q1_bound(1); Q1_bound(2);Q1_bound(3);Q3_bound(1);Q3_bound(2);Q3_bound(3);grid_size];
no_grid = ((Q1_bound(1)-Q3_bound(1))/grid_size)*((Q1_bound(2)-Q3_bound(2))/grid_size)*((Q1_bound(3)-Q3_bound(3))/grid_size);
ALL_grid(no_grid).point = [];
FOV = 90;
viewingDist = offset_distance/cos(deg2rad(FOV/2));
Vist_radius3D = offset_distance*tan(deg2rad(FOV/2))*0.5;
viewingAngle_threshold = 60;


%%%%%%%%%%%%%%%%%%% Compute centroid for each triangle mesh %%%%%%%%%%%%%%%%%%%
CCentroid = zeros(size(connectivity,1),2);
CCentroid_3D = zeros(size(connectivity,1),3);
for i = 1:size(connectivity,1)
    vertex_centroid = [];
    DDD_centroid = [];
    for ii = 1:size(connectivity,2)
        vertex_centroid = [vertex_centroid; UV(connectivity(i,ii),:)];
        DDD_centroid = [DDD_centroid; Vertex(connectivity(i,ii),:)];
    end
    mean_u = mean(vertex_centroid(:,1));
    mean_v = mean(vertex_centroid(:,2));
    CCentroid(i,:) = [mean_u mean_v];
    [gridID,~] = getGridID(CCentroid(i,:), grid_size_2D, u_lower, u_upper, v_lower, v_upper, [], [],2, false);
    ALL_2Dgrid(gridID).point = [ALL_2Dgrid(gridID).point; CCentroid(i,:)];
    mean_x = mean(DDD_centroid(:,1));
    mean_y = mean(DDD_centroid(:,2));
    mean_z = mean(DDD_centroid(:,3));
    CCentroid_3D(i,:) = [mean_x mean_y mean_z];
    [gridID,~] = getGridID(CCentroid_3D(i,:), grid_size, Q3_bound(1), Q1_bound(1), Q3_bound(2), Q1_bound(2), Q3_bound(3), Q1_bound(3),3,false);
    ALL_grid(gridID).point = [ALL_grid(gridID).point; CCentroid_3D(i,:) i];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uv23Dratio = ratiocompute(CCentroid, CCentroid_3D);
disp("UV to 3D unit ratio = ")
disp(uv23Dratio)


uv_radius = Vist_radius3D*uv23Dratio;

%%%%%%%%%%%%%%%%%%%%%%%%% REWORKED UV-to-3D & Visibility %%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1)
traj_3D = zeros(size(test_traj,1),3);
traj_FaceID = zeros(size(traj_3D,1),1);
Visibility_UV(1).Visibility = [];
for i = 1:size(test_traj,1)
    uv_traj = test_traj(i,:);
    [gridID,neighour] = getGridID(uv_traj, grid_size_2D, u_lower, u_upper, v_lower, v_upper, [], [],2,true);
    list_of_2Dcentroid = ALL_2Dgrid(gridID).point;
    visible = [];
    for j = 1:size(neighour,1)
        list_of_2Dcentroid = [list_of_2Dcentroid; ALL_2Dgrid(neighour(j)).point];
    end
    %%%%%%%%
    list_of_2Dcentroid(all(list_of_2Dcentroid== 0 ,2),:) = [];
    %%%%%%%%
    for ii = 1:size(list_of_2Dcentroid,1)
        if ii == 1
            best_dist = norm(uv_traj-list_of_2Dcentroid(ii,:));
            closest = list_of_2Dcentroid(ii,:);
            continue
        end
        temp_dist = norm(uv_traj-list_of_2Dcentroid(ii,:));
        if temp_dist < best_dist
            best_dist = temp_dist;
            closest = list_of_2Dcentroid(ii,:);
        end
        if temp_dist < uv_radius
            temp_id = find(all(CCentroid == list_of_2Dcentroid(ii,:),2));
            visible = [visible; temp_id];
        end

    end
    traj_FaceID(i) = find(all(CCentroid == closest,2));
    Visibility_UV(i).Visibility = visible;
    %%%%%%%%%%% TO DO  (22/03/2024)  %%%%%%%%%%%%%%%%
    %       FIX THIS PART
    barycentric_UV = cartesianToBarycentric(UV_map,traj_FaceID(i),uv_traj);
    mesh_point = barycentricToCartesian(model_triangulation, traj_FaceID(i), barycentric_UV);
    traj_3D(i,:) = mesh_point;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UV_Vist = vertcat(Visibility_UV.Visibility);
UV_Vist = unique(UV_Vist);
trisurf(connectivity(UV_Vist,:), uv_x, uv_y, zeros(size(uv_x,1),1),'facecolor','g');
plot(test_traj(:,1), test_traj(:,2),'r','LineWidth',2.5);
subplot(1,2,2)
plot3(traj_3D(:,1),traj_3D(:,2),traj_3D(:,3),'r');
figure
hold on
axis on
axis equal
grid on
trisurf(connectivity, Vertex(:,1), Vertex(:,2), Vertex(:,3),'facecolor',[.7 .7 .7]);
% plot3(traj_3D(:,1),traj_3D(:,2),traj_3D(:,3),'b');
view([35 25])
no_segment = 1;
prev_end = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Offset Part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(traj_FaceID,1)-1
    current = traj_FaceID(i);
    next = traj_FaceID(i+1);
    if current~= next
        OG_segment(no_segment).segment = traj_3D(prev_end:i,:);
        OG_segment(no_segment).FaceID = current;
        no_segment = no_segment + 1;
        prev_end = i+1;
    elseif i == size(traj_FaceID,1)-1
        OG_segment(no_segment).segment = traj_3D(prev_end:end,:);
        OG_segment(no_segment).FaceID = current;
    end
end

untrim_segment = untrim_offset3D(OG_segment, offset_distance, Normal, CCentroid_3D);
unclip_traj = [];
for i = 1:size(untrim_segment,2)
    unclip_traj = [unclip_traj; untrim_segment(i).segment];
    temp_seg = untrim_segment(i).segment;
    % plot3(temp_seg(:,1), temp_seg(:,2), temp_seg(:,3),'g')
end
plot3(unclip_traj(:,1), unclip_traj(:,2), unclip_traj(:,3),'r','LineWidth',2.5);
drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Clipping part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
untrim_segment(1).collisionID = [];
untrim_segment(1).GroundClear = [];

[untrim_segment, replan_segment, replan_segmentEND] = clipping(untrim_segment,offset_distance, CCentroid_3D, unclip_traj,connectivity, Vertex, Normal,ALL_grid, Grid_INFO);
for i = 1:size(untrim_segment,2)
    if isempty(untrim_segment(i).collisionID) && isempty(untrim_segment(i).GroundClear)
        plot3(untrim_segment(i).segment(:,1),untrim_segment(i).segment(:,2),untrim_segment(i).segment(:,3),'g');
    end
end

replan_segment_HEADs = zeros(size(replan_segmentEND,2),1);
replan_segment_TAILs = zeros(size(replan_segmentEND,2),1);
for i = 1:size(replan_segmentEND,2)
    replan_segment_HEADs(i) = replan_segmentEND(i).Head;
    replan_segment_TAILs(i) = replan_segmentEND(i).Tail;
end

subplot(1,2,1)
hold on
axis on
axis equal
grid on
trisurf(connectivity, Vertex(:,1), Vertex(:,2), Vertex(:,3),'facecolor',[.7 .7 .7]);
plot3(traj_3D(:,1),traj_3D(:,2),traj_3D(:,3),'b');
view([35 25])
unclip_traj = [];
for i = 1:size(untrim_segment,2)
    unclip_traj = [unclip_traj; untrim_segment(i).segment];
    temp_seg = untrim_segment(i).segment;
    plot3(temp_seg(:,1), temp_seg(:,2), temp_seg(:,3),'r')
end

for i = 1:size(replan_segmentEND,2)
    % replan_traj = [replan_traj;replan_segmentEND(i).segment];
    temp_seg = replan_segmentEND(i).segment;
    plot3(temp_seg(:,1), temp_seg(:,2), temp_seg(:,3),'b')
end

subplot(1,2,2)
hold on
axis on
axis equal
grid on
trisurf(connectivity, Vertex(:,1), Vertex(:,2), Vertex(:,3),'facecolor',[.7 .7 .7]);
% plot3(traj_3D(:,1),traj_3D(:,2),traj_3D(:,3),'b');
view([35 25])
replan_traj = [];
i = 1;
while i < size(untrim_segment,2)
    if any(replan_segment_HEADs == i) || any(replan_segment_TAILs == i)
        if any(replan_segment_HEADs == i) % Replan Segment starts
            if ~isempty(untrim_segment(i).segment) % Replan Segment starts within this segment
                replan_ID_list = find(replan_segment_HEADs == i);
                if size(replan_ID_list,2) > 1
                    GIGA = 1;
                end
                for j = 1:size(replan_ID_list,2) 
                    replan_ID = replan_ID_list(j);
                    if replan_segmentEND(replan_ID).Head == replan_segmentEND(replan_ID).Tail % IF this is a local segment
                        start = replan_segmentEND(replan_ID).insertion_pt;
                        end_point = replan_segmentEND(replan_ID).insertion_end;
                        if start == 0 % IF local segment start at first node
                            replan_traj = [replan_traj; replan_segmentEND(replan_ID).segment];
                            replan_traj = [replan_traj; untrim_segment(i).segment(end_point:end,:)];
                        else % IF local segment starts at middle
                            replan_traj = [replan_traj; untrim_segment(i).segment(1:start,:); replan_segmentEND(replan_ID).segment];
                            replan_traj = [replan_traj; untrim_segment(i).segment(end_point:end,:)];
                        end
                    else % IF this is NOT a local segment: head ~= tail
                        start = replan_segmentEND(replan_ID).insertion_pt;
                        if start == 0 % This should be the same as replan segment replaced the original segment --> untrim is empty
                            replan_traj = [replan_traj; replan_segmentEND(replan_ID).segment];
                            replan_traj = [replan_traj; untrim_segment(i).segment(end_point:end,:)];
                        else
                            replan_traj = [replan_traj; untrim_segment(i).segment(1:start,:); replan_segmentEND(replan_ID).segment];
                        end
                    end
                end
            else % Replan Segment replaced the original segment
                replan_ID = find(replan_segment_HEADs == i);
                replan_traj = [replan_traj; replan_segmentEND(replan_ID).segment];
            end
        else % Replan Segment END
            % if ~isempty(untrim_segment(i).segment) % Replan Segment end within this segment
            %     replan_ID = find(replan_segment_TAILS == i);
            %     if size(replan_ID,2)>1
            %         GIGA = 1;
            %     end
            %     replan_ID = replan_ID(1);
            %     replan_traj = [replan_traj; replan_segmentEND(replan_ID).segment];
            %     replan_traj = [replan_traj; untrim_segment(i).segment];
            % else % Replan Segment replaced the original segment
            %     replan_ID = find(replan_segment_TAILs == i);
            %     replan_traj = [replan_traj; replan_segmentEND(replan_ID).segment];
            % end
        end
    else
        replan_traj = [replan_traj; untrim_segment(i).segment];
    end
    i = i+1;
end

plot3(replan_traj(:,1), replan_traj(:,2), replan_traj(:,3),'r');
disp("Replan Finished")

start_point = replan_traj(1,:);
end_point = replan_traj(end,:);


down_sampled_traj = replan_traj;
% disp("Starting Two-opt")
% down_sampled_traj = two_opt(down_sampled_traj,100);

% Mean distance between points
% constant lower --> more detailed (more points); constant higher --> less
% detailed (less points)
down_sample_goal = 3*(max(max(replan_traj))*0.01);


% Use Distance to down sample
% Calculate point to point distance -> Use median/ mean as reference
while true
    distance = [];
    for i = 1:size(down_sampled_traj,1)
        if i ~= size(down_sampled_traj,1)
            temp_dist = norm(down_sampled_traj(i+1,:)-down_sampled_traj(i,:));
            distance = [distance; temp_dist];
        end
    end
    % The downsample threshold is the mean distance between points %
    threshold = mean(distance);
    i = 1;
    down_sampled_ID = [];
    temp_dist = 0;
    while i<size(down_sampled_traj,1)
        temp_dist = temp_dist + distance(i);
        if temp_dist > threshold
            down_sampled_ID = [down_sampled_ID; i];
            temp_dist = 0;
        end
        i = i +1;
    end
    down_sampled_traj = down_sampled_traj(down_sampled_ID,:);
    if threshold >= down_sample_goal
        break
    end
end
down_sampled_traj = [start_point; down_sampled_traj; end_point];
disp("Starting Two-opt")
down_sampled_traj = two_opt(down_sampled_traj,20);

figure
hold on
axis on
axis equal
grid on
trisurf(connectivity, Vertex(:,1), Vertex(:,2), Vertex(:,3),'facecolor',[.7 .7 .7]);
view([35 25])
plot3(down_sampled_traj(:,1), down_sampled_traj(:,2), down_sampled_traj(:,3),'r');
plot3(Q1_bound(1),Q1_bound(2),Q1_bound(3),'m*','MarkerSize',30)
plot3(Q3_bound(1),Q3_bound(2),Q3_bound(3),'m*','MarkerSize',30)

dist2surf = zeros(size(down_sampled_traj,1),1);
target_point = zeros(size(down_sampled_traj,1),3);
spraying_dir = zeros(size(down_sampled_traj,1),3);

for i = 1:size(down_sampled_traj,1)
    temp_point = down_sampled_traj(i,:);
    [gridID,neighbour] = getGridID(temp_point, grid_size, Q3_bound(1), Q1_bound(1), Q3_bound(2), Q1_bound(2), Q3_bound(3), Q1_bound(3),3,true);
    Part1 = ALL_grid(gridID).point;
    PointID_P1 = [];
    PointID_P2 = [];
    if ~isempty(Part1)
        PointID_P1 = ALL_grid(gridID).point(:,end);
        Part1 = Part1(:,1:3);
    end
    Part2 = [];
    for k = 1:size(neighbour,1)
        Part2 = [Part2; ALL_grid(neighbour(k)).point];
    end

    if ~isempty(Part2)
        PointID_P2 = Part2(:,end);
        Part2(:,end) = [];
        PointID = [PointID_P1; PointID_P2];
    else
        PointID = PointID_P1;
    end
    neighbour_centroid = [Part1; Part2];
    giga = 1;
    closest_point = [];
    current_shortest_dist = [];
    for j = 1:size(neighbour_centroid,1)
        temp_centroid = neighbour_centroid(j,:);
        temp_dist = norm(temp_point-temp_centroid);
        if isempty(current_shortest_dist)
            current_shortest_dist = temp_dist;
            closest_point = temp_centroid;
        else
            if temp_dist < current_shortest_dist
                current_shortest_dist = temp_dist;
                closest_point = temp_centroid;
            end

        end
    end
    vectorToSurface = temp_point - closest_point;
    temp_normal = Normal(all(CCentroid_3D == closest_point,2),:);
    perpendicularDistance = abs(dot(vectorToSurface, temp_normal));
    dist2surf(i) = perpendicularDistance;
    target_point(i,:) = closest_point;
    spraying_dir(i,:) = -temp_normal;
end

spraying_dir(spraying_dir(:,3)>0,:) = zeros(size(spraying_dir(spraying_dir(:,3)>0,:)));
remove_id = [];
for i = 1:size(down_sampled_traj,1)
    temp_point = down_sampled_traj(i,:);
    [gridID,neighbour] = getGridID(temp_point, grid_size, Q3_bound(1), Q1_bound(1), Q3_bound(2), Q1_bound(2), Q3_bound(3), Q1_bound(3),3,true);
    if gridID < 0 || isnan(gridID) || isempty(gridID)
        GIGA_ERROR = 1;
    end
    Part1 = ALL_grid(gridID).point;
    PointID_P1 = [];
    PointID_P2 = [];
    if ~isempty(Part1)
        PointID_P1 = ALL_grid(gridID).point(:,end);
        Part1 = Part1(:,1:3);
    else
        GIGA = 1;
    end
    Part2 = [];
    for k = 1:size(neighbour,1)
        Part2 = [Part2; ALL_grid(neighbour(k)).point];
    end

    if ~isempty(Part2)
        PointID_P2 = Part2(:,end);
        Part2(:,end) = [];
        PointID = [PointID_P1; PointID_P2];
    else
        PointID = PointID_P1;
    end
    neighbour_centroid = [Part1; Part2];
    Collided_point = [];
    dist_list = [];
    for j = 1:size(neighbour_centroid,1)
        temp_centroid = neighbour_centroid(j,:);
        temp_dist = norm(temp_point-temp_centroid);
        % Dont need to find closest, only check for collision
        if round(temp_dist - offset_distance,1) < 0
            Collided_point = [Collided_point; temp_centroid];
            dist_list = [dist_list; temp_dist];
        end
    end
    if ~isempty(Collided_point)
        remove_id = [remove_id;i];
    end
end

down_sampled_traj(remove_id,:) = [];
spraying_dir(remove_id,:) = [];


disp('Mean viewing distance = ');
disp(mean(dist2surf));
disp('Median viewing distance = ');
disp(median(dist2surf));
disp('Max viewing distance = ');
disp(max(dist2surf));
disp('Min viewing distance = ');
disp(min(dist2surf));
disp('Standard Deviation = ');
disp(std(dist2surf));

% test = sort(dist2surf,'ascend');
% figure
% hold on
% plot(1:size(dist2surf,1),test);
% axis tight

figure
h = histfit(dist2surf);
title('Distribution of viewing distance across the inspection path')
xlabel('Viewing Distance (m)')
ylabel('Count of waypoints')
CurveX = h(2).XData;
CurveY = h(2).YData;

figure
hold on
histogram(dist2surf, 'BinWidth',0.025);
plot(CurveX, CurveY);
title('Distribution of viewing distance across the inspection path')
xlabel('Viewing Distance (m)')
ylabel('Count of waypoints')



figure
hold on
axis off
axis equal
triplot(UV_map,'k')
trisurf(connectivity(UV_Vist,:), uv_x, uv_y, zeros(size(uv_x,1),1),'facecolor','g');
plot(test_traj(:,1), test_traj(:,2),'r','LineWidth',2);
legend('','Expected coverage','Planned 2D sweep path')
title('UV map')

result = fitdist(dist2surf,'Normal');
disp(result)


Visibility_3D = calculateVist3D(down_sampled_traj,viewingDist, spraying_dir, FOV, CCentroid_3D, Normal, viewingAngle_threshold, ALL_grid, Grid_INFO);
Vist3D = vertcat(Visibility_3D.Visibility);
Vist3D = unique(Vist3D);

figure
hold on
axis on
axis equal
grid on
trisurf(connectivity, Vertex(:,1), Vertex(:,2), Vertex(:,3),'facecolor',[.9 .9 .9],'EdgeColor',[.4 .4 .4]);
view([35 25])
Hermite_curve = Hermite_plot(down_sampled_traj);
hold on
quiver3(down_sampled_traj(:,1),down_sampled_traj(:,2), down_sampled_traj(:,3), ...
    spraying_dir(:,1), spraying_dir(:,2), spraying_dir(:,3),0.5,'Color',[0.6350 0.0780 0.1840])
title('The planned 3D inspection path')

trisurf(connectivity(Vist3D,:), Vertex(:,1), Vertex(:,2), Vertex(:,3),'facecolor','g','EdgeColor',[0 0 0]);
legend('','Waypoints','Inspection Path','Viewing Direction','Covered surface')

coverage_percent = sum(ismember(UV_Vist, Vist3D))/size(UV_Vist,1);
disp("Coverage Percentage = ")
disp(coverage_percent)


% Portal_reduced offset: x:0; y:-18; z:-4.5
% A380 offset : x:+40; y:-35.9; z:0
% Jet offset : x: +15; y: 0; z: 0
output_path = [-(down_sampled_traj(:,1)) -(down_sampled_traj(:,2)) down_sampled_traj(:,3)];
output_viewing_dir = [-spraying_dir(:,2) -spraying_dir(:,1) spraying_dir(:,3)];
output_cameraCMD = zeros(size(output_viewing_dir,1),3);
for i = 1:size(output_viewing_dir,1)
    temp_vector = output_viewing_dir(i,:);
    temp_vector = temp_vector/norm(temp_vector);
    yaw = atan2(temp_vector(2),temp_vector(1));
    if yaw < 0
        yaw = yaw+ 2*pi;
    end
    pitch = asin(temp_vector(3));
    yaw = rad2deg(yaw);
    if yaw > 180
        yaw = yaw - 360;
    end
    pitch = rad2deg(pitch);
    if (isnan(yaw) && isnan(pitch) && i~= 1) || (yaw == 0 && pitch == 0 && i~= 1) % If looking up only
        output_cameraCMD(i,:) = output_cameraCMD(i-1,:); % maintain previous camera angle
    else
        output_cameraCMD(i,:) = [pitch  0 yaw];
    end
end

writematrix(output_path,'pillar_path.txt','Delimiter','tab')
writematrix(output_cameraCMD,'pillar_CMD.txt','Delimiter','tab')

toc

% Compute the viewing angle of each waypoints
% Point to closest mesh
