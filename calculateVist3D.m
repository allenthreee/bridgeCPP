function Visibility_3D = calculateVist3D(traj,viewingDist, viewingAngle, FOV, centroid3D, normal, threshold, ALL_grid, grid_INFO)
Q1_bound = [grid_INFO(1) grid_INFO(2) grid_INFO(3)];
Q3_bound = [grid_INFO(4) grid_INFO(5) grid_INFO(6)];
grid_size = grid_INFO(end);

for i = 1:size(traj,1)-1
    point1 = traj(i,:);
    point2 = traj(i+1,:);
    segment = [linspace(point1(1),point2(1),5)' linspace(point1(2),point2(2),5)' linspace(point1(3),point2(3),5)'];
    vist = [];
    for j = 1:size(segment,1)
        temp_point = segment(j,:);
        angle = viewingAngle(i,:);
        [gridID,neighbour] = getGridID(temp_point, grid_size, Q3_bound(1), Q1_bound(1), Q3_bound(2), ...
            Q1_bound(2), Q3_bound(3), Q1_bound(3),3,true);
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
        for k = 1:size(neighbour_centroid,1)
            temp_centroid = neighbour_centroid(k,:);
            temp_dist = norm(temp_point-temp_centroid);
            id = all(centroid3D == temp_centroid,2);
            temp_normal = normal(id,:);
            if size(temp_normal,1) > 1
                temp_normal = temp_normal(1,:);
            end
            isAngleOfIncidenceLarge = detectLargeAngleOfIncidence(temp_normal, angle, threshold);
            if temp_dist > viewingDist || isAngleOfIncidenceLarge
                continue
            else
                isWithinViewCone = CalPointVist(temp_point, angle, temp_centroid, FOV);
                if isWithinViewCone
                    id = find(id);
                    vist = [vist; id];
                end
            end
        end
    end
    vist = unique(vist);
    Visibility_3D(i).Visibility = vist;
end
end