function [untrim_segment, replan_segment, replan_segmentEND] = clipping(untrim_segment,offset_distance, ccentroid, unclip_traj,connectivity, Vertex, normal, ALL_grid, grid_INFO)
% figure
% hold on
% axis on
% axis equal
% grid on
% trisurf(connectivity, Vertex(:,1), Vertex(:,2), Vertex(:,3),'facecolor',[.7 .7 .7]);
% plot3(unclip_traj(:,1), unclip_traj(:,2), unclip_traj(:,3),'g');
% view([35 25])
array_size = [round(size(ccentroid,1)/10) 1];
Q1_bound = [grid_INFO(1) grid_INFO(2) grid_INFO(3)];
Q3_bound = [grid_INFO(4) grid_INFO(5) grid_INFO(6)];
grid_size = grid_INFO(end);

ground_lvl = 1;
for i = 1:size(untrim_segment,2)
    temp_segment = untrim_segment(i).segment;
    for ii = 1:size(temp_segment,1)
        ColliFaceID = 1;
        temp_point = temp_segment(ii,:);
        if i == 6
            giga = 1;
        end
        untrim_segment(i).closest(ii).id = zeros(array_size);
        untrim_segment(i).closest(ii).dist = Inf(array_size);
        temp_dist_array = Inf(array_size);
        if temp_point(3) < ground_lvl
            untrim_segment(i).GroundClear = [untrim_segment(i).GroundClear; ii];
            % plot3(temp_point(1),temp_point(2),temp_point(3),'r.');
        end
        
        % REWORK COLLISION CHECK %%%%%%%%%%%%%%%%
        % for j = 1:size(ccentroid,1) % THIS IS COLLISION CHECK
        %     temp_centroid = ccentroid(j,:);
        %     temp_dist = norm(temp_point-temp_centroid);
        %     if any(temp_dist < untrim_segment(i).closest(ii).dist) || any(isinf(untrim_segment(i).closest(ii).dist))
        %         if any(isinf(untrim_segment(i).closest(ii).dist))
        %             id_list = find(isinf(untrim_segment(i).closest(ii).dist));
        %             untrim_segment(i).closest(ii).dist(id_list(1)) = temp_dist;
        %             untrim_segment(i).closest(ii).id(id_list(1)) = j;
        %         else
        %             id_list = find(temp_dist < untrim_segment(i).closest(ii).dist);
        %             untrim_segment(i).closest(ii).dist(id_list(1)) = temp_dist;
        %             untrim_segment(i).closest(ii).id(id_list(1)) = j;
        %         end
        %     end
        %     if temp_dist < offset_distance
        %         untrim_segment(i).collisionID = [untrim_segment(i).collisionID; ii];
        %         untrim_segment(i).collisionFace(ii).id(ColliFaceID) = j;
        %         ColliFaceID = ColliFaceID+1;
        %         % plot3(temp_point(1),temp_point(2),temp_point(3),'b.');
        %     end
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%% REWORKED COLLISION CHECK %%%%%%%%%
        [gridID,neighbour] = getGridID(temp_point, grid_size, Q3_bound(1), Q1_bound(1), Q3_bound(2), Q1_bound(2), Q3_bound(3), Q1_bound(3),3,true);
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
        giga = 1;
        for j = 1:size(neighbour_centroid,1)
            temp_centroid = neighbour_centroid(j,:);
            temp_dist = norm(temp_point-temp_centroid);
            if any(temp_dist < untrim_segment(i).closest(ii).dist) || any(isinf(untrim_segment(i).closest(ii).dist))
                if any(isinf(untrim_segment(i).closest(ii).dist))
                    id_list = find(isinf(untrim_segment(i).closest(ii).dist));
                    untrim_segment(i).closest(ii).dist(id_list(1)) = temp_dist;
                    centroid_ID = PointID(j);
                    untrim_segment(i).closest(ii).id(id_list(1)) = centroid_ID;
                else
                    id_list = find(temp_dist < untrim_segment(i).closest(ii).dist);
                    untrim_segment(i).closest(ii).dist(id_list(1)) = temp_dist;
                    centroid_ID = PointID(j);
                    untrim_segment(i).closest(ii).id(id_list(1)) = centroid_ID;
                end
            end
            if temp_dist < offset_distance
                untrim_segment(i).collisionID = [untrim_segment(i).collisionID; ii];
                centroid_ID = PointID(j);
                untrim_segment(i).collisionFace(ii).id(ColliFaceID) = centroid_ID;
                ColliFaceID = ColliFaceID+1;
                % plot3(temp_point(1),temp_point(2),temp_point(3),'b.');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        untrim_segment(i).GroundClear = unique(untrim_segment(i).GroundClear);
        untrim_segment(i).collisionID = unique(untrim_segment(i).collisionID);
        delete_list = (untrim_segment(i).closest(ii).id == 0);
        untrim_segment(i).closest(ii).id(delete_list) = [];
        untrim_segment(i).closest(ii).dist(delete_list) = [];
        if any(untrim_segment(i).closest(ii).id == 0)
            giga = 1;
        end
        % drawnow
    end
end
% untrim_segment data meaning:
% untrim_segment(1).segment: points of the corresponding segment
% untrim_segment(1).FaceID: The corresponding Face to be inspected by this segment
% untrim_segment(1).collisionID: List of segment points with the risk of collision (needs replan)
% untrim_segment(1).GroundClear: List of segment points with the risk of collision with ground (needs replan)
% untrim_segment(1).closest(1).id: Top 10% closest Face with segment points. E.g. untrim_segment(1).closest(1).id = List of FaceID closest to untrim_segment(1).segment(1)
replan_segID = 1;
seg_cont = false;
is_start = false;
temp_replanSeg = [];
temp_replanEND = [];
GroundClear_IDlist = [];
Collision_IDlist = [];
segment_IDlist = [];
replan_segment = [];
replan_segmentEND = [];
for i = 1:size(untrim_segment,2)
    disp('Current Segment = ')
    disp(i)
    if i == 833
        GIGA = 1;
    end
    if ~isempty(untrim_segment(i).collisionID) || ~isempty(untrim_segment(i).GroundClear)
        if ~isempty(untrim_segment(i).GroundClear)
            if i == 833
                GIGA = 1;
            end
            if i ~= 1
                if ~seg_cont
                    is_start = true; % If previous segment dont have collision --> is_start = true
                else
                    is_start = false;
                end
            else
                is_start = true;
            end
            if untrim_segment(i).GroundClear(end) == size(untrim_segment(i).segment,1) % Collision til the end of the current segment
                if i ~= size(untrim_segment,2)
                    if isempty(untrim_segment(i+1).GroundClear)
                        end_here = true; % If next segment dont have collision --> end_here = true
                    else
                        if untrim_segment(i+1).GroundClear(1) == 1 % If next segment collision start at first point --> seg_cont = true; end_here = false
                            seg_cont = true;
                            end_here = false;
                        else
                            seg_cont = false; % If not, current segment = end --> seg_cont = false; end_here = true;
                            end_here = true;
                        end
                    end
                else % If current the the last segment --> seg_cont = false; end_here = true
                    seg_cont = false;
                    end_here = true;
                end
            else
                % if seg_cont % If current is the continue of the previous segment AND segment do not continue to next segment
                    % seg_cont = false;
                    % is_start = true;
                    % end_here = true;
                % else % If current is NOT the continue of the previous AND segment do not conitnue to next
                    seg_cont = false;
                    % is_start = true;
                    end_here = true;
                % end
            end
            % disp('is_start = '); disp(is_start);
            % disp('seg_cont = '); disp(seg_cont);
            % disp('end_here = ');disp(end_here);
            giga = 1;

            for ii = 1:size(untrim_segment(i).GroundClear,1)
                point_id = untrim_segment(i).GroundClear(ii);
                if ii ~= size(untrim_segment(i).GroundClear,1)
                    if untrim_segment(i).GroundClear(ii+1) ~= point_id+1
                        % This untrim_segment contain 2 replan segment
                        is_start = true;
                        seg_cont = false;
                        end_here = true;
                    end
                else

                end
                if point_id == 51
                    GIGA = 1;
                end
                plot3(untrim_segment(i).segment(point_id,1),untrim_segment(i).segment(point_id,2),untrim_segment(i).segment(point_id,3),'g.');
                vector_list = zeros(size(untrim_segment(i).closest(point_id).id,1),3);
                old_point = untrim_segment(i).segment(point_id,:);
                collision_ID = [];
                for j = 1:size(vector_list,1)
                    face_id = untrim_segment(i).closest(point_id).id(j);
                    temp_dist = norm(old_point- ccentroid(face_id,:));
                    if temp_dist < offset_distance
                        collision_ID = [collision_ID; j];
                    end
                    temp_vector = normal(untrim_segment(i).closest(point_id).id(j),:);
                    temp_vector = temp_vector*(1/untrim_segment(i).closest(point_id).dist(j));
                    vector_list(j,:) = temp_vector;
                end
                if ~isempty(collision_ID)
                    vector_list = vector_list(collision_ID,:);
                end
                temp_vector = [0 0 1];
                if size(vector_list,1) < 2
                    test_sum = vector_list/norm(vector_list);
                else
                    test_sum = sum(vector_list)/norm(sum(vector_list));
                end
                test_sum(3) = abs(test_sum(3));
                % sum_vector = temp_vector + test_sum;
                sum_vector = test_sum;
                sum_vector = sum_vector/norm(sum_vector);
                if sum_vector(3) <= 0
                    giga = 1; % Error detection, waypoint should go up, so if vector(3) <0 --> Error
                end
                closest_pointID = untrim_segment(i).closest(point_id).id;
                closest_centroidLIST = ccentroid(closest_pointID,:);
                closest_vertexLIST = connectivity(closest_pointID,:);
                closest_vertexLIST = reshape(closest_vertexLIST,[],1);
                closest_vertexLIST = Vertex(closest_vertexLIST,:);
                closest_pointlist = [closest_centroidLIST; closest_vertexLIST];
                
                colli_type = 1;
                giga_vector = [temp_vector; vector_list];
                % giga_vector = [sum_vector; vector_list];
                if i == 80 && point_id == 28
                    pause = 1;
                end
                replanned_point = pointwise_replan(old_point, giga_vector, offset_distance, ground_lvl, closest_pointlist, colli_type,Vertex, ...
                    connectivity, ALL_grid, grid_INFO, ccentroid, normal);
                temp_replanSeg = [temp_replanSeg; replanned_point];
                % untrim_segment(i).segment(point_id,:) = [];
                % untrim_segment(i).GroundClear(ii) = [];
                % untrim_segment(i).closest(point_id) = [];
                GroundClear_IDlist = [GroundClear_IDlist; ii];
                segment_IDlist = [segment_IDlist;point_id];
                plot3(replanned_point(1),replanned_point(2),replanned_point(3),'b.');
                giga = 1;

            end
            temp_replanEND = [temp_replanEND; replanned_point];
            if is_start
                Prev_seg = i;
                % ADD INSERTION POINT
                Insertion_pt = segment_IDlist(1)-1;
                if end_here
                    Insertion_END = segment_IDlist(1);
                end
            end
            if end_here
                Next_seg = i;
                replan_segment(replan_segID).segment = temp_replanSeg;
                replan_segment(replan_segID).Head = Prev_seg;
                replan_segment(replan_segID).Tail = Next_seg;
                replan_segmentEND(replan_segID).segment = temp_replanEND;
                replan_segmentEND(replan_segID).Head = Prev_seg;
                replan_segmentEND(replan_segID).Tail = Next_seg;
                if ~is_start
                    Insertion_END = segment_IDlist(1);
                end
                replan_segmentEND(replan_segID).insertion_pt = [];
                replan_segmentEND(replan_segID).insertion_end = [];
                temp_replanSeg = [];
                temp_replanEND = [];
                replan_segID = replan_segID + 1;
            end
            plot3(replanned_point(1),replanned_point(2),replanned_point(3),'r.','MarkerSize',20);
            untrim_segment(i).GroundClear(GroundClear_IDlist) = [];
            GroundClear_IDlist = [];
            untrim_segment(i).segment(segment_IDlist,:) = [];
            untrim_segment(i).closest(segment_IDlist) = [];
            segment_IDlist = [];
            giga =1;
        else
            %%%%%%% COLLISION ONLY %%%%%%%
        % if ~isempty(untrim_segment(i).collisionID)
            if i == 833
                GIGA = 1;
            end
            if i ~= 1
                if ~seg_cont
                    is_start = true; % If previous segment dont have collision --> is_start = true
                else
                    is_start = false;
                end
            else
                is_start = true;
            end
            if untrim_segment(i).collisionID(end) == size(untrim_segment(i).segment,1) % Collision til the end of the current segment
                if i ~= size(untrim_segment,2)
                    if isempty(untrim_segment(i+1).collisionID)
                        end_here = true; % If next segment dont have collision --> end_here = true
                        seg_cont = false;
                    else
                        if untrim_segment(i+1).collisionID(1) == 1 % If next segment collision start at first point --> seg_cont = true; end_here = false
                            seg_cont = true;
                            end_here = false;
                        else
                            seg_cont = false; % If not, current segment = end --> seg_cont = false; end_here = true;
                            end_here = true;
                        end
                    end
                else % If current the the last segment --> seg_cont = false; end_here = true
                    seg_cont = false;
                    end_here = true;
                end
            else % Collision END before the end of the current segment
                    seg_cont = false;
                    end_here = true;
            end
            % disp('is_start = '); disp(is_start);
            % disp('seg_cont = '); disp(seg_cont);
            % disp('end_here = ');disp(end_here);
            giga = 1;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % WORK: HOW TO FIND DISCONTINUE ELEMENT IN AN ARRAY
            % OUTPUT: NO. OF SUB-SEGMENTS
            % https://www.mathworks.com/matlabcentral/answers/216921-need-to-remove-repeated-adjacent-elements-in-an-array
            test_array = untrim_segment(i).collisionID(1):untrim_segment(i).collisionID(end);
            temp_collisionID = untrim_segment(i).collisionID;
            dis_cont = diff(temp_collisionID);
            no_local_seg = sum(dis_cont>1)+1;
            if no_local_seg > 1
                GIGA = 1;
            end
            breakpoint_list = find(dis_cont > 1);
            if isempty(breakpoint_list)
                local_seg(1).segment = untrim_segment(i).collisionID;
            else
                local_seg_id = 1;
                temp_pointer = 1;
                for ii = 1:no_local_seg-1
                    breakpoint = breakpoint_list(ii);
                    local_seg(local_seg_id).segment = untrim_segment(i).collisionID(temp_pointer:breakpoint);
                    local_seg_id = local_seg_id + 1;
                    temp_pointer = breakpoint+1;
                    if ii == no_local_seg - 1
                        local_seg(local_seg_id).segment = untrim_segment(i).collisionID(temp_pointer:end);
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for id = 1:size(local_seg,2)
                if 1 < id && id < size(local_seg,2)
                    is_start = true;
                    seg_cont = false;
                    end_here = true;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % TO DO: Need to redeclare start_here, seg_cont, and
                    % end_here for the second/ n>1 th local segment
                    % IDEA: Should follow previous rules, e.g.: If end ==
                    % size AND next.collisionID(1) == 1...
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                elseif id == size(local_seg,2) && 1 < size(local_seg,2)
                    if untrim_segment(i).collisionID(end) == size(untrim_segment(i).segment,1) % Collision til the end of the current segment
                        if i ~= size(untrim_segment,2)
                            if isempty(untrim_segment(i+1).collisionID)
                                is_start = true;
                                end_here = true; % If next segment dont have collision --> end_here = true
                                seg_cont = false;
                            else
                                if untrim_segment(i+1).collisionID(1) == 1 % If next segment collision start at first point --> seg_cont = true; end_here = false
                                    is_start = true;
                                    seg_cont = true;
                                    end_here = false;
                                else
                                    is_start = true;
                                    seg_cont = false; % If not, current segment = end --> seg_cont = false; end_here = true;
                                    end_here = true;
                                end
                            end
                        else % If current the the last segment --> seg_cont = false; end_here = true
                            is_start = true;
                            seg_cont = false;
                            end_here = true;
                        end
                    else % Collision END before the end of the current segment
                        is_start = true;
                        seg_cont = false;
                        end_here = true;
                    end
                elseif id == 1 && 1 < size(local_seg,2)
                    seg_cont = false;
                    end_here = true;
                end
            for ii = 1:size(local_seg(id).segment,1)
                point_id = local_seg(id).segment(ii);

                plot3(untrim_segment(i).segment(point_id,1),untrim_segment(i).segment(point_id,2),untrim_segment(i).segment(point_id,3),'g.');

                vector_list = zeros(size(untrim_segment(i).closest(point_id).id,1),3);
                % Calculate vectors for the closest 10% surface [vector_list]
                for j = 1:size(vector_list,1)
                    face_id = untrim_segment(i).closest(point_id).id(j);
                    temp_vector = normal(face_id,:);
                    temp_vector = temp_vector*(1/untrim_segment(i).closest(point_id).dist(j));
                    vector_list(j,:) = temp_vector;
                end
                temp_vector = [];
                closest_collide_pt = [];
                collided_point = [];
                escaping_normal = [];
                old_point = untrim_segment(i).segment(point_id,:);
                % temp_vector : Collision occurs
                for j = 1:size(untrim_segment(i).collisionFace(point_id).id,2)
                    face_id = untrim_segment(i).collisionFace(point_id).id(j);
                    temp_centroid = [ccentroid(face_id,1) ccentroid(face_id,2) ccentroid(face_id,3)];
                    % closest_collide_pt = [closest_collide_pt; temp_centroid];
                    temp_normal = (old_point - temp_centroid)/norm(old_point - temp_centroid);
                    escaping_normal = [escaping_normal; temp_normal];
                    temp_vector = [temp_vector; normal(face_id,:)];
                    temp_dist = norm(old_point - temp_centroid);
                    if isempty(closest_collide_pt)
                        closest_collide_pt = face_id;
                        closest_collide_dist = temp_dist;
                    elseif temp_dist < closest_collide_dist
                        closest_collide_pt = face_id;
                        closest_collide_dist = temp_dist;
                    end
                    if i == 372
                        % plot3(ccentroid(face_id,1),ccentroid(face_id,2),ccentroid(face_id,3),'y.','MarkerSize',10);
                    end
                end
                closest_collide_normal = normal(closest_collide_pt,:);
                opposite_side_id = [];
                for j = 1:size(temp_vector,1)
                    temp_normal = temp_vector(j,:);
                    dot_product = dot(closest_collide_normal, temp_normal);
                    dot_product = acosd(dot_product/(norm(closest_collide_normal)*norm(temp_normal)));
                    if dot_product > 135
                        opposite_side_id = [opposite_side_id; j];
                    end
                end
                temp_vector(opposite_side_id,:) = [];
                escaping_normal(opposite_side_id,:) = [];
                temp_vector = sum(temp_vector,1)/norm(sum(temp_vector,1));
                escaping_normal = sum(escaping_normal,1)/norm(sum(escaping_normal,1));
                % 1) Find the closest point
                % 2) For all collide point, compute the dot product
                % 3) If dot product > 135 --> Opposite side --> discarded

                closest_pointID = untrim_segment(i).closest(point_id).id;
                closest_centroidLIST = ccentroid(closest_pointID,:);
                closest_vertexLIST = connectivity(closest_pointID,:);
                closest_vertexLIST = reshape(closest_vertexLIST,[],1);
                closest_vertexLIST = Vertex(closest_vertexLIST,:);
                closest_pointlist = [closest_centroidLIST; closest_vertexLIST]; % Include centroid and vertex of the closest 10% points
                % vector_list = [vector_list;vector_list;vector_list;vector_list];
                
                colli_type = 2;
                giga_vector = [escaping_normal; vector_list];
                replanned_point = pointwise_replan(old_point, giga_vector, offset_distance, ground_lvl, closest_pointlist, colli_type,Vertex, ...
                    connectivity, ALL_grid, grid_INFO, ccentroid, normal);

                temp_replanSeg = [temp_replanSeg; replanned_point];
                Collision_IDlist = [Collision_IDlist; ii];
                segment_IDlist = [segment_IDlist;point_id];
                plot3(replanned_point(1),replanned_point(2),replanned_point(3),'b.');

                giga = 1;
            end
            temp_replanEND = [temp_replanEND; replanned_point];
            % disp('is_start = '); disp(is_start);
            % disp('seg_cont = '); disp(seg_cont);
            % disp('end_here = ');disp(end_here);
            giga = 1;
            if is_start
                Prev_seg = i;
                % ADD INSERTION POINT
                Insertion_pt = segment_IDlist(1)-1;
                if end_here
                    Insertion_END = segment_IDlist(1);
                end
            end
            if end_here
                Next_seg = i;
                replan_segment(replan_segID).segment = temp_replanSeg;
                replan_segment(replan_segID).Head = Prev_seg;
                replan_segment(replan_segID).Tail = Next_seg;
                replan_segmentEND(replan_segID).segment = temp_replanEND;
                replan_segmentEND(replan_segID).Head = Prev_seg;
                replan_segmentEND(replan_segID).Tail = Next_seg;
                if ~is_start
                    Insertion_END = segment_IDlist(1);
                end
                replan_segmentEND(replan_segID).insertion_pt = Insertion_pt;
                replan_segmentEND(replan_segID).insertion_end = Insertion_END;
                temp_replanSeg = [];
                temp_replanEND = [];
                replan_segID = replan_segID + 1;
            end
            plot3(replanned_point(1),replanned_point(2),replanned_point(3),'r.','MarkerSize',20);
            untrim_segment(i).collisionID(Collision_IDlist) = [];
            if id ~= size(local_seg,2) && 1 < size(local_seg,2)
                % if segment_IDlist(1) ~= 1
                    % GIGA = 1;
                    % local_seg(id+1).segment = local_seg(id+1).segment - max(segment_IDlist)+ min(segment_IDlist)-1; % - size(segment_IDlist,1)
                % else
                    if abs(- max(segment_IDlist)+ min(segment_IDlist)-1) ~= size(segment_IDlist,1)
                        GIGA = 1;
                    end
                    for kk = 1:size(local_seg,2)
                        local_seg(kk).segment = local_seg(kk).segment - size(segment_IDlist,1);
                    end
                    
                    untrim_segment(i).collisionID = untrim_segment(i).collisionID - size(segment_IDlist,1);
                % end
            end
            Collision_IDlist = [];
            untrim_segment(i).segment(segment_IDlist,:) = [];
            untrim_segment(i).closest(segment_IDlist) = [];
            segment_IDlist = [];
            giga =1;
            end
            clear local_seg
        end

    end
    if i == size(untrim_segment,2)
        giga = 1;
    end
end
end