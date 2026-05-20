function update_point = pointwise_replan(old_point, vector, offset_distance, ground_lvl, point_list,colli_type, Vertex, ...
    connectivity, ALL_grid, grid_INFO, ccentroid, normal)

prim_dir = vector(1,:); % vector(1,:) = temp_vector --> Collision occur, escaping vector
second_dir = sum(vector(2:end,:));
second_dir = second_dir/norm(second_dir); % second_dir --> closest 10% points
vector_list = vector(2:end,:);

Q1_bound = [grid_INFO(1) grid_INFO(2) grid_INFO(3)];
Q3_bound = [grid_INFO(4) grid_INFO(5) grid_INFO(6)];
grid_size = grid_INFO(end);

if colli_type == 1

    closest_point = [old_point(1) old_point(2) 0];
    % move along  positive Z, calculate new coordinate
    [update_point, ~] = calculateNewCoordinate(old_point, prim_dir, closest_point, 0, ground_lvl);
    no_of_iter = 1;

    % Collision check of replanned point
    while true
    disp("No of Iteration in Replanning = ")
    disp(no_of_iter)
       if update_point(3) < ground_lvl
            GIGA = 1;
        end

    [gridID,neighbour] = getGridID(update_point, grid_size, Q3_bound(1), Q1_bound(1), Q3_bound(2), Q1_bound(2), Q3_bound(3), Q1_bound(3),3,true);
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
        temp_dist = norm(update_point-temp_centroid);
        % Dont need to find closest, only check for collision
        if round(temp_dist - offset_distance,1) < 0
            Collided_point = [Collided_point; temp_centroid];
            dist_list = [dist_list; temp_dist];
        end
    end
    [~,closest] = min(dist_list);


    if ~isempty(closest) %|| update_point(3) < ground_lvl
        % if ~isempty(closest) && ~(update_point(3) < ground_lvl) % Closest is NOT empty, but ground clearance ok
            % escape_dir = normal(all(ccentroid == Collided_point(closest,:),2),:);
            % escape_dir = (update_point-Collided_point(closest,:))/norm(update_point-Collided_point(closest,:));
            % if escape_dir(3) < 0
            %     escape_dir(3) = 0;
            % end
            need_replan = true;
    else
        need_replan = false;
    end

    if no_of_iter > 20
        need_replan = false;
        if ~isempty(Collided_point)
            ERROR = 1;
        end
    end

    if need_replan == false
        if update_point(3) < ground_lvl
            ERROR = 1;
        end
        break
    end


    normal_list = [];
    weighted_normal = [];
    elite_weighted_normal = [];
    weight = [];
    elite = floor(size(Collided_point,1)*0.1);
    if elite < 1
        elite = 1;
    end

        GIGA = 1;
        normal_list = [];
        for jj = 1:size(Collided_point,1)
            % id = all(ccentroid == Collided_point(jj,:),2);
            % % normal_list = [normal_list; normal(id,:)]; % Normal should should be obtain between the RELATIVE POSITION BETWEEN THE POINTS,
            % % NOT SURFACE NORMAL
            % temp_normal = (update_point - Collided_point(jj,:))/norm(update_point - Collided_point(jj,:));
            % normal_list = [normal_list; temp_normal];

            % id = all(ccentroid == Collided_point(jj,:),2);
            temp_normal = (update_point - Collided_point(jj,:))/norm(update_point - Collided_point(jj,:));
            normal_list = [normal_list; temp_normal];
            % normal_list = [normal_list; normal(id,:)];
            % total_weight = 1./(dist_list);
            temp_ratio = dist_list(jj)/sum(dist_list);
            temp_weight = 1/temp_ratio;
            weight = [weight; temp_weight];
            temp_weighted_normal = temp_weight*temp_normal;
            weighted_normal = [weighted_normal; temp_weighted_normal];
        end
        [sorted_weight, sorted_id] = sort(weight,'descend');
        sorted_normal = normal_list(sorted_id,:);
        for jj = 1:elite
            temp_weighted_normal = sorted_weight(jj)*sorted_normal(jj,:);
            elite_weighted_normal = [elite_weighted_normal; temp_weighted_normal];
        end
        elite_weighted_normal = sum(elite_weighted_normal,1)/norm(sum(elite_weighted_normal,1));
        pause = 1;

        escape_dir = sum(normal_list,1)/norm(sum(normal_list,1));
        if escape_dir(3) < 0
            escape_dir(3) = 0;
        end
        pause = 1;
        escape_dir = elite_weighted_normal;

    %%%%% NEED to redo collision detection after every replanning %%%%%%
    if inpolyhedron(connectivity, Vertex, update_point)
        GIGA = 1;
    end
    

    % second_dir should be obtained from replanned point, NOT og
    % The same for point_list

    %%%%%%%%%%% Original %%%%%%%%%%%
    % move along xy plane, calculate new coordinate
    % second_dir(3) = 0;
    % [update_point, ~] = calculateNewCoordinate(update_point, second_dir, point_list, offset_distance, ground_lvl);
    %%%%%%%%%%% End of Original %%%%%%%%%%%
    prev_replan_point = update_point;
    [update_point, ~] = calculateNewCoordinate(update_point, escape_dir, neighbour_centroid, offset_distance, ground_lvl);
    no_of_iter = no_of_iter + 1;
    end
    if update_point(3) < ground_lvl
        GIGA = 1;
    end

    giga = 1;


elseif colli_type == 2
    [update_point, ~] = calculateNewCoordinate(old_point, prim_dir, point_list, offset_distance, ground_lvl);
    newDistances = vecnorm(point_list - update_point, 2, 2);
    % while still have collision
    % Retrive the colliding point and corresponding normal
    % If collision with vertex --> ignore
    % Replan based on those point and normal only
    no_of_iter = 1;
    while true
        disp("No of Iteration in Replanning = ")
        disp(no_of_iter)


        [gridID,neighbour] = getGridID(update_point, grid_size, Q3_bound(1), Q1_bound(1), Q3_bound(2), Q1_bound(2), Q3_bound(3), Q1_bound(3),3,true);
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
            temp_dist = norm(update_point-temp_centroid);
            % Dont need to find closest, only check for collision
            if round(temp_dist - offset_distance,1) < 0
                Collided_point = [Collided_point; temp_centroid];
                dist_list = [dist_list; temp_dist];
            end
        end
        [~,closest] = min(dist_list);
    
        if ~isempty(closest)
            % escape_dir = normal(all(ccentroid == Collided_point(closest,:),2),:);
            escape_dir = (update_point-Collided_point(closest,:))/norm(update_point-Collided_point(closest,:));

            need_replan = true;
        else
            need_replan = false;
        end

        if no_of_iter > 30
            need_replan = false;
            if ~isempty(Collided_point)
                ERROR = 1;
            end
        end
    
        if need_replan == false
            if update_point(3) < ground_lvl
                ERROR = 1;
            end
            break
        end

        if no_of_iter > 1
            % if all(round(prev_replan_point,4) == round(update_point,4))
                GIGA = 1;
                normal_list = [];
                weighted_normal = [];
                elite_weighted_normal = [];
                weight = [];
                elite = floor(size(Collided_point,1)*(floor(no_of_iter/5)+1)*0.1);
                if elite < 1
                    elite = 1;
                end
                for jj = 1:size(Collided_point,1)
                    % id = all(ccentroid == Collided_point(jj,:),2);
                    temp_normal = (update_point - Collided_point(jj,:))/norm(update_point - Collided_point(jj,:));
                    normal_list = [normal_list; temp_normal];
                    % normal_list = [normal_list; normal(id,:)];
                    % total_weight = 1./(dist_list);
                    temp_ratio = dist_list(jj)/sum(dist_list);
                    temp_weight = 1/temp_ratio;
                    weight = [weight; temp_weight];
                    temp_weighted_normal = temp_weight*temp_normal;
                    weighted_normal = [weighted_normal; temp_weighted_normal];
                end
                [sorted_weight, sorted_id] = sort(weight,'descend');
                sorted_normal = normal_list(sorted_id,:);
                if size(sorted_weight,1) == 0 
                    ERROR = 1;
                end
                for jj = 1:elite
                    temp_weighted_normal = sorted_weight(jj)*sorted_normal(jj,:);
                    elite_weighted_normal = [elite_weighted_normal; temp_weighted_normal];
                end
                elite_weighted_normal = sum(elite_weighted_normal,1)/norm(sum(elite_weighted_normal,1));
                pause = 1;
                % if no_of_iter > 10
                %     escape_dir = sum(normal_list,1)/norm(sum(normal_list,1));
                %     % escape_dir(3) = 0;
                % else
                %     escape_dir = sum(normal_list,1)/norm(sum(normal_list,1));
                % end
                escape_dir = elite_weighted_normal;
                if size(escape_dir,2) < 3
                    ERROR = 1;
                end
                if update_point(3) < ground_lvl
                    escape_dir(3) = abs(escape_dir(3));
                end
    
            % end
        end


        %%%%% NEED to redo collision detection after every replanning %%%%%%
        if inpolyhedron(connectivity, Vertex, update_point)
            GIGA = 1;
        end
        
    
        % second_dir should be obtained from replanned point, NOT og
        % The same for point_list
    
        prev_replan_point = update_point;
        [update_point, ~] = calculateNewCoordinate(update_point, escape_dir, neighbour_centroid, offset_distance, ground_lvl);
        no_of_iter = no_of_iter + 1;
    end
    % while any(newDistances - offset_distance < 0)
    %     ID = find(newDistances - offset_distance < 0);
    %     new_vector_list = [];
    %     for i = 1:size(ID,1)
    %         temp_vector = vector_list(i,:);
    %         new_vector_list = [new_vector_list; temp_vector];
    %     end
    %     new_vector = sum(new_vector_list,1)/norm(sum(new_vector_list,1));
    %     [update_point, ~] = calculateNewCoordinate(update_point, new_vector, point_list, offset_distance, ground_lvl);
    %     if ~inpolyhedron(connectivity, Vertex, update_point)
    %         GIGA = 1;
    %     end
    %     newDistances = vecnorm(point_list - update_point, 2, 2);
    %     iter = iter + 1;
    %     % if iter > 10
    %     %     giga = 1;
    %     %     break
    %     % end
    % end
end
% Calculate new coordinate
% [update_point, stepSize] = calculateNewCoordinate(old_point, vector, point_list, offset_distance);

% Calculate the distances between the current coordinate and the point_list
currentDistances = vecnorm(point_list - old_point, 2, 2);

% Calculate the distances between the updated coordinate and the point_list
newDistances = vecnorm(point_list - update_point, 2, 2);

% Calculate the current collision distance and the updated distance
currentCollisionDistance = min(currentDistances) - offset_distance;
updatedDistance = min(newDistances) - offset_distance;
giga = 1;

    function [newCoordinate, stepSize] = calculateNewCoordinate(old_point, vector, point_list, offset_distance, Ground_lvl)
    % Calculate the new coordinate based on the given direction vector and optimize the step size
    % while satisfying the collision constraint
    
    % Define the objective function
    fun = @(stepSize) distanceChanged(old_point, vector, stepSize);
    
    % Set the initial guess for the step size
    stepSize0 = offset_distance; % You can adjust the initial guess as needed
    
    % Define the lower and upper bounds for the step size
    lb = 0; % Lower bound of step size (non-negative)
    ub = Inf; % Upper bound of step size (unbounded)
    
    % Define the nonlinear constraint function
    nonlcon = @(stepSize) collisionConstraint(old_point, vector, point_list, stepSize, offset_distance, Ground_lvl);
    
    % Set the optimization options
    options = optimoptions('fmincon','MaxIterations',100, 'Display','off');
    
    % Perform the optimization to find the step size that minimizes the distance changed
    [stepSize,~,exitflag,~] = fmincon(fun, stepSize0, [], [], [], [], lb, ub, nonlcon, options);
    if exitflag == -1 || exitflag == -2
        ERROR = 1;
    end
    
    % Calculate the displacement vector by multiplying the direction vector with the optimized step size
    displacementVector = vector / norm(vector) * stepSize;
    
    % Calculate the new coordinate by adding the displacement vector to the current coordinate
    newCoordinate = old_point + displacementVector;
end

function distance = distanceChanged(old_point, vector, stepSize)
    % Calculate the distance changed by moving in the given direction with the specified step size
    
    % Normalize the direction vector
    vector = vector / norm(vector);
    
    % Calculate the displacement vector by multiplying the direction vector with the step size
    displacementVector = vector * stepSize;
    
    % Calculate the new coordinate by adding the displacement vector to the current coordinate
    newCoordinate = old_point + displacementVector;
    
    % Calculate the distance changed
    distance = norm(newCoordinate - old_point);
end

    function [c, ceq] = collisionConstraint(old_point, vector, point_list, stepSize, offset_distance, Ground_lvl)
    % Collision constraint: Ensure the new distance between the updated coordinates and point_list is greater than offset_distance
    
    % Calculate the displacement vector by multiplying the direction vector with the step size
    displacementVector = vector / norm(vector) * stepSize;
    
    % Calculate the new coordinate by adding the displacement vector to the current coordinate
    newCoordinate = old_point + displacementVector;
    
    % Calculate the distances between the new coordinate and the point_list
    % v = newCoordinate - point_list;
    % dotProduct = dot(v, vector_list,2);

    distances1 = vecnorm(point_list - newCoordinate, 2, 2);
    distances2 = newCoordinate(3);
    
    % Calculate the constraint violation
    violation1 = offset_distance - min(distances1);
    violation2 = Ground_lvl - distances2;
    violation = max([violation1 violation2]);
    
    % Set the constraint and equality values
    c = violation; % Inequality constraint: c <= 0
    ceq = []; % No equality constraints
end
end