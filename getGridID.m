function [gridID, neighour_list] = getGridID(point, grid_size, x_lower, x_upper, y_lower, y_upper, z_lower, z_upper,dimension,find_neighbor)
    % point: 3D point coordinates [x, y, z]
    % grid_size: size of each grid in the discretized 3D space
    % x_lower, x_upper: lower and upper bounds for the x dimension
    % y_lower, y_upper: lower and upper bounds for the y dimension
    % z_lower, z_upper: lower and upper bounds for the z dimension

    neighour_list = [];


    if dimension == 3
        if point(3) < 0
            point(3) = 0;
        end
    % Calculate grid indices in each dimension
    x_size = (x_upper - x_lower) / grid_size;
    y_size = (y_upper - y_lower) / grid_size;
    z_size = (z_upper - z_lower) / grid_size;
    no_of_grid = x_size*y_size*z_size;
    grid_x = floor((point(1) - x_lower) / grid_size);
    grid_y = floor((point(2) - y_lower) / grid_size);
    grid_z = floor((point(3) - z_lower) / grid_size);

    % Calculate single grid ID
    gridID = grid_x + grid_y * ((x_upper - x_lower) / grid_size) + ...
        grid_z * ((x_upper - x_lower) / grid_size ) * ((y_upper - y_lower) / grid_size) + 1;

    if gridID < 0
        GIGA_ERROR = 1;
    end
    if find_neighbor
        neighour_list = [gridID+1; gridID-1; gridID+x_size; gridID+x_size+1; gridID+x_size-1; gridID-x_size; gridID-x_size+1; ... 
            gridID-x_size-1; gridID+(x_size*y_size); gridID+(x_size*y_size)+1; gridID+(x_size*y_size)-1; ...
            gridID+(x_size*y_size)+x_size; gridID+(x_size*y_size)+x_size+1; gridID+(x_size*y_size)+x_size-1; ... 
            gridID+(x_size*y_size)-x_size; gridID+(x_size*y_size)-x_size+1; gridID+(x_size*y_size)-x_size-1;... 
            gridID-(x_size*y_size); gridID-(x_size*y_size)+1; gridID-(x_size*y_size)-1; ... 
            gridID-(x_size*y_size)+x_size; gridID-(x_size*y_size)+x_size+1; gridID-(x_size*y_size)+x_size-1; ... 
            gridID-(x_size*y_size)-x_size; gridID-(x_size*y_size)-x_size+1; gridID-(x_size*y_size)-x_size-1];
        min_check = neighour_list < 0;
        max_check = neighour_list > no_of_grid;
        check = (min_check+max_check)>0;
        neighour_list(check) = [];
    end
    else
        grid_x = floor((point(1) - x_lower) / grid_size);
        grid_y = floor((point(2) - y_lower) / grid_size);
    
        % Calculate single grid ID
        gridID = grid_x + grid_y * ((x_upper - x_lower) / grid_size) + 1;
        if find_neighbor
            if grid_x == 0 % No Left neighbour
                neighour_list = [neighour_list; gridID+1];
            elseif grid_x == ((x_upper - x_lower) / grid_size)-1 % No Right neighbour
                neighour_list = [neighour_list; gridID-1];
            else % HAVE BOTH left and right neighbour
                neighour_list = [neighour_list; gridID-1; gridID+1];
            end

            if grid_y == 0 % NO down neighbour
                tempID = grid_x + (grid_y+1) * ((x_upper - x_lower) / grid_size) + 1;
                neighour_list = [neighour_list; tempID];
                    if grid_x == 0 % No Left neighbour
                        neighour_list = [neighour_list; tempID+1];
                    elseif grid_x == ((x_upper - x_lower) / grid_size)-1 % No Right neighbour
                        neighour_list = [neighour_list; tempID-1];
                    else % HAVE BOTH left and right neighbour
                        neighour_list = [neighour_list; tempID-1; tempID+1];
                    end
            elseif grid_y == ((y_upper - y_lower) / grid_size)-1 % NO up neighbour
                tempID = grid_x + (grid_y-1) * ((x_upper - x_lower) / grid_size) + 1;
                neighour_list = [neighour_list; tempID];
                    if grid_x == 0 % No Left neighbour
                        neighour_list = [neighour_list; tempID+1];
                    elseif grid_x == ((x_upper - x_lower) / grid_size)-1 % No Right neighbour
                        neighour_list = [neighour_list; tempID-1];
                    else % HAVE BOTH left and right neighbour
                        neighour_list = [neighour_list; tempID-1; tempID+1];
                    end
            else
                tempID1 = grid_x + (grid_y-1) * ((x_upper - x_lower) / grid_size) + 1; % UP neighbour
                tempID2 = grid_x + (grid_y+1) * ((x_upper - x_lower) / grid_size) + 1; % DOWN neighour
                neighour_list = [neighour_list; tempID1; tempID2];
                    if grid_x == 0 % No Left neighbour
                        neighour_list = [neighour_list; tempID1+1; tempID2+1];
                    elseif grid_x == ((x_upper - x_lower) / grid_size)-1 % No Right neighbour
                        neighour_list = [neighour_list; tempID1-1; tempID2-1];
                    else % HAVE BOTH left and right neighbour
                        neighour_list = [neighour_list; tempID1-1; tempID1+1; tempID2-1; tempID2+1];
                    end
            end
        end
    end
end