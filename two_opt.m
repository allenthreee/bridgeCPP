function optimizedPath = two_opt(path, n)
    numPoints = size(path, 1);  % Number of waypoints
    bestPath = path;    % Initialize the best path as the given path
    improved = true;    % Flag to indicate if improvement is made

    while improved
        improved = false;
        
        for i = 1:numPoints-2
            for j = i+2:min(i+n, numPoints-1)
                % Swap two edges and check if it improves the path length
                newPath = twoOptSwap(bestPath, i, j);
                
                if pathLength(newPath) < pathLength(bestPath)
                    bestPath = newPath;
                    improved = true;
                end
            end
        end
    end
    
    optimizedPath = bestPath;
end

function newPath = twoOptSwap(path, i, j)
    newPath = [path(1:i,:); path(j:-1:i+1,:); path(j+1:end,:)];
end

function length = pathLength(path)
    length = 0;
    numPoints = size(path, 1);
    
    for i = 1:numPoints-1
        length = length + norm(path(i,:) - path(i+1,:));
    end
    
    % Add the distance between the last and first waypoints
    length = length + norm(path(numPoints,:) - path(1,:));
end