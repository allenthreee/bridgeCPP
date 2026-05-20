function isAngleOfIncidenceLarge = detectLargeAngleOfIncidence(surfaceNormal, viewingDirection, threshold)
    % Normalize the surface normal and viewing direction
    surfaceNormal = surfaceNormal / norm(surfaceNormal);
    viewingDirection = viewingDirection / norm(viewingDirection);

    % Calculate the angle of incidence
    if size(surfaceNormal,1) ~= size(viewingDirection,1) || size(surfaceNormal,2) ~= size(viewingDirection,2)
        ERROR = 1;
    end

    dot_product = dot(surfaceNormal, viewingDirection);
    if dot_product >=0 
        isAngleOfIncidenceLarge = true;
    else
        angleOfIncidence = acosd(-dot_product);
        
        % Check if the angle of incidence is larger than the threshold
        isAngleOfIncidenceLarge = angleOfIncidence > threshold;
    end
end