function isWithinViewCone = CalPointVist(viewpointLocation, viewingAngle, pointLocation, fovDegrees)
    % Calculate the vector from the viewpoint to the point
    vectorToPoint = pointLocation - viewpointLocation;
    
    % Normalize the vectors
    vectorToPoint = vectorToPoint / norm(vectorToPoint);
    viewingAngle = viewingAngle / norm(viewingAngle);
    
    % Calculate the dot product
    dotProduct = dot(vectorToPoint, viewingAngle);
    
    % Calculate the cosine of half the FOV
    cosHalfFOV = cosd(fovDegrees / 2);
    
    % Check if the dot product is greater than or equal to the cosine of half the FOV
    isWithinViewCone = dotProduct >= cosHalfFOV;
    
end