function uv23dRatio = ratiocompute(centroid_uv, centroid_3D)
tic
if size(centroid_uv,1) ~= size(centroid_3D,1)
    ERROR = 1;
end
C2C_dist_uv = nan(size(centroid_uv,1),size(centroid_uv,1));
C2C_dist_3d = nan(size(centroid_uv,1),size(centroid_uv,1));
for i = 1:size(centroid_uv,1)
    for j = 1:size(centroid_uv,1)
        temp_dist = norm(centroid_uv(i,:)-centroid_uv(j,:));
        C2C_dist_uv(i,j) = temp_dist;
        if temp_dist ~= 0
            temp_dist = norm(centroid_3D(i,:)-centroid_3D(j,:));
            C2C_dist_3d(i,j) = temp_dist;
        end
    end
end
ratio = C2C_dist_uv./C2C_dist_3d;
uv23dRatio = mean(mean(ratio,"omitnan"));
toc
end