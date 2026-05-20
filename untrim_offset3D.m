function  [untrim_segments] = untrim_offset3D(OG_segment, offset_distance, Normal, ccentroid)
N = 50;
circr = @(radius,rad_ang1, rad_ang2)  [radius*cos(rad_ang1).*cos(rad_ang2);  radius*sin(rad_ang1).*cos(rad_ang2); radius*sin(rad_ang2)]; 
for i = 1:size(OG_segment,2)
    temp_faceID = OG_segment(i).FaceID;
    temp_traj = OG_segment(i).segment;
    offset_direction = Normal(temp_faceID,:);
    xy_angle = atan2(offset_direction(2), offset_direction(1));
    z_angle = atan2(offset_direction(3),sqrt(offset_direction(1)^2+offset_direction(2)^2));
    temp_offset_traj = zeros(size(temp_traj));
    for ii = 1:size(temp_traj,1)
        temp_x = temp_traj(ii,1) + (offset_distance*cos(z_angle))*cos(xy_angle);
        temp_y = temp_traj(ii,2) + (offset_distance*cos(z_angle))*sin(xy_angle);
        temp_z = temp_traj(ii,3) + offset_distance*sin(z_angle);
        temp_offset_traj(ii,:) = [temp_x temp_y temp_z];
    end
    untrim_segments(i).segment = temp_offset_traj;
    untrim_segments(i).FaceID = temp_faceID;
    if i == 394
        giga = 1;
    end

    if i ~= 1
        segment1 = untrim_segments(i-1).segment;
        segment2 = untrim_segments(i).segment;
        % plot3(segment1(:,1),segment1(:,2), segment1(:,3),'g');
        % plot3(segment2(:,1), segment2(:,2),segment2(:,3),'m');
        % drawnow
        intersection = [];
        for j = 1:size(segment1,1)-1
            line1_startx = segment1(j,1);
            line1_starty = segment1(j,2);
            line1_startz = segment1(j,3);
            line1_endx = segment1(j+1,1);
            line1_endy = segment1(j+1,2);
            line1_endz = segment1(j+1,3);
            dist1 = sqrt((line1_endx-line1_startx)^2+(line1_endy-line1_starty)^2+(line1_endz-line1_startz)^2);
            for jj = 1:size(segment2,1)-1
                line2_startx = segment2(jj,1);
                line2_starty = segment2(jj,2);
                line2_startz = segment2(jj,3);
                line2_endx = segment2(jj+1,1);
                line2_endy = segment2(jj+1,2);
                line2_endz = segment2(jj+1,3);
                dist2 = sqrt((line2_endx-line2_startx)^2+(line2_endy-line2_starty)^2+(line2_endz-line2_startz)^2);
                startpoints = [line1_startx line1_starty line1_startz; line2_startx line2_starty line2_startz];
                endpoints = [line1_endx line1_endy line1_endz; line2_endx line2_endy line2_endz];
                [point, dist] = lineIntersect3D(startpoints, endpoints);
                if ~((any(dist > dist1) || any(dist > dist2))) && any((~isnan(dist)))
                    intersection = [intersection; point];
                end
            end
        end
        if isempty(intersection) % Two adjacent line segment does not intersect --> connect with arc
            flag = 0;
            % quiver3(ccentroid(untrim_segments(i).FaceID,1),ccentroid(untrim_segments(i).FaceID,2),ccentroid(untrim_segments(i).FaceID,3),Normal(untrim_segments(i).FaceID,1),Normal(untrim_segments(i).FaceID,2),Normal(untrim_segments(i).FaceID,3));
            % quiver3(ccentroid(untrim_segments(i-1).FaceID,1),ccentroid(untrim_segments(i-1).FaceID,2),ccentroid(untrim_segments(i-1).FaceID,3),Normal(untrim_segments(i-1).FaceID,1),Normal(untrim_segments(i-1).FaceID,2),Normal(untrim_segments(i-1).FaceID,3));
            mean_normal = [mean([Normal(untrim_segments(i).FaceID,1) Normal(untrim_segments(i-1).FaceID,1)]) mean([Normal(untrim_segments(i).FaceID,2) Normal(untrim_segments(i-1).FaceID,2)]) mean([Normal(untrim_segments(i).FaceID,3) Normal(untrim_segments(i-1).FaceID,3)])];
            segment_endx = OG_segment(i-1).segment(end,1);
            segment_endy = OG_segment(i-1).segment(end,2);
            segment_endz = OG_segment(i-1).segment(end,3);
            % quiver3(segment_endx, segment_endy, segment_endz, mean_normal(1), mean_normal(2), mean_normal(3));
            % plot3(segment_endx, segment_endy, segment_endz,'r.');
            % drawnow

            direction1 = Normal(untrim_segments(i-1).FaceID,:);
            angle_xy = atan2(direction1(2), direction1(1));
            angle_z = atan2(direction1(3),sqrt(direction1(1)^2+direction1(2)^2));
            point1_x = segment_endx + (offset_distance*cos(angle_z))*cos(angle_xy);
            point1_y = segment_endy + (offset_distance*cos(angle_z))*sin(angle_xy);
            point1_z = segment_endz + offset_distance*sin(angle_z);
            % plot3(point1_x, point1_y, point1_z,'g.');

            
            direction2 = mean_normal;
            angle_xy = atan2(direction2(2), direction2(1));
            angle_z = atan2(direction2(3),sqrt(direction2(1)^2+direction2(2)^2));
            point2_x = segment_endx + (offset_distance*cos(angle_z))*cos(angle_xy);
            point2_y = segment_endy + (offset_distance*cos(angle_z))*sin(angle_xy);
            point2_z = segment_endz + offset_distance*sin(angle_z);
            % plot3(point2_x, point2_y, point2_z,'g.');

            direction3 = Normal(untrim_segments(i).FaceID,:);
            angle_xy = atan2(direction3(2), direction3(1));
            angle_z = atan2(direction3(3),sqrt(direction3(1)^2+direction3(2)^2));
            point3_x = segment_endx + (offset_distance*cos(angle_z))*cos(angle_xy);
            point3_y = segment_endy + (offset_distance*cos(angle_z))*sin(angle_xy);
            point3_z = segment_endz + offset_distance*sin(angle_z);
            % plot3(point3_x, point3_y, point3_z,'g.');

            % drawnow
            % giga = 1;
            arc = [point1_x point1_y point1_z];
            if ~all([point1_x point1_y point1_z] == [point2_x point2_y point2_z])
                temp_seg = [linspace(point1_x, point2_x,20)' linspace(point1_y, point2_y,20)' linspace(point1_z, point2_z,20)'];
                arc(end,:) = [];
                arc = [arc; temp_seg];
            end

            if ~all([point1_x point1_y point1_z] == [point3_x point3_y point3_z]) && ~all([point2_x point2_y point2_z] == [point3_x point3_y point3_z])
                temp_seg = [linspace(point2_x, point3_x,20)' linspace(point2_y, point3_y,20)' linspace(point2_z, point3_z,20)'];
                arc(end,:) = [];
                arc = [arc; temp_seg];
            end
            % OG General Version %
            % angle1_xy = atan2(segment1(end,2)-segment_endy, segment1(end,1)-segment_endx);
            % angle1_z = atan2(segment1(end,3)-segment_endz, sqrt((segment1(end,2)-segment_endy)^2+(segment1(end,1)-segment_endx)^2));
            % angle2_xy = atan2(segment2(1,2)-segment_endy, segment2(1,1)-segment_endx);
            % angle2_z = atan2(segment2(1,3)-segment_endz,sqrt((segment2(1,2)-segment_endy)^2+(segment2(1,1)-segment_endx)^2));
            % if (angle1_xy > 0 && angle2_xy <0) || (angle1_xy < 0 && angle2_xy > 0)
            %     direction1 = linspace(angle1_xy, angle2_xy, N);
            %     if angle1_xy < 0
            %         angle1_xy_alt = angle1_xy +2*pi;
            %         direction2 = linspace(angle1_xy_alt, angle2_xy, N);
            %     elseif angle2_xy < 0
            %         angle2_xy_alt = angle2_xy + 2*pi;
            %         direction2 = linspace(angle1_xy, angle2_xy_alt, N);
            %     end
            %     change1 = abs(direction1(2)-direction1(1));
            %     change2 = abs(direction2(2)-direction2(1));
            %     if change1 <= change2 % Original direction is shortest
            %         angle_xy = direction1;
            %     elseif change2 < change1 % Alternative direction is shortest
            %         angle_xy = direction2;
            %     end
            % else
            %     angle_xy = linspace(angle1_xy, angle2_xy,N);
            % end
            % angle_z = linspace(angle1_z, angle2_z,N);

            % Too narrow separation of different cases % 
            % if sum(segment1(end,:) == [segment_endx segment_endy segment_endz])==2 || sum(segment2(1,:) == [segment_endx segment_endy segment_endz])==2 % offset segment and OG segment have the same xy-coordinates
            %      if all(segment1(end,1:2) == [segment_endx segment_endy]) || all(segment2(1,1:2) == [segment_endx segment_endy])
            %          if round(segment1(end,3),4) == round(segment2(1,3),4)
            %              flag = 2;
            %              x = linspace(segment1(end,1), segment2(1,1),N);
            %              y = linspace(segment1(end,2), segment2(1,2),N);
            %              z = linspace(segment1(end,3), segment1(end,3),N);
            %          else
            %              if angle1_xy == 0
            %                  angle_xy = linspace(angle2_xy, angle2_xy,N);
            %                  angle_z = linspace(angle1_z, angle2_z,N);
            %              elseif angle2_xy == 0
            %                  angle_xy = linspace(angle1_xy, angle1_xy,N);
            %                  angle_z = linspace(angle1_z, angle2_z,N);
            %              end
            %          end
            %      else
            %          flag = 1;
            %          if round(segment1(end,3),4) == round(segment2(1,3),4) && round(segment1(end,3),4) == round(segment_endz,4)
            %              flag = 1.5;
            %             angle_z = linspace(angle1_z, angle2_z,N);
            %             if (angle1_xy > 0 && angle2_xy <0) || (angle1_xy < 0 && angle2_xy > 0)
            %                 direction1 = linspace(angle1_xy, angle2_xy, N);
            %                 if angle1_xy < 0
            %                     angle1_xy_alt = angle1_xy +2*pi;
            %                     direction2 = linspace(angle1_xy_alt, angle2_xy, N);
            %                 elseif angle2_xy < 0
            %                     angle2_xy_alt = angle2_xy + 2*pi;
            %                     direction2 = linspace(angle1_xy, angle2_xy_alt, N);
            %                 end
            %                 change1 = abs(direction1(2)-direction1(1));
            %                 change2 = abs(direction2(2)-direction2(1));
            %                 if change1 <= change2 % Original direction is shortest
            %                     angle_xy = direction1;
            %                 elseif change2 < change1 % Alternative direction is shortest
            %                     angle_xy = direction2;
            %                 end
            %             else
            %                 angle_xy = linspace(angle1_xy, angle2_xy,N);
            %             end
            %          else
            %              angle_z = linspace(angle1_z, angle2_z,N);
            %              x = linspace(segment1(end,1), segment2(1,1),N);
            %              y = linspace(segment1(end,2), segment2(1,2),N);
            %          end
            %      end
            % else
            %     if (angle1_xy > 0 && angle2_xy <0) || (angle1_xy < 0 && angle2_xy > 0)
            %         direction1 = linspace(angle1_xy, angle2_xy, N);
            %         if angle1_xy < 0
            %             angle1_xy_alt = angle1_xy +2*pi;
            %             direction2 = linspace(angle1_xy_alt, angle2_xy, N);
            %         elseif angle2_xy < 0
            %             angle2_xy_alt = angle2_xy + 2*pi;
            %             direction2 = linspace(angle1_xy, angle2_xy_alt, N);
            %         end
            %         change1 = abs(direction1(2)-direction1(1));
            %         change2 = abs(direction2(2)-direction2(1));
            %         if change1 <= change2 % Original direction is shortest
            %             angle_xy = direction1;
            %         elseif change2 < change1 % Alternative direction is shortest
            %             angle_xy = direction2;
            %         end
            %     else
            %         angle_xy = linspace(angle1_xy, angle2_xy,N);
            %         angle_z = linspace(angle1_z, angle2_z,N);
            %     end
            % end
            
            % OG General create arc %
            % if flag == 1
            %     arc = circr(offset_distance, angle_xy, angle_z);
            %     arc(1,:) = arc(1,:) + segment_endx;
            %     arc(2,:) = arc(2,:) + segment_endy;
            %     arc(3,:) = arc(3,:) + segment_endz;
            %     arc(1,:) = x;
            %     arc(2,:) = y;
            % elseif flag == 1.5
            %     arc = circr(offset_distance, angle_xy, angle_z);
            %     arc(1,:) = arc(1,:) + segment_endx;
            %     arc(2,:) = arc(2,:) + segment_endy;
            %     arc(3,:) = arc(3,:) + segment_endz;
            % elseif flag == 2
            %     arc(1,:) = x;
            %     arc(2,:) = y;
            %     arc(3,:) = z;
            % elseif flag == 0
            %     arc = circr(offset_distance, angle_xy, angle_z);
            %     arc(1,:) = arc(1,:) + segment_endx;
            %     arc(2,:) = arc(2,:) + segment_endy;
            %     arc(3,:) = arc(3,:) + segment_endz;
            % end
            % % plot3(arc(1,:), arc(2,:), arc(3,:),'r')
            % arc(:,1) = [];
            % arc(:,end) = [];
            % 
            untrim_segments(i).segment = [arc; untrim_segments(i).segment];
            % drawnow
            % if i == 73
            %     giga = 1;
            % end
        end
    end
end
end