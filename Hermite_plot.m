function interpPoints = Hermite_plot(points)
    % Number of points
    n = size(points, 1);

    % Calculate tangent vectors using finite differencing
    tangents = zeros(n, 3);

    % Calculate tangent vectors for interior points
    for i = 2:n-1
        tangents(i, :) = (points(i+1, :) - points(i-1, :)) / norm(points(i+1, :) - points(i-1, :));
    end

    % Estimate tangent vectors for the first and last points
    tangents(1, :) = (points(2, :) - points(1, :)) / norm(points(2, :) - points(1, :));
    tangents(n, :) = (points(n, :) - points(n-1, :)) / norm(points(n, :) - points(n-1, :));

    % Create parameter t ranging from 0 to 1
    t = linspace(0, 1, n);

    % Initialize arrays to store interpolated points
    interpPoints = zeros(5000, 3);  % Increase 100 if you want more interpolated points

    % Generate Hermite spline for each segment
    for i = 1:n-1
        % Create a parameter t within the range of the segment
        ti = linspace(t(i), t(i+1), 100);  % Increase 100 for more interpolated points

        % Calculate Hermite blending functions
        h00 = 2*ti.^3 - 3*ti.^2 + 1;
        h01 = ti.^3 - 2*ti.^2 + ti;
        h10 = -2*ti.^3 + 3*ti.^2;
        h11 = ti.^3 - ti.^2;

        % Calculate interpolated points using Hermite spline formula
        interpPoints((i-1)*100+1:i*100, :) = h00' * points(i, :) + h01' * tangents(i, :) + h10' * points(i+1, :) + h11' * tangents(i+1, :);
    end

    % Plot the original points and the interpolated points
    hold on;
    plot3(points(:, 1), points(:, 2), points(:, 3), '.','Color',[0.8500 0.3250 0.0980],'MarkerSize',15);  % Original points
    plot3(interpPoints(:, 1), interpPoints(:, 2), interpPoints(:, 3), 'color',[0 0.4470 0.7410],'LineWidth',2);  % Interpolated points
    hold off;
    grid on;
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    % legend('Original Points', 'Interpolated Points');
end