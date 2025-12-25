function vars = defectOrient1(orientation, xP, yP, xM, yM, gridSpacing, vars)
% defect coordinates and charge with both plus and minus defect
x_defect = [xP; xM];
y_defect = [yP; yM];
charge = [ones(length(xP), 1) * (1/2); ones(length(xM), 1) * (-1/2)];

radius =  2 * gridSpacing;      % Define the radius of the small circle around the defect
markerSize = 0.5 * gridSpacing;   % Define size of the plus/minus marker

% Create a grid for the orientation field
[X, Y] = meshgrid(1:size(orientation, 1), 1:size(orientation, 2));

defect_orientation = zeros(size(x_defect)); % Initialize the array to store defect orientations

% Loop over each defect to compute its orientation and visualize it
for i = 1:size(x_defect, 1)
    % Compute the distance from each grid point to the defect's position
    dist = sqrt((X - x_defect(i)).^2 + (Y - y_defect(i)).^2);
    
    % Create a mask for points within the radius around the defect
    mask = dist <= radius;
    
    % Get the orientations within the circle
    orientations_in_circle = orientation(mask);
    
    % Get the coordinates of the points inside the mask
    [circle_y, circle_x] = find(mask);
    
    % Compute the angle of each point relative to the defect's position
    phi = atan2(circle_y - y_defect(i), circle_x - x_defect(i)); % Angle in radians, relative to the defect
    
    % Compute theta0, defect orientation, using complex phase analysis
    m = charge(i);
    complex_terms = exp(1i * 2 * (m * phi - orientations_in_circle));
    avg_complex = mean(complex_terms);
    theta0 = angle(avg_complex) / 2;
    
    % Store orientation
    defect_orientation(i) = theta0 / (charge(i) - 1);
end

% Save outputs
vars.defectCharge = [vars.defectCharge; charge];
vars.defectX = [vars.defectX; x_defect];
vars.defectY = [vars.defectY; y_defect];
vars.defectOrientation = [vars.defectOrientation; defect_orientation];

hold on;


% Use one-sided red lines to show orientation of plus defects
% Plot orientation for plus defects: red dot + outward line
for i = 1:length(xP)
    % Dot at the center
    plot(xP(i), yP(i), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3);

    % Outward line
    x_end = xP(i) + markerSize * cos(defect_orientation(i));
    y_end = yP(i) + markerSize * sin(defect_orientation(i));
    plot([xP(i), x_end], [yP(i), y_end], 'r-', 'LineWidth', 2);
end

% Use custom 3-segment marker for minus defects
% Plot orientation for minus defects: blue dot + 3-spoke lines
angles_deg = [0, 120, 240];
angles_rad = deg2rad(angles_deg);
base_marker = [cos(angles_rad); sin(angles_rad)] * markerSize;

for i = 1:length(xM)
    idx = i + length(xP); % index in orientation array

    % Dot at the center
    plot(xM(i), yM(i), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 3);

    % Rotate and plot 3 spokes
    R = [cos(defect_orientation(idx)), -sin(defect_orientation(idx));
         sin(defect_orientation(idx)),  cos(defect_orientation(idx))];
    rotated_segments = R * base_marker;

    for k = 1:3
        x1 = xM(i);
        y1 = yM(i);
        x2 = x1 + rotated_segments(1, k);
        y2 = y1 + rotated_segments(2, k);
        line([x1, x2], [y1, y2], 'Color', 'b', 'LineWidth', 2);
    end
end

%axis equal;
hold on;

end