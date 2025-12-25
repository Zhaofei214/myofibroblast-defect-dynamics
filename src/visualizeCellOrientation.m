 function [directors, S] = visualizeCellOrientation(inputImage, params)
    % % Ensure figure exists or create new one
    % figure('Name', 'Cell Orientation Visualization', ...
    %        'Position', [400 100 800 600]);


    % Process image with Gaussian smoothing
    [Gx, Gy] = imgradientxy(inputImage,"sobel");
    
    % Structure tensor gradients
    Axx = Gx.^2;
    Ayy = Gy.^2;
    Axy = Gx .* Gy;

    % Smooth the tensor component
    Axx = imgaussfilt(Axx,params.smoothingSigma);
    Ayy = imgaussfilt(Ayy,params.smoothingSigma);
    Axy = imgaussfilt(Axy,params.smoothingSigma);
    
    % Compute the eigenvalues of the structure tensor
    lambda1 = (Axx + Ayy)/2 + sqrt((Axx - Ayy).^2/4 + Axy.^2);
    lambda2 = (Axx + Ayy)/2 - sqrt((Axx - Ayy).^2/4 + Axy.^2);
    
    % Calculate orientation field
    orientation = atan2(2 * Axy, Axx - Ayy) / 2 + pi/2;
    %confidence = (lambda1 - lambda2) ./ (lambda1 + lambda2 + eps);
    
    % Adjust orientation to director to fit Prof Kevin Mitchell's code(from UC Merced)
    [row, col] = size(orientation);
    directors = zeros(row, col, 2);
    directors(:,:,1) = cos(orientation);
    directors(:,:,2) = sin(orientation);
 
    % Compute the order parameter S = sqrt(<cos 2theta>^2 + <sin 2theta>^2)
    S = sqrt(mean(cos(2* orientation(:))).^2 + mean(sin(2*orientation(:))).^2);

    % Generate sampling grid
    [rows, cols] = size(inputImage);
    [X, Y] = meshgrid(1:params.gridSpacing:cols, 1:params.gridSpacing:rows);
    % Interpolate orientation and confidence fields
    orientationInterp = interp2(orientation, X, Y, 'linear');
    % confidenceInterp = interp2(confidence, X, Y, 'linear');

    % Display the image
    imshow(inputImage, 'Colormap', gray, 'DisplayRange', [0 0.5]);
    axis ij;
    hold on;
    
    % Draw orientation lines
    for i = 1:size(X, 1)
        for j = 1:size(X, 2)
            %if confidenceInterp(i, j) > params.confidenceThreshold add end
                angle = orientationInterp(i, j);
                dx = params.lineLength * cos(angle);
                dy = params.lineLength * sin(angle);

                x1 = X(i, j) - dx / 2;
                y1 = Y(i, j) - dy / 2;
                x2 = X(i, j) + dx / 2;
                y2 = Y(i, j) + dy / 2;

                line([x1 x2], [y1 y2], 'Color', params.lineColor, ...
                     'LineWidth', params.lineWidth);
        end
    end
end


