function [defect_points_m, defect_points_p,defect_points_2m, defect_points_2p,loops] = get_all_defects(directors,S)


% Set hperiodic = 1 or vperiodic = 1 to impose periodic boundary conditions on the grid in
% the horizontal or vertical direction.

hperiodic = 1;
vperiodic = 1;
gridSpacing = 80;
anglethreshold = pi/4;

% Given the director field, this code computes all of the defects
%
% This code is self contained and requires no additional files.
%
% Determines patch_size directly from directors.

% Version date 10/13/2020

% Things that still need doing
% hperiodic = 0 is not properly implemented.
% need to update use of patch_size.  Currently, only works for patch_size =
% 0.
% Code needs to be verified on Amanda's data.

% Computes whether there is a band of zero-length directors framing the
% image.  This happens when the directors cannot be determined too close to
% the boundary.
frame_size = find(directors(:,round(end/2),1) ~= 0,1)-1;

% Delete frame
directors = directors(1+frame_size:end-frame_size,1+frame_size:end-frame_size,:);


[m,n] = size(directors(:,:,1));
% First entry of directors is y coordinate.

% % Normalize directors  
% In our(Yimin'b) code, input is unit vectoralready
% dnorms = sqrt(directors(:,:,1).^2 + directors(:,:,2).^2);
% directors(:,:,1) = directors(:,:,1) ./ dnorms;
% directors(:,:,2) = directors(:,:,2) ./ dnorms;

% Find the angles between adjacent grid points in horizontal and vertical
% directions.  Uses periodic boundary conditions.
hangle = asin( circshift(directors(:,:,1),-1,2).*directors(:,:,2) ...
    - circshift(directors(:,:,2),-1,2).*directors(:,:,1) );
vangle = asin(circshift(directors(:,:,1),-1,1).*directors(:,:,2) ...
    - circshift(directors(:,:,2),-1,1).*directors(:,:,1));

% 
% % Construct connection and cell graphs between points and cells.
% cellG is a graph connecting cells with a large-angle-difference edge
% between them.
% pointG is a graph connecting points with a small angle difference between
% them.

% Find adjacent horizontal points with small angle difference
indices = abs(hangle)<anglethreshold;
hpconnect = find(indices(:));

% Find points to the right.
hpconnectp1 = hpconnect + m;
% Impose periodic boundary conditions.
indices2 = find(hpconnectp1 > m*n);
hpconnectp1(indices2) = hpconnectp1(indices2) - m*n;

% Find adjacent vertical cells with large angle difference
vcconnect = find(~indices(:));

% Find cells above.
vcconnectm1 = vcconnect - 1;
% Impose periodic boundary conditions.
indices2 = find(mod(vcconnectm1,m) == 0);
vcconnectm1(indices2) = vcconnectm1(indices2) + m;

% Find adjacent vertical points with small angle difference
indices = abs(vangle)<anglethreshold;
vpconnect = find(indices(:));

% Find points below.
vpconnectp1 = vpconnect + 1;
% Impose periodic boundary conditions.
indices2 = find(mod(vpconnectp1,m) == 1);
vpconnectp1(indices2) = vpconnectp1(indices2) - m;

% Find adjacent horizontal cells with large angle difference
hcconnect = find(~indices(:));

hcconnectm1 = hcconnect - m;
indices2 = find(hcconnectm1 < 1);
hcconnectm1(indices2) = hcconnectm1(indices2) + m*n;


pointG = graph([vpconnect;hpconnect],[vpconnectp1;hpconnectp1],1,m*n);
cellG = graph([vcconnect;hcconnect],[vcconnectm1;hcconnectm1],1,m*n);


% %         
% Impose boundary conditions; default construction above uses periodic
% boundary conditions

% if ~vperiodic
%     % Eliminate periodic boundary conditions.
%     pointconnectv(1:m:end,m:m:end) = 0;
%     pointconnectv(m:m:end,1:m:end) = 0;
%     cellconnectv(:,m:m:end) = 0;
%     cellconnectv(m:m:end,:) = 0;
%     cellconnecth(:,m:m:end) = 0;
%     cellconnecth(m:m:end,:) = 0;
% end
% 
% if ~hperiodic
%     % Eliminate periodic boundary conditions.
%     
%     pointconnecth(1:m,end-m+1:end) = 0;
%     pointconnecth(end-m+1:end,1:m) = 0;
%     cellconnectv(:,end-m+1:end) = 0;
%     cellconnectv(end-m+1:end,:) = 0;
%     cellconnecth(:,end-m+1:end) = 0;
%     cellconnecth(end-m+1:end,:) = 0;
% end

Xcoords = (1:(n))+1/2;
Ycoords = (1:(m))+1/2;
[Xmatrix,Ymatrix] = meshgrid(Xcoords,Ycoords);
% plot(G,'Xdata',Xmatrix(:),'Ydata',Ymatrix(:),'Linewidth',10)
% hold on

% 
%Gpoints = graph(pointconnect);
% 
%  Xcoords2 = 1:n;
%  Ycoords2 = 1:m;
%  [Xmatrices2,Ymatrices2] = meshgrid(Xcoords2,Ycoords2);
%  plot(Gpoints,'Xdata',Xmatrices2(:),'Ydata',Ymatrices2(:),'Linewidth',2)

% Find connected components of cells with more than 1 cell 
[bins,binsizes] = conncomp(cellG);
bigbins = find(binsizes ~= 1);
loops = {};
% Each loop is a list of points in an nx2 array.  The x-coordinate is in
% the first column.

defect_points = zeros(0,2);
for bin = bigbins
    blobcells = find(bins == bin);
    % Convert blobcells indices to ij pairs of indices
    
    %plot(Xmatrix(blobcells),Ymatrix(blobcells),'.','markersize',10);
    inds1 = mod(blobcells-1,m)+1;
    inds2 = floor((blobcells-1)/m) + 1;
    
    loop = walkaroundblob(inds1',inds2',m,n,pointG);
    
    %if trace(pointconnect(loop(1:end-1),loop(2:end))) == length(loop)-1
    if ~isempty(loop)
        % Only keep loop if it is confirmed that all connections between
        % vertices have small angles.  This is only a possibility for loops that touch the boundary.
                       
        % Translate to i j indexing
        inds2 = floor((loop-1)/m) + 1;
        inds1 = mod(loop-1,m)+1;
        loop = [inds2',inds1'];
                
        loops{end+1} = loop;
        
        % Find representative point for defect
        if exist('S')
            [~,index] = min(max(S(blobcells)));
            index = blobcells(index);
            defect_points(end+1,:) = [Xmatrix(index),Ymatrix(index)];            
        else
            D = distances(cellG,blobcells,blobcells);
            [~,index] = min(max(D));
            index = blobcells(index);
            defect_points(end+1,:) = [Xmatrix(index),Ymatrix(index)];
        end
        
    end
end

% Add back in frame_size

defect_points(:,1) = defect_points(:,1) + frame_size;
defect_points(:,2) = defect_points(:,2) + frame_size;

% Post processing

% Determine charge of points.

defect_points_p = zeros(0,2);
defect_points_m = zeros(0,2);
defect_points_2p = zeros(0,2);
defect_points_2m = zeros(0,2);

charges = [];
for ind = 1:length(loops)
    loop = loops{ind};
    charges(ind) = compute_charge(loop,directors);
end

% Keep only +1/2 and -1/2 charges

% plot_loops(loops,directors)

pindices = find(charges == 0.5);
mindices = find(charges == -0.5);

defect_points_p =defect_points(pindices,:);
defect_points_m =defect_points(mindices,:);

p2indices = find(charges == 1);
m2indices = find(charges == -1);

defect_points_2p =defect_points(p2indices,:);
defect_points_2m =defect_points(m2indices,:);

% Remove defects that are too close to opposite charges (pairwise removal)
% Initialize logical arrays to track removed defects
removedP = false(size(defect_points_p, 1), 1);
removedM = false(size(defect_points_m, 1), 1);

numP = size(defect_points_p, 1);
numM = size(defect_points_m, 1);

for i = 1:numP
    if removedP(i)  % Skip if already marked for removal
        continue;
    end
    for j = 1:numM
        if removedM(j)  % Skip if already marked for removal
            continue;
        end
        
        dist = sqrt((defect_points_p(i,1) - defect_points_m(j,1))^2 + ...
                    (defect_points_p(i,2) - defect_points_m(j,2))^2); % Calculate distance
        
        if dist < gridSpacing  % If defects are within a specified grid spacing
            removedP(i) = true;  % Mark positive defect for removal
            removedM(j) = true;  % Mark negative defect for removal
            break;  % Ensure only one pair is removed
        end
    end
end

% Remove marked defects
defect_points_p(removedP, :) = [];
defect_points_m(removedM, :) = [];

return


function loop = walkaroundblob(inds1,inds2,m,n,pointG)

% inds1 and inds2 are lists of cells.

% Loop here is a list of points in 1xn vector.

% Next step is to fill in the interior of cells




% Find extreme upper right point in cells 

ind1 = inds1(1);
ind2 = inds2(1);

done = 0;
while ~done
    
    ind11 = mod(ind1,m)+1;
    ind21 = mod(ind2,n)+1;
    
    if ismember([ind1,ind21],[inds1,inds2],'rows')
        ind2 = ind21;
    elseif ismember([ind11,ind2],[inds1,inds2],'rows')
        ind1 = ind11;
    else
        done = 1;
    end
end

cellind1 = ind1;
cellind2 = ind2;


% Find topmost rightmost point in loop,

ind1 = mod(cellind1,m)+1;
ind2 = mod(cellind2,n)+1;
ind = ind1 + (ind2-1)*m; 

ind1 = mod(cellind1,m)+1;
ind2 = cellind2;
prevind = ind1 + (ind2-1)*m; 


% Find all points that are vertices of cells.

% Bottom left point
pointinds1 = inds1;
pointinds2 = inds2;

% Bottom right point
pointinds1 = [pointinds1;inds1];
pointinds2 = [pointinds2;mod(inds2,n)+1];

% Top left point
pointinds1 = [pointinds1;mod(inds1,m)+1];
pointinds2 = [pointinds2;inds2];

% Top right point
pointinds1 = [pointinds1;mod(inds1,m)+1];
pointinds2 = [pointinds2;mod(inds2,n)+1];

pointinds = [pointinds1,pointinds2];
pointinds = unique(pointinds,'rows');

points = pointinds(:,1) + (pointinds(:,2)-1)*m;

Gloop = subgraph(pointG,points);

ind1 = find(points == ind);
prevind1 = find(points == prevind);

[bins,binsizes] = conncomp(Gloop);
loopbin = find(binsizes ~= 1);
if length(loopbin) ~= 1
    loop = [];
    return
end
looppoints = points(find(bins == loopbin));

Gloop = rmedge(Gloop,prevind1,ind1);
Gloop = rmedge(Gloop,ind1,prevind1);

loop = points(shortestpath(Gloop,ind1,prevind1))';
if ~isempty(loop)
    loop(end+1) = loop(1);
end

% Reverse loop so that it is CCW.

loop = loop(end:-1:1);


function charge = compute_charge(points,directors)

if length(points) == 0
    charge = NaN;
    return
end

new_director = squeeze(directors(points(1,2),points(1,1),:));
new_director=new_director/norm(new_director);

theta = 0;
for ind = 2:size(points,1)
    old_director = new_director;
    new_director=squeeze(directors(points(ind,2),points(ind,1),:));
    new_director=new_director/norm(new_director);
    if new_director'*old_director <0
        new_director = -new_director;
    end
    theta = theta + asin(old_director(1)*new_director(2)-old_director(2)*new_director(1));
end

charge = theta/pi/2;
charge = round(2*charge)/2;


