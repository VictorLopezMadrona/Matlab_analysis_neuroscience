function B = custom_bwboundaries(mask_max)
% This function finds the boundaries of objects in a binary image.
% Input:
%   - mask_max: A binary image (1 for objects, 0 for background)
% Output:
%   - B: A cell array where each cell contains the coordinates (row, col)
%         of the boundary pixels for each connected object in the image

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Aug. 2025; Last revision: 07-Apr-2025   

% Add a border to the image
mask_aux = zeros(size(mask_max,1)+2,size(mask_max,2)+2);
mask_aux(2:end-1,2:end-1) = mask_max;
mask_max = mask_aux;    

% Initialize the cell array to store boundary coordinates
B = {};

% Initialize visited matrix to keep track of processed pixels
visited = false(size(mask_max));

% Directions for 8-connected neighborhood (clockwise starting from the top-left)
directions = [
    0,  1;  % East
    1,  1;  % SouthEast
    1,  0;  % South
    1, -1;  % SouthWest
    0, -1;  % West
    -1, -1;  % NorthWest
    -1,  0;  % North
    -1,  1   % NorthEast
    ];

% Loop through each pixel in the image (excluding borders)
for r = 2:size(mask_max, 1)-1
    for c = 2:size(mask_max, 2)-1
        % If it's a foreground pixel and hasn't been visited yet
        if mask_max(r, c) == 1 && ~visited(r, c)
            % Start tracing the boundary from this pixel
            [boundary, visited] = trace_boundary(mask_max, r, c, visited, directions);
            % Add the found boundary to the cell array
            B{end+1} = boundary;
            
            % Mark the entire region as visited (flood-fill)
            visited = flood_fill(visited, B{end});  % Mark the entire object as visited
        else
            visited(r,c) = 1;
        end
    end
end

%Remove 1 from each value, which corresponds to the border
for bi=1:length(B)
    B{bi}=B{bi}-1;
end
end

function [boundary, visited] = trace_boundary(mask_max, r_start, c_start, visited, directions)
    % This function traces the boundary of an object starting from a given pixel.
    % Input:
    %   - mask_max: The binary image
    %   - r_start, c_start: Starting row and column for boundary tracing
    %   - visited: Matrix of visited pixels
    %   - directions: The 8-connected direction array
    % Output:
    %   - boundary: A list of boundary pixels' coordinates [row, col]
    
    boundary = [];  % Initialize boundary list
    r = r_start;
    c = c_start;
    
    boundary = [boundary; r, c];  % Add the starting point to the boundary
    first_pixel = [r, c];  % Save the starting pixel
    
    % Mark the starting pixel as visited
    %visited(r, c) = true;
    
    % Starting direction: top-left (North-West)
    dir = 0;  % Start direction: North (index 1 in the directions array)
    
    flag = true;
    while flag
        % Try the next direction in a clockwise manner (8-connected)
        for i = 1:8
            % Move to the next direction
            dir = dir+1;
            dir = mod(dir-1, 8) + 1;  % Ensure direction cycles through 1 to 8
            dr = directions(dir, 1);  % Row change
            dc = directions(dir, 2);  % Column change
            
            nr = r + dr;
            nc = c + dc;
            
            % Check if the new position is within bounds and is part of the object
            if nr > 0 && nr <= size(mask_max, 1) && nc > 0 && nc <= size(mask_max, 2)
                if mask_max(nr, nc) == 1 && ~visited(nr, nc)  % It's part of the object and not visited
                    % Add it to the boundary and mark it as visited
                    boundary = [boundary; nr, nc];
                    visited(nr, nc) = true;  % Mark this pixel as visited
                    r = nr;
                    c = nc;
                    % Go to opposite direction and restart from there
                    dir=dir+4;
                    break;
                end
            else
                visited(nr, nc) = 1;
            end
            % There is no ending point, like a line
            if i==8
                flag = false;
                boundary = [boundary; first_pixel(1), first_pixel(2)];
            end
        end
        
        % If we are back at the starting pixel, we have closed the boundary
        if r == first_pixel(1) && c == first_pixel(2)
            boundary = [boundary; r, c];
            break;
        end
    end
end

function visited = flood_fill(visited, B)
% This function fills the entire object starting from (r, c), marking it as visited.
% It turns all pixels of the object into visited (true) pixels.

% Check all  points within a square around the are to check if they are
% inside the boundaries or not
for r = min(B(:,1)):max(B(:,1))
    for c = min(B(:,2)):max(B(:,2))
        pix_aux = find(B(:,1)==r);
        if sum(B(pix_aux,2)>c)>0 && sum(B(pix_aux,2)<c)>0
            visited(r,c)=1;
        end
    end
end
end




