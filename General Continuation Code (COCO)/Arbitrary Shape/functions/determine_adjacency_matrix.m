function [adjacency_matrix] = determine_adjacency_matrix(points)
%DETERMINE_ADJACENCY Determines the adjacency matrix from a matrix of points
%   INPUTS
%   ===================================================
%   points(:,1) = x positions
%   points(:,2) = y positions
%   Any points that are within 1 unit of one another will be connected
%
%   OUTPUTS
%   ===================================================
%   adjacency_matrix = adjacency matrix where the (i,j) element is 1 if the
%   i and j elements are connect and 0 if the elements are not connected
%%
x = points(:,1);
y = points(:,2);

adjacency_matrix = zeros(length(x),length(x)) ;

for i = 1:height(points)
    % [x,y] position of the ith point
    xi = x(i);
    yi = y(i);

    % Determine the distance from the ith point to all other points
    distances = sqrt(sum(([xi yi] - [x y]).^2,2));
    
    % Find all the points that are distance 1 from the ith point
    neighbors = find(round(distances,1)==1);
    
    i_sized_neighbors = i*ones(size(neighbors));

    adjacency_matrix(i_sized_neighbors,neighbors) = 1;
    adjacency_matrix(neighbors,i_sized_neighbors) = 1;
end

end