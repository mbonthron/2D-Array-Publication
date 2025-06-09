function [data] = add_periodicity(data)
%ADD_PERIODICITY Adds the required arches for the periodic unit cell
%% Load Data
points = data.points;

%% Add Periodicity
% Determine which points are "ground nodes" that begin the supercell
[root_node_x, root_node_idx] = min(points(:,1));
ground_nodes_idx = find(points(:,1) - root_node_x < 1);
ground_nodes_points = points(ground_nodes_idx,:);

% Add another column of ground nodes on the right side temporarily
% Determine minimum distance for L = 1 on same y level
% ASSUMES horizontal arch in periodic 
temp_points = zeros(size(ground_nodes_points,1),2);
min_dist = -9999;
for GN_point_idx = 1:size(ground_nodes_points,1)
    x = ground_nodes_points(GN_point_idx,1);
    y = ground_nodes_points(GN_point_idx,2);
    
    [min_dist_curr, min_dist_idx] = max(points(points(:,2) == y,1)-points(GN_point_idx,1));
    min_dist = max([min_dist, min_dist_curr-x]);
end

% Check if the x points are spaced by '1' or by sqrt(3)/2
xpoints = points(:,1) / (sqrt(3)/2);
if all(round(xpoints,0) == round(xpoints,2)) 
    extra_hor_offset = sqrt(3);
else
    extra_hor_offset = 1;
end

temp_points = ground_nodes_points + [ones(size(ground_nodes_points,1),1)*(min_dist+extra_hor_offset), zeros(size(ground_nodes_points,1),1)];
temp_adjacency_matrix = data.adjacency_matrix;
temp_adjacency_matrix2 = [data.adjacency_matrix zeros(size(data.adjacency_matrix,1), size(ground_nodes_points,1)); zeros(size(ground_nodes_points,1), size(data.adjacency_matrix,1)+size(ground_nodes_points,1))];

% Find which points are 1 unit away from these points
% This can be vectorized using matrix math, do later
x = points(:,1);
y = points(:,2);
for i = 1:size(temp_points,1)
    % i is the index of the ground node and the new "temp" point
    % Determine the distance from the ith point to all other points
    xi = temp_points(i,1);
    yi = temp_points(i,2);
    distances = sqrt(sum(([xi yi] - [x y]).^2,2));
    
    % Find all the points that are distance 1 from the ith point
    neighbors = find(round(distances,1)==1);
    
    i_sized_neighbors = i*ones(size(neighbors));

    temp_adjacency_matrix(i_sized_neighbors,neighbors) = 1;
    temp_adjacency_matrix(neighbors,i_sized_neighbors) = 1;

    temp_adjacency_matrix2(i_sized_neighbors+size(points,1),neighbors) = 1;
    temp_adjacency_matrix2(neighbors,i_sized_neighbors+size(points,1)) = 1;
end
% Add connections between the rightmost points to the leftmost points

%% Store into data
data.adjacency_matrix        = temp_adjacency_matrix;
data.adjacency_matrix_finite = temp_adjacency_matrix2;

data.N = sum(triu(temp_adjacency_matrix,1) ==1,'all');

data.points_finite = [points;temp_points];
data.b_vector = zeros(data.N,1);
[data] = initialize_elastic_deformation(zeros(data.N,1),zeros(data.V,1),data);

% Map the periodic vertecies to the finite structure
data.vertex_map_p2f = [ground_nodes_idx size(points,1)*ones(size(ground_nodes_idx,1),1)+ground_nodes_idx];

% Determine the expanded points for time integration
points_time_integration = points;
for i = 1:(data.N_cells-1)
    add_in = points + i*[ones(size(points,1),1)*(min_dist+extra_hor_offset), zeros(size(points,1),1)];
    points_time_integration = [points_time_integration; add_in];
end

% Add in the last ground points
% Check if the ground points all have the same x value
% if all(ground_nodes_points(:,1) == ground_nodes_points(1,1))
%     % Ground nodes have the same x values -> add in all
%     add_in = ground_nodes_points + data.N_cells*[ones(size(ground_nodes_points,1),1)*(min_dist+extra_hor_offset), zeros(size(ground_nodes_points,1),1)];
% 
% else
%     % Ground nodes have different x values -> add in smallest x
%     [val, idx] = min(ground_nodes_points(:,1));
%     ground_nodes_to_add = ground_nodes_points(ground_nodes_points(:,1) == val,:);
%     add_in = ground_nodes_to_add + data.N_cells*[ones(size(ground_nodes_to_add,1),1)*(min_dist+extra_hor_offset), zeros(size(ground_nodes_to_add,1),1)];
% end
add_in = ground_nodes_points + data.N_cells*[ones(size(ground_nodes_points,1),1)*(min_dist+extra_hor_offset), zeros(size(ground_nodes_points,1),1)];

points_time_integration = [points_time_integration ; add_in];

data.points_time_integration = points_time_integration;
data.L_super_cell = min_dist+extra_hor_offset;

end

