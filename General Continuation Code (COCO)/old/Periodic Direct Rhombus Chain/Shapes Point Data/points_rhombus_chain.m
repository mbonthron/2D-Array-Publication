yy = sqrt(3)/2;
number_of_elements = 1;

points0 = [0 0;
         0.5 yy;
         1 0;
         1.5 yy;
         2 0;
         2.5 yy;
          3 0;
          3.5 yy
          4 0
          4.5 yy];

points = points0;

for i = 2:number_of_elements
    points_to_add = points0(3:end,:) + [4*(i-1)*ones(8,1) zeros(8,1)];
    points = [points ; points_to_add];
end


% Determine the adjacency matrix
adjacency_matrix = determine_adjacency_matrix(points);

% Visualize the point and connection between the nodes
plot_grid(adjacency_matrix,points,1);

N = sum(triu(adjacency_matrix,1) ==1,'all');
unit_cell_count = (N - 1) / 16;
