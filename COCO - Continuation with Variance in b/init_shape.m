function [data] = init_shape(shapeNum, data)
if shapeNum == 1
    % Rhombus
    run('rhombus.m')
elseif shapeNum == 2
    % Double Triangle Chain
    run('double_triangle.m')

elseif shapeNum == 3
    % Alternating Triangle Chain
    run('alternating_traingle.m' )
elseif shapeNum == 4
    % Hexagon
    run('hexagon.m')

elseif shapeNum == 5
    % Diamond Chain
    run('diamond_chain.m')

elseif shapeNum == 6
    % Diagonal Rhombus
    run('diagonal_rhombus.m')

elseif shapeNum == 7
    % Diagonal Rhombus with Arch
    run('diagonal_rhombus_with_arch.m')

else
    error("SHAPE DOES NOT EXIST (please fix)")
end

%% Prepare 'data' for COCO
% Remove any unneeded vertecies
data = determine_adjacency_matrix(data);

data = remove_node(data, nodes_to_remove);
[~] = plot_grid(data, 1); %debugging, ask Michael for his nice plotting code

% Determine the adjacency_matrix (assuming ALL connections)
data = determine_adjacency_matrix(data);

% Determine what is needed to make structure periodic
% Differentiates between points, points_finite, and points_time_integration
data = add_periodicity(data);
[~] = plot_grid(data, 1); %debugging, ask Michael for his nice plotting code


% Create the adjacency matrix for the time integration
data2.points = data.points_time_integration;
data2 = determine_adjacency_matrix(data2);

% Save the time integration adjacency matrix into data
data.adjacency_matrix_time_integration = data2.adjacency_matrix;

data2.points = data.points_time_integration;
data2.adjacency_matrix = data.adjacency_matrix_time_integration;
[~] = plot_grid(data2, 1);
data2 = remove_node(data2,nodes_to_remove2);
data.points_time_integration = data2.points;
data.N_time_integration = data2.N;
data.adjacency_matrix_time_integration = data2.adjacency_matrix;

% If needed to help remove connections
[~] = plot_grid(data, 1); %debugging, ask Michael for his nice plotting code


data = remove_connection(data,connections_to_remove);

[~] = plot_grid(data, 1); %debugging, ask Michael for his nice plotting code

data = determine_per_to_finite(data);
data = determine_time_to_periodic(data);

[~] = plot_grid(data, 1); %debugging, ask Michael for his nice plotting code
data.b_vector = zeros(data.N,1);

% Nodes to hold stationary
data.nodes_to_hold = nodes_to_hold;
data.arches_to_displace = arches_to_displace;
data.nodes_to_rotate = nodes_to_rotate;
data.arches_to_force_positive = arches_to_force_positive;
data.arches_to_force_negative = arches_to_force_negative;

%% Start with elastic deformation
[data] = initialize_elastic_deformation(zeros(data.N,1),zeros(data.V,1),data);

%Consider what's actually necessary since this is going into COCO
data.e_vector = 0*ones(data.N,1);

% Determine the coefficient matrix and number of constraints of the system
data = determine_coefficient_matrix(data);
data = determine_modes_to_skip(data);

COCO_plot_system_once(zeros(data.N*data.N_modes),data);

end