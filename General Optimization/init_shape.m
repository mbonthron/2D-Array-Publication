function [data] = init_shape(shapeNum, data)
run('points_chain_direct')

data = determine_adjacency_matrix(data);
% plot_grid(data, 1) %debugging, ask Michael for his nice plotting code

if shapeNum == 1
    %Rhombus
    %% === Load a points data for the system
    %run('points_rhombus_direct')
    nodes_to_remove = [3,6,9,12];
    connections_to_remove = []; %THIS IS AFTER NODES ARE REMOVED
    data.shape_name = 'Rhombus';

elseif shapeNum == 2
    % Double Triangle Chain

    nodes_to_remove = [];

    connections_to_remove = [2 5
                            5 8
                            11 8
                            11 2
                            3 6
                            6 9
                            9 12
                            12 3];
    data.shape_name = 'Double Triangle';

elseif shapeNum == 3
    % Alternating Triangle Chain

    nodes_to_remove = [3, 5, 9, 11];

    connections_to_remove = [];
    data.shape_name = 'Alternating Triangle';

elseif shapeNum == 4

elseif shapeNum == 5

elseif shapeNum == 6


end

data = remove_node(data, nodes_to_remove);
data = determine_adjacency_matrix(data);

data = add_periodicity(data);

% If needed to help remove connections
% plot_grid(data, 1) %debugging, ask Michael for his nice plotting code

data = remove_connection(data,connections_to_remove);

data = determine_per_to_finite(data);

plot_grid(data, 1) %debugging, ask Michael for his nice plotting code
COCO_plot_system_once(zeros(data.N*data.N_modes),data)
data.b_vector = zeros(data.N,1);


%% Start with elastic deformation
[data] = initialize_elastic_deformation(zeros(data.N,1),zeros(data.V,1),data);

%Consider what's actually necessary since this is going into COCO
data.e_vector = 0*ones(data.N,1);
data.t_vector = 0.01*pi*ones(data.N,1);

% Determine the coefficient matrix and number of constraints of the system
data = determine_coefficient_matrix(data);
data = determine_modes_to_skip(data);
end