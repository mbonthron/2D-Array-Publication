function [data] = remove_connection(data,connections_to_remove)
%ADD_PERIODICITY Adds the required arches for the periodic unit cell
%   This only applies for the unit cell!!!!

%% Remove the values for the PERIODIC adjacency_matrix
adjacency_matrix = data.adjacency_matrix;

for i = 1:length(connections_to_remove)
    node1 = connections_to_remove(i,1);
    node2 = connections_to_remove(i,2);
    
    adjacency_matrix(node1,node2) = 0;
    adjacency_matrix(node2,node1) = 0;

end

data.adjacency_matrix = adjacency_matrix;
data.N = sum(triu(adjacency_matrix,1) ==1,'all');   

%% Remove the values for the FINITE adjacency_matrix
adjacency_matrix_finite = data.adjacency_matrix_finite;
vertex_map_p2f = data.vertex_map_p2f;

for i = 1:length(connections_to_remove)
    node1 = connections_to_remove(i,1);
    node2 = connections_to_remove(i,2);
    
    adjacency_matrix_finite(node1,node2) = 0;
    adjacency_matrix_finite(node2,node1) = 0;

    % Check if either node is mapped to a different the finite structure
    if any(ismember([node1 node2],vertex_map_p2f(:,1)))
        % Find which node needs mapping
        mode_needing_mapping =  ismember([node1 node2],vertex_map_p2f(:,1));
        
        old_nodes = [node1 node2];
        new_nodes = old_nodes;

        % Find which of the nodes is mapped to a new one
        idx_to_change = find(mode_needing_mapping == 1);

        for k = idx_to_change
            row_from_map_p2f = find(vertex_map_p2f(:,1) == old_nodes(k));
            new_nodes(k) = vertex_map_p2f(row_from_map_p2f,2);
        end

        % If they are mapped, switch out the numbering
        node1 = new_nodes(1);
        node2 = new_nodes(2);

        adjacency_matrix_finite(node1,node2) = 0;
        adjacency_matrix_finite(node2,node1) = 0;

    end
end

% Remove any "floating" points (i.e. unconnected points)
floating_points = find(sum(adjacency_matrix_finite) == 0);

% Remove floating points from points
data.points_finite(floating_points,:) = [];

% Remove floating rows and columns from adjacency matrix
adjacency_matrix_finite(floating_points,:) = [];
adjacency_matrix_finite(:,floating_points) = [];

data.adjacency_matrix_finite = adjacency_matrix_finite;

%% = Remove the Connections for the Time Integration Adj. Matrix
N_cells = data.N_cells;
V = data.V;
adjacency_matrix_time_integration = data.adjacency_matrix_time_integration;

V_time_integration = size(data.points_time_integration,1);

% Loop over each super cell and remove connections
for i = 1:N_cells
    connections_to_remove_cell = connections_to_remove + (i-1)*V;
    vertex_map_p2f_cell = vertex_map_p2f + (i-1)*V;
    
    % Any points in connection_to_remove_cell exceed points,
    any_row_matches = any(vertex_map_p2f_cell > V_time_integration,2);
    vertex_map_p2f_cell(any_row_matches,:) = [];

    for k = 1:length(connections_to_remove_cell)
        node1 = connections_to_remove_cell(k,1);
        node2 = connections_to_remove_cell(k,2);
        
        adjacency_matrix_time_integration(node1,node2) = 0;
        adjacency_matrix_time_integration(node2,node1) = 0;

        % Check if either node is mapped to a different the finite structure
        if any(ismember([node1 node2],vertex_map_p2f_cell(:,1)))
            % Find which node needs mapping
            mode_needing_mapping =  ismember([node1 node2],vertex_map_p2f_cell(:,1));

            old_nodes = [node1 node2];
            new_nodes = old_nodes;

            % Find which of the nodes is mapped to a new one
            idx_to_change = find(mode_needing_mapping == 1);

            for j = idx_to_change
                row_from_map_p2f = find(vertex_map_p2f_cell(:,1) == old_nodes(j));
                new_nodes(j) = vertex_map_p2f_cell(row_from_map_p2f,2);
            end

            % If they are mapped, switch out the numbering
            node1 = new_nodes(1);
            node2 = new_nodes(2);

            adjacency_matrix_time_integration(node1,node2) = 0;
            adjacency_matrix_time_integration(node2,node1) = 0;

        end

    end
end

data.adjacency_matrix_time_integration = adjacency_matrix_time_integration;
data.N_time_integration = sum(triu(adjacency_matrix_time_integration,1) ==1,'all');

end

