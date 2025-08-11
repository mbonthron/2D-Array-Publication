function [data] = determine_vertex_map_p2f(data)
%ADD_PERIODICITY Adds the required arches for the periodic unit cell
%% Load Data
points        = data.points;
points_finite = data.points_finite;

L_super_cell  = data.L_super_cell;

vertex_map_p2f = [];
for i = 1:length(points_finite)
    % Check if the point exactly appears in points_finite
    x_pos = points_finite(i,1);
    y_pos = points_finite(i,2);
    
    matched = find(ismember(points,[x_pos y_pos],'rows') == 1);

    %% Try to rematch the ground nodes
    ground_nodes_matched = find(ismember(data.ground_nodes_points, points_finite(i,:),'rows')==1, 1);
    
    % Renumber the ground nodes
    if ~isempty(ground_nodes_matched) && ~isempty(matched)
        data.ground_nodes_points(ground_nodes_matched,:) = data.points(matched,:);
        data.ground_nodes_idx(ground_nodes_matched,:) = matched;
        
    end

    if matched == i
        % The periodic hinge is exactly the finite hinge
        % Do not need to include any mapping
    elseif ~isempty(matched)
        % The periodic hinge is renumber (this can happen if our numbering
        % convention got changed)
        vertex_map_p2f = [vertex_map_p2f ; 
                          matched i];

    else
        % Check if off by length of super cell
        x_pos = points_finite(i,1) - L_super_cell;
        matched = find(ismember(points,[x_pos y_pos],'rows') == 1);
        
        % If we need to remap, add this to the matrix vertex_map_p2f
        % This states that in the finite the ith hinge -> matched hinge on
        % periodic
        vertex_map_p2f = [vertex_map_p2f ; 
                          matched i];

        data.ground_nodes_points(ground_nodes_matched,:) = [];
        data.ground_nodes_idx(ground_nodes_matched,:) = [];
    end

  


end

data.vertex_map_p2f = vertex_map_p2f;

end

