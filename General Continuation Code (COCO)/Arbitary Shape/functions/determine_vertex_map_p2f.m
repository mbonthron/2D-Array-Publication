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

    if ~isempty(matched)
        % The periodic hinge is exactly the finite hinge
    else
        % Check if off by length of super cell
        x_pos = points_finite(i,1) - L_super_cell;
        matched = find(ismember(points,[x_pos y_pos],'rows') == 1);
        
        % If we need to remap, add this to the matrix vertex_map_p2f
        % This states that in the finite the ith hinge -> matched hinge on
        % periodic
        vertex_map_p2f = [vertex_map_p2f ; 
                          matched i];
    end
end

data.vertex_map_p2f = vertex_map_p2f;

end

