function [data] = remove_node(data,nodes_to_remove)
%ADD_PERIODICITY Adds the required arches for the periodic unit cell
%   This only applies for the unit cell!!!!

points = data.points;

if isfield(data, 'adjacency_matrix')
   adjacency_matrix = data.adjacency_matrix; %Might have issues when using finite
end

for i = 1:length(nodes_to_remove)
    node1 = nodes_to_remove(i);
    
    points(node1,1) = NaN;
    points(node1,2) = NaN;
    
    if isfield(data, 'adjacency_matrix')
        adjacency_matrix(node1,:) = NaN;
        adjacency_matrix(:,node1) = NaN;
    end
    
end



points_cleaned = points(~any(isnan(points), 2), :);
data.points = points_cleaned;

if isfield(data, 'adjacency_matrix')
    adjacency_matrix_cleaned = adjacency_matrix(~all(isnan(adjacency_matrix), 2), ~all(isnan(adjacency_matrix), 1));
    data.adjacency_matrix = adjacency_matrix_cleaned;
    data.N = sum(triu(adjacency_matrix,1) ==1,'all');   
end

end

