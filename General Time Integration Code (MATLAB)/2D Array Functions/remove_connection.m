function [data] = remove_connection(data,connections_to_remove)
%ADD_PERIODICITY Adds the required arches for the periodic unit cell
%   This only applies for the unit cell!!!!

adjacency_matrix = data.adjacency_matrix;
adjacency_matrix_finite = data.adjacency_matrix_finite;

for i = 1:length(connections_to_remove)
    node1 = connections_to_remove(i,1);
    node2 = connections_to_remove(i,2);
    
    adjacency_matrix(node1,node2) = 0;
    adjacency_matrix(node2,node1) = 0;

    adjacency_matrix_finite(node1,node2) = 0;
    adjacency_matrix_finite(node2,node1) = 0;
end

data.adjacency_matrix = adjacency_matrix;
data.adjacency_matrix_finite = adjacency_matrix_finite;

data.N = sum(triu(adjacency_matrix,1) ==1,'all');   

end

