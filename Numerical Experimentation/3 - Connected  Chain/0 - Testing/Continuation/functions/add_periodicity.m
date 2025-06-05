function [data] = add_periodicity(data)
%ADD_PERIODICITY Adds the required arches for the periodic unit cell
%   This only applies for the unit cell!!!!
connections_to_add = [1 11;
                      1 12;
                      1 10];

adjacency_matrix = data.adjacency_matrix;

for i = 1:length(connections_to_add)
    node1 = connections_to_add(i,1);
    node2 = connections_to_add(i,2);
    
    adjacency_matrix(node1,node2) = 1;
    adjacency_matrix(node2,node1) = 1;
end

data.adjacency_matrix = adjacency_matrix;
data.N = sum(triu(adjacency_matrix,1) ==1,'all');   

end

