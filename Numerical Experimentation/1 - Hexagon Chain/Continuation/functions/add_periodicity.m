function [adjacency_matrix] = add_periodicity(adjacency_matrix)
%ADD_PERIODICITY Adds the required arches for the periodic unit cell
%   This only applies for the unit cell!!!!
connections_to_add = [1 9;
                      2 9;
                      2 10;
                      3 10];

for i = 1:length(connections_to_add)
    node1 = connections_to_add(i,1);
    node2 = connections_to_add(i,2);
    
    adjacency_matrix(node1,node2) = 1;
    adjacency_matrix(node2,node1) = 1;

end
end

