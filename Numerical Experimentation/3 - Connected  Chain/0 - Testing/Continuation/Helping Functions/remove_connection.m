function [data] = remove_connection(i,j,data)
%ADD_PERIODICITY Adds the required arches for the periodic unit cell
%   This only applies for the unit cell!!!!

adjacency_matrix = data.adjacency_matrix;

adjacency_matrix(i,j) = 0;
adjacency_matrix(j,i) = 0;

data.adjacency_matrix = adjacency_matrix;

end

