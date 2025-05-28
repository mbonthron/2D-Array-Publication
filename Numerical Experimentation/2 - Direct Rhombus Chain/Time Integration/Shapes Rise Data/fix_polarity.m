function [data] = fix_polarity(data)
%FIX_POLARITY Given an e_vector of constant
%   rise, determines which elements are neighbors and alternates the
%   direction of plastic deformation so that a wave can travel through the
%   structure
%
%   INPUTS
%   ===================================================
%   e_vector: Nx1 vector with all positive values describing the magnitude
%       of plastic deformation. In this function we determine the direction
%   adjacency_matrix: Describes if the i^th node is connected to the j^th
%       nodes
%
%   OUTPUTS
%   ===================================================
%   e_vector_fixed: Nx1 vector with correct directions (local coordinate)
%       describing the plastic deformation of each element
%
%%
b_vector = data.b_vector;
e_vector = data.e_vector;
adjacency_matrix = data.adjacency_matrix;
points = data.points;

%%
% Make the grid
plot_grid(data,true)

x = points(:,1);
y = points(:,2);

% Create a vector defining the "spin" of the point
%   0: No value has been written yet
%   1: Spin "up" - out of the page moment 
%   -1: Spin "down" - into the page moment
point_count = length(x);
polarity = zeros(point_count,1);

% Add Arch number label
polarity(1) = 1;

% Use adjacency matrix to assign opposite polarity to neighboring points
while any(polarity == 0)
    for i = 1:point_count
        if polarity(i) ~= 0
            neighboring_points = adjacency_matrix(i,:) == 1;
            polarity(neighboring_points) = -1 * polarity(i);
        end
    end
end


up_adjac = triu(adjacency_matrix,1);
N = sum(up_adjac,'all');

[left, right] = find(up_adjac == 1);

% Go arch by arch and look at the polarity of each element's left node
% if polarity of left element is positive -> e > 0
% if polarity of left element is negative -> e < 0
e_vector_fixed = zeros(size(e_vector));
b_vector_fixed = zeros(size(b_vector));

e_vector_fixed( polarity(left) > 0 ) =    e_vector( polarity(left) > 0 );
e_vector_fixed( polarity(left) < 0 ) = -1*e_vector( polarity(left) < 0 );

b_vector_fixed( polarity(left) > 0 ) =    b_vector( polarity(left) > 0 );
b_vector_fixed( polarity(left) < 0 ) = -1*b_vector( polarity(left) < 0 );


data.b_vector = b_vector_fixed;
data.e_vector = e_vector_fixed;
end

