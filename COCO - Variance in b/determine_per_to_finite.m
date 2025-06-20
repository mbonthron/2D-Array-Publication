function [data] = determine_per_to_finite(data)
%DETERMINE_PER_TO_FINITE Summary of this function goes here
%   Detailed explanation goes here
%%
N = data.N;     % Number of arches (periodic)

points = data.points;   % Points of periodic
points_finite  = data.points_finite; % Point of finite

adjacency_matrix = data.adjacency_matrix;
adjacency_matrix_finite = data.adjacency_matrix_finite;

% Create empty vector 'expand'
% the index refers to the arch number of the periodic system
% the value in the index maps to the finite system
% e.g. [1 3 -2] means
%   1 periodic -> 1 finite
%   2 periodic -> 3 finite
%   3 periodic -> 2 finite with a flipped sign convention
expand = zeros(1,N);

%% Finite the arch numbering convention for periodic and finite structure
up_adjac = triu(adjacency_matrix,1);
[left, right] = find(up_adjac == 1);

up_adjac_finite = triu(adjacency_matrix_finite,1);
[left_finite, right_finite] = find(up_adjac_finite == 1);

%% Change left_finite and right_finite to hold the mapped values from points
vertex_map_p2f = data.vertex_map_p2f;

[m,n] = size(vertex_map_p2f);

for i = 1:m
    periodic_vertex = vertex_map_p2f(i,1);
    finite_vertex   = vertex_map_p2f(i,2);
    
    left_finite(left_finite==finite_vertex) = periodic_vertex;
    right_finite(right_finite==finite_vertex) = periodic_vertex;
end

%% Go through the periodic arch and try to find an element that matches left and right
leftright = [left_finite right_finite];
for i = 1:N
    % i is the arch number
    left_hinge = left(i);
    right_hinge = right(i);

    % Find the row in finite structure that shares the hinges
    % from the periodic structure
    matched = find(sum(ismember(leftright,[left_hinge right_hinge]),2)==2);

    % Check if the left and right hinge match
    if left_hinge == left_finite(matched) && right_hinge == right_finite(matched)
        % Left and Right Convention remains
        expand(i) = matched;
    else
        % Need to flip left and right convention
        expand(i) = -matched;
    end
end

data.expand = expand;
end