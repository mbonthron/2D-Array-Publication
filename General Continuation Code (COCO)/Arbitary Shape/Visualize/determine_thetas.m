function [theta_vector] = determine_thetas(A,data)
%DETERMINE_THETAS Summary of this function goes here
%   Determines the theta's of the hinges
%% Load from data
adjacency_matrix        = data.adjacency_matrix;

V  = data.V;    % Number of Hinges
N  = data.N;    % Number of Arches
N_modes = data.N_modes;     % Number of Modes
%% Determine Numbering convention of arches and hinges
up_adjac = triu(adjacency_matrix,1);
[left, right] = find(up_adjac == 1);

left_and_right = [left, right];

% Vector which describes all the arch numbers
arches = 1:N;

% Set up empty theta_vector to be populated
theta_vector = zeros(V,1);

% Iterate over each hinge
for i = 1:V
    % Find all arches that intersect the hinge
    all_arches_at_i = arches(any(left_and_right == i,2));

    % Take the lowest number arch that is at the hinge
    lowest_arch = min(all_arches_at_i);

    % Find the mode coefficients associated with this arch
    Apart = A(N_modes*lowest_arch-(N_modes - 1):N_modes*lowest_arch)';

    % Find if the lowest arch is denoted as being 'left' or 'right' 
    % at the hinge
    if left(lowest_arch) == i
        % Then the arch counts the ith hinge as its left hinge
        slope = sum((1:N_modes).*Apart);
    else
        % Then the arch counts the ith hinge as its right hinge
        slope = sum((1:N_modes).*((-1).^(1:N_modes)).*Apart);
    end

    % Save the arctan of the slope
    theta_vector(i) = atand(slope);
end

end

