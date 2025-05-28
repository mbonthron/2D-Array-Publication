function [f] = plot_system_once(A,N,N_modes,e_vector,adjacency_matrix,points)
%PLOT_SYSTEM_ONCE Plots the system for a single point
%
%   INPUTS
%   ===================================================
%   A: 2*N row vector corresponding to the magnitude of each mode;
%   N: Number of arches
%   N_modes: Number of modes describing each arch
%   e_vector: Nx1 vector containing the magnitude of plastic deformation for each arch
%   adjacency_matrix: Describing the connection of the 
%   points: Nx2 vector with the x and y position of each node in 1st and 2nd column (respectively)

%% Plot Styles
arch_color = 'k';
arch_linewidth = 4;

%% Make the grid without any of the deformed arches plotted
f = plot_grid(adjacency_matrix,points,false);

%% Add the shape of each arch
up_adjac = triu(adjacency_matrix,1);
[left, right] = find(up_adjac == 1);

% Load the x and y position of the endpoints
x = points(:,1);
y = points(:,2);



for i = 1:N
    % Determine the left and right coordinates of the ends of the arch
    x_left = x(left(i));
    x_right = x(right(i));
    y_left = y(left(i));
    y_right = y(right(i));

    % Determine the horizontal length between the supports
    horiz_length = sqrt((y_right-y_left)^2+(x_right-x_left)^2);

    % Take the portion of A that corresponds to the ith arch
    Apart = A(N_modes*i-(N_modes - 1):N_modes*i);

    % Find x, w(x) for the data
    [xi,wi] = determine_shape_from_modes(Apart,e_vector(i),horiz_length);

    % Find the angle between the two endpoints
    angle = atan2(y_right-y_left,x_right-x_left);
    rotationMatrix = [cos(angle), -sin(angle); sin(angle), cos(angle)];
    
    rotated_xw = rotationMatrix * [xi;wi];
    
    figure(f); hold on
    line = plot(rotated_xw(1,:)+x_left,rotated_xw(2,:)+y_left,"linewidth",arch_linewidth,"color",arch_color,"LineStyle",'-');
end


end

