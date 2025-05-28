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


addpath("Shapes Point Data\")
% Loads points2 to describe the full system
run("points_hexagon2")
adjacency_matrix2 = determine_adjacency_matrix(points2);

% The ith index of expand points to the arch number of the expanded system
expand = [1 2 3 4 5 6 7 8 9 10 11 12 13 -19 -20 14 15 -21 -23 16 17 18];
%% Make the grid without any of the deformed arches plotted
f = plot_grid(adjacency_matrix2,points2,false);

%% Add the shape of each arch
up_adjac = triu(adjacency_matrix,1);
[left, right] = find(up_adjac == 1);

up_adjac2 = triu(adjacency_matrix2,1);
[left2, right2] = find(up_adjac2 == 1);

% Load the x and y position of the endpoints
% Need to load from the full system
x = points2(:,1);
y = points2(:,2);


% Go arch by arch
for i = 1:N
    i2 = expand(i); % In the expanded system, we really want to plot the i2th arch

    if i2 < 0
            % Determine the left and right coordinates of the ends of the arch
        i2 = -1*i2;
        x_left = x(right2(i2));
        x_right = x(left2(i2));
        y_left = y(right2(i2));
        y_right = y(left2(i2));
    else
        % Determine the left and right coordinates of the ends of the arch
        x_left = x(left2(i2));
        x_right = x(right2(i2));
        y_left = y(left2(i2));
        y_right = y(right2(i2));
    end



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

