function [f] = COCO_plot_system_once(A,data)
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

%%
N = data.N;
N_modes = data.N_modes;
e_vector = zeros(data.N,1);
adjacency_matrix = data.adjacency_matrix;
adjacency_matrix_finite = data.adjacency_matrix_finite;

points = data.points;
points_finite = data.points_finite;

%% Plot Styles
arch_color = 'k';
arch_linewidth = 4;


% The ith index of expand points to the arch number of the expanded system
expand = data.expand;

%% Make the grid without any of the deformed arches plotted
data2.points = points_finite;
data2.adjacency_matrix = adjacency_matrix_finite;

f = plot_grid(data2,true);

%% Add the shape of each arch
up_adjac = triu(adjacency_matrix,1);
[left, right] = find(up_adjac == 1);

up_adjac2 = triu(adjacency_matrix_finite,1);
[left2, right2] = find(up_adjac2 == 1);

% Load the x and y position of the endpoints
% Need to load from the full system
x = points_finite(:,1);
y = points_finite(:,2);


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

    wi = (1/pi)*wi;

    % Find the angle between the two endpoints
    angle = atan2(y_right-y_left,x_right-x_left);
    rotationMatrix = [cos(angle), -sin(angle); sin(angle), cos(angle)];
    
    rotated_xw = rotationMatrix * [xi;wi];
    
    figure(f); hold on
    line = plot(rotated_xw(1,:)+x_left,rotated_xw(2,:)+y_left,"linewidth",arch_linewidth,"color",arch_color,"LineStyle",'-');
end

end

