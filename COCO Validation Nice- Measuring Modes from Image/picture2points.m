function [picture_points_height, picture_points_position] = picture2points(data, file_name)
%PICTURE2POINTS Summary of this function goes here
%   Detailed explanation goes here
% N = 11; %Number of arches
% arch_center_distance_mm = 120;
% number_of_points = 5;

N = data.N;
arch_center_distance_mm = data.arch_center_distance_mm;
number_of_points = data.number_of_points;



objectFrame = imread(file_name);

f = figure(1); %clf;
f.Position = [10 10 900 600];
for arch_num=1:N
    %% Grab Values for the hinges and the rise
    
    
    t = tiledlayout(1,2);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    nexttile;
    figure(1);
    plot_grid(data,1,1);
    % hold off;
    % figure(1);
    nexttile;
    imshow(objectFrame); title("Select Left Hinge of Arch " + num2str(arch_num));
    drawn_circle = drawcircle;
    pause()
    hinge_L = drawn_circle.Position;
    hinge_L_radius = drawn_circle.Radius;

    figure(1); imshow(objectFrame); title("Select Right Hinge " + num2str(arch_num));
    drawn_circle = drawcircle;
    pause()
    hinge_R = drawn_circle.Position;
    hinge_R_radius = drawn_circle.Radius;

    % Determine axis vector as difference between hinges
    axis_vector = hinge_R - hinge_L;

    % Convert pixels to mm using distance between hinges (pixels) and known
    % distance between hinges (mm)
    mm_per_pixel = arch_center_distance_mm / norm(axis_vector);

    % Normalize axis_vector for projection calculations
    axis_vector = axis_vector/norm(axis_vector);

    midpoint = 0.5*(hinge_L + hinge_R);

    % Length of perpendicular line (you can change this)
    len = 400;

    v_perp = [-axis_vector(2), axis_vector(1)];
    v_perp = v_perp / norm(v_perp);

    % Set perpendicular line length equal to original
    perp_1 = midpoint + v_perp * len/2;
    perp_2 = midpoint - v_perp * len/2;


    %% Markup the first frame
    figure(1);
    startingFrame = objectFrame;

    % Plot the two hinges
    startingFrame = insertShape(startingFrame,'Circle',[hinge_L hinge_L_radius],'Color','green','linewidth',4);
    startingFrame = insertShape(startingFrame,'Circle',[hinge_R hinge_R_radius],'Color','green','linewidth',4);

    % Plot the x axis
    startingFrame = insertShape(startingFrame,'Line',[hinge_R  hinge_L],'Color','green','linewidth',4);

    % Plot the midpoint of the axis
    startingFrame = insertShape(startingFrame,'Circle',[midpoint 10],'Color','green','linewidth',4);

    % Add perpendicular line to the axis
    %startingFrame = insertShape(startingFrame,'Line',[midpoint  x_perp 1920],'Color','green','linewidth',4);
    %startingFrame = insertShape(startingFrame,'Line',[midpoint  x_perp y_perp],'Color','green','linewidth',4);
    startingFrame = insertShape(startingFrame,'Line',[perp_1 perp_2],'Color','green','linewidth',4);


    % Add perpendicular lines to know where the grab the points from
    for i = 1:number_of_points
        % Evenly divide axis to number_of_nodes+2 points and choose i+1th point
        point_on_axis = hinge_L + (i)/(number_of_points+1)*(hinge_R-hinge_L);
        % Find the x point perpendicular to the axis
        %perpendicular_point = -point_on_axis(2)/perpindicular + point_on_axis(1);
        perp_1 = point_on_axis + v_perp * len/2;
        perp_2 = point_on_axis - v_perp * len/2;

        % Add perpendicular line to the axis
        %startingFrame = insertShape(startingFrame,'Line',[point_on_axis  perpendicular_point 0],'Color','red','linewidth',4);
        %startingFrame = insertShape(startingFrame,'Line',[point_on_axis  perpendicular_point 1920],'Color','red','linewidth',4);
        startingFrame = insertShape(startingFrame,'Line',[perp_1  perp_2],'Color','red','linewidth',4);
    end

    imshow(startingFrame)
    %% Grab the Nodes
    points_vector = {};

    for i = 1:number_of_points
        figure(1); imshow(startingFrame); title("Select Point "+string(i));
        drawn_point=drawpoint;
        pause()
        approx_center =drawn_point.Position;
        
        % Not Perpendicular 
%         node_height_pixels(counter,i) = norm(approx_center - node0_intersection{i})*sign(node0_intersection{i}(2) - approx_center(2));
        
        % Perpendicular
        node_vector = approx_center - hinge_L;

        % Project the node vector onto the axis_vector
        node_projection = dot(node_vector,axis_vector);

        % Determine the intersect of the node down to the axis
        node_intersect = hinge_L + node_projection*axis_vector;

        % Store the 'y' value in pixels
        sign_of_height = -1*cross([node_vector-node_projection*axis_vector, 0],[node_projection*axis_vector,0]); % Get sign from cross product
        picture_points_height(arch_num,i) = mm_per_pixel * norm(approx_center - node_intersect)*sign(sign_of_height(3)) % Accounts for both possibilities of directionality 
        
        % Store the 'x' value of the point in UL parameters
        node_distance_UL = node_projection / norm(hinge_R - hinge_L);
        picture_points_position(arch_num,i) = node_distance_UL;

        % % Add a line from the center of node projected to axis
        % frame = insertShape(frame,'Line',[approx_center  node_intersect],'Color','green','linewidth',4);
        % 
        % % Add the points that have been tracked
        % frame = insertMarker(frame,points{i}(validity{i}, :),'+');
        % 
        % % Add a Circle at the center of each node
        % frame = insertShape(frame,'Circle',[approx_center 10],'Color','blue','linewidth',4);
    end


end
end

