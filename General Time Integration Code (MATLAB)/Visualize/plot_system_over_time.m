function [] = plot_system_over_time(t,A,data)
%PLOT_SYSTEM_OVER_TIME For a time vector and number of frames, interpolates
%   the matrix of data A and makes a little video of the frames
%
%   INPUTS
%   ===================================================
%   t:  time vector of the data
%   A:  matrix of the modal responses in state space form
%   frames: number of frames for the video 
%   file_name: file name to save the video as
%   N:  number of arches in the system
%   N_modes: Number of modes to describe each arch
%   e_vector: vector of values associated with plastic deformation of the
%       arch
%   adjacency_matrix: adjacency matrix describing connections of the system
%   points: x,y coordinate pairs of the nodes of the system
%% Load Data
frames = data.frames;
file_name = data.file_name;
N = data.N;
N_modes = data.N_modes;
e_vector = data.e_vector;
adjacency_matrix = data.adjacency_matrix;
points = data.points;

%% Plot Styles
arch_color = 'k';
arch_linewidth = 2;

%%  Interpolate the data evenly for the desired number of frames
tinterp = linspace(t(1),t(end),frames);
Ainterp = interp1(t,A,tinterp);

if ~exist("Videos\", 'dir')
    mkdir("Videos\");
end

file_name = "Videos\"+file_name;

%  Initialize the video write and open the video
fprintf("CREATING VIDEO\n")
v = VideoWriter(file_name);
open(v)

% Make the empty grid
f = plot_grid(data,false);
f.Position(3:4) = [2500 1000];


up_adjac = triu(adjacency_matrix,1);
[left, right] = find(up_adjac == 1);

% Load the x and y position of the endpoints
x = points(:,1);
y = points(:,2);



for j = 1:frames
    A = Ainterp(j,:);
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
    
        wi = (1/pi)*wi;

        % Find the angle between the two endpoints
        angle = atan2(y_right-y_left,x_right-x_left);
        rotationMatrix = [cos(angle), -sin(angle); sin(angle), cos(angle)];
        
        rotated_xw = rotationMatrix * [xi;wi];
        
        figure(f); hold on
        line = plot(rotated_xw(1,:)+x_left,rotated_xw(2,:)+y_left,"linewidth",arch_linewidth,"color",arch_color,"LineStyle",'-');
        line_cell{i} = line;
    end

    writeVideo(v, getframe(f));
    
    if j == frames
        exportgraphics(gcf,file_name + " - Last Frame.png","Resolution",600)            
    end

    for i = 1:N 
        delete(line_cell{i})
    end
    

end

%% ========================================================================
%  If making a video, close the video file
fprintf("CLOSING VIDEO\n")
close(v)
close(f);
end
