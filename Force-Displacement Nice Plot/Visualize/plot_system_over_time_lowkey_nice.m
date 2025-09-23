%% ==============
function [] = plot_system_over_time_lowkey_nice(t,A,data)

%% Load Data
frames      = data.frames;
file_name   = data.file_name;
N           = data.N;
N_modes     = data.N_modes;
e_vector    = data.e_vector;
adjacency_matrix    = data.adjacency_matrix;
points              = data.points;

%% Plot Styles
arch_linewidth = 10;

%%  Interpolate the data evenly for the desired number of frames
tinterp = linspace(t(1),t(end),frames);
Ainterp = interp1(t,A,tinterp);

% Determine the starting and ending position for coloring purposes
for i = 1:N
    Apart = Ainterp(:,N_modes*i-(N_modes - 1):N_modes*i);
    [~,w0] = determine_shape_from_modes(Apart(1,:),0,1);
    [~,w1] = determine_shape_from_modes(Apart(end,:),0,1);

    starting_positions{i} = w0;
    ending_position{i} = w1;
end

%%
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
axis off

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
        Apart    = A(N_modes*i-(N_modes - 1):N_modes*i);
        Apartvel = A(N_modes*N+N_modes*i-(N_modes - 1):N_modes*N+N_modes*i);
    
        % Find x, w(x) for the data
        [xi,wi] = determine_shape_from_modes(Apart,e_vector(i),horiz_length);
    
        % Determine what color to make the arch
        arch_color = find_color(wi,starting_positions{i},ending_position{i});
        if anynan(arch_color)
            arch_color = [0 0 0 ];
        end
        % arch_color = find_color_velocity(Apartvel);

        % Scale the arch (if needed)
        wi = (1/pi)*wi*(.15*pi/data.b_vector(1));

        % Find the angle between the two endpoints
        angle = atan2(y_right-y_left,x_right-x_left);
        rotationMatrix = [cos(angle), -sin(angle); sin(angle), cos(angle)];
        
        rotated_xw = rotationMatrix * [xi;wi];
        
        figure(f); hold on
        line = plot(rotated_xw(1,:)+x_left,rotated_xw(2,:)+y_left,"linewidth",arch_linewidth,"color",arch_color,"LineStyle",'-');
        line_cell{i} = line;
    end

    scatter_atop = just_scatter(data);

    writeVideo(v, getframe(f));
    
    if j == frames
        exportgraphics(gcf,file_name + " - Last Frame.png","Resolution",600)            
    end

    for i = 1:N 
        delete(line_cell{i})
    end

    for i = 1:2
        delete(scatter_atop{i})
    end
    

end

%% ========================================================================
%  If making a video, close the video file
fprintf("CLOSING VIDEO\n")
close(v)
close(f);
end

function [color] = find_color(w,w0,w1)
    % Given starting position w0 and ending position w1
    % determine what the color ought to be based on how close
    % it is between the two
    max_diff = sum((w1 - w0).^2);
    cur_diff = sum((w1 - w ).^2);

    color1 = [229 75 75]/255;   % Red representing unstable
    color2 = [58 83 164]/255;    % Blue representing stable
    if cur_diff > max_diff
        cur_diff = max_diff;
    end
    color = (cur_diff/max_diff)*color1 + (1 - cur_diff/max_diff)*color2;
end


function [color] = find_color_velocity(Adotpart)
    % Given starting position w0 and ending position w1
    % determine what the color ought to be based on how close
    % it is between the two
    velocity = norm(Adotpart);
    maxvelocity  = 1e-3;

    color1 = [229 75 75]/255;   % Red representing unstable
    % color2 = [58 83 164]/255;    % Blue representing stable
    color2 = 0.25*[1 1 1];
    if velocity > maxvelocity
        velocity = maxvelocity;
    end
    color = (velocity/maxvelocity)*color1 + (1 - velocity/maxvelocity)*color2;
end


function [f] = just_scatter(data)
node_big_circle_color = 0.75*[1 1 1];
node_big_circle_size  = 400;

node_little_circle_color = 'k';
node_little_circle_size = 100;

x = data.points(:,1);
y = data.points(:,2);


s1 = scatter(x,y,node_big_circle_size, "MarkerFaceColor",node_big_circle_color,"MarkerEdgeColor","k");
s2 = scatter(x,y,node_little_circle_size, "MarkerFaceColor",node_little_circle_color,"MarkerEdgeColor","k");

f = {s1 s2};

end