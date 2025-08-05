%%
arch_center_distance_mm = 120;
number_of_nodes = 5;
video_file_name = dir('*.mp4').name;

%%
close all
videoReader = VideoReader(video_file_name);
videoPlayer = vision.VideoPlayer('Position',[100,100,680,520]);
objectFrame = readFrame(videoReader);

%% Grab Values for the hinges and the rise
figure(1);  imshow(objectFrame); title("Select Left Hinge");
drawn_circle = drawcircle;
pause()
hinge_L = drawn_circle.Position;
hinge_L_radius = drawn_circle.Radius;
close(figure(1))

figure(1); imshow(objectFrame); title("Select Right Hinge");
drawn_circle = drawcircle;
pause()
hinge_R = drawn_circle.Position;
hinge_R_radius = drawn_circle.Radius;
close(figure(1))

% Determine axis vector as difference between hinges
axis_vector = hinge_R - hinge_L;

% Convert pixels to mm using distance between hinges (pixels) and known
% distance between hinges (mm)
mm_per_pixel = arch_center_distance_mm / norm(axis_vector);

% Normalize axis_vector for projection calculations
axis_vector = axis_vector/norm(axis_vector);

midpoint = 0.5*(hinge_L + hinge_R);
slope = (hinge_R(2) - hinge_L(2))/(hinge_R(1) - hinge_L(1));
perpindicular = - 1 / slope;

x_perp = -midpoint(2)/perpindicular + midpoint(1);


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
startingFrame = insertShape(startingFrame,'Line',[midpoint  x_perp 0],'Color','green','linewidth',4);
startingFrame = insertShape(startingFrame,'Line',[midpoint  x_perp 1920],'Color','green','linewidth',4);


% Add perpendicular lines to know where the grab the points from
for i = 1:number_of_nodes
    % Evenly divide axis to number_of_nodes+2 points and choose i+1th point
    point_on_axis = hinge_L + (i)/(number_of_nodes+1)*(hinge_R-hinge_L);
    % Find the x point perpendicular to the axis
    perpendicular_point = -point_on_axis(2)/perpindicular + point_on_axis(1);
    
    % Add perpendicular line to the axis
    startingFrame = insertShape(startingFrame,'Line',[point_on_axis  perpendicular_point 0],'Color','red','linewidth',4);
    startingFrame = insertShape(startingFrame,'Line',[point_on_axis  perpendicular_point 1920],'Color','red','linewidth',4);
end

imshow(startingFrame)
%% Grab the Nodes
points_vector = {};

for i = 1:number_of_nodes
    figure(1); imshow(startingFrame); title("Select Mode "+string(i));
    drawn_point=drawpoint;
    pause()
    points_vector{i} =drawn_point.Position;
    close all
end

%%
% Plot the Rectangles 
node0_intersection = {};
for i = 1:number_of_nodes
    point = points_vector{i};
    center_x = point(1);
    center_y = point(2);
    
    %Calculate the center of the "region" found from the rectangles
    node_center = [center_x center_y];
    
    % Determine the vector from the left hinge to the node center
    node_vector = node_center - hinge_L;
    
    % Project the node vector onto the axis_vector
    node_projection = dot(node_vector,axis_vector);
    
    % Determine the intersect of the node down to the axis
    node_intersect = hinge_L + node_projection*axis_vector;
    node0_intersection{i} = node_intersect;
    
    % Add the Line to the starting frame
    startingFrame = insertShape(startingFrame,'Line',[node_center  node_intersect],'Color','green','linewidth',4);
    
    % Add the text to denote the proportion
    node_distance_UL = node_projection / norm(hinge_R - hinge_L);
    node_distance_string = sprintf("%.3f",node_distance_UL);
    startingFrame = insertText(startingFrame,node_intersect,node_distance_string,"FontSize",32);

end

% Find the initial rise of the arch by prompting the user to select a point
imshow(startingFrame);
drawn_point = drawpoint;
pause()

rise_location = drawn_point.Position;
rise_distance_pixels = norm(rise_location - midpoint);
rise_distance_mm = rise_distance_pixels*mm_per_pixel;
rise_distance_mm_string = sprintf("%.3f",rise_distance_mm);

startingFrame = insertShape(startingFrame,'Circle',[rise_location 10],'Color','green','linewidth',4);
startingFrame = insertText(startingFrame,rise_location-[0 100],rise_distance_mm_string+" mm","FontSize",32);

imshow(startingFrame);
exportgraphics(gcf,"Starting Frame Marked Up.png")
pause()


%% Create a tracker object.
tracker = {};
for i = 1:number_of_nodes
    tracker{i} =vision.PointTracker('MaxBidirectionalError',1);
end

%% Initialize the tracker.
for i = 1:number_of_nodes
    initialize(tracker{i},points_vector{i},objectFrame);
end


%% Read, track, display points, and results in each video frame.
node_height_pixels = zeros(1,number_of_nodes);
counter = 1;

% Create vector to store the 'x' positions as they could 'drift'
x_position_UL      = zeros(1,number_of_nodes);

v = VideoWriter("Tracked Video.avi");
open(v)

while hasFrame(videoReader)
    frame = readFrame(videoReader);
    
    points = {};
    validity = {};
    
    for i = 1:number_of_nodes
        [points{i},validity{i}] = tracker{i}(frame);
    end
    
    validity_array = zeros(1,number_of_nodes);
    for i = 1:number_of_nodes
        validity_array(i) = mean(validity{i});
    end
    
    %% If we lose validity of half the points, regrab the values
    if any(validity_array <= 0.5)
        regrab = 1:number_of_nodes;                        % At snap through grab all the points again
        frame_regrab = frame;
        frame_regrab = insertShape(frame_regrab,'Line',[hinge_R  hinge_L],'Color','green','linewidth',4);

        % Plot the vertical node lines
        for k = 1:number_of_nodes
            % Evenly divide axis to number_of_nodes+2 points and choose i+1th point
%             point_on_axis = hinge_L + (k)/(number_of_nodes+1)*(hinge_R-hinge_L);          

            % Use the previosuly calculated 'x' values
            point_on_axis = hinge_L + x_position_UL(counter-1,k)*(hinge_R-hinge_L);
            
            % Find the x point perpendicular to the axis
            perpendicular_point = -point_on_axis(2)/perpindicular + point_on_axis(1);
                     
            % Add perpendicular line to the axis
            frame_regrab = insertShape(frame_regrab,'Line',[point_on_axis  perpendicular_point 0],'Color','red','linewidth',4);
            frame_regrab = insertShape(frame_regrab,'Line',[point_on_axis  perpendicular_point 1920],'Color','red','linewidth',4);
        end
        
        for i = regrab
            % Release the tracker
            release(tracker{i});
            
            % Reshow the figure and grab a new rectangle
            figure(2); imshow(frame_regrab);  
            title('Regrab Node ' +  string(i),"FontSize",14);
            drawn_point=drawpoint;
            pause()
            points_vector{i} =drawn_point.Position;
            
            close(figure(2))          
            
            % Reinitilize the Tracker
            tracker{i} =vision.PointTracker('MaxBidirectionalError',1);

            initialize(tracker{i},points_vector{i},frame);
            
            % Regrab the points
            [points{i},validity{i}] = tracker{i}(frame);
            
            validity_array(i) = sum(validity{i});
        end
    end
        
    for i = 1:number_of_nodes
        % Measure the mean of points distance to the node0_intersection
        approx_center = points{i};
        
        % Not Perpendicular 
%         node_height_pixels(counter,i) = norm(approx_center - node0_intersection{i})*sign(node0_intersection{i}(2) - approx_center(2));
        
        % Perpendicular
        node_vector = approx_center - hinge_L;

        % Project the node vector onto the axis_vector
        node_projection = dot(node_vector,axis_vector);

        % Determine the intersect of the node down to the axis
        node_intersect = hinge_L + node_projection*axis_vector;

        % Store the 'y' value in pixels
        node_height_pixels(counter,i) = norm(approx_center - node_intersect)*sign(node_intersect(2) - approx_center(2));
        
        % Store the 'x' value of the point in UL parameters
        node_distance_UL = node_projection / norm(hinge_R - hinge_L);
        x_position_UL(counter,i) = node_distance_UL;

        % Add a line from the center of node projected to axis
        frame = insertShape(frame,'Line',[approx_center  node_intersect],'Color','green','linewidth',4);
        
        % Add the points that have been tracked
        frame = insertMarker(frame,points{i}(validity{i}, :),'+');
        
        % Add a Circle at the center of each node
        frame = insertShape(frame,'Circle',[approx_center 10],'Color','blue','linewidth',4);
    end
    

    
    % Add some nice annotations to the frame
    frame = insertShape(frame,'Line',[hinge_R  hinge_L],'Color','green','linewidth',4);   
    
    % Add perpendicular line to the axis
    frame = insertShape(frame,'Line',[midpoint  x_perp 0],'Color','green','linewidth',4);
    frame = insertShape(frame,'Line',[midpoint  x_perp 1920],'Color','green','linewidth',4);
    
    videoPlayer(frame);
    writeVideo(v,frame)
    counter = counter + 1;
end
%% 
% Release the video player.
release(videoPlayer);
close(v)

%% Plot the Results
figure(1)
clf, hold on
for i = 1:number_of_nodes
    plot(node_height_pixels(:,i))

end
xlabel("Time")

save("Tracked Video.mat","node_height_pixels","mm_per_pixel","x_position_UL")