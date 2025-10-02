global file_name
close all; clc
%% Measure Distances in an Image
% Read image into the workspace. 
file_name = '20250922_090753';
im = imread(file_name+".jpg");


mm_measured = 10;      % [mm]
pixels_measured = 489;

fprintf("1.) CLICK AND DRAG TO MEASURE A DISTANCE\n")
fprintf("2.) DOUBLE LEFT CLICK TO SET AS KNOWN DISTANCE\n")
fprintf("3.) RIGHT CLICK ON MEASUREMENT TO SAVE FIGURE\n")

%%
% Gather data about the image, such as its size, and store the data in a
% structure that you can pass to callback functions.
sz = size(im);
myData.Units = 'mm';
myData.MaxValue = hypot(sz(1),sz(2));
myData.Colormap = hot;

% mm_per_pixel = mm_measured / pixels_measured;

% myData.ScaleFactor = mm_per_pixel;
myData.ScaleFactor = 1;

%% ADD THE AXIS CIRCLES AND PLOT TO CONNECT THE DOTS
% Display the image in an axes.
im = imrotate(im,90);

[y,x,z] = size(im);

% Grab left hinge
hinge_count = 3;
hinge_position = [];

for i = 1:hinge_count
    f1 = figure(1); hIm = imshow(im); title("Select Hinge");
    drawn_circle = drawpoint;
    pause()
    hinge_position(i,:) = drawn_circle.Position;
    hinge_radius(i,:)   = drawn_circle.Radius;
    close(figure(1))

end


%%
im = imread(file_name+".jpg");
im = imrotate(im,90);


adjacency_matrix = [0 1 1 ; 1 0 1 ; 1 1 0];
up_adjac = triu(adjacency_matrix,1);
[left, right] = find(up_adjac == 1);

% Add Hinges
for i = 1:length(adjacency_matrix)
    % Plot the two hinges
    im = insertShape(im,'Circle',[hinge_position(i,:) hinge_radius(i,:)],'Color','green','linewidth',4);
end


% Add Centerlines
axis_length = [];
for i = 1:length(left)
    hinge_R = hinge_position(left(i),:);
    hinge_L = hinge_position(right(i),:);

    axis_vector = hinge_R - hinge_L;

    axis_length(i) = norm(axis_vector);

    % Plot the x axis
    im = insertShape(im,'Line',[hinge_R  hinge_L],'Color','green','linewidth',4);

    midpoint = 0.5*(hinge_L + hinge_R);
    im = insertText(im,midpoint,[num2str(axis_length(i),'%.0f') ' pixels'], ...
        FontSize = 48,TextBoxColor="w");

end


hIm = imshow(im);
% Determine axis vector as difference between hinges
% axis_vector = hinge_R - hinge_L;
% 
% midpoint = 0.5*(hinge_L + hinge_R);
% slope = (hinge_R(2) - hinge_L(2))/(hinge_R(1) - hinge_L(1));
% perpindicular = - 1 / slope;
% 
% x_perp1 = -midpoint(2)/perpindicular + midpoint(1);
% x_perp2 = (1920-midpoint(2))/perpindicular + midpoint(1);



%%
% Specify a callback function for the |ButtonDownFcn| callback on the
% image. Pass the |myData| structure to the callback function. This
% callback function creates the line objects and starts drawing the ROIs.
hIm.ButtonDownFcn = @(~,~) startDrawing(hIm.Parent,myData);

%% Create Callback Function to Start Drawing ROIs
% Create the function used with the |ButtonDownFcn| callback to create line
% ROIs. This function: 
%
% 1. Instantiates a line ROI object. 
%
% 2. Sets up listeners to react to clicks and movement of the ROI. 
%
% 3. Adds a custom context menu to the ROIs that includes a 'Delete All'
%    option. 
%
% 4. Begins drawing the ROI, using the point clicked in the image as the
%    starting point.
% 
function startDrawing(hAx,myData)

% Create a line ROI object. Specify the initial color of the line and
% store the |myData| structure in the |UserData| property of the ROI.
h = images.roi.Line('Color',[0, 0, 0.5625],'UserData',myData);

% Set up a listener for movement of the line ROI. When the line ROI moves,
% the |updateLabel| callback updates the text in the line ROI label and
% changes the color of the line, based on its length.
addlistener(h,'MovingROI',@updateLabel);

% Set up a listener for clicks on the line ROI. When you click on the line
% ROI, the |updateUnits| callback opens a GUI that lets you specify the
% known distance in real-world units, such as, meters or feet. 
addlistener(h,'ROIClicked',@updateUnits);


% Get the current mouse location from the |CurrentPoint| property of the
% axes and extract the _x_ and _y_ coordinates.
cp = hAx.CurrentPoint;
cp = [cp(1,1) cp(1,2)];

% Begin drawing the ROI from the current mouse location. Using the
% |beginDrawingFromPoint| method, you can draw multiple ROIs.
h.beginDrawingFromPoint(cp);

% Add a custom option to the line ROI context menu to delete all existing
% line ROIs.
c = h.UIContextMenu;
uimenu(c,'Label','Delete All','Callback',@deleteAll);
uimenu(c,'Label','Save','Callback',@saveAll);
end

%%
%
% <<../drawROIs.gif>>
%

%% Create Callback Function to Update ROI Label and Color
% Create the function that is called whenever the line ROI is moving, that
% is, when the |'MovingROI'| event occurs. This function updates the ROI
% label with the length of the line and changes the color of the line based
% on its length.
% 
% This function is called repeatedly when the ROI moves. If you want to
% update the ROI only when the movement has finished, listen for the
% |'ROIMoved'| event instead.
function updateLabel(src,evt)

% Get the current line position.
pos = evt.Source.Position;

% Determine the length of the line.
diffPos = diff(pos);
mag = hypot(diffPos(1),diffPos(2));

% Choose a color from the colormap based on the length of the line. The
% line changes color as it gets longer or shorter.

% Apply the scale factor to line length to calibrate the measurements.
mag = mag*src.UserData.ScaleFactor;
color = [1 1 1];
% Update the label.
set(src,'Label',[num2str(mag,'%30.1f') ' ' src.UserData.Units],'Color',color);

end

%% Create Callback Function to Update Measurement Units
% Create the function that is called whenever you double-click the ROI
% label. This function opens a popup dialog box in which you can enter
% information about the real-world distance and units. 
% 
% This function listens for the |'ROIClicked'| event, using event data to
% check the type of click and the part of the ROI that was clicked.
%
% The popup dialog box prompts you to enter the known distance and
% units for this measurement. With this information, you can calibrate all
% the ROI measurements to real world units. 
function updateUnits(src,evt)

% When you double-click the ROI label, the example opens a popup dialog box
% to get information about the actual distance. Use this information to
% scale all line ROI measurements.
if strcmp(evt.SelectionType,'double') && strcmp(evt.SelectedPart,'label')
    
    % Display the popup dialog box.
    answer = inputdlg({'Known distance','Distance units'},...
        'Specify known distance',[1 20],{'10','meters'});
    
    % Determine the scale factor based on the inputs.
    num = str2double(answer{1});
    
    % Get the length of the current line ROI.
    pos = src.Position;
    diffPos = diff(pos);
    mag = hypot(diffPos(1),diffPos(2));
    
    % Calculate the scale factor by dividing the known length value by the
    % current length, measured in pixels.
    scale = num/mag;
    
    % Store the scale factor and the units information in the |myData|
    % structure.
    myData.Units = answer{2};
    myData.MaxValue = src.UserData.MaxValue;
    myData.Colormap = src.UserData.Colormap;
    myData.ScaleFactor = scale;
    
    % Reset the data stored in the |UserData| property of all existing line
    % ROI objects. Use |findobj| to find all line ROI objects in the axes.
    hAx = src.Parent;
    hROIs = findobj(hAx,'Type','images.roi.Line');
    set(hROIs,'UserData',myData);
    
    % Update the label in each line ROI object, based on the information
    % collected in the input dialog.
    for i = 1:numel(hROIs)
        
        pos = hROIs(i).Position;
        diffPos = diff(pos);
        mag = hypot(diffPos(1),diffPos(2));


        
        set(hROIs(i),'Label',[num2str(mag*scale,'%30.1f') ' ' answer{2}]);
        
    end
    
    % Reset the |ButtonDownFcn| callback function with the current |myData|
    % value.
    hIm = findobj(hAx,'Type','image');
    hIm.ButtonDownFcn = @(~,~) startDrawing(hAx,myData);
    
end

end

%%
%
% <<../setUnits.gif>>
%

%% Create Callback Function to Delete All ROIs
% Create the function to delete all ROIs. You added a custom context menu
% item to each line ROI in the |startDrawing| callback function. This is
% the callback associated with that custom context menu. This callback uses
% the |findobj| function to search for the ROI Type and deletes any found
% ROIs.
function deleteAll(src,~)

hFig = ancestor(src,'figure');
hROIs = findobj(hFig,'Type','images.roi.Line');
delete(hROIs)

end


function saveAll(src,~)
global file_name
fig1 = gcf;
% fig1.Units = "inches";
% fig1.Position(3:4) = [4 4];
set(fig1,"Color","white")

save_name = file_name + " Measured.png";
xlabel(upper(save_name))

exportgraphics(gcf,save_name,"Resolution",300)


end

%%
%
% <<../deleteAll.gif>>
%