%% Clear Everything so there are no stragglers
clear; clc; close all

%% Add the Paths to the Required Functions
%restoredefaultpath
startup
addpath('..\General Time Integration Code (MATLAB)\Visualize')
addpath('..\General Optimization\')
addpath('..\General Optimization\Shapes Point Data\')
addpath('..\General Time Integration Code (MATLAB)\2D Array Functions')
addpath('..\General Continuation Code (COCO)\Arbitary Shape\functions\')
%addpath('..\General Continuation Code (COCO)\Arbitary Shape\Visualize\')
%addpath('..\General Continuation Code (COCO)\Arbitary Shape\')
%addpath('Shapes Point Data/')
%addpath("Debugging\")
%% Create Empty Data Structure to be Populated
data = struct();
data.N_modes = 3;   % Number of modes used to describe the system
data.N_cells = 1;
data.plot_grids = 1;

%% Run Continuation to Get Stable Configurations at each b
% Choose which shape
shapeNum = 8;
data.static_only = 1;
data = init_shape(shapeNum, data);

% Get points for 2D picture
data.arch_center_distance_mm = 120;
data.number_of_points = 1;
[picture_points_height, picture_points_position] = picture2points(data);

% Run COCO and find stable configurations   
    % Do we care which configuration? All of them?
    % If all of them, then load picture for each configuration?
    % If multiple, would be nice to have a mapping from physical picture number to
    % COCO number
    % Autodetection of mapping based on coarsely autodetecting which are
    % most similar would be good (ML?)

t = tvals(t_idx);
% Run COCO
data = raw_data;
data.t = t;
data.t_vector = t*ones(data.N,1);
bpoints = OG_bpoints;

runs_exist = 1;
troubleshooting_flag = 0;
for b = bpoints
    if ~isfile(data.shape_name + " b = "+ num2str(round(b/pi,4)) +" pi t = "+num2str(round(t/pi,4)) +" pi.mat")
        runs_exist = 0;
    end
end
if ~runs_exist || troubleshooting_flag
    [data,run_max_E_per_b,bpoints] = general_COCO(data, bpoints);
end

% Show plot of configuration with left and right side of each arch labeled
% To make labeling the photo easier

% Load picture and manually select points of configuration (In new figure)
% Maybe ideally subfigure to display both photo and labels

    % Option 1?:
    % For each arch label left and right corresponding to COCO
    % configuration, then points
    % At this stage we have all the positional photo data

    % Option 2?:
    % For first arch label left and right corresponding to COCO
    % Autocomplete each remaining left right, based on point data
    % Give user the option to change circles around to better fit picture
    % For each arch match point
    % At this stage we have all the positional photo data

    % Option 3?:
    % For first arch label left and right corresponding to COCO
    % Suggest points to pick based on COCO values (converted) with user
    % finetuning
    % For each following arch, suggest (based on COCO) left right, then points 
    % with user finetuning
    % At this stage we have all the positional photo data
   
% Convert into same units
% Plot Percentage error per arch

% Overlay Real vs COCO Arches on the same plot to see differences