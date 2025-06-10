%% Clear Everything so there are no stragglers
clear; clc; close all

%% Add the Paths to the Required Functions
restoredefaultpath
startup
addpath('..\..\General Time Integration Code (MATLAB)\Visualize')
addpath('..\..\General Time Integration Code (MATLAB)\2D Array Functions')
addpath('..\..\General Continuation Code (COCO)\Arbitary Shape\functions\')
addpath('..\..\General Continuation Code (COCO)\Arbitary Shape\Visualize\')
addpath('..\..\General Continuation Code (COCO)\Arbitary Shape\')
addpath('..\Shapes Point Data/')
addpath('..\')
data = struct();
data.N_modes = 3;   % Number of modes used to describe the system
data.N_cells = 2;
data.plot_grids = 1;
for i=1:44
    data = init_shape(4,data); data = initialize_time_integration(0.2*pi,data,i); plot_system_once(data.A0, data);
end
