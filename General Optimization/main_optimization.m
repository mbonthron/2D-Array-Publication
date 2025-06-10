%% Clear Everything so there are no stragglers
clear; clc; close all

%% Add the Paths to the Required Functions
% addpath('2D Array Functions')
% addpath('COCO Continuation/Shapes Point Data/')
% addpath('Shapes Rise Data')
% addpath('Visualize')
% addpath('COCO Continuation/Visualize')
% addpath('COCO Continuation/functions')
% addpath('COCO Continuation')
% addpath('Time Integration')
restoredefaultpath
startup
addpath('..\General Time Integration Code (MATLAB)\Visualize')
addpath('..\General Time Integration Code (MATLAB)\2D Array Functions')
addpath('..\General Continuation Code (COCO)\Arbitary Shape\functions\')
addpath('..\General Continuation Code (COCO)\Arbitary Shape\Visualize\')
addpath('..\General Continuation Code (COCO)\Arbitary Shape\')
addpath('Shapes Point Data/')
%% Create Empty Data Structure to be Populated
data = struct();
data.N_modes = 3;   % Number of modes used to describe the system
data.N_cells = 2;
data.plot_grids = 1;

% initialize parameters to sweep over (b, maybe t)
% bpoints = [0.05:0.01:0.20]; %times pi
% bpoints = [.05 0.1 0.15 .2]*pi;
% betavals = [.00002 .0025  .005 .0075];
bpoints = [.05 .1 0.15 .2]*pi;
betavals = [.01 .05 .1];

%% Run Continuation to Get Stable Configurations at each b
% Choose which shape
shapeNum = 4;
data = init_shape(shapeNum, data);

%%
% Run COCO
[data,run_max_E_per_b,bpoints] = general_COCO(data, bpoints);
 
%% Run Optimization
% Get the current date and time
nowTime = datetime('now');

% Format the datetime as a string (e.g., '2025-06-08_14-30-15')
data.timeStr = string(datestr(nowTime, 'yyyy-mm-dd_HH-MM-SS'));


optimize(data, bpoints, betavals);