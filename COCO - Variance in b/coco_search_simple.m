%% Clear Everything so there are no stragglers
restoredefaultpath;
startup
addpath('..\General Time Integration Code (MATLAB)\Visualize')
addpath('..\General Time Integration Code (MATLAB)\2D Array Functions')
% addpath('..\General Continuation Code (COCO)\Arbitary Shape\functions\')
% addpath('..\General Continuation Code (COCO)\Arbitary Shape\Visualize\')
% addpath('..\General Continuation Code (COCO)\Arbitary Shape\')
addpath('Shapes Point Data/')

clear; clc; close all
%% Create Empty Data Structure to be Populated
data = struct();
data.N_modes = 4;   % Number of modes used to describe the system
data.N_cells = 2;
data.plot_grids = 1;

%% Initialize the shape and the 'data' structure
% Choose which shape
% shapeNum = 1;   % Rhombus
% shapeNum = 2;   % Double Triangle 
% shapeNum = 3;   % Alternating Triangle 
% shapeNum = 4;   % Hexagon
% shapeNum = 5;   % Diamond Chain
% shapeNum = 6;   % Diagonal Rhombus
% shapeNum = 7;   % Diagonal Rhombus w/ Arch

% data = init_shape(shapeNum, data);

% run('super_simple_triangle.m')
run('super_simple_hexagon.m')
%%
data.t = 0.01*pi;
data.t_vector = data.t*ones(data.N,1);

bpoints = 0.15*pi;


%% Add variance
data.mu = 1;
data.sigma = 0.025;

% rng('default')
% data.variance = normrnd(data.mu,data.sigma,[1,data.N])';
data.variance = ones(data.N,1);
% data.variance = [1.05 ; 1 ; 1];

%% Run COCO with variance
variance_COCO(data,bpoints)