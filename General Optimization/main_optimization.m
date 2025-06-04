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
addpath('Shapes Point Data/')
%% Create Empty Data Structure to be Populated
data = struct();
data.N_modes = 3;   % Number of modes used to describe the system

% initialize parameters to sweep over (b, maybe t)
bpoints = [0.05:0.01:0.20]; %times pi

%% Run Continuation to Get Stable Configurations at each b
% Choose which shape
shapeNum = 1;
data = init_shape(shapeNum, data);
data = general_COCO(data, bpoints);

%% Determine High Energy State(s) at each b
%Assumes highest energy state goes to another ideal state
data = get_mode_shape_from_coco(data);
% Inside for loop for each b


%% Run time integration for each b
% Pattern periodic into long chain

% Time integration and mitigate edge effects

% Determine which arch to force/displace %MICHAEL QUESTION

% We have tiled state with edge effects, now need to compare transition


% for loop for different beta values


% Determine if a transition occurred and save info (boolean? or distance of wave?)
% Within each unit cell, calc potential energy, see which ones went from
% high to low (within say 10% and within 3 unit/super cells)





% data = rhombus_direct(data, bpoints);
% 
% data = get_mode_shape_from_coco(data);
% 
% for b=bpoints
%     data = system_elastic_from_chiral(data,b);
% end