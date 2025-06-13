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
addpath("Debugging\")
%% Create Empty Data Structure to be Populated
data = struct();
data.N_modes = 3;   % Number of modes used to describe the system
data.N_cells = 5;
data.plot_grids = 1;

% initialize parameters to sweep over (b, maybe t)
% bpoints = [0.05:0.01:0.20]; %times pi
% bpoints = [.05 0.1 0.15 .2]*pi;
% betavals = [.00002 .0025  .005 .0075];
%bpoints = [.08 .1 .12 .15 .2]*pi;
%bpoints = [.25];
%betavals = [.01 .05 .1];
%tvals = [.07 .08 .09 .1]*pi;

bpoints = [.25 .2*pi];

betavals = [.01 .1];
tvals = [.01 .1]*pi;

%% Run Continuation to Get Stable Configurations at each b
% Choose which shape
shapeNum = 4;
data = init_shape(shapeNum, data);

%%
% Get the current date and time
nowTime = datetime('now');

% Format the datetime as a string (e.g., '2025-06-08_14-30-15')
data.timeStr = string(datestr(nowTime, 'yyyy-mm-dd_HH-MM-SS'));
raw_data = deepCopyStruct(data);
trans_percent_tensor = zeros(length(bpoints),length(betavals), length(tvals));
for t_idx = 1:length(tvals)
    t = tvals(t_idx);
    % Run COCO
    data = raw_data;
    data.t = t;
    data.t_vector = t*ones(data.N,1);

    [data,run_max_E_per_b,bpoints] = general_COCO(data, bpoints);
    
    %% Run Optimization
    trans_percent_tensor(:,:,t_idx) = optimize(data, bpoints, betavals);
end

for beta_idx = 1:length(betavals)
    beta = betavals(beta_idx);
    results_cell = [ {'b/t'}, num2cell(tvals/pi); ...
             num2cell(bpoints/pi)', squeeze(num2cell(trans_percent_tensor(:,beta_idx,:)))];
    
    if ~exist("TransitionExcel\"+ data.timeStr + "\", 'dir')
        mkdir("TransitionExcel\"+data.timeStr + "\");
    end
    data.file_name_trans = data.timeStr + "\"+ data.shape_name + " beta = " + num2str(beta) + " NumCells = "+ num2str(data.N_cells);

    
    writecell(results_cell, "TransitionExcel\"+data.file_name_trans +".xlsx");
end