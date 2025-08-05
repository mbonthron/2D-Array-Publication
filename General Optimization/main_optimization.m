%% Clear Everything so there are no stragglers
clear; clc; close all

%profile -memory on;
%setpref('profiler','showJitLines',1);
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
%bpoints = [.07 .08 .09 .1 .11 .12 .13 .14 .15 .175 .2 .225 .25 .275 .3 .325 .35]*pi;
%bpoints = [.25];

% bpoints = [.070 .075 .080 .085 .090 .095 .10 .105 .110 .115 .120 .125 .13]*pi;
% betavals = [.1:.1:1 2:1:5];
% tvals = [.01 .025 .05 .06 .065 .070 .075 .08 .085 .09 .095 .1 .105]*pi;

bpoints = [.09]*pi;
betavals = [.1];
tvals = [.1]*pi;
sigmavals = [.0001];

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
tensor_true = 1;
OG_bpoints = bpoints;
%% Run COCO
for t_idx = 1:length(tvals)
    t = tvals(t_idx);
    % Run COCO
    data = raw_data;
    data.t = t;
    data.t_vector = t*ones(data.N,1);
    bpoints = OG_bpoints;
    
    runs_exist = 1;
    troubleshooting_flag = 1;
    for b = bpoints
        if ~isfile(data.shape_name + " b = "+ num2str(round(b/pi,4)) +" pi t = "+num2str(round(t/pi,4)) +" pi.mat")
            runs_exist = 0;
        end
    end
    if ~runs_exist || troubleshooting_flag
        [data,run_max_E_per_b,bpoints] = general_COCO(data, bpoints,-1); %close all;
    
    end
end

%% Run Optimization
myCluster = parcluster('Processes');
delete(myCluster.Jobs);
myCluster.NumWorkers = 4;
saveProfile(myCluster);
for t_idx = 1:length(tvals)
    t = tvals(t_idx);
    % Run COCO
    data = raw_data;
    data.t = t;
    data.t_vector = t*ones(data.N,1);
    bpoints = OG_bpoints;
    %% Run Optimization
    if ~exist("TransitionExcel\"+ data.timeStr + "\", 'dir')
            mkdir("TransitionExcel\"+data.timeStr + "\");
    end
    trans_percent_tensor(:,:,t_idx) = optimize(data, bpoints, betavals, sigmavals);
    close all hidden;
    % if tensor_true % Check for errors if bpoints change size in COCO
    %     if length(bpoints) == length(OG_bpoints)
    %         trans_percent_tensor(:,:,t_idx) = optimize(data, bpoints, betavals);
    %         save("temp", "trans_percent_tensor");
    %     else
    %         if any(trans_percent_tensor) % Check if data has been written
    %             for i=1:(t_idx-1)
    %                 trans_percent_cell{i,1} = trans_percent_tensor(:,:,i);
    %                 trans_percent_cell{i,2} = OG_bpoints;
    %             end
    %         end
    %         trans_percent_cell{t_idx,1} = optimize(data, bpoints, betavals);
    %         trans_percent_cell{t_idx,2} = bpoints;
    %         tensor_true = 0;
    %     end
    % else
    %     trans_percent_cell{t_idx,1} = optimize(data, bpoints, betavals);
    %     trans_percent_cell{t_idx,2} = bpoints;
    %     save("temp", "trans_percent_cell");
    % end
end

data = raw_data;
bpoints = OG_bpoints;

if tensor_true
    for beta_idx = 1:length(betavals)
        beta = betavals(beta_idx);
        results_cell = [ {'b/t'}, num2cell(tvals/pi); ...
                 num2cell(bpoints/pi)', squeeze(num2cell(trans_percent_tensor(:,beta_idx,:)))];

        if ~exist("TransitionExcel\"+ data.timeStr + "\", 'dir')
            mkdir("TransitionExcel\"+data.timeStr + "\");
        end
        data.file_name_trans = data.timeStr + "\"+ data.shape_name + " beta = " + num2str(beta) + " NumCells = "+ num2str(data.N_cells);


        writecell(results_cell, "TransitionExcel\"+data.file_name_trans +".xlsx");

        figure(77); clf; hold on;
        title("Transition Percent " + data.shape_name + " beta = " + num2str(beta) + " NumCells = "+ num2str(data.N_cells));
        xlabel("t");
        ylabel("b");
        clim([0 100]);
        a=colorbar;
        a.Label.String = 'Transition Percent %';
        imagesc(tvals/pi, bpoints/pi, squeeze(trans_percent_tensor(:,beta_idx,:)))
        print(gcf, "TransitionExcel\"+data.file_name_trans + " - Visual.png", '-dpng', '-r600');
    end
else
    for t_idx = 1:length(tvals)
        t = tvals(t_idx);
        bvals = trans_percent_cell{t_idx,2};
        results_cell = [ {'b/beta'}, num2cell(betavals); ...
                 num2cell(bvals/pi)', squeeze(num2cell(trans_percent_cell{t_idx,1}))];

        if ~exist("TransitionExcel\"+ data.timeStr + "\", 'dir')
            mkdir("TransitionExcel\"+data.timeStr + "\");
        end
        data.file_name_trans = data.timeStr + "\"+ data.shape_name + " t = " + num2str(t) + " NumCells = "+ num2str(data.N_cells);


        writecell(results_cell, "TransitionExcel\"+data.file_name_trans +".xlsx");
    end
end
