%% === Define a path to the functions and clear screen
addpath('functions')
addpath('Visualize')
addpath('Shapes Point Data')
clear; clc; close all

%% Set the Number of Modes to Consider
data = struct();
data.N_modes = 3; % Number of modes to describe each triangle

%% === Load a points data for the system
run('points_chain_direct')
run('points_chain_direct_finite')

% Determine the adjacency matrix and number of arches
data = determine_adjacency_matrix(data);
data = remove_connection(data);
data = add_periodicity(data);

% Find the periodic to constant 
data.vertex_map_p2f = [1 13];
data = determine_per_to_finite(data);

% Visualize the point and connection between the nodes
plot_grid(data,true);

% Determin the coefficient matrix and number of constraints of the system
data = determine_coefficient_matrix(data);
data = determine_modes_to_skip(data);

%% General Stuff for Coco Continuation
run_number = 1;
run_name1 = ['coco_run_' sprintf('%.0f',run_number)];


Ahat0 = zeros(2*(data.N*(data.N_modes)-data.constraint_count),1);

% Define the function as ode_triangle
f = @(x,p) COCO_arbitrary_grid_ODE(x,p,data);

parameter_names = {'b' 't'};            % Names our two parameters 'b' and 't'
initial_parameter_value = [0;0.01*pi];      % Starting values of b and t

active_continuation_parameter = 'b';    % Which parameter do we want to vary
computational_domain = [0 0.7];         % What is the domain of b to explore
% UZpoints = [0.05 0.10 0.15 0.20]*pi;                       % Values to explicitly call out to be used to plot
UZpoints = [0.20]*pi;                       % Values to explicitly call out to be used to plot


%% Continuation Constants
iterations_max = 5000;                  % Maximum number of iterations before continuation terminates
hmin = 0.0005;                           % Minimum step size of the continuation
hmax = 0.001;                            % Maximum step size of the continuation

%% ===============
%  INITIAL RUN FROM ZERO
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Ahat0,parameter_names,initial_parameter_value);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',hmax,'h_min',hmin);

coco(prob,run_name1,[],1,parameter_names,computational_domain)

clc
% Load all the HBs as UZR Points and Run Again
bd = coco_bd_read(run_name1);
HBlbls = coco_bd_labs(run_name1, 'HB');

bcrits = zeros(1,length(HBlbls));
acrits = zeros(2*(data.N*(data.N_modes)-data.constraint_count),length(HBlbls));

for k = 1:length(HBlbls)
    bcrits(k) = coco_bd_val(bd,HBlbls(k),'b');
    acrits(:,k) = coco_bd_val(bd,HBlbls(k),'x');
end

prob = coco_add_event(prob,'UZ','b',bcrits);
prob = coco_add_event(prob,'UZ','b',UZpoints);


fprintf("Run %.0f =========================================",1)
coco(prob,run_name1,[],1,parameter_names,computational_domain)

run_number = run_number + 1;

%% ===============
% Continue from UZ points
UZ = coco_bd_labs(run_name1, 'UZ'); 

for i = 1:2
    run_name = ['coco_run_' sprintf('%.0f',run_number)];
    prob = coco_prob();
    prob = ode_ep2ep(prob,'',run_name1,UZ(i));
    prob = coco_set(prob,'cont','branch','switch');
    prob = coco_set(prob,'cont','ItMX', iterations_max);
    prob = coco_set(prob,'cont','NPR',0);
    prob = coco_set(prob,'cont','h_max',hmax,'h_min',hmin);
    prob = coco_add_event(prob,'UZ','b',UZpoints);

    fprintf("Run %.0f =========================================",run_number)
    coco(prob,run_name,[],1,parameter_names,computational_domain)
    run_number = run_number + 1;
end


%% ===============
% Continue from the branch points of the second run
run_name_start_from = ['coco_run_' sprintf('%.0f',1)];
BP2 = coco_bd_labs(run_name_start_from, 'BP'); % labels for BP points in run1

for i = 1:2
    run_name = ['coco_run_' sprintf('%.0f',run_number)];
    prob = coco_prob();
    prob = ode_ep2ep(prob,'',run_name_start_from,BP2(i));
    prob = coco_set(prob,'cont','branch','switch');
    prob = coco_set(prob,'cont','ItMX', iterations_max);
    prob = coco_set(prob,'cont','NPR',0);
    prob = coco_set(prob,'cont','h_max',hmax,'h_min',hmin);
    prob = coco_add_event(prob,'UZ','b',UZpoints);
    coco(prob,run_name,[],1,parameter_names,computational_domain)  
    run_number = run_number + 1;
end


%% Plot the Results from Coco
close all
cd('data');
rhombus_runs = dir('coco_run*');
cd ..

theme1 = struct('special', {{'EP','FP','HB','BP'}});
figure(9899); clf; hold on

for i = 1:length(rhombus_runs)
    coco_plot_bd(theme1, rhombus_runs(i).name, 'x',1,'x',2,'b')
    if i ~= 1
        plot_shape_from_COCO(rhombus_runs(i).name,data)
    end
end


figure(9899)
axis tight; grid on; view(3)
zlim([0 0.5]); xlim([-1 1]); ylim(0.5*[-1 1])


%% Coco Plotting Trouble Shoot
theme1 = struct('special', {{'EP','FP','HB','BP'}});

figure(9900); clf; hold on; view(3); zlim([0 0.1]); xlim(0.3*[-1 1]); ylim(0.2*[-1 1])
grid()
coco_plot_bd(theme1, 'coco_run_1', 'x',1,'x',2,'b')
coco_plot_bd(theme1, 'coco_run_2', 'x',1,'x',2,'b')
coco_plot_bd(theme1, 'coco_run_3', 'x',1,'x',2,'b')