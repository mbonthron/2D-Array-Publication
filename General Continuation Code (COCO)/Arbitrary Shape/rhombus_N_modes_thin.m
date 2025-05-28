%% === Define a path to the functions and clear screen
addpath('functions')
addpath('Visualize')
addpath('Shapes Point Data')
clear; clc; close all

%% Set the Number of Modes to Consider
N_modes = 5; % Number of modes to describe each triangle

%% === Load a points data for the system
run('points_rhombus')

% Determine the adjacency matrix and number of arches
adjacency_matrix = determine_adjacency_matrix(points);
N = sum(triu(adjacency_matrix,1) ==1,'all');

% Visualize the point and connection between the nodes
plot_grid(adjacency_matrix,points,1);

% Determin the coefficient matrix and number of constraints of the system
[coeff_matrix,constraint_count] = determine_coefficient_matrix(adjacency_matrix,0,N,N_modes);
modes_to_skip = determine_modes_to_skip(coeff_matrix,constraint_count,N,N_modes);

%% General Stuff for Coco Continuation
run_name1 = 'thin_rhombus_run1';

A0 = zeros(2*(N*(N_modes)-constraint_count),1);

% Define the function as ode_triangle
f = @(x,p) arbitrary_grid_ODE(x,p,coeff_matrix,N,N_modes,modes_to_skip);

parameter_names = {'b' 't'};            % Names our two parameters 'b' and 't'
initial_parameter_value = [0;0.01*pi];      % Starting values of b and t

active_continuation_parameter = 'b';    % Which parameter do we want to vary
computational_domain = [0 0.4];         % What is the domain of b to explore
UZpoints = [0.1];                       % Values to explicitly call out to be used to plot

%% Continuation Constants
iterations_max = 5000;                  % Maximum number of iterations before continuation terminates
hmin = 0.001;                           % Minimum step size of the continuation
hmax = 0.01;                            % Maximum step size of the continuation

%% ===============
%  INITIAL RUN FROM ZERO
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,A0,parameter_names,initial_parameter_value);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',hmax,'h_min',hmin);
prob = coco_add_event(prob,'UZ','b',UZpoints);

coco(prob,run_name1,[],1,parameter_names,computational_domain)

% Load all the HBs as UZR Points and Run Again
bd = coco_bd_read(run_name1);
HBlbls = coco_bd_labs(run_name1, 'HB');

bcrits = zeros(1,length(HBlbls));
acrits = zeros(2*(N*(N_modes)-constraint_count),length(HBlbls));

for k = 1:length(HBlbls)
    bcrits(k) = coco_bd_val(bd,HBlbls(k),'b');
    acrits(:,k) = coco_bd_val(bd,HBlbls(k),'x');
end

prob = coco_add_event(prob,'UZ','b',bcrits);

coco(prob,run_name1,[],1,parameter_names,computational_domain)

%% ===============
% Continue from the first branch point
BP = coco_bd_labs(run_name1, 'BP'); % labels for BP points in run1

for i = 1:length(BP)
    run_name = ['thin_rhombus_run' sprintf('%.0f',1+i)];
    prob = coco_prob();
    prob = ode_ep2ep(prob,'',run_name1,BP(i));
    prob = coco_set(prob,'cont','branch','switch');
    prob = coco_set(prob,'cont','ItMX', iterations_max);
    prob = coco_set(prob,'cont','NPR',0);
    prob = coco_set(prob,'cont','h_max',hmax,'h_min',hmin);
    prob = coco_add_event(prob,'UZ','b',UZpoints);
    coco(prob,run_name,[],1,parameter_names,computational_domain)  
end


%% ===============
% Continue from UZ points
UZ = coco_bd_labs(run_name1, 'UZ'); 

for i = 1:length(UZ)
    run_name = ['thin_rhombus_run' sprintf('%.0f',1+i+length(BP))];
    prob = coco_prob();
    prob = ode_ep2ep(prob,'',run_name1,UZ(i));
    prob = coco_set(prob,'cont','branch','switch');
    prob = coco_set(prob,'cont','ItMX', iterations_max);
    prob = coco_set(prob,'cont','NPR',0);
    prob = coco_set(prob,'cont','h_max',hmax,'h_min',hmin);
    prob = coco_add_event(prob,'UZ','b',UZpoints);

    coco(prob,run_name,[],1,parameter_names,computational_domain)
end


%% Plot the Results from Coco
close all;
cd('data');
rhombus_runs = dir('thin_rhombus*');
cd ..

theme1 = struct('special', {{'EP','FP','HB','BP'}});
figure(9899); clf; hold on

for i = 1:length(rhombus_runs)
% for i = [1 2 3]
    coco_plot_bd(theme1, rhombus_runs(i).name, 'x',1,'x',2,'b')
    if i ~= 1
        plot_shape_from_COCO(rhombus_runs(i).name,N,N_modes,adjacency_matrix,points,coeff_matrix,modes_to_skip)
    end
end

%
figure(9899)
axis tight; grid on; view(3)
zlim([0 0.15]); xlim(0.125*[-1 1]); ylim(0.05*[-1 1])

