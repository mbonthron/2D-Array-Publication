function [data] = general_COCO(data, bpoints)
if shapeNum == 1
    %% === Define a path to the functions and clear screen
    
    %% === Load a points data for the system
    run('points_rhombus_direct')
    
    % Determine the adjacency matrix and number of arches
    data = determine_adjacency_matrix(data);
    run('rise_rhombus_constant')
    %% Start with elastic deformation
    [data] = initialize_elastic_deformation(zeros(data.N,1),zeros(data.V,1),data);
    
    % Take a look at the initial condition
    %plot_system_once(data.A0,data)
    %%
    data = add_periodicity(data);
    
    % Visualize the point and connection between the nodes
    plot_grid(data,1);
    
    % Determin the coefficient matrix and number of constraints of the system
    data = determine_coefficient_matrix(data);
    modes_to_skip = determine_modes_to_skip(data);
    
    %% General Stuff for Coco Continuation
    run_number = 1;
    run_name1 = ['rhombus_direct_run' sprintf('%.0f',run_number)];

end

N = data.N;
N_modes = data.N_modes;
constraint_count = 2*N-data.V; %double check line


A0 = zeros(2*(N*(N_modes)-constraint_count),1);

% Define the function as ode_triangle
f = @(x,p) arbitrary_grid_ODE_COCO(x,p,data,modes_to_skip);

parameter_names = {'b' 't'};            % Names our two parameters 'b' and 't'
initial_parameter_value = [0;0.01*pi];      % Starting values of b and t

active_continuation_parameter = 'b';    % Which parameter do we want to vary
computational_domain = [0 pi*0.2];         % What is the domain of b to explore
%UZpoints = [0.05*pi:0.01*pi:0.20*pi];                       % Values to explicitly call out to be used to plot
UZpoints = bpoints;
%% Continuation Constants
iterations_max = 5000;                  % Maximum number of iterations before continuation terminates
hmin = 0.0005;                           % Minimum step size of the continuation
hmax = 0.001;                            % Maximum step size of the continuation

%% ===============
%  INITIAL RUN FROM ZERO
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,A0,parameter_names,initial_parameter_value);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',hmax,'h_min',hmin);

coco(prob,run_name1,[],1,parameter_names,computational_domain)

clc
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
prob = coco_add_event(prob,'UZ','b',UZpoints);


fprintf("Run %.0f =========================================",1)
coco(prob,run_name1,[],1,parameter_names,computational_domain)

run_number = run_number + 1;

%% ===============
% Continue from UZ points
UZ = coco_bd_labs(run_name1, 'UZ'); 

for i = 1:2
    run_name = ['rhombus_direct_run' sprintf('%.0f',run_number)];
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
% run_name_start_from = ['rhombus_direct_run' sprintf('%.0f',9)];
% BP2 = coco_bd_labs(run_name_start_from, 'BP'); % labels for BP points in run1
% 
% for i = 1:length(BP2)
%     run_name = ['rhombus_direct_run' sprintf('%.0f',run_number)];
%     prob = coco_prob();
%     prob = ode_ep2ep(prob,'',run_name_start_from,BP2(i));
%     prob = coco_set(prob,'cont','branch','switch');
%     prob = coco_set(prob,'cont','ItMX', iterations_max);
%     prob = coco_set(prob,'cont','NPR',0);
%     prob = coco_set(prob,'cont','h_max',hmax,'h_min',hmin);
%     prob = coco_add_event(prob,'UZ','b',UZpoints);
%     coco(prob,run_name,[],1,parameter_names,computational_domain)  
%     run_number = run_number + 1;
% 
% end
% 
% run_name_start_from = ['rhombus_direct_run' sprintf('%.0f',14)];
% BP2 = coco_bd_labs(run_name_start_from, 'BP'); % labels for BP points in run1
% 
% for i = 1:length(BP2)
%     run_name = ['rhombus_direct_run' sprintf('%.0f',run_number)];
%     prob = coco_prob();
%     prob = ode_ep2ep(prob,'',run_name_start_from,BP2(i));
%     prob = coco_set(prob,'cont','branch','switch');
%     prob = coco_set(prob,'cont','ItMX', iterations_max);
%     prob = coco_set(prob,'cont','NPR',0);
%     prob = coco_set(prob,'cont','h_max',hmax,'h_min',hmin);
%     prob = coco_add_event(prob,'UZ','b',UZpoints);
%     coco(prob,run_name,[],1,parameter_names,computational_domain)  
%     run_number = run_number + 1;
% 
% end

%% Plot the Results from Coco
close all
cd('data');
rhombus_runs = dir('rhombus_direct*');
cd ..

theme1 = struct('special', {{'EP','FP','HB','BP'}});
figure(9899); clf; hold on

% 5 - metric
% 6 - Asymmetric
for i = 1:length(rhombus_runs)
% for i = [11]
    coco_plot_bd(theme1, rhombus_runs(i).name, 'x',1,'x',2,'b')
    if i ~= 1
        plot_shape_from_COCO(rhombus_runs(i).name,data,modes_to_skip)
    end
end


figure(9899)
axis tight; grid on; view(3)
zlim([0 0.5]); xlim([-1 1]); ylim(0.5*[-1 1])


%% Coco Plotting Trouble Shoot
theme1 = struct('special', {{'EP','FP','HB','BP'}});

figure(9900); clf; hold on; view(3); zlim([0 0.1]); xlim(0.3*[-1 1]); ylim(0.2*[-1 1])
grid()
coco_plot_bd(theme1, 'rhombus_direct_run1', 'x',1,'x',2,'b')
coco_plot_bd(theme1, 'rhombus_direct_run2', 'x',1,'x',2,'b')
coco_plot_bd(theme1, 'rhombus_direct_run3', 'x',1,'x',2,'b')
% coco_plot_bd(theme1, 'rhombus_direct_run4', 'x',1,'x',2,'b')
% coco_plot_bd(theme1, 'rhombus_direct_run5', 'x',1,'x',2,'b')
% coco_plot_bd(theme1, 'rhombus_direct_run6', 'x',1,'x',2,'b')
% coco_plot_bd(theme1, 'rhombus_direct_run7', 'x',1,'x',2,'b')
% coco_plot_bd(theme1, 'rhombus_direct_run8', 'x',1,'x',2,'b')




