clear; clc
%% 
addpath("functions")
addpath("visualize")

%% Initialize 'data' for the specific problem
run('initialize_single_arch.m')
t_val = 0.01*pi;
data.t_vector = t_val*[1];
data.e_vector = 0*[1];

rng('default')
data.mu = 1;
data.sigma = 0;
data.variance = normrnd(data.mu,data.sigma,[1,data.N])';
data.variance = [1 1]';

%% Define the function as the arbitrary grid ODE
f = @(x,p) COCO_arbitrary_grid_ODE(x,p,data);

% Set up for COCO specific things
data.parameter_names = {'b' 't'};
data.initial_parameter_values = [0;t_val];

data.computational_domain = [-0.01 0.2*pi];
data.UZpoint = [0.1]*pi;

data.iterations_max = 5000;
data.h_min = 0.0005;
data.h_max = 0.001;

% Initial Guess for Continuation
Ahat0 = zeros(2*(data.N*(data.N_modes)-data.constraint_count),1);

%% Run the Initial Continuation Problem
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Ahat0,data.parameter_names,data.initial_parameter_values);
prob = coco_set(prob,'cont','ItMX', data.iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',data.h_max,'h_min',data.h_min);
    prob = coco_add_event(prob,'UZ','b',data.UZpoint);

coco(prob,'single_arch1',[],1,data.parameter_names,data.computational_domain)

%% BP1 
continue_from_BP('single_arch1',1,'single_arch1-1',data)

%% BP2
continue_from_BP('single_arch1',2,'single_arch1-2',data)

%% BP3
continue_from_BP('single_arch1',3,'single_arch1-3',data)

%% BP4
continue_from_BP('single_arch1',4,'single_arch1-4',data)

%% BP5
continue_from_BP('single_arch1',5,'single_arch1-5',data)

%% BP6
continue_from_BP('single_arch1',6,'single_arch1-6',data)

%% BP7
continue_from_BP('single_arch1',7,'single_arch1-7',data)

%% Plot all the Results COCO
idx1 = 1;
idx2 = 2;

figure(9899); clf; hold on
theme1 = struct('special', {{'HB','BP','EP'}});
coco_plot_bd(theme1, 'single_arch1'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'single_arch1-1'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'single_arch1-2'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'single_arch1-3'      ,'x',idx1,'x',idx2,'b')

view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.1])


%% Make the bifurcation diagram for the figure
[b1,x1stab,x1unst] = load_data_from_coco('single_arch1');
[b2,x2stab,x2unst] = load_data_from_coco('single_arch1-1');
[b3,x3stab,x3unst] = load_data_from_coco('single_arch1-2');
[b4,x4stab,x4unst] = load_data_from_coco('single_arch1-3');
[b5,x5stab,x5unst] = load_data_from_coco('single_arch1-4');
[b6,x6stab,x6unst] = load_data_from_coco('single_arch1-5');
[b7,x7stab,x7unst] = load_data_from_coco('single_arch1-6');
[b8,x8stab,x8unst] = load_data_from_coco('single_arch1-7');



f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = [2.5 1.5];
set(gca,"FontSize",10);
xlabel("b")
ylabel("$||u||$")

plot(b1,x1stab,"Color","k")
plot(b1,x1unst,"Color","r")

plot(b2,x2stab,"Color","k")
plot(b2,x2unst,"Color","r")

plot(b3,x3stab,"Color","k")
plot(b3,x3unst,"Color","r")

plot(b4,x4stab,"Color","k")
plot(b4,x4unst,"Color","r")

plot(b5,x5stab,"Color","k")
plot(b5,x5unst,"Color","r")

plot(b6,x6stab,"Color","k")
plot(b6,x6unst,"Color","r")

plot(b7,x7stab,"Color","k")
plot(b7,x7unst,"Color","r")

plot(b8,x8stab,"Color","k")
plot(b8,x8unst,"Color","r")