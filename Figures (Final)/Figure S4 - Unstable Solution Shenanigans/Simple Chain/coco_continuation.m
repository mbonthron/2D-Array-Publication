clear; clc
%% 
addpath("functions")
addpath("visualize")

%% Initialize 'data' for the specific problem
run('initialize_simple_chain.m')
t_val = 0.01*pi;
data.t_vector = t_val*[1 1];
data.e_vector = 0*[1 1];

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

coco(prob,'simple_chain1',[],1,data.parameter_names,data.computational_domain)

%% BP1 
continue_from_BP('simple_chain1',1,'simple_chain1-1',data)

%% BP2
continue_from_BP('simple_chain1',2,'simple_chain1-2',data)

continue_from_BP('simple_chain1-2',1,'simple_chain1-2-1',data)
continue_from_BP('simple_chain1-2',3,'simple_chain1-2-3',data)


%% BP3
continue_from_BP('simple_chain1',3,'simple_chain1-3',data)

%% BP3
continue_from_BP('simple_chain1',4,'simple_chain1-4',data)

continue_from_BP('simple_chain1-4',1,'simple_chain1-4-1',data)
continue_from_BP('simple_chain1-4',3,'simple_chain1-4-3',data)


%% BP3
continue_from_BP('simple_chain1',5,'simple_chain1-5',data)

continue_from_BP('simple_chain1-5',1,'simple_chain1-5-1',data)
continue_from_BP('simple_chain1-5',3,'simple_chain1-5-3',data)

%% More BP
continue_from_BP('simple_chain1',6,'simple_chain1-6',data)

continue_from_BP('simple_chain1',7,'simple_chain1-7',data)

continue_from_BP('simple_chain1-7',1,'simple_chain1-7-1',data)
continue_from_BP('simple_chain1-7',3,'simple_chain1-7-3',data)

continue_from_BP('simple_chain1',8,'simple_chain1-8',data)

continue_from_BP('simple_chain1',9,'simple_chain1-9',data)
continue_from_BP('simple_chain1-9',1,'simple_chain1-9-1',data)
continue_from_BP('simple_chain1-9',3,'simple_chain1-9-3',data)

continue_from_BP('simple_chain1',10,'simple_chain1-10',data)

continue_from_BP('simple_chain1',11,'simple_chain1-11',data)
continue_from_BP('simple_chain1-11',1,'simple_chain1-11-1',data)
continue_from_BP('simple_chain1-11',3,'simple_chain1-11-3',data)

continue_from_BP('simple_chain1',12,'simple_chain1-12',data)

continue_from_BP('simple_chain1',13,'simple_chain1-13',data)
continue_from_BP('simple_chain1-13',1,'simple_chain1-13-1',data)
continue_from_BP('simple_chain1-13',3,'simple_chain1-13-3',data)



%% Plot all the Results COCO
idx1 = 1;
idx2 = 2;

figure(9899); clf; hold on
theme1 = struct('special', {{'HB','BP','EP'}});
coco_plot_bd(theme1, 'simple_chain1'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'simple_chain1-1'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'simple_chain1-2'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'simple_chain1-3'      ,'x',idx1,'x',idx2,'b')

view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.1])

%% Load the Data
[b0,x0stab,x0unst] = load_data_from_coco('simple_chain1');
[b1,x1stab,x1unst] = load_data_from_coco('simple_chain1-1');
[b2,x2stab,x2unst] = load_data_from_coco('simple_chain1-2');
[b3,x3stab,x3unst] = load_data_from_coco('simple_chain1-3');
[b4,x4stab,x4unst] = load_data_from_coco('simple_chain1-4');
[b5,x5stab,x5unst] = load_data_from_coco('simple_chain1-5');
[b6,x6stab,x6unst] = load_data_from_coco('simple_chain1-6');
[b7,x7stab,x7unst] = load_data_from_coco('simple_chain1-7');
[b8,x8stab,x8unst] = load_data_from_coco('simple_chain1-8');
[b9,x9stab,x9unst] = load_data_from_coco('simple_chain1-9');
[b10,x10stab,x10unst] = load_data_from_coco('simple_chain1-10');
[b11,x11stab,x11unst] = load_data_from_coco('simple_chain1-11');
[b12,x12stab,x12unst] = load_data_from_coco('simple_chain1-12');
[b13,x13stab,x13unst] = load_data_from_coco('simple_chain1-13');



f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = [2.5 1.5];
set(gca,"FontSize",10);
xlabel("b")
ylabel("$||u||$")

plot(b0,x0stab,"Color","k")
plot(b0,x0unst,"Color","r")

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

plot(b6,x6stab,"Color","k"); plot(b6,x6unst,"Color","r")
plot(b7,x7stab,"Color","k"); plot(b7,x7unst,"Color","r")
plot(b8,x8stab,"Color","k"); plot(b8,x8unst,"Color","r")
plot(b9,x9stab,"Color","k"); plot(b9,x9unst,"Color","r")
plot(b10,x10stab,"Color","k"); plot(b10,x10unst,"Color","r")
plot(b11,x11stab,"Color","k"); plot(b11,x11unst,"Color","r")
plot(b12,x12stab,"Color","k"); plot(b12,x12unst,"Color","r")
%%
plot_shape_from_COCO('simple_chain1',data);
plot_shape_from_COCO('simple_chain1-1',data);
plot_shape_from_COCO('simple_chain1-2',data);
plot_shape_from_COCO('simple_chain1-3',data);
plot_shape_from_COCO('simple_chain1-4',data);
plot_shape_from_COCO('simple_chain1-5',data);
plot_shape_from_COCO('simple_chain1-6',data);
plot_shape_from_COCO('simple_chain1-7',data);
plot_shape_from_COCO('simple_chain1-8',data);
plot_shape_from_COCO('simple_chain1-9',data);
plot_shape_from_COCO('simple_chain1-10',data);
plot_shape_from_COCO('simple_chain1-11',data);
plot_shape_from_COCO('simple_chain1-12',data);
plot_shape_from_COCO('simple_chain1-13',data);

%%
[b2_1,x2_1stab,x2_1unst] = load_data_from_coco('simple_chain1-2-1');
[b2_3,x2_3stab,x2_3unst] = load_data_from_coco('simple_chain1-2-3');

[b4_1,x4_1stab,x4_1unst] = load_data_from_coco('simple_chain1-4-1');
[b4_3,x4_3stab,x4_3unst] = load_data_from_coco('simple_chain1-4-3');

[b5_1,x5_1stab,x5_1unst] = load_data_from_coco('simple_chain1-5-1');
[b5_3,x5_3stab,x5_3unst] = load_data_from_coco('simple_chain1-5-3');

[b7_1,x7_1stab,x7_1unst] = load_data_from_coco('simple_chain1-7-1');
[b7_3,x7_3stab,x7_3unst] = load_data_from_coco('simple_chain1-7-3');

[b9_1,x9_1stab,x9_1unst] = load_data_from_coco('simple_chain1-9-1');
[b9_3,x9_3stab,x9_3unst] = load_data_from_coco('simple_chain1-9-3');

[b11_1,x11_1stab,x11_1unst] = load_data_from_coco('simple_chain1-11-1');
[b11_3,x11_3stab,x11_3unst] = load_data_from_coco('simple_chain1-11-3');

[b13_1,x13_1stab,x13_1unst] = load_data_from_coco('simple_chain1-13-1');
[b13_3,x13_3stab,x13_3unst] = load_data_from_coco('simple_chain1-13-3');


%%
f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = [2.5 1.5];
set(gca,"FontSize",10);
xlabel("b")
ylabel("$||u||$")

plot(b0,x0stab,"Color","k"); plot(b0,x0unst,"Color","r")
plot(b1,x1stab,"Color","k"); plot(b1,x1unst,"Color","r")
plot(b2,x2stab,"Color","k"); plot(b2,x2unst,"Color","r")
plot(b3,x3stab,"Color","k"); plot(b3,x3unst,"Color","r")
plot(b4,x4stab,"Color","k"); plot(b4,x4unst,"Color","r")
plot(b5,x5stab,"Color","k"); plot(b5,x5unst,"Color","r")
plot(b6,x6stab,"Color","k"); plot(b6,x6unst,"Color","r")
plot(b7,x7stab,"Color","k"); plot(b7,x7unst,"Color","r")
plot(b8,x8stab,"Color","k"); plot(b8,x8unst,"Color","r")
plot(b9,x9stab,"Color","k"); plot(b9,x9unst,"Color","r")
plot(b10,x10stab,"Color","k"); plot(b10,x10unst,"Color","r")
plot(b11,x11stab,"Color","k"); plot(b11,x11unst,"Color","r")
plot(b12,x12stab,"Color","k"); plot(b12,x12unst,"Color","r")

plot(b2_1,x2_1stab,"Color","k"); plot(b2_1,x2_1unst,"Color","r")
plot(b2_3,x2_3stab,"Color","k"); plot(b2_3,x2_3unst,"Color","r")

plot(b4_1,x4_1stab,"Color","k"); plot(b4_1,x4_1unst,"Color","r")
plot(b4_3,x4_3stab,"Color","k"); plot(b4_3,x4_3unst,"Color","r")

plot(b5_1,x5_1stab,"Color","k"); plot(b5_1,x5_1unst,"Color","r")
plot(b5_3,x5_3stab,"Color","k"); plot(b5_3,x5_3unst,"Color","r")

plot(b7_1,x7_1stab,"Color","k"); plot(b7_1,x7_1unst,"Color","r")
plot(b7_3,x7_3stab,"Color","k"); plot(b7_3,x7_3unst,"Color","r")

plot(b9_1,x9_1stab,"Color","k"); plot(b9_1,x9_1unst,"Color","r")
plot(b9_3,x9_3stab,"Color","k"); plot(b9_3,x9_3unst,"Color","r")

plot(b11_1,x11_1stab,"Color","k"); plot(b11_1,x11_1unst,"Color","r")
plot(b11_3,x11_3stab,"Color","k"); plot(b11_3,x11_3unst,"Color","r")

plot(b13_1,x13_1stab,"Color","k"); plot(b13_1,x13_1unst,"Color","r")
plot(b13_3,x13_3stab,"Color","k"); plot(b13_3,x13_3unst,"Color","r")