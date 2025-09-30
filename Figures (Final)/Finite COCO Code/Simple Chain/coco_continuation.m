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
data.UZpoint = [0.025 0.050 0.075 0.1]*pi;

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

%% BP3
continue_from_BP('simple_chain1',3,'simple_chain1-3',data)

%% BP3
continue_from_BP('simple_chain1',4,'simple_chain1-4',data)

%% BP3
continue_from_BP('simple_chain1',5,'simple_chain1-5',data)

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
[a1stab1,a1unst1,a3stab1,a3unst1,b1] = get_data_from_coco('simple_chain1');
[a1stab2,a1unst2,a3stab2,a3unst2,b2] = get_data_from_coco('simple_chain1-1');
[a1stab3,a1unst3,a3stab3,a3unst3,b3] = get_data_from_coco('simple_chain1-2');


deltaL1 = (b1 / 2).^2;
deltaL2 = (b2 / 2).^2;
deltaL3 = (b3 / 2).^2;
%%

L_dimensional = 100;

f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = [5/3 1.05];
% f.Position(3:4) = [1.25 1.1229];
plot(deltaL1*100/pi,a1stab1*100/pi,"k-","DisplayName","Stable")
plot(deltaL1*100/pi,a1unst1*100/pi,"r-","DisplayName","Unstable")

plot(deltaL2*100/pi,a1stab2*100/pi,"k-","HandleVisibility","off")
plot(deltaL2*100/pi,a1unst2*100/pi,"r-","HandleVisibility","off")

plot(deltaL3*100/pi,a1stab3*100/pi,"k-","HandleVisibility","off")
plot(deltaL3*100/pi,a1unst3*100/pi,"r-","HandleVisibility","off")

x1 = 0;
x2 = 0.05;
y1 = -2.457;
y2 = 2.457;

plot(4*[x1 x2 x2 x1 x1],3*[y1 y1 y2 y2 y1],"k-")

% xticks([])
% yticks([])

% legend()

xlabel("$\Delta L$ - [mm]")
ylabel("$b$ - [mm]")
set(gca,'FontSize',10)