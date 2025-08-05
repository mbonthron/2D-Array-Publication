%% Clear Everything so there are no stragglers
clear; clc; close all

%% Add the Paths to the Required Functions
addpath('2D Array Functions')
addpath('Visualize')

%%
data = struct();
data.N = 1;
data.N_modes = 5;

data.b = 0.1*pi;
data.e = 0;
data.t = 0.01*pi;

data.points = [0 0 ; 1 0];
data.adjacency_matrix = [0 1 ; 1 0];

data.A0 = zeros(2*data.N_modes,1);
data.A0(1) = sqrt(12*(data.b/2)^2 - data.t^2)/sqrt(3);

data.eta = 0.50;
data.alpha = 1.9432e-07;

% Determine the coefficient matrix
data = determine_coefficient_matrix(data);

% Determine Ahat from A
[A0hat] = determine_Ahat_from_A(data.A0,data);

data.A0hat =  A0hat;

%% Prepare for Time Integration
data.beta = 1;
data = determine_starting_vals(data);

%% Start with coco
data.imposed = true;
f = @(x,p) single_arch_ODE(p,x,data);

toolbox_family = 'ode';                   %
initial_point = 'isol';                   %
branch_type   = 'ep';                     %

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', false);
prob = ode_isol2ep(prob,'',f,A0hat,{'t'},[0]);


prob = coco_set(prob,'cont','ItMX', 50000);
prob = coco_set(prob,'cont','NPR',1000);
prob = coco_set(prob,'cont','h_max',1000,'h_min',100);

coco(prob,'run1lowb',[],1,{'t'},[-0.1 1.5/data.alpha])

%%
bd = coco_bd_read('run1lowb');


col_idx = find(strcmp(bd(1,:), 't'));
t = cell2mat(bd(2:end,col_idx));

% Find the column index with header 'x'
col_idx = find(strcmp(bd(1,:), 'x'));
data_cells = bd(2:end, col_idx);
Ahat = reshape(cell2mat(cellfun(@(x) x(:), data_cells, 'UniformOutput', false)),8,[])';

% Recover A
A = determine_A_from_Ahat(t,Ahat',data)';

%
Q = determine_force(A,data);

%%
figure(1); clf; hold on
thm = struct('special', {{'SN' 'BP' 'EP'}});
coco_plot_bd(thm, 'run1lowb', 't', 'x',2)

coeff = sin(pi*data.eta*(1:data.N_modes));

figure(2); clf; hold on
plot(data.alpha*t,sum(coeff.*A(:,1:5),2))

figure(3); clf; hold on
plot(data.alpha*t,Q)
