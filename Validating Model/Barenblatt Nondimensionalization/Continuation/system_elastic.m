%% Clear Everything so there are no stragglers
clear; clc; close all

%% Add the Paths to the Required Functions
addpath('2D Array Functions')
addpath('Visualize')

%% Create Empty Data Structure to be Populated
data = struct();
data.N = 1;
data.N_modes = 5;   % Number of modes used to describe the system

data.b = 0.1*pi;
data.e = 0;
data.t = 0.01*pi;

data.points = [0 0 ; 1 0];
data.adjacency_matrix = [0 1 ; 1 0];


data.A0 = zeros(2*data.N_modes,1);
data.A0(1) = sqrt(12*(data.b/2)^2 - data.t^2)/sqrt(3);
data.A0(2) = 1e-2;

% Imposed displacement values
data.eta = 0.5;
data.alpha = 1.9432e-07;
% data.alpha = 0.001;


% Determine the coefficient matrix
data = determine_coefficient_matrix(data);

% Take a look at the initial condition
plot_system_once(data.A0,data)

% Determine Ahat from A
[A0hat] = determine_Ahat_from_A(data.A0,data);

data.A0hat =  A0hat;

%% Prepare for Time Integration
data.beta = 1;
data = determine_starting_vals(data);

%% Run Initial Time Integration
data.imposed = true;

T_end = 1.5/data.alpha;
tic
[t1,Ahat1] = ode45(@(t,A) single_arch_ODE(t,A,data),linspace(0,T_end,50000),data.A0hat);
toc

%% Recover A
A1 = determine_A_from_Ahat(t1,Ahat1',data)';

%% Determine the forces
Q = determine_force(A1,data);

first_positive = find(Q > 0,1);
Q = Q(first_positive:end);
t1 = t1(first_positive:end);
A1 = A1(first_positive:end,:);

last_positive = find(Q > 0,1,'last');

A0 = A1(last_positive+1,:)';
t0 = t1(last_positive+1);

%% Switch to Non Imposed Displacement
data.imposed = false;

tic
[t2,A2] = ode45(@(t,A) single_arch_ODE(t,A,data),[0 0.1*T_end],A0);
toc

%% Combine the Results
A = [A1(1:last_positive-1,:) ; A2];
t = [t1(1:last_positive-1) ; t2+t1(last_positive)];

%%
figure(1); clf; hold on
plot(t1,Q)
scatter(t1(last_positive),Q(last_positive),100,'filled')


figure(2); clf; hold on
plot(t,A(:,1:3))
plot(t1(last_positive)*[1 1],ylim,"k--","LineWidth",2)

%% Save the Results
file_string = "b = "+data.b/pi + "pi t = "+data.t/pi+"pi.mat";
save(file_string)