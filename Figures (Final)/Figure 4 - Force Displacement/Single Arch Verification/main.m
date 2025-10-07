%% =========================================================
%  Physical Constants
%  =========================================================
clear; clc
run('physical_constants.m')
run('initialize_single_arch.m')
data.variance = ([1])';

%% =========================================================
%  Run coco w/o variance
%  =========================================================
run('run_coco_no_variance.m')
% 
% % Stable State
% [Ahat,b] = COCO_grab_UZ('singlearch1-1',11);
% A = determine_A_from_Ahat(Ahat',data)';
% plot_system_once(A,data);
% save("Starting Shape.mat","Ahat","b")

b  = data.rise_ND;
A0 = sqrt(3*data.rise_ND^2 - data.t_ND^2) / sqrt(3);
A = zeros(2*data.N_modes,1);
A(1) = A0;
Ahat = determine_Ahat_from_A(A,data);
save("Starting Shape.mat","Ahat","b")

%% =========================================================
%  Run Time Integration for given eta value
%  =========================================================
data.beta = 0.05;
data.alpha_D = .5/1000;  % m/s displacement
data.interp_count = 3*50;

% data.eta  = 0.500;
% data.mms_to_run = 10/1000;
% run_time_integration(data,'Single Arch eta 3')

data.eta  = 0.575;
% data.eta  = 0.714;
data.mms_to_run = 3*10/1000;
run_time_integration(data,'Single Arch eta 5')

run('process_data.m')