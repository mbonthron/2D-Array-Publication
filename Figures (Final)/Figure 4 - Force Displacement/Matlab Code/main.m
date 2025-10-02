%% =========================================================
%  Physical Constants
%  =========================================================
clear; clc
run('physical_constants.m')
data.variance = ([8.7 8.8 8.7 9.2 8.9]/8.7)';

%% =========================================================
%  Run coco w/o variance
%  =========================================================
run('run_coco_no_variance.m')

% Save the two 'S' shapes to be used as initial guesses in the variance run
% Negative a2
[Ahat,b] = COCO_grab_UZ('dogbone1-2',6);
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data);
save("Frustrated State Negative a2.mat","Ahat","b")

% Positive a2
[Ahat,b] = COCO_grab_UZ('dogbone1-2',15);
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data);
save("Frustrated State Positive a2.mat","Ahat","b")

%% =========================================================
%  Run coco w/o variance
%  =========================================================
run('run_coco_variance.m')

% Save the positive 'C' shape w/ variance as starting point for time
% integration
[Ahat,b] = COCO_grab_UZ('vdogbone1-1',4);
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data);
save("C State.mat","Ahat","b")

% Save the positive 'S' shape w/ variance as starting point for time
% integration
[Ahat,b] = COCO_grab_UZ('vdogbone1-8',2);
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data);
save("S State.mat","Ahat","b")


%% =========================================================
%  Run Time Integration for given eta value
%  =========================================================
data.beta = 0.1;
data.alpha_D = .1/1000;  % m/s displacement

data.eta  = 0.50;
data.mms_to_run = 15/1000;
run_time_integration('C',data,'C eta 50')

data.eta  = 0.633;
data.mms_to_run = 15/1000;
run_time_integration('C',data,'C eta 75')

data.eta  = 0.362;
% data.eta  = 0.232;
data.mms_to_run = 5/1000;
run_time_integration('S',data,'S eta 25')

data.eta  = 0.5;
% data.eta  = 0.464;
data.mms_to_run = 10/1000;
run_time_integration('S',data,'S eta 50')

data.eta  = 0.635;
% data.eta  = 0.8;
data.mms_to_run = 15/1000;
run_time_integration('S',data,'S eta 75')