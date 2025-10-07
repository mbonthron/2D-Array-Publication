%% =========================================================
%  Physical Constants
%  =========================================================
clear; clc
run('physical_constants.m')

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
% data.variance = ([8.7 8.8 8.7 9.2 8.9]/8.7)';
data.variance = ([10 11 11 11 11 ]/10)';

run('run_coco_variance.m')

% Save the positive 'C' shape w/ variance as starting point for time
% integration
[Ahat,b] = COCO_grab_UZ('vdogbone1-1',11);
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data);
if A(1)<0
    [Ahat,b] = COCO_grab_UZ('vdogbone1-1',3);
    A = determine_A_from_Ahat(Ahat',data)';
    plot_system_once(A,data);
end

save("C State.mat","Ahat","b")

% Save the positive 'S' shape w/ variance as starting point for time
% integration
[Ahat,b] = COCO_grab_UZ('vdogbone1-9',3);
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data);
save("S State.mat","Ahat","b")


%% =========================================================
%  Run Time Integration for given eta value
%  =========================================================
% data.beta = 0.05;
% data.alpha_D = .5/1000;  % m/s displacement
% data.interp_count = 50;
% data.mms_to_run = 10/1000;
% eta_vector = 0.1:0.05:0.9;
% 
% for i = 1:length(eta_vector)
%     data.eta = eta_vector(i);
%     % run_name = 'C eta ' + string(i);
%     % run_time_integration('C',data,run_name)
% 
%     run_name = 'S eta ' + string(i);
%     run_time_integration('S',data,run_name)
% 
% end


%%
data.beta = 0.05;
data.alpha_D = .5/1000;  % m/s displacement
data.interp_count = 50;

data.eta  = 0.259;
data.mms_to_run = 10/1000;
run_time_integration('C',data,'C eta 1')

data.eta  = 0.354;
data.mms_to_run = 10/1000;
run_time_integration('C',data,'C eta 2')

data.eta  = 0.5;
data.mms_to_run = 10/1000;
run_time_integration('C',data,'C eta 3')

data.eta  = 0.52;
data.mms_to_run = 10/1000;
run_time_integration('C',data,'C eta 4')

data.eta  = 0.701;
data.mms_to_run = 10/1000;
run_time_integration('C',data,'C eta 5')

%%
% data.eta  = 0.25;
% data.mms_to_run = 10/1000;    
% run_time_integration('S',data,'S eta 1')
% 
% data.eta  = 0.33;
% data.mms_to_run = 10/1000;
% run_time_integration('S',data,'S eta 2')
% 
% data.eta  = 0.6;
% data.mms_to_run = 10/1000;
% run_time_integration('S',data,'S eta 3')
% 
% data.eta  = 0.66;
% data.mms_to_run = 5/1000;
% run_time_integration('S',data,'S eta 4')
% 
% data.eta  = 0.75;
% data.mms_to_run = 2/1000;
% run_time_integration('S',data,'S eta 5')
