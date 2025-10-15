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
UZ = coco_bd_labs('dogbone1-2', 'UZ');
[Ahat,b] = COCO_grab_UZ('dogbone1-2',UZ(2));
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data);
save("Frustrated State Negative a2.mat","Ahat","b")

% Positive a2
[Ahat,b] = COCO_grab_UZ('dogbone1-2',UZ(1));
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data);
save("Frustrated State Positive a2.mat","Ahat","b")

%% =========================================================
%  Run coco w/o variance
%  =========================================================
% data.variance = ([8.7 9.2 8.9 8.5 8.5]/8.7)';
data.variance = ([10 10.5 10.5 10.5 10.5 ]/10)';
% data.variance = ([1 1.01 1.02 1.02 0.99])';

run('run_coco_variance.m')

% Save the positive 'C' shape w/ variance as starting point for time
% integration
UZ = coco_bd_labs('vdogbone1-1', 'UZ');
[Ahat,b] = COCO_grab_UZ('vdogbone1-1',UZ(1));
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data);
if A(1)<0
    [Ahat,b] = COCO_grab_UZ('vdogbone1-1',UZ(2));
    A = determine_A_from_Ahat(Ahat',data)';
    plot_system_once(A,data);
end

save("C State.mat","Ahat","b")

% Save the positive 'S' shape w/ variance as starting point for time
% integration
UZ = coco_bd_labs('vdogbone1-9', 'UZ');
[Ahat,b] = COCO_grab_UZ('vdogbone1-9',UZ(1));
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data);
save("S State.mat","Ahat","b")


%% =========================================================
%  Run Time Integration for given eta value
%  =========================================================
data.beta = 0.1;
data.alpha_D = .403/1000;  % m/s displacement
data.interp_count = 500;
data.mms_to_run = 10/1000;

eta_vector = 0.375:0.125:0.625;
for i = 1:length(eta_vector)
    data.eta = eta_vector(i);
    run_name = 'C eta ' + string(data.eta);
    run_time_integration('C',data,run_name)

    run_name = 'S eta ' + string(data.eta);
    run_time_integration('S',data,run_name)
end

%% =========================================================
%  Make Videos of displacement
%  =========================================================
for i = 1:length(eta_vector)
    data.eta = eta_vector(i);
    run_name = 'C eta ' + string(data.eta);
    load(run_name)
    plot_system_over_time(t,A,data,"C")
end

for i = 1:length(eta_vector)
    data.eta = eta_vector(i);
    run_name = 'S eta ' + string(data.eta);
    load(run_name)
    plot_system_over_time(t,A,data,"S")
    
end

%% =========================================================
%  Make Videos of force displacement
%  =========================================================
plot_force_displacement_curves('C - eta 3 - v = 0.5 mms','C eta 0.375.mat',"C")
plot_force_displacement_curves('C - eta 4 - v = 0.5 mms','C eta 0.5.mat',"C")
plot_force_displacement_curves('C - eta 5 - v = 0.5 mms','C eta 0.625.mat',"C")

plot_force_displacement_curves('S - eta 3 - v = 0.5 mms','S eta 0.375.mat',"S")
plot_force_displacement_curves('S - eta 4 - v = 0.5 mms','S eta 0.5.mat',"S")
plot_force_displacement_curves('S - eta 5 - v = 0.5 mms','S eta 0.625.mat',"S")