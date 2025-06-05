function [data,run_max_E_per_b] = general_COCO(data, bpoints)
%% Visualize the point and connection between the nodes
plot_grid(data,1);

%% Determine the coefficient matrix and the modes to skip
% (Remove after a push / pull)
data = determine_coefficient_matrix(data);
data = determine_modes_to_skip(data);

%% General Stuff for Coco Continuation
run_number = 1;
run_name1 = [data.shape_name '_run' sprintf('%.0f',run_number)];

Ahat0 = zeros(2*(data.N*(data.N_modes)-data.constraint_count),1);

% Define the function as the arbitrary gride ODE
f = @(x,p) COCO_arbitrary_grid_ODE(x,p,data);

parameter_names = {'b' 't'};                % Names our two parameters 'b' and 't'
initial_parameter_value = [0;0.01*pi];      % Starting values of b and t

active_continuation_parameter = 'b';        % Which parameter do we want to vary
computational_domain = [0 pi*0.2];          % What is the domain of b to explore
UZpoints = bpoints;

%% Continuation Constants
iterations_max = 5000;                  % Maximum number of iterations before continuation terminates
hmin = 0.0005;                           % Minimum step size of the continuation
hmax = 0.001;                            % Maximum step size of the continuation

%% ========================================================================
%  INITIAL RUN FROM ZERO
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Ahat0,parameter_names,initial_parameter_value);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',hmax,'h_min',hmin);

coco(prob,run_name1,[],1,parameter_names,computational_domain)


%% ========================================================================
%  RERUN ANY HB POINTS AS UZ POINTS TO DO CONTINUATION FROM
bd = coco_bd_read(run_name1);
HBlbls = coco_bd_labs(run_name1, 'HB');

bcrits = zeros(1,length(HBlbls));
acrits = zeros(2*(data.N*(data.N_modes)-data.constraint_count),length(HBlbls));

for k = 1:length(HBlbls)
    bcrits(k) = coco_bd_val(bd,HBlbls(k),'b');
    acrits(:,k) = coco_bd_val(bd,HBlbls(k),'x');
end

prob = coco_add_event(prob,'UZ','b',bcrits);
prob = coco_add_event(prob,'UZ','b',UZpoints);

fprintf("Run %.0f =========================================",1)
coco(prob,run_name1,[],1,parameter_names,computational_domain)

run_number = run_number + 1;

%% ========================================================================
% Continue from UZ points
UZ = coco_bd_labs(run_name1, 'UZ'); 

for i = 1:2
    run_name = [data.shape_name '_run' sprintf('%.0f',run_number)];
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

%% ========================================================================
% CONTINUE FROM BP POINTS
run_name_start_from = run_name1;
BP2 = coco_bd_labs(run_name_start_from, 'BP'); % labels for BP points in run1

for i = 1:2
    run_name = [data.shape_name '_run' sprintf('%.0f',run_number)];
    prob = coco_prob();
    prob = ode_ep2ep(prob,'',run_name_start_from,BP2(i));
    prob = coco_set(prob,'cont','branch','switch');
    prob = coco_set(prob,'cont','ItMX', iterations_max);
    prob = coco_set(prob,'cont','NPR',0);
    prob = coco_set(prob,'cont','h_max',hmax,'h_min',hmin);
    prob = coco_add_event(prob,'UZ','b',UZpoints);

    fprintf("Run %.0f =========================================",run_number)
    coco(prob,run_name,[],1,parameter_names,computational_domain)  
    run_number = run_number + 1;
end


%% Plot the Results from Coco
close all
cd('data');
coco_runs = dir([data.shape_name '*']);
cd ..

theme1 = struct('special', {{'EP','FP','HB','BP'}});
figure(9899); clf; hold on


run_max_E_per_b = [bpoints' zeros(length(bpoints),2)];
for i = 1:length(coco_runs)
    coco_plot_bd(theme1, coco_runs(i).name, 'x',1,'x',2,'b')
    
    b_V_matrix = plot_shape_from_COCO(coco_runs(i).name,data);

    % for bpnt_idx = 1:length(bpoints)
    %     bpnt = bpoints(bpnt_idx);
    %     matrix_b = b_V_matrix(round(b_V_matrix(:,1),7) == round(bpnt,7),:);
    %     [max_E, max_E_idx] = max(matrix_b(:,2));
    %     if max_E > run_max_E_per_b(bpnt_idx,2)
    %         run_max_E_per_b(bpnt_idx, 2) = i;
    %         run_max_E_per_b(bpnt_idx, 3) = matrix_b(max_E_idx,3);
    %     end
    % end
end
% Step 1: max over j (2nd dim) → gives size [I x K]
% max_over_j = squeeze(max(b_V_matrix, [], 2));  % size [I x K]
% 
% % Step 2: find index of max over k (3rd dim) for each i
% [~, run_max_E_per_b] = max(max_over_j, [], 2);  % size [I x 1]

figure(9899)
axis tight; grid on; view(3)
zlim([0 0.5]); xlim([-1 1]); ylim(0.5*[-1 1])


%% Coco Plotting Trouble Shoot
theme1 = struct('special', {{'EP','FP','HB','BP'}});

figure(9900); clf; hold on; view(3); zlim([0 0.1]); xlim(0.3*[-1 1]); ylim(0.2*[-1 1])
grid()
coco_plot_bd(theme1, [data.shape_name '_run' sprintf('%.0f',1)], 'x',1,'x',2,'b')
coco_plot_bd(theme1, [data.shape_name '_run' sprintf('%.0f',2)], 'x',1,'x',2,'b')
coco_plot_bd(theme1, [data.shape_name '_run' sprintf('%.0f',3)], 'x',1,'x',2,'b')
coco_plot_bd(theme1, [data.shape_name '_run' sprintf('%.0f',4)], 'x',1,'x',2,'b')
coco_plot_bd(theme1, [data.shape_name '_run' sprintf('%.0f',5)], 'x',1,'x',2,'b')

end