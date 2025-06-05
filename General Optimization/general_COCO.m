function [data,run_max_E_per_b] = general_COCO(data, bpoints)
%% Visualize the point and connection between the nodes
plot_grid(data,1);

%% General Stuff for Coco Continuation
run_number = 1;
run_name1 = [data.shape_name '_run' sprintf('%.0f',run_number)];

Ahat0 = zeros(2*(data.N*(data.N_modes)-data.constraint_count),1);

% Define the function as the arbitrary gride ODE
f = @(x,p) COCO_arbitrary_grid_ODE(x,p,data);

parameter_names = {'b' 't'};                % Names our two parameters 'b' and 't'
initial_parameter_value = [0;0.01*pi];      % Starting values of b and t

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

fprintf("\nRun %.0f =========================================",1)
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

    fprintf("\nRun %.0f =========================================",run_number)
    coco(prob,run_name,[],1,parameter_names,computational_domain)
    run_number = run_number + 1;
end

%% ========================================================================
% CONTINUE FROM BP POINTS
% run_name_start_from = run_name1;
% BP2 = coco_bd_labs(run_name_start_from, 'BP'); % labels for BP points in run1
% 
% for i = 1:2
%     run_name = [data.shape_name '_run' sprintf('%.0f',run_number)];
%     prob = coco_prob();
%     prob = ode_ep2ep(prob,'',run_name_start_from,BP2(i));
%     prob = coco_set(prob,'cont','branch','switch');
%     prob = coco_set(prob,'cont','ItMX', iterations_max);
%     prob = coco_set(prob,'cont','NPR',0);
%     prob = coco_set(prob,'cont','h_max',hmax,'h_min',hmin);
%     prob = coco_add_event(prob,'UZ','b',UZpoints);
% 
%     fprintf("\nRun %.0f =========================================",run_number)
%     coco(prob,run_name,[],1,parameter_names,computational_domain)  
%     run_number = run_number + 1;
% end


%% Plot the Results from Coco
close all
cd('data');
coco_runs = dir([data.shape_name '*']);
cd ..

theme1 = struct('special', {{'EP','FP','HB','BP'}});
figure(9899); clf; hold on

% Initilize matrix for pertinant values from coco
UZ_data_from_coco = zeros(1,5);

run_max_E_per_b = [bpoints' zeros(length(bpoints),2)];
for i = 1:length(coco_runs)
    coco_plot_bd(theme1, coco_runs(i).name, 'x',1,'x',2,'b')
    
    % Plot shape and find b_V_matrix
    % b_V_matrix = [run_number bcrit energy UZ stability]
    b_V_matrix = plot_shape_from_COCO(coco_runs(i).name,data);
    b_V_matrix(:,1) = i;

    UZ_data_from_coco = [UZ_data_from_coco ; b_V_matrix];
end

figure(9899)
axis tight; grid on; view(3)
zlim([0 0.5]); xlim([-1 1]); ylim(0.5*[-1 1])

%% Find time integration starting points
% Round 'b' values from numerical inaccuracy
UZ_data_from_coco(:,2) = round(UZ_data_from_coco(:,2),5);
UZpoints = round(UZpoints,5);


for i = 1:length(UZpoints)
    % Find branches that had values matching the UZ points
    % (i.e. values that match the current b value)
    b_val = UZpoints(i);
    possible_branches = UZ_data_from_coco(UZ_data_from_coco(:,2)==b_val,:);

    % Ignore all the unstable branches
    possible_branches = possible_branches(possible_branches(:,5)==-1,:);

    % Find the row with the highest energy state
    [~,row] = max(possible_branches(:,3));

    % Take the Run Number and the UZ index for the higher energy stable
    % state
    run_number = possible_branches(row,1);
    uz_idx     = possible_branches(row,4);

    % Have COCO grab the A0hatp from reading the solution
    run_name_to_grab = [data.shape_name '_run' sprintf('%.0f',run_number)];
    A0hatp = COCO_grab_UZ(run_name_to_grab,uz_idx);

    % Save the A0hatp
    save("b = "+ num2str(b_val/pi) +" pi.mat","A0hatp")
    
end



%% Coco Plotting Trouble Shoot
theme1 = struct('special', {{'EP','FP','HB','BP'}});

figure(9900); clf; hold on; view(3); zlim([0 0.1]); xlim(0.3*[-1 1]); ylim(0.2*[-1 1])
grid()
coco_plot_bd(theme1, [data.shape_name '_run' sprintf('%.0f',1)], 'x',1,'x',2,'b')
coco_plot_bd(theme1, [data.shape_name '_run' sprintf('%.0f',2)], 'x',1,'x',2,'b')
coco_plot_bd(theme1, [data.shape_name '_run' sprintf('%.0f',3)], 'x',1,'x',2,'b')
% coco_plot_bd(theme1, [data.shape_name '_run' sprintf('%.0f',4)], 'x',1,'x',2,'b')
% coco_plot_bd(theme1, [data.shape_name '_run' sprintf('%.0f',5)], 'x',1,'x',2,'b')

end