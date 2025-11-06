function [position_interp,F_mean,lower_error,upper_error,stiffness_mean] = load_data_from_folder(folder_name,mm_max)
% Go into directory
cd("Trial 6")
cd(folder_name)

%%
file_name = dir('*.xlsx');
trial_count = length(file_name);

displacement_matrix = zeros(750,trial_count);
force_matrix        = zeros(750,trial_count);
energy_matrix       = zeros(750,trial_count);
stiffness_matrix    = zeros(1,trial_count);
velocity_vector     = [];

for i = 1:trial_count
    % Load the data
    A = readtable(file_name(i).name,"ReadVariableNames",false);
    time        = A.Var1;
    position    = A.Var2; 
    force       = A.Var3;   % Force offset  0.0843

    startidx = find(force<-0.025,1,"last");

    time     = time(startidx:end);
    position = position(startidx:end);
    force = force(startidx:end);

    position = position - position(1);
    time     = seconds(time - time(1));

    p = polyfit(time,position,1);
    velocity_vector = [velocity_vector p(1)];

    % Remove any non-unqiue values from position
    [position,idx] = unique(position);
    force = force(idx);

    % Apply Interpolation
    position_interp = linspace(0,mm_max,750)';
    force_interp    = interp1(position,force,position_interp);
    % force_interp    = smooth(force_interp);

    % Save the data
    displacement_matrix(:,i) = position_interp;
    force_matrix(:,i) = force_interp;

    energy_matrix(:,i) = cumtrapz(position_interp,force_interp);
   


end

velocity_mean = mean(velocity_vector);

%% Error Bar Stuff
F_mean = mean(force_matrix,2,'omitnan');
F_min  = min(force_matrix,[],2,'omitnan');
F_max  = max(force_matrix,[],2,'omitnan');

lower_error = F_mean - F_min;
upper_error = F_max - F_mean;

energy_mean = mean(energy_matrix,2,'omitnan');
e_min = min(energy_matrix,[],2,'omitnan');
e_max = max(energy_matrix,[],2,'omitnan');

e_lower_error = e_max - energy_mean;
e_upper_error = energy_mean - e_min;

p = polyfit(position_interp(1:100),F_mean(1:100), 1);


stiffness_mean = p(1);


%%
% Return to 
cd ..
cd ..
end

