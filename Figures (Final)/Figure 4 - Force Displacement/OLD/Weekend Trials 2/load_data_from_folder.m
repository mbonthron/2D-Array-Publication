function [displacement_cell,force_cell,energy_cell,velocity_mean] = load_data_from_folder(folder_name,mm_max)
% Go into directory
cd ('Trials 4 - Michael')
cd(folder_name)

%%
file_name = dir('*.xlsx');
trial_count = length(file_name);

displacement_cell = {};
force_cell = {};
energy_cell = {};
velocity_vector = [];

for i = 1:trial_count
    % Load the data
    A = readtable(file_name(i).name,ReadVariableNames=false);
    time        = A.Var1;
    position    = A.Var2; 
    force       = A.Var3;

    startidx = find(force>0.01,1,"first");

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
    position_interp = linspace(0,mm_max,50);
    force_interp    = interp1(position,force,position_interp);


    % Save the data
    displacement_cell{i} = position_interp;
    force_cell{i} = force_interp;

    energy_cell{i} = cumsum(position.*force);

end

velocity_mean = mean(velocity_vector);

%%
% Return to 
cd ..
cd ..
end

