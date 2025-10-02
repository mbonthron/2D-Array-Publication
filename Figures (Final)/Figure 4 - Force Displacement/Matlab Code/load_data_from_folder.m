function [displacement_cell,force_cell,energy_cell] = load_data_from_folder(folder_name,trial_name)
% Go into directory
cd(trial_name)
cd(folder_name)

%%
file_name = dir('*.xlsx');
trial_count = length(file_name);

displacement_cell = {};
force_cell = {};
energy_cell = {};

for i = 1:trial_count
    % Load the data
    A = readtable(file_name(i).name,ReadVariableNames=false);
    position    = A.Var2; 
    force       = abs(A.Var3);

    startidx = find(force>0.05,1,"first");

    position = position(startidx:end);
    force = force(startidx:end);

    position = position - position(1);

    displacement_cell{i} = position;
    force_cell{i} = force;

    energy_cell{i} = cumsum(position.*force);

end

%%
% Return to 
cd ..
cd ..
end

