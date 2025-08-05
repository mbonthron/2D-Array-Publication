clear; clc;
folders = dir('b =*');

%% t = 0.05
cd(folders(1).name)
trial_files = dir('*.xlsx');

for k = 1:5
    A = readtable(trial_files(k).name,ReadVariableNames=false);
    position = A.Var2; 
    force = A.Var3;

    startidx = find(position > 128,1,'first');
    endidx = find(position > 138,1,'first');
    
    position = position(startidx:endidx);
    force = force(startidx:endidx);

    position_t_05{k} = position;
    force_t_05{k} = force;

    energy_t_05{k} = cumtrapz(position,force);
end

cd ..

%% t = 0.10
cd(folders(2).name)
trial_files = dir('*.xlsx');

for k = 1:5
    A = readtable(trial_files(k).name,ReadVariableNames=false);
    position = A.Var2; 
    force = A.Var3;

    startidx = find(position > 126,1,'first');
    endidx = find(position > 136,1,'first');
    
    position = position(startidx:endidx);
    force = force(startidx:endidx);

    figure(1); clf; hold on
    plot(position,force)

    position_t_10{k} = position;
    force_t_10{k} = force;

    energy_t_10{k} = cumtrapz(position,force);
end

cd ..

%% t = 0.15
cd(folders(3).name)
trial_files = dir('*.xlsx');

for k = 1:5
    A = readtable(trial_files(k).name,ReadVariableNames=false);
    position = A.Var2; 
    force = A.Var3;

    figure(1); clf; hold on
    plot(position,force)


    startidx = find(position > 126,1,'first');
    endidx = find(position > 135,1,'first');
    
    position = position(startidx:endidx);
    force = force(startidx:endidx);


    position_t_15{k} = position;
    force_t_15{k} = force;

    energy_t_15{k} = cumtrapz(position,force);
end

cd ..


%%
save('Loaded Data.mat')