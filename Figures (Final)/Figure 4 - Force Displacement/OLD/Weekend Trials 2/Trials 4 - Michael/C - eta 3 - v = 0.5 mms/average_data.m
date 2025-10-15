files = dir('*.xlsx')

disp_max = 10;

figure(1); clf; hold on

for i = 1:length(files)
    file_name = files(i).name;
    
    A = readtable(file_name,'PreserveVariableNames',true);
    
    position = A.Untitled;
    force    = A.('Untitled 1');
    
    % Find first positive force
    idx = find(force>0,1,'first');
    position = position(idx:end);
    force    = force(idx:end);
    
    displacement = position - position(1);
    
    plot(displacement,force)
end