addpath('PRL')
file_name = dir('*.xlsx');

figure(4); clf; hold on

look_at = setdiff(1:18,[1 4 9 12 13 14 15 17]);
counter = 1;

for i = look_at
    A = readtable(file_name(i).name,ReadVariableNames=false);
    position = A.Var2; 
    force = A.Var3;

    startidx = find(force>0.05);

    position = position(startidx:end);
    force = force(startidx:end);

    plot(position-position(1),force,'color',[0 0 0 0.5],"DisplayName","Trial "+counter,'LineWidth',2)

    counter = counter + 1;
end

cd('PRL')
matfiles = dir('*.mat');
for i = 1:length(matfiles)
    load(matfiles(i).name)

    plot(disp_dimensional*1000,Rxn_dimensional(1:snap_through_index),"LineWidth",3,'DisplayName',matfiles(i).name(1:end-4))
end
cd ..

legend()
xlabel('Position')
ylabel('Force')

