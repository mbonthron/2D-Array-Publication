N = 11;
color1 = [0 0 1 0.5];
color2 = [1 0 0 0.5];
colormat = color1 + (color2 - color1).*linspace(0,1,N).';


%%
figure(1); clf; hold on

matfiles = dir('*.mat');
vecc = fliplr(1:length(matfiles));

for i = vecc
    load(matfiles(i).name)

    plot(disp_dimensional*1000,Rxn_dimensional(1:snap_through_index),"LineWidth",3,'DisplayName',matfiles(i).name(1:end-4),"Color",colormat(i,:))
end
xlabel("Position")
ylabel("Force")
legend()