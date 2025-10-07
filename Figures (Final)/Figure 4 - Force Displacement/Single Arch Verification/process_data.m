[disp_C_eta_3 ,force_C_eta_3,energy_C_eta_3] = load_data_from_folder('Single Arch - eta 3 - v = 0.5 mms',15);
[disp_C_eta_5 ,force_C_eta_5,energy_C_eta_5] = load_data_from_folder('Single Arch - eta 5 - v = 0.5 mms',15);


alphaval = 0.5;
color1 = [0.894, 0.102, 0.110, alphaval];  % red
color2 = [0.216, 0.494, 0.722, alphaval];  % blue
color3 = [0.302, 0.686, 0.290, alphaval];  % green
color4 = [0.596, 0.306, 0.639, alphaval];  % purple
color5 = [1.000, 0.498, 0.000, alphaval];  % orange
% color6 = [1.000, 1.000, 0.200, alphaval];  % yellow
% color7 = [0.651, 0.337, 0.157, alphaval];  % brown
% color8 = [0.969, 0.506, 0.749, alphaval];  % pink
% color9 = [0.600, 0.600, 0.600, alphaval];  % gray
% color10 = [0.121, 0.470, 0.705, alphaval]; % steel blue

color6 = color1;
color7 = color2;
color8 = color3;
color9 = color4;
color10 = color5;

%%
f = figure(1); clf; hold on
plot(disp_C_eta_3{1},force_C_eta_3{1},"Color",color3,"DisplayName","C $\eta = 3$","LineWidth",2)
plot(disp_C_eta_5{1},force_C_eta_5{1},"Color",color5,"DisplayName","C $\eta = 5$","LineWidth",2)




for i = 2:5
    plot(disp_C_eta_3{i},force_C_eta_3{i},"Color",color3,"HandleVisibility","off","LineWidth",2)
    plot(disp_C_eta_5{i},force_C_eta_5{i},"Color",color5,"HandleVisibility","off","LineWidth",2)
end


%%

load("Single Arch eta 3.mat")
scatter(displacement*1000,Q,50,"MarkerEdgeColor","k","MarkerFaceColor",color3(1:3),"HandleVisibility","off")

load("Single Arch eta 5.mat")
scatter(displacement*1000,Q,50,"MarkerEdgeColor","k","MarkerFaceColor",color5(1:3),"HandleVisibility","off")

% load("Single Arch eta 3 - N = 11.mat")
% scatter(displacement*1000,Q,50,"Marker","+","MarkerEdgeColor","k","MarkerFaceColor",color3(1:3),"HandleVisibility","off")
% 
% load("Single Arch eta 5 - N = 11.mat")
% scatter(displacement*1000,Q,50,"Marker","+","MarkerEdgeColor","k","MarkerFaceColor",color5(1:3),"HandleVisibility","off")
% 
% load("Single Arch eta 3 - N = 21.mat")
% scatter(displacement*1000,Q,50,"Marker","x","MarkerEdgeColor","k","MarkerFaceColor",color3(1:3),"HandleVisibility","off")
% 
% load("Single Arch eta 5 - N = 21.mat")
% scatter(displacement*1000,Q,50,"Marker","x","MarkerEdgeColor","k","MarkerFaceColor",color5(1:3),"HandleVisibility","off")


f.Units = "inches";
f.Position(3:4) = 3*[3.4 1.5];
set(gca,"FontSize",10)
% legend("Location","northwest")
set(gca,"FontSize",12)
xlabel("Displacement - [mm]")
ylabel("Force - [N]")
ylim([0 3])