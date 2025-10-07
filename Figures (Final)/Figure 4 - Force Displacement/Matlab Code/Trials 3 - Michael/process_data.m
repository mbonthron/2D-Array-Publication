[disp_C_eta_1 ,force_C_eta_1,energy_C_eta_1] = load_data_from_folder('C - eta 1 - v = 0.5 mms');
[disp_C_eta_2 ,force_C_eta_2,energy_C_eta_2] = load_data_from_folder('C - eta 2 - v = 0.5 mms');
[disp_C_eta_3 ,force_C_eta_3,energy_C_eta_3] = load_data_from_folder('C - eta 3 - v = 0.5 mms');
[disp_C_eta_4 ,force_C_eta_4,energy_C_eta_4] = load_data_from_folder('C - eta 4 - v = 0.5 mms');
[disp_C_eta_5 ,force_C_eta_5,energy_C_eta_5] = load_data_from_folder('C - eta 5 - v = 0.5 mms');

[disp_S_eta_1 ,force_S_eta_1,energy_S_eta_1] = load_data_from_folder('S - eta 1 - v = 0.5 mms');
[disp_S_eta_2 ,force_S_eta_2,energy_S_eta_2] = load_data_from_folder('S - eta 2 - v = 0.5 mms');
[disp_S_eta_3 ,force_S_eta_3,energy_S_eta_3] = load_data_from_folder('S - eta 3 - v = 0.5 mms');
[disp_S_eta_4 ,force_S_eta_4,energy_S_eta_4] = load_data_from_folder('S - eta 4 - v = 0.5 mms');
[disp_S_eta_5 ,force_S_eta_5,energy_S_eta_5] = load_data_from_folder('S - eta 5 - v = 0.5 mms');

alphaval = 0.5;
color1 = [0.894, 0.102, 0.110, alphaval];  % red
color2 = [0.216, 0.494, 0.722, alphaval];  % blue
color3 = [0.302, 0.686, 0.290, alphaval];  % green
color4 = [0.596, 0.306, 0.639, alphaval];  % purple
color5 = [1.000, 0.498, 0.000, alphaval];  % orange
color6 = [1.000, 1.000, 0.200, alphaval];  % yellow
color7 = [0.651, 0.337, 0.157, alphaval];  % brown
color8 = [0.969, 0.506, 0.749, alphaval];  % pink
color9 = [0.600, 0.600, 0.600, alphaval];  % gray
color10 = [0.121, 0.470, 0.705, alphaval]; % steel blue

%%
figure(1); clf; hold on
plot(disp_C_eta_1{1},force_C_eta_1{1},"Color",color1,"DisplayName","C \eta = 1","LineWidth",2)
plot(disp_C_eta_2{1},force_C_eta_2{1},"Color",color2,"DisplayName","C \eta = 2","LineWidth",2)
plot(disp_C_eta_3{1},force_C_eta_3{1},"Color",color3,"DisplayName","C \eta = 3","LineWidth",2)
plot(disp_C_eta_4{1},force_C_eta_4{1},"Color",color4,"DisplayName","C \eta = 4","LineWidth",2)
plot(disp_C_eta_5{1},force_C_eta_5{1},"Color",color5,"DisplayName","C \eta = 5","LineWidth",2)

% plot(disp_S_eta_1{1},force_S_eta_1{1},"Color",color6,"DisplayName","S \eta = 1","LineWidth",2)
% plot(disp_S_eta_2{1},force_S_eta_2{1},"Color",color7,"DisplayName","S \eta = 2","LineWidth",2)
% plot(disp_S_eta_3{1},force_S_eta_3{1},"Color",color8,"DisplayName","S \eta = 3","LineWidth",2)
% plot(disp_S_eta_4{1},force_S_eta_4{1},"Color",color9,"DisplayName","S \eta = 4","LineWidth",2)
% plot(disp_S_eta_5{1},force_S_eta_5{1},"Color",color10,"DisplayName","S \eta = 5","LineWidth",2)



for i = 2:5
    plot(disp_C_eta_1{i},force_C_eta_1{i},"Color",color1,"HandleVisibility","off","LineWidth",2)
    plot(disp_C_eta_2{i},force_C_eta_2{i},"Color",color2,"HandleVisibility","off","LineWidth",2)
    plot(disp_C_eta_3{i},force_C_eta_3{i},"Color",color3,"HandleVisibility","off","LineWidth",2)
    plot(disp_C_eta_4{i},force_C_eta_4{i},"Color",color4,"HandleVisibility","off","LineWidth",2)
    plot(disp_C_eta_5{i},force_C_eta_5{i},"Color",color5,"HandleVisibility","off","LineWidth",2)
    
    % plot(disp_S_eta_1{i},force_S_eta_1{i},"Color",color6,"HandleVisibility","off","LineWidth",2)
    % plot(disp_S_eta_2{i},force_S_eta_2{i},"Color",color7,"HandleVisibility","off","LineWidth",2)
    % plot(disp_S_eta_3{i},force_S_eta_3{i},"Color",color8,"HandleVisibility","off","LineWidth",2)
    % plot(disp_S_eta_4{i},force_S_eta_4{i},"Color",color9,"HandleVisibility","off","LineWidth",2)
    % plot(disp_S_eta_5{i},force_S_eta_5{i},"Color",color10,"HandleVisibility","off","LineWidth",2)
end

cd ..
load("C eta 1.mat")
plot(displacement*1000,Q,"linewidth",2,"Color",color1(1:3),"HandleVisibility","off")

load("C eta 2.mat")
plot(displacement*1000,Q,"linewidth",2,"Color",color2(1:3),"HandleVisibility","off")

load("C eta 3.mat")
plot(displacement*1000,Q,"linewidth",2,"Color",color3(1:3),"HandleVisibility","off")

load("C eta 4.mat")
plot(displacement*1000,Q,"linewidth",2,"Color",color4(1:3),"HandleVisibility","off")

load("C eta 5.mat")
plot(displacement*1000,Q,"linewidth",2,"Color",color5(1:3),"HandleVisibility","off")

legend("Location","northwest")
set(gca,"FontSize",12)
xlabel("Displacement - [mm]")
ylabel("Force - [N]")
ylim([0 6])