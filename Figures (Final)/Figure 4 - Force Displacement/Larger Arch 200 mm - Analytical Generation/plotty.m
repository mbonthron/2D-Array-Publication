[disp_C_eta_3 ,force_C_eta_3,force_low_C_eta_3,force_high_C_eta_3,energy_C_eta_3,e_C3_low,e_C3_hig,C_velocity_3] = load_data_from_folder('C - eta 3 - v = 0.5 mms',10);
[disp_C_eta_4 ,force_C_eta_4,force_low_C_eta_4,force_high_C_eta_4,energy_C_eta_4,e_C4_low,e_C4_hig,C_velocity_4] = load_data_from_folder('C - eta 4 - v = 0.5 mms',10);
[disp_C_eta_5 ,force_C_eta_5,force_low_C_eta_5,force_high_C_eta_5,energy_C_eta_5,e_C5_low,e_C5_hig,C_velocity_5] = load_data_from_folder('C - eta 5 - v = 0.5 mms',10);

[disp_S_eta_3 ,force_S_eta_3,force_low_S_eta_3,force_high_S_eta_3,energy_S_eta_3,e_S3_low,e_S3_hig,C_velocity_3] = load_data_from_folder('S - eta 3 - v = 0.5 mms',10);
[disp_S_eta_4 ,force_S_eta_4,force_low_S_eta_4,force_high_S_eta_4,energy_S_eta_4,e_S4_low,e_S4_hig,C_velocity_4] = load_data_from_folder('S - eta 4 - v = 0.5 mms',10);
[disp_S_eta_5 ,force_S_eta_5,force_low_S_eta_5,force_high_S_eta_5,energy_S_eta_5,e_S5_low,e_S5_hig,C_velocity_5] = load_data_from_folder('S - eta 5 - v = 0.5 mms',10);



%%
C_color = [52 93 168]/255;
S_color = [196 43 43]/255;

fig_size = [1.65 1.25];
xlimss = [0 10];
%% Eta = 3
f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = fig_size;
set(gca,"FontSize",10)
% legend("Location","northwest")
% xlabel("Displacement - [mm]")
% ylabel("Force - [N]")
ylim([0 1])
xlim(xlimss)

fill_in(disp_C_eta_3,force_C_eta_3,force_low_C_eta_3,force_high_C_eta_3,"C")
plot(disp_C_eta_3,force_C_eta_3,"-","LineWidth",1.5,"Color",C_color)

run_name = "C eta 0.375.mat";
load(run_name)
plot(displacement*1000,Q,"--","linewidth",1.5,"Color",C_color,"HandleVisibility","off")

fill_in(disp_S_eta_3,force_S_eta_3,force_low_S_eta_3,force_high_S_eta_3,"S")
plot(disp_S_eta_3,force_S_eta_3,"LineWidth",1.5,"Color",S_color)

run_name = "S eta 0.375.mat";
load(run_name)
plot(displacement*1000,Q,"--","linewidth",1.5,"Color",S_color,"HandleVisibility","off")

%%
f = figure(2); clf; hold on
f.Units = "inches";
f.Position(3:4) = fig_size;
set(gca,"FontSize",10)
% legend("Location","northwest")
% xlabel("Displacement - [mm]")
% ylabel("Force - [N]")
ylim([0 1])
xlim(xlimss)


fill_in(disp_C_eta_4,force_C_eta_4,force_low_C_eta_4,force_high_C_eta_4,"C")
plot(disp_C_eta_4,force_C_eta_4,"LineWidth",1.5,"Color",C_color)

run_name = "C eta 0.5.mat";
load(run_name)
plot(displacement*1000,Q,"--","linewidth",1.5,"Color",C_color,"HandleVisibility","off")

fill_in(disp_S_eta_4,force_S_eta_4,force_low_S_eta_4,force_high_S_eta_4,"S")
plot(disp_S_eta_4,force_S_eta_4,"LineWidth",1.5,"Color",S_color)

run_name = "S eta 0.5.mat";
load(run_name)
plot(displacement*1000,Q,"--","linewidth",1.5,"Color",S_color,"HandleVisibility","off")

%%
f = figure(3); clf; hold on
f.Units = "inches";
f.Position(3:4) = fig_size;
set(gca,"FontSize",10)
% legend("Location","northwest")
% xlabel("Displacement - [mm]")
% ylabel("Force - [N]")
ylim([0 1])
xlim(xlimss)

fill_in(disp_C_eta_5,force_C_eta_5,force_low_C_eta_5,force_high_C_eta_5,"C")
plot(disp_C_eta_5,force_C_eta_5,"LineWidth",1.5,"Color",C_color)

run_name = "C eta 0.625.mat";
load(run_name)
plot(displacement*1000,Q,"--","linewidth",1.5,"Color",C_color,"HandleVisibility","off")

fill_in(disp_S_eta_5,force_S_eta_5,force_low_S_eta_5,force_high_S_eta_5,"S")
plot(disp_S_eta_5,force_S_eta_5,"LineWidth",1.5,"Color",S_color)

run_name = "S eta 0.625.mat";
load(run_name)
plot(displacement*1000,Q,"--","linewidth",1.5,"Color",S_color,"HandleVisibility","off")


%% Eta = 3 - ENERGY
f = figure(4); clf; hold on
f.Units = "inches";
f.Position(3:4) = fig_size;
set(gca,"FontSize",10)
% legend("Location","northwest")
% xlabel("Displacement - [mm]")
% ylabel("Force - [N]")
ylim([0 5])
xlim(xlimss)

fill_in(disp_C_eta_3,energy_C_eta_3,e_C3_low,e_C3_hig,"C")
plot(disp_C_eta_3,energy_C_eta_3,"LineWidth",1.5,"Color",C_color)

run_name = "C eta 0.375.mat";
load(run_name)
plot(displacement*1000,energy*1000,"--","linewidth",1.5,"Color",C_color,"HandleVisibility","off")

fill_in(disp_S_eta_3,energy_S_eta_3,e_S3_low,e_S3_hig,"S")
plot(disp_S_eta_3,energy_S_eta_3,"LineWidth",1.5,"Color",S_color)

run_name = "S eta 0.375.mat";
load(run_name)
plot(displacement*1000,energy*1000,"--","linewidth",1.5,"Color",S_color,"HandleVisibility","off")


%%
f = figure(5); clf; hold on
f.Units = "inches";
f.Position(3:4) = fig_size;
set(gca,"FontSize",10)
% legend("Location","northwest")
% xlabel("Displacement - [mm]")
% ylabel("Force - [N]")
ylim([0 5])
xlim(xlimss)

fill_in(disp_C_eta_4,energy_C_eta_4,e_C4_low,e_C4_hig,"C")
plot(disp_C_eta_4,energy_C_eta_4,"LineWidth",1.5,"Color",C_color)

run_name = "C eta 0.5.mat";
load(run_name)
plot(displacement*1000,energy*1000,"--","linewidth",1.5,"Color",C_color,"HandleVisibility","off")

fill_in(disp_S_eta_4,energy_S_eta_4,e_S4_low,e_S4_hig,"S")
plot(disp_S_eta_4,energy_S_eta_4,"LineWidth",1.5,"Color",S_color)

run_name = "S eta 0.5.mat";
load(run_name)
plot(displacement*1000,energy*1000,"--","linewidth",1.5,"Color",S_color,"HandleVisibility","off")

%%
f = figure(6); clf; hold on
f.Units = "inches";
f.Position(3:4) = fig_size;
set(gca,"FontSize",10)
% legend("Location","northwest")
% xlabel("Displacement - [mm]")
% ylabel("Force - [N]")
ylim([0 5])
xlim(xlimss)


fill_in(disp_C_eta_5,energy_C_eta_5,e_C4_low,e_C5_hig,"C")
plot(disp_C_eta_5,energy_C_eta_5,"LineWidth",1.5,"Color",C_color)

run_name = "C eta 0.625.mat";
load(run_name)
plot(displacement*1000,energy*1000,"--","linewidth",1.5,"Color",C_color,"HandleVisibility","off")

fill_in(disp_S_eta_5,energy_S_eta_5,e_S5_low,e_S5_hig,"S")
plot(disp_S_eta_5,energy_S_eta_5,"LineWidth",1.5,"Color",S_color)

run_name = "S eta 0.625.mat";
load(run_name)
plot(displacement*1000,energy*1000,"--","linewidth",1.5,"Color",S_color,"HandleVisibility","off")