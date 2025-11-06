[disp_C_eta_2 ,force_C_eta_2,force_low_C_eta_2,force_high_C_eta_2,stiffness_C_eta_2] = load_data_from_folder('C - eta 2 - v = 0.5 mms',10);
[disp_C_eta_3 ,force_C_eta_3,force_low_C_eta_3,force_high_C_eta_3,stiffness_C_eta_3] = load_data_from_folder('C - eta 3 - v = 0.5 mms',10);
[disp_C_eta_4 ,force_C_eta_4,force_low_C_eta_4,force_high_C_eta_4,stiffness_C_eta_4] = load_data_from_folder('C - eta 4 - v = 0.5 mms',10);
[disp_C_eta_5 ,force_C_eta_5,force_low_C_eta_5,force_high_C_eta_5,stiffness_C_eta_5] = load_data_from_folder('C - eta 5 - v = 0.5 mms',10);
[disp_C_eta_6 ,force_C_eta_6,force_low_C_eta_6,force_high_C_eta_6,stiffness_C_eta_6] = load_data_from_folder('C - eta 6 - v = 0.5 mms',10);

[disp_S_eta_2 ,force_S_eta_2,force_low_S_eta_2,force_high_S_eta_2,stiffness_S_eta_2] = load_data_from_folder('S - eta 2 - v = 0.5 mms',10);
[disp_S_eta_3 ,force_S_eta_3,force_low_S_eta_3,force_high_S_eta_3,stiffness_S_eta_3] = load_data_from_folder('S - eta 3 - v = 0.5 mms',10);
[disp_S_eta_4 ,force_S_eta_4,force_low_S_eta_4,force_high_S_eta_4,stiffness_S_eta_4] = load_data_from_folder('S - eta 4 - v = 0.5 mms',10);
[disp_S_eta_5 ,force_S_eta_5,force_low_S_eta_5,force_high_S_eta_5,stiffness_S_eta_5] = load_data_from_folder('S - eta 5 - v = 0.5 mms',10);
[disp_S_eta_6 ,force_S_eta_6,force_low_S_eta_6,force_high_S_eta_6,stiffness_S_eta_6] = load_data_from_folder('S - eta 6 - v = 0.5 mms',10);



%%
C_color = [52 93 168]/255;
S_color = [196 43 43]/255;

fig_size = [1 0.85];
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
plot([0 4],stiffness_C_eta_3*[0 4],'k-','LineWidth',1.5)

run_name = "C eta 0.375.mat";
load(run_name)
plot(displacement*1000,Q,"-","linewidth",1.5,"Color",C_color,"HandleVisibility","off")

fill_in(disp_S_eta_3,force_S_eta_3,force_low_S_eta_3,force_high_S_eta_3,"S")
plot(disp_S_eta_3,force_S_eta_3,"LineWidth",1.5,"Color",S_color)
plot([0 4],stiffness_S_eta_3*[0 4],'k-','LineWidth',1.5)

run_name = "S eta 0.375.mat";
load(run_name)
plot(displacement*1000,Q,"-","linewidth",1.5,"Color",S_color,"HandleVisibility","off")

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
plot([0 4],stiffness_C_eta_4*[0 4],'k-','LineWidth',1.5)

run_name = "C eta 0.5.mat";
load(run_name)
plot(displacement*1000,Q,"-","linewidth",1.5,"Color",C_color,"HandleVisibility","off")

fill_in(disp_S_eta_4,force_S_eta_4,force_low_S_eta_4,force_high_S_eta_4,"S")
plot(disp_S_eta_4,force_S_eta_4,"LineWidth",1.5,"Color",S_color)
plot([0 4],stiffness_S_eta_4*[0 4],'k-','LineWidth',1.5)

run_name = "S eta 0.5.mat";
load(run_name)
plot(displacement*1000,Q,"-","linewidth",1.5,"Color",S_color,"HandleVisibility","off")

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
plot(displacement*1000,Q,"-","linewidth",1.5,"Color",C_color,"HandleVisibility","off")
plot([0 4],stiffness_C_eta_5*[0 4],'k-','LineWidth',1.5)

fill_in(disp_S_eta_5,force_S_eta_5,force_low_S_eta_5,force_high_S_eta_5,"S")
plot(disp_S_eta_5,force_S_eta_5,"LineWidth",1.5,"Color",S_color)
plot([0 4],stiffness_S_eta_5*[0 4],'k-','LineWidth',1.5)

run_name = "S eta 0.625.mat";
load(run_name)
plot(displacement*1000,Q,"-","linewidth",1.5,"Color",S_color,"HandleVisibility","off")
