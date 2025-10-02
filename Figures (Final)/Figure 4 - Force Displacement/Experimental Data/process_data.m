[disp_C_eta_050 ,force_C_eta_050,energy_C_eta_050] = load_data_from_folder('C Shape, eta = .5, v = .1');
[disp_C_eta_075 ,force_C_eta_075,energy_C_eta_075] = load_data_from_folder('C Shape, eta = .75, v = .1');


[disp_S_eta_075 ,force_S_eta_075,energy_S_eta_075] = load_data_from_folder('S Shape, eta = .75, v = .1');
[disp_S_eta_050 ,force_S_eta_050,energy_S_eta_050] = load_data_from_folder('S Shape, eta = .5, v = .1');
[disp_S_eta_025 ,force_S_eta_025,energy_S_eta_025] = load_data_from_folder('S Shape, eta = .25, v = .1');

alphaval = 0.5;
color1 = [0.0660    0.4430    0.7450 alphaval];
color2 = [0.8660    0.3290    0.0000 alphaval];
color3 = [0.9290    0.6940    0.1250 alphaval];
color4 = [0.5210    0.0860    0.8190 alphaval];
color5 = [0.2310    0.6660    0.1960 alphaval];

%%
figure(1); clf; hold on
plot(disp_C_eta_050{1},force_C_eta_050{1},"Color",color2,"DisplayName","C $\eta = 0.50$")
plot(disp_C_eta_075{1},force_C_eta_075{1},"Color",color1,"DisplayName","C $\eta = 0.75$")
plot(disp_S_eta_025{1},force_S_eta_025{1},"Color",color3,"DisplayName","S $\eta = 0.25$")
plot(disp_S_eta_050{1},force_S_eta_050{1},"Color",color4,"DisplayName","S $\eta = 0.50$")
plot(disp_S_eta_075{1},force_S_eta_075{1},"Color",color5,"DisplayName","S $\eta = 0.75$")
for i = 2:5
    plot(disp_C_eta_050{i},force_C_eta_050{i},"Color",color2,"HandleVisibility","off")
    plot(disp_C_eta_075{i},force_C_eta_075{i},"Color",color1,"HandleVisibility","off")
    plot(disp_S_eta_025{i},force_S_eta_025{i},"Color",color3,"HandleVisibility","off")
    plot(disp_S_eta_050{i},force_S_eta_050{i},"Color",color4,"HandleVisibility","off")
    plot(disp_S_eta_075{i},force_S_eta_075{i},"Color",color5,"HandleVisibility","off")
end

%
load('C eta 50.mat')
plot(displacement*1000,Q,"linewidth",2,"Color",color2(1:3),"HandleVisibility","off")

load('C eta 75.mat')
plot(displacement*1000,Q,"linewidth",2,"Color",color1(1:3),"HandleVisibility","off")

load('S eta 25.mat')
plot(displacement*1000,Q,"linewidth",2,"Color",color3(1:3),"HandleVisibility","off")

load('S eta 50.mat')
plot(displacement*1000,Q,"linewidth",2,"Color",color4(1:3),"HandleVisibility","off")

load('S eta 75.mat')
plot(displacement*1000,Q,"linewidth",2,"Color",color5(1:3),"HandleVisibility","off")

legend()
set(gca,"FontSize",12)
xlabel("Displacement - [mm]")
ylabel("Force - [N]")
ylim([0 6])
%%
% figure(2); clf; hold on
% plot(disp_C_eta_050{1},energy_C_eta_050{1},"Color",color2,"DisplayName","C $\eta = 0.50$")
% plot(disp_C_eta_075{1},energy_C_eta_075{1},"Color",color1,"DisplayName","C $\eta = 0.75$")
% plot(disp_S_eta_025{1},energy_S_eta_025{1},"Color",color3,"DisplayName","S $\eta = 0.25$")
% plot(disp_S_eta_050{1},energy_S_eta_050{1},"Color",color4,"DisplayName","S $\eta = 0.50$")
% plot(disp_S_eta_075{1},energy_S_eta_075{1},"Color",color5,"DisplayName","S $\eta = 0.75$")
% for i = 2:5
%     plot(disp_C_eta_050{i},energy_C_eta_050{i},"Color",color2,"HandleVisibility","off")
%     plot(disp_C_eta_075{i},energy_C_eta_075{i},"Color",color1,"HandleVisibility","off")
%     plot(disp_S_eta_025{i},energy_S_eta_025{i},"Color",color3,"HandleVisibility","off")
%     plot(disp_S_eta_050{i},energy_S_eta_050{i},"Color",color4,"HandleVisibility","off")
%     plot(disp_S_eta_075{i},energy_S_eta_075{i},"Color",color5,"HandleVisibility","off")
% end
% 
% legend()
% set(gca,"FontSize",12)
% xlabel("Displacement -[mm]")
% ylabel("Energy - [Nmm]")