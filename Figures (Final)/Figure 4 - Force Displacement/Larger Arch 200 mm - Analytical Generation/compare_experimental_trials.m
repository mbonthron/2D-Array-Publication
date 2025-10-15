[disp_C_eta_050 ,force_C_eta_050,energy_C_eta_050] = load_data_from_folder('C - eta = 0.50 - v = 0.1 mms','Trials 2 - Michael');
[disp_C_eta_025 ,force_C_eta_025,energy_C_eta_025] = load_data_from_folder('C - eta = 0.25 - v = 0.1 mms','Trials 2 - Michael');

[disp_S_eta_075 ,force_S_eta_075,energy_S_eta_075] = load_data_from_folder('S - eta = 0.75 - v = 0.1 mms','Trials 2 - Michael');
[disp_S_eta_050 ,force_S_eta_050,energy_S_eta_050] = load_data_from_folder('S - eta = 0.50 - v = 0.1 mms','Trials 2 - Michael');
[disp_S_eta_025 ,force_S_eta_025,energy_S_eta_025] = load_data_from_folder('S - eta = 0.25 - v = 0.1 mms','Trials 2 - Michael');


figure(1); clf; hold on
plot(disp_C_eta_050{1},force_C_eta_050{1},"Color",color2,"DisplayName","C \eta = 0.50")
plot(disp_C_eta_025{1},force_C_eta_025{1},"Color",color1,"DisplayName","C \eta = 0.75")
plot(disp_S_eta_025{1},force_S_eta_025{1},"Color",color3,"DisplayName","S \eta = 0.25")
plot(disp_S_eta_050{1},force_S_eta_050{1},"Color",color4,"DisplayName","S \eta = 0.50")
plot(disp_S_eta_075{1},force_S_eta_075{1},"Color",color5,"DisplayName","S \eta = 0.75")
for i = 2:5
    plot(disp_C_eta_050{i},force_C_eta_050{i},"Color",color2,"HandleVisibility","off")
    plot(disp_C_eta_025{i},force_C_eta_025{i},"Color",color1,"HandleVisibility","off")
    plot(disp_S_eta_025{i},force_S_eta_025{i},"Color",color3,"HandleVisibility","off")
    plot(disp_S_eta_050{i},force_S_eta_050{i},"Color",color4,"HandleVisibility","off")
    plot(disp_S_eta_075{i},force_S_eta_075{i},"Color",color5,"HandleVisibility","off")
end

