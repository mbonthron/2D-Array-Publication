f = figure(1); clf; hold on

color1 = "r";
color6 = "b";

for i = 1:5
    % plot(disp_C_eta_1{i},force_C_eta_1{i},"Color",color1,"HandleVisibility","off","LineWidth",2)
    % plot(disp_C_eta_2{i},force_C_eta_2{i},"Color",color1,"HandleVisibility","off","LineWidth",2)
    % plot(disp_C_eta_3{i},force_C_eta_3{i},"Color",color1,"HandleVisibility","off","LineWidth",2)
    % plot(disp_C_eta_4{i},force_C_eta_4{i},"Color",color1,"HandleVisibility","off","LineWidth",2)
    plot(disp_C_eta_5{i},force_C_eta_5{i},"Color",color1,"HandleVisibility","off","LineWidth",2)

    % plot(disp_S_eta_1{i},force_S_eta_1{i},"Color",color6,"HandleVisibility","off","LineWidth",2)
    % plot(disp_S_eta_2{i},force_S_eta_2{i},"Color",color6,"HandleVisibility","off","LineWidth",2)
    % plot(disp_S_eta_3{i},force_S_eta_3{i},"Color",color6,"HandleVisibility","off","LineWidth",2)
    % plot(disp_S_eta_4{i},force_S_eta_4{i},"Color",color6,"HandleVisibility","off","LineWidth",2)
    plot(disp_S_eta_5{i},force_S_eta_5{i},"Color",color6,"HandleVisibility","off","LineWidth",2)
end

f.Units = "inches";
f.Position(3:4) = [3.4 3];
set(gca,"FontSize",10)
% legend("Location","northwest")
set(gca,"FontSize",12)
xlabel("Displacement - [mm]")
ylabel("Force - [N]")
ylim([0 4])