%%
C_color = [52 93 168]/255;
S_color = [196 43 43]/255;

fig_size = 2*[1.65 0.85];
xlimss = [0 10];

%%
color1 = [52 93 168]/255;
color2 = [52 235 94]/255;
color3 = [212 188 55]/255;
color4 = [196 43 43]/255;
color5 = [255 160 97]/255;
color6 = [247 92 206]/255;

%%
f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = fig_size;
set(gca,"FontSize",10)
% legend("Location","northwest")
% xlabel("Displacement - [mm]")
% ylabel("Force - [N]")
ylim([0 1])
xlim(xlimss)

fill_in(disp_C_eta_3,force_C_eta_3,force_low_C_eta_3,force_high_C_eta_3,color1)
plot(disp_C_eta_3,force_C_eta_3,"-","LineWidth",1.5,"Color",color1)
run_name = "C eta 0.375.mat";
load(run_name)
plot(displacement*1000,Q,"--","linewidth",1.5,"Color",color1,"HandleVisibility","off")

fill_in(disp_C_eta_4,force_C_eta_4,force_low_C_eta_4,force_high_C_eta_4,color2)
plot(disp_C_eta_4,force_C_eta_4,"-","LineWidth",1.5,"Color",color2)
run_name = "C eta 0.5.mat";
load(run_name)
plot(displacement*1000,Q,"--","linewidth",1.5,"Color",color2,"HandleVisibility","off")

fill_in(disp_C_eta_5,force_C_eta_5,force_low_C_eta_5,force_high_C_eta_5,color3)
plot(disp_C_eta_5,force_C_eta_5,"-","LineWidth",1.5,"Color",color3)
run_name = "C eta 0.625.mat";
load(run_name)
plot(displacement*1000,Q,"--","linewidth",1.5,"Color",color3,"HandleVisibility","off")



fill_in(disp_S_eta_3,force_S_eta_3,force_low_S_eta_3,force_high_S_eta_3,color4)
plot(disp_S_eta_3,force_S_eta_3,"LineWidth",1.5,"Color",color4)
run_name = "S eta 0.375.mat";
load(run_name)
plot(displacement*1000,Q,"--","linewidth",1.5,"Color",color4,"HandleVisibility","off")

fill_in(disp_S_eta_4,force_S_eta_4,force_low_S_eta_4,force_high_S_eta_4,color5)
plot(disp_S_eta_4,force_S_eta_4,"-","LineWidth",1.5,"Color",color5)
run_name = "S eta 0.5.mat";
load(run_name)
plot(displacement*1000,Q,"--","linewidth",1.5,"Color",color5,"HandleVisibility","off")

fill_in(disp_S_eta_5,force_S_eta_5,force_low_S_eta_5,force_high_S_eta_5,color6)
plot(disp_S_eta_5,force_S_eta_5,"-","LineWidth",1.5,"Color",color6)
run_name = "S eta 0.625.mat";
load(run_name)
plot(displacement*1000,Q,"--","linewidth",1.5,"Color",color6,"HandleVisibility","off")