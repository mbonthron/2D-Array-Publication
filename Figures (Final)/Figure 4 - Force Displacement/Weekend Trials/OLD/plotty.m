addpath("Data from")
cd("Data from")
matfiles = dir("*.mat");
cd ..

%%
[disp1_025,Q1_025] = load_from_mat(matfiles(1).name,0.01);
[disp1_050,Q1_050] = load_from_mat(matfiles(2).name,0.01);
[disp1_075,Q1_075] = load_from_mat(matfiles(3).name,0.01);


%%
[disp2_025,Q2_025] = load_from_mat(matfiles(4).name,0.01);
[disp2_050,Q2_050] = load_from_mat(matfiles(5).name,0.01);
[disp2_075,Q2_075] = load_from_mat(matfiles(6).name,0.01);

%%
figure(100); clf; hold on; 
plot(disp1_025*1000,Q1_025,"LineWidth",2,"DisplayName","State 1 $\eta = 0.25$")
plot(disp1_050*1000,Q1_050,"LineWidth",2,"DisplayName","State 1 $\eta = 0.50$")

plot(disp2_025*1000,Q2_025,"LineWidth",2,"DisplayName","State 2 $\eta = 0.25$")
plot(disp2_050*1000,Q2_050,"LineWidth",2,"DisplayName","State 2 $\eta = 0.50$")
plot(disp2_075*1000,Q2_075,"LineWidth",2,"DisplayName","State 2 $\eta = 0.75$")


ylim([0 10])
xlim([0 10])

leg = legend();
leg.Location = "eastoutside";

grid('on')
xlabel("Displacement - [mm]")
ylabel("Force -[N]")
set(gca,'FontSize',18)