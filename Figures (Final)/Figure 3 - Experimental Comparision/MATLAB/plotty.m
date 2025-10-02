clear; clc
load("Triad Bifurcation.mat")

L1      = (b1/2).^2;
L1_1    = (b1_1/2).^2;

f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = [5/3 1.05];

plot(b1*100/pi,midptstab1*100/pi,"k","LineWidth",2)
plot(b1*100/pi,midptunst1*100/pi,"r","LineWidth",2)

plot(b1_1*100/pi,midptstab1_1*100/pi,"k","LineWidth",2)
plot(b1_1*100/pi,midptunst1_1*100/pi,"r","LineWidth",2)

xlim([0 10])
xlabel("$b$ [mm]")
ylabel("$w(L/2)$ [mm]")
set(gca,'FontSize',9.5)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc
load("Square Bifurcation.mat")

f = figure(2); clf; hold on
f.Units = "inches";
f.Position(3:4) = [5/3 1.05];

plot(b1*100/pi,midptstab1*100/pi,"k","LineWidth",2)
plot(b1*100/pi,midptunst1*100/pi,"r","LineWidth",2)

plot(b1_1*100/pi,midptstab1_1*100/pi,"k","LineWidth",2)
% plot(b1_1*100/pi,midptunst1_1*100/pi,"r","LineWidth",2)

plot(b1_7*100/pi,midptstab1_7*100/pi,"k","LineWidth",2)
% plot(b1_7*100/pi,midptunst1_7*100/pi,"r","LineWidth",2)

plot(b1_8*100/pi,midptstab1_8*100/pi,"k","LineWidth",2)
% plot(b1_8*100/pi,midptunst1_8*100/pi,"r","LineWidth",2)



xlim([0 20])
xlabel("$b$ [mm]")
ylabel("$w(L/2)$ [mm]")
set(gca,'FontSize',9.5)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc
load("Dogbone Bifurcation.mat")

L1      = (b1/2).^2;
L1_1    = (b1_1/2).^2;
L1_8    = (b1_8/2).^2;
L1_9    = (b1_9/2).^2;

f = figure(2); clf; hold on
f.Units = "inches";
f.Position(3:4) = [5/3 1.05];

plot(b1*100/pi,midptstab1*100/pi,"k","LineWidth",2)
plot(b1*100/pi,midptunst1*100/pi,"r","LineWidth",2)

plot(b1_1*100/pi,midptstab1_1*100/pi,"k","LineWidth",2)
plot(b1_1*100/pi,midptunst1_1*100/pi,"r","LineWidth",2)

plot(b1_8*100/pi,midptstab1_8*100/pi,"k","LineWidth",2)
% plot(b1_8*100/pi,midptunst1_8*100/pi,"r","LineWidth",2)

plot(b1_9*100/pi,midptstab1_9*100/pi,"k","LineWidth",2)
% plot(b1_9*100/pi,midptunst1_9*100/pi,"r","LineWidth",2)

xlim([0 10])
xlabel("$b$ [mm]")
ylabel("$w(L/2)$ [mm]")
set(gca,'FontSize',9.5)