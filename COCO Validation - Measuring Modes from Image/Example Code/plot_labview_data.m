excel_file_name = dir('*.xlsx').name;
T = readtable(excel_file_name);

t_raw = datetime(T.Time);
x_raw = T.Untitled;
F_raw = T.Untitled1;

%%
t = t_raw(1:end);
x = x_raw(1:end);
F = F_raw(1:end);

t_sec = seconds(t-t(1));
save("LabView Data.mat","t_sec","x","F")

%%
close all
f=figure(1);
f.Color = 'white';
f.Position(3) = f.Position(3)*2;
clf, hold on
tiledlayout(1,2)
nexttile
plot(t_sec,x,'linewidth',3)
speed = (x(end)-x(1)) / ((t_sec(end) - t_sec(1)));
speed_string = sprintf("%.3f",speed);

xlabel("Time - [sec]")
ylabel("Position - [mm]")
title("Position vs. Time - "+speed_string + " mm/s")
set(gca,"FontSize",18)
grid()

nexttile
plot(t_sec,F,'linewidth',3)

xlabel("Time - [sec]")
ylabel("Force - [N]")
title("Force vs. Time")
set(gca,"FontSize",18)
grid()

exportgraphics(gcf,"LabView Data.png")