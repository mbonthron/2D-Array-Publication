function [outputArg1,outputArg2] = plot_force_displacement_curves(folder_name,run_name,shape)
% Check if 'C' or 'S' Shape
if shape == "C"
    color = [52 93 168]/255;
else
    color = [196 43 43]/255;
end

%% Load the Data from the experimental Trials
[disp,force,force_low,force_high,~,~,~,~] = load_data_from_folder(folder_name,10);

load(run_name);
force_exp_interp = interp1(displacement*1000,Q,disp);

%%
frames = data.frames;
file_name = data.file_name + " - Force Displacement";

if ~exist("Videos\", 'dir')
    mkdir("Videos\");
end

file_name = "Videos\"+file_name;

%% Initialize the Video
fprintf("CREATING VIDEO\n")
v = VideoWriter(file_name, 'MPEG-4');
v.FrameRate = 30;
v.Quality = 100;
open(v)

%% Initialize the Frame
f = figure(); clf; hold on;
f.Position(3:4) = [1000 1000];
f.Renderer = "painters";
f.Color = "k";
xlim([0 10]);
ylim([0 1]);
ylabel("Force [N]");
xlabel("Displacement [mm]")

ax = gca;
ax.Color = 'k';                 % axes background
ax.XColor = 'w';                % x-axis color
ax.YColor = 'w';                % y-axis color
set(gca,'FontSize',32)

lw = 4;
%% Iterate over the frames
for i = 1:frames
    % Plot the Experimental Data
    f1 = fill_in(disp(1:i),force(1:i),force_low(1:i),force_high(1:i),shape);
    l1 = plot(disp(1:i),force(1:i),"-","LineWidth",lw,"Color",color);

    % Plot the Analytical Prediction
    l2 = plot(disp(1:i),force_exp_interp(1:i),"--","linewidth",lw,"Color",color,"HandleVisibility","off");

    % Write the frame to the file
    writeVideo(v, getframe(f));

    % Delete the plotted values
    delete(f1)
    delete(l1);
    delete(l2);
end

%% Close the Video and the Figure
fprintf("CLOSING VIDEO\n")
close(v)
close(f);
end

