%% Open connection to webcam
clear cam
fprintf("IF CAMERA DOES NOT OPEN, CHANGE LINE webcam(i) VALUE\n")
cam = webcam(0);

%% Open the webcamera feed
fig = figure('NumberTitle','off','MenuBar','none');
fig.Name = 'My Camera';

% Create an axes the size of the webcamera
ax = axes(fig); 
frame = snapshot(cam); 
im = image(ax,zeros(size(frame),'uint8')); 
axis(ax,'image');

% Preview the Camera feed
preview(cam,im)


%%
hold on
frame_size = size(frame);
for i = 1:10
    plot(frame_size(2)*[i i]/10,[0 frame_size(1)],":",'color','white','linewidth',2)
    plot([0 frame_size(2)],frame_size(1)*[i i]/10,":",'color','white','linewidth',2)
end

for i = 0.25:0.25:10
    if ~ismember(i,1:10)
        plot(frame_size(2)*[i i]/10,[0 frame_size(1)],":",'color','white','linewidth',0.5)
        plot([0 frame_size(2)],frame_size(1)*[i i]/10,":",'color','white','linewidth',0.5)
    end
end
