% Script to find side lengths of arbitrary polygon from an image with calibration
% -------------------------------------------------------------------

% Load image
file_name = 'DSCN2896';
im_OG = imread(file_name+".jpg");


% Check if white or black background
if size(im_OG,3) == 3
    grayImg = rgb2gray(im_OG);
else
    grayImg = im_OG;
end
avgBrightness = mean(grayImg(:)) / 255;
if avgBrightness < 0.5
    invert = true; % black for bright images
else
    invert = false; % white for dark images
end

if invert
    img = imcomplement(im_OG);
else
    img = im_OG;
end


figure; imshow(img); hold on;

%% --- Calibration step ---
title('Click two points defining a known distance for calibration');
[x1,y1] = ginput(1);
plot(x1,y1,'ro','MarkerSize',10,'LineWidth',2);
[x2,y2] = ginput(1);
plot(x2,y2,'ro','MarkerSize',10,'LineWidth',2);

calibDistPixels = norm([x1-x2, y1-y2]);
prompt = 'Enter the known real-world distance (in your preferred units): ';
calibDistReal = input(prompt);

scaleFactor = calibDistReal / calibDistPixels;  % units per pixel
disp(['Calibration factor: ', num2str(scaleFactor), ' units/pixel']);

%% --- Polygon step ---
numSides = input('Enter the number of sides of the polygon: ');

title(['Now click two points along each side (total ', num2str(2*numSides), ' points)']);
points = zeros(2*numSides, 2);
for i = 1:2*numSides
    [x,y] = ginput(1);
    plot(x, y, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    points(i,:) = [x,y];
end

% Line equations (ax + by + c = 0)
lines = zeros(numSides, 3);
for i = 1:numSides
    p1 = points(2*i-1, :);
    p2 = points(2*i, :);
    a = p1(2) - p2(2);
    b = p2(1) - p1(1);
    c = p1(1)*p2(2) - p2(1)*p1(2);
    lines(i,:) = [a b c];
end

% Find intersections (polygon vertices)
verts = zeros(numSides, 2);
for i = 1:numSides
    L1 = lines(i,:);
    L2 = lines(mod(i,numSides)+1, :);
    A = [L1(1) L1(2); L2(1) L2(2)];
    B = -[L1(3); L2(3)];
    xy = A\B;
    verts(i,:) = xy';
    plot(xy(1), xy(2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
end

% Compute side lengths
sideLengthsPixels = zeros(numSides,1);
sideLengthsReal = zeros(numSides,1);
for i = 1:numSides
    p1 = verts(i,:);
    p2 = verts(mod(i,numSides)+1,:);
    dPix = norm(p1 - p2);
    sideLengthsPixels(i) = dPix;
    sideLengthsReal(i) = dPix * scaleFactor;
end

% Display results
disp('Polygon side lengths:');
T = table((1:numSides)', sideLengthsPixels, sideLengthsReal, ...
          'VariableNames', {'Side','Length_in_Pixels','Length_in_Units'});
disp(T);

%%
% Draw polygon overlay on original image
figure; imshow(im_OG); hold on;
for i = 1:numSides
    p1 = verts(i,:);
    p2 = verts(mod(i,numSides)+1,:);
    plot([p1(1) p2(1)], [p1(2) p2(2)], 'g-', 'LineWidth', 2);
end


% for i = 1:numSides
%     p1 = verts(i,:);
%     p2 = verts(mod(i,numSides)+1,:);
%     mid = (p1+p2)/2; % midpoint of side
%     text(mid(1), mid(2), sprintf('%.2f', sideLengthsReal(i)), ...
%          'Color','y','FontSize',10,'FontWeight','bold', ...
%          'HorizontalAlignment','center');
% end