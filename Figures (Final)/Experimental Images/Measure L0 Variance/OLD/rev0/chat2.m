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
[x1,y1] = ginput(1); plot(x1,y1,'ro','MarkerSize',10,'LineWidth',2);
[x2,y2] = ginput(1); plot(x2,y2,'ro','MarkerSize',10,'LineWidth',2);

calibDistPixels = norm([x1-x2, y1-y2]);
calibDistReal   = input('Enter the known real-world distance: ');
scaleFactor     = calibDistReal / calibDistPixels;
disp(['Calibration factor = ', num2str(scaleFactor), ' units/pixel']);

%% --- Edge definition ---
numEdges = input('Enter the number of edges (lines): ');

lines = zeros(numEdges,3); % store ax + by + c = 0
for e = 1:numEdges
    title(['Click two points along edge ', num2str(e)]);
    [xA,yA] = ginput(1); plot(xA,yA,'rx','MarkerSize',10,'LineWidth',2);
    [xB,yB] = ginput(1); plot(xB,yB,'rx','MarkerSize',10,'LineWidth',2);
    
    % Line equation
    a = yA - yB;
    b = xB - xA;
    c = xA*yB - xB*yA;
    lines(e,:) = [a b c];
end

%% --- Find intersections (nodes) ---
nodes = [];
for i = 1:numEdges
    for j = i+1:numEdges
        L1 = lines(i,:); L2 = lines(j,:);
        A = [L1(1) L1(2); L2(1) L2(2)];
        B = -[L1(3); L2(3)];
        if abs(det(A)) > 1e-9
            xy = A\B;
            nodes = [nodes; xy'];
        end
    end
end

% Remove duplicate nodes (within tolerance)
tol = 5; % pixels
keep = true(size(nodes,1),1);
for i = 1:size(nodes,1)
    if ~keep(i), continue; end
    d = sqrt(sum((nodes - nodes(i,:)).^2,2));
    dupIdx = find(d<tol & (1:size(nodes,1))'~=i);
    keep(dupIdx) = false;
end
nodes = nodes(keep,:);

%% --- Build edges from intersections ---
figure; imshow(im_OG); hold on;

edges = [];
edgeLengthsPix = [];
for e = 1:numEdges
    % Find intersections with this line
    pts = [];
    for n = 1:size(nodes,1)
        if abs(lines(e,1)*nodes(n,1) + lines(e,2)*nodes(n,2) + lines(e,3)) < tol
            pts = [pts; nodes(n,:)]; %#ok<AGROW>
        end
    end
    
    % If line has >=2 intersection nodes, connect them
    if size(pts,1) >= 2
        % Sort points along the line (project onto line direction)
        dirVec = null(lines(e,1:2)); % direction vector of line
        proj = pts*dirVec;
        [~,order] = sort(proj);
        pts = pts(order,:);
        
        for k = 1:size(pts,1)-1
            p1 = pts(k,:); p2 = pts(k+1,:);
            edges = [edges; p1 p2];
            edgeLengthsPix(end+1,1) = norm(p1-p2); %#ok<AGROW>
            plot([p1(1) p2(1)], [p1(2) p2(2)], 'g-', 'LineWidth', 2);
        end
    end
end

%% --- Display nodes ---

for i = 1:size(nodes,1)
    plot(nodes(i,1), nodes(i,2), 'bo', 'MarkerSize',8,'LineWidth',2);
    text(nodes(i,1)+5, nodes(i,2)+5, num2str(i), ...
        'Color','w','FontSize',50,'FontWeight','bold');
end

%% --- Output results ---
edgeLengthsReal = edgeLengthsPix * scaleFactor;
T = table((1:numel(edgeLengthsPix))', edgeLengthsPix, edgeLengthsReal, ...
    'VariableNames', {'Edge_ID','Length_in_Pixels','Length_in_Units'});
disp(T);