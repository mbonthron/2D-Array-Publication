% Script to find side lengths of arbitrary polygon from an image with calibration
% -------------------------------------------------------------------

% Load image
file_name = 'DSCN2896';
im_OG = imread(file_name+".jpg");





figure; imshow(img); hold on;

%% --- Calibration step ---
title('Click two points defining a known distance for calibration');
[x1,y1] = ginput(1); plot(x1,y1,'ro','MarkerSize',10,'LineWidth',2);
[x2,y2] = ginput(1); plot(x2,y2,'ro','MarkerSize',10,'LineWidth',2);

calibDistPixels = norm([x1-x2, y1-y2]);
calibDistReal   = input('Enter the known real-world distance: ');
scaleFactor     = calibDistReal / calibDistPixels;
disp(['Calibration factor = ', num2str(scaleFactor), ' units/pixel']);

%% --- Vertex definition ---
numVerts = input('Enter the number of vertices: ');

verts = zeros(numVerts,2);

for v = 1:numVerts
    title(['Vertex ',num2str(v),': select 2 points along edge 1 meeting at this vertex']);
    [xA,yA] = ginput(1); plot(xA,yA,'rx','MarkerSize',10,'LineWidth',2);
    [xB,yB] = ginput(1); plot(xB,yB,'rx','MarkerSize',10,'LineWidth',2);

    title(['Vertex ',num2str(v),': select 2 points along edge 2 meeting at this vertex']);
    [xC,yC] = ginput(1); plot(xC,yC,'rx','MarkerSize',10,'LineWidth',2);
    [xD,yD] = ginput(1); plot(xD,yD,'rx','MarkerSize',10,'LineWidth',2);

    % Line equations
    a1 = yA - yB; b1 = xB - xA; c1 = xA*yB - xB*yA;
    a2 = yC - yD; b2 = xD - xC; c2 = xC*yD - xD*yC;

    % Intersection
    A = [a1 b1; a2 b2];
    B = -[c1; c2];
    if abs(det(A)) < 1e-9
        warning('Lines for vertex %d are parallel or nearly so', v);
        verts(v,:) = [NaN,NaN];
    else
        xy = A\B;
        verts(v,:) = xy';
        plot(xy(1), xy(2), 'bo','MarkerSize',8,'LineWidth',2);
        text(xy(1)+5, xy(2)+5, sprintf('V%d',v), ...
            'Color','c','FontSize',10,'FontWeight','bold');
    end
end

%% --- Compute all pairwise distances ---
numVerts = size(verts,1);
distPix = zeros(numVerts);
distReal = zeros(numVerts);
for i = 1:numVerts
    for j = i+1:numVerts
        dPix = norm(verts(i,:) - verts(j,:));
        distPix(i,j) = dPix;
        distPix(j,i) = dPix;
        distReal(i,j) = dPix*scaleFactor;
        distReal(j,i) = distReal(i,j);
    end
end

% Display results
disp('Pairwise distances between vertices:');
ca = cell(numVerts,1);
for i = 1:numVerts, ca{i} = ['V',num2str(i)]; end
T = array2table(distReal, 'VariableNames', ca, 'RowNames', ca);
disp(T);
