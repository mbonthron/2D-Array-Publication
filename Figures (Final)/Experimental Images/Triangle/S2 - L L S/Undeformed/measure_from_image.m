% Load image
file_name = 'DSCN2999';
im_OG = imread(file_name+".jpg");


run("initialize_triangle.m")

%% --- Check if the image should be inverted ---
if check_if_dark(im_OG)
    img = imcomplement(im_OG);  % Invert the image if dark
else
    img = im_OG;    % Otherwise don't invert the image
end

figure; imshow(img); hold on;

%% --- Calibration step ---
title('Click two points defining a known distance for calibration');
[x1,y1] = ginput(1); plot(x1,y1,'r+','MarkerSize',10,'LineWidth',2);
[x2,y2] = ginput(1); plot(x2,y2,'r+','MarkerSize',10,'LineWidth',2);

calibDistPixels = norm([x1-x2, y1-y2]);
calibDistReal   = input('Enter the known real-world distance: ');
scaleFactor     = calibDistReal / calibDistPixels;
disp(['Calibration factor = ', num2str(scaleFactor), ' units/pixel']);

%% --- Vertex definition ---
numVerts = data.N;

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

%% --- Adjacency Matrix
figure(1); clf; imshow(im_OG); hold on;


up_adjac = triu(data.adjacency_matrix,1);
[left, right] = find(up_adjac == 1);

numSides = length(left);
sideLengthsPixels = zeros(numSides,1);
sideLengthsReal = zeros(numSides,1);

% Add in the vertex points
for i = 1:length(verts)
    scatter(verts(i,1),verts(i,2),100,"MarkerEdgeColor","r","LineWidth",2,"Marker","+")
end


% Go side by side - Adding to the plot
for i = 1:length(left)
    p1 = verts(left(i),:);
    p2 = verts(right(i),:);

    dPix = norm(p1 - p2);

    sideLengthsPixels(i) = dPix;
    sideLengthsReal(i) = dPix * scaleFactor;

    % Plot the edge axis
    % plot([p1(1) p2(1)], [p1(2) p2(2)], 'g-', 'LineWidth', 2);

    mid = (p1+p2)/2; % midpoint of side

    % Add in the distances
    text(mid(1), mid(2), sprintf('%.2f mm', sideLengthsReal(i)), ...
     'Color','w','FontSize',14,'FontWeight','bold', ...
     'HorizontalAlignment','center');

end

T = table((1:numSides)', sideLengthsPixels, sideLengthsReal, ...
          'VariableNames', {'Side','Length_in_Pixels','Length_in_Units'});
disp(T);

exportgraphics(gcf,file_name+" - Annotated.png")

save(file_name+" - Measurement.mat","T","verts")