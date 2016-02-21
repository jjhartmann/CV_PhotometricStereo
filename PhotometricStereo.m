%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assingment 2 - Photostereo Imaging

% Calab Sphere Array
data = ['Photostereo_RealImages/sphere-lamp1.tif'; 'Photostereo_RealImages/sphere-lamp2.tif'; 'Photostereo_RealImages/sphere-lamp3.tif'];
sphereArray = cellstr(data);

%% Init calabration Sphere. 
imgSphere = imread(sphereArray{3});

%Convert to gray
imgSphereGray = rgb2gray(imgSphere);
imshow(imgSphereGray)
[cx, cy] = ginput(1);
[px, py] = ginput(1);

% Find the radius of circle
radius = floor(sqrt((cx - px)^2 + (cy - py)^2));

% Draw circle on image
theta = 0 : (2 * pi /10000) : (2 * pi);
pline_x = radius * cos(theta) + cx;
pline_y = radius * sin(theta) + cy;
hold on; 
plot(pline_x, pline_y)
hold off;

% Crop the image based on circle. 
dim = [(cx - radius - 10) (cy - radius - 10) (2*radius + 10) (2*radius + 10)];
imgCropSphere3 = imcrop(imgSphereGray, dim);
imgCropSphere2 = rgb2gray(imcrop(imread(sphereArray{2}), dim));
imgCropSphere1 = rgb2gray(imcrop(imread(sphereArray{1}), dim));
figure(1)
imshow(imgCropSphere3)

% Get new center
[h, w] = size(imgCropSphere3);
ccx = w/2;
ccy = h/2;

%% Collect fg value pairs and E1, E2, and E3. 

% Sphere 1: Iterate over inter image x and y coordinates and build f, g and E values
fge = double(zeros(w * h, 7));
indRow = [];
index = 1;
for x = 1:w
   for y = 1:h
      
       xx = ceil(x - ccx);
       yy = ceil(y - ccy);
       % check to make sure x,y is in cicle. 
       rtmp = ceil(sqrt(xx^2 + yy^2));
       if (rtmp < radius)
          % Use sterographic projection to detect f and g
           zz = ceil(sqrt(radius^2 - (xx^2 + yy^2)));
           
           nxx = xx/radius;
           nyy = yy/radius;
           nzz = zz/radius;
           
           f = nxx/(1 - nzz);
           g = nyy/(1 - nzz);
           E1 = imgCropSphere1(x, y);
           E2 = imgCropSphere2(x, y);
           E3 = imgCropSphere3(x, y);

           %Add to matrix
           % fge = [fge; double(x) double(y) double(f) double(g) double(E1) double(E2) double(E3)];
           fge((x - 1) * (h - 1) + y, 1) = double(x);
           fge((x - 1) * (h - 1) + y, 2) = double(y);
           fge((x - 1) * (h - 1) + y, 3) = double(f);
           fge((x - 1) * (h - 1) + y, 4) = double(g);
           fge((x - 1) * (h - 1) + y, 5) = double(E1);
           fge((x - 1) * (h - 1) + y, 6) = double(E2);
           fge((x - 1) * (h - 1) + y, 7) = double(E3);
       else
           % build indeces of empty rows
            indRow(index) = (x - 1) * (h - 1) + y;
            index = index + 1;
       end
   end
end
fge = removerows(fge, indRow);

%% Build Lookup Table indexed by E1/E2 E2/E3
[fgesize, w]  = size(fge);

epsilon = 20;
BinScale = 100;
xsize = 30;
ysize = 30;
LookUpTable = [];
LookUpTable(xsize * BinScale, ysize * BinScale).f = 0;
LookUpTable(xsize * BinScale, ysize * BinScale).g = 0;
AvgFG = nan(xsize * BinScale, ysize * BinScale);


E1E2Vec = [];
E2E3Vec = [];
VVec = [];
for y = 1:fgesize
    E1E2 = ceil((fge(y,5) + 1)/(fge(y,6) + 1) * 10);
    E2E3 = ceil((fge(y,6) + 1)/(fge(y,7) + 1) * 10);
    
    f = fge(y,3);
    g = fge(y,4);
    
    curr = LookUpTable(E1E2, E2E3);
    if (isempty(curr.f))
        % Populate Container
        LookUpTable(E1E2, E2E3).f = f;
        LookUpTable(E1E2, E2E3).g = g;
        AvgFG(E1E2, E2E3) = double((f + g)/2);
        
        % build sample data
        E1E2Vec = [E1E2Vec; E1E2];
        E2E3Vec = [E2E3Vec; E2E3];
        VVec = [VVec;  double((f + g)/2)];
    else
        % Check values and averge or new spot
        N = 33;
    end
    
end




