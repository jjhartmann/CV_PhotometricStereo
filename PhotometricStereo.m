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
[cx, cy] = ginput(1)
[px, py] = ginput(1)

% Find the radius of circle
radius = sqrt((cx - px)^2 + (cy - py)^2);

% Draw circle on image
theta = 0 : (2 * pi /10000) : (2 * pi);
pline_x = radius * cos(theta) + cx;
pline_y = radius * sin(theta) + cy;
hold on; 
plot(pline_x, pline_y)
hold off;

% Crop the image based on circle. 
dim = [(cx - radius - 10) (cy - radius - 10) (2*radius + 10) (2*radius + 10)];
imgCropSphere = imcrop(imgSphereGray, dim);
figure
imshow(imgCropSphere)

% Get new center
[h, w] = size(imgCropSphere);
ccx = w/2;
ccy = h/2;
%% Collect fg value pairs and E1, E2, and E3. 


