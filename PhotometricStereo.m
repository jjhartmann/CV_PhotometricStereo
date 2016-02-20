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
fge1 = double([]);
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
           E = imgCropSphere1(x, y);
           tuple = [];
           
           %Add to matrix
           fge1 = [fge1; double(x) double(y) double(f) double(g) double(E)];
       end
   end
end


% Sphere 2: Iterate over inter image x and y coordinates and build f, g and E values
fge2 = double([]);
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
           E = imgCropSphere2(x, y);
           tuple = [];
           
           %Add to matrix
           fge2 = [fge2; double(x) double(y) double(f) double(g) double(E)];
       end
   end
end


% Sphere 3: Iterate over inter image x and y coordinates and build f, g and E values
fge3 = double([]);
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
           E = imgCropSphere3(x, y);
           tuple = [];
           
           %Add to matrix
           fge3 = [fge3; double(x) double(y) double(f) double(g) double(E)];
       end
   end
end

%% Build Lookup Table indexed by E1/E2 E2/E3
[fgesize w]  = size(fge1);

epsilon = 20;
BinScale = 10;
LookUpTable(255 * BinScale, 255 * BinScale).f = 0;
LookUpTable(255 * BinScale, 255 * BinScale).g = 0;
AvgFG = double(zeros(255 * BinScale, 255 * BinScale));

for y = 1:fgesize
    E1E2 = ceil((fge1(y,5) + 1)/(fge2(y,5) + 1) * BinScale);
    E2E3 = ceil((fge2(y,5) + 1)/(fge3(y,5) + 1) * BinScale);
    
    f = fge1(y,3);
    g = fge1(y,4);
    
    curr = LookUpTable(E1E2, E2E3);
    if (isempty(curr.f))
        % Populate Container
        LookUpTable(E1E2, E2E3).f = f;
        LookUpTable(E1E2, E2E3).g = g;
        AvgFG(E1E2, E2E3) = double((f + g)/2);
    else
        % Check values and averge or new spot
        N = 33;
    end
    
end




