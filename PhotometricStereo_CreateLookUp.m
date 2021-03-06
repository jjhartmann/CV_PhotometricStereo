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
testmat = zeros(h, w);
for x = 1:w
   for y = 1:h
      
       xx = ceil(x - ccx);
       yy = ceil(y - ccy);
       % check to make sure x,y is in cicle. 
       rtmp = ceil(sqrt(xx^2 + yy^2));
       if (rtmp < radius)
%            if (rtmp < 10)
%                 tt = 4;
%            end
          % Use sterographic projection to detect f and g
           zz = ceil(sqrt(radius^2 - (xx^2 + yy^2)));
           testmat(y, x) = zz;
           
           nxx = xx/radius;
           nyy = yy/radius;
           nzz = zz/radius;
           
           f = nxx/(1.000001 - nzz);
           g = nyy/(1.000001 - nzz);
           E1 = imgCropSphere1(y, x);
           E2 = imgCropSphere2(y, x);
           E3 = imgCropSphere3(y, x);

           % Check for infinit
           if (f == inf || isnan(f) || f == -inf)
              continue; 
           end
           if (g == inf || isnan(g) || g == -inf)
              continue; 
           end
           
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


% TODO: Change look-up table to use log. Map the values between -5 to 5. 
epsilon = 20;
BinScale = 100;
xsize = 10;
ysize = 10;
LookUpTable = [];
LookUpTable(ysize * BinScale, xsize * BinScale).f = [0 0];
LookUpTable(ysize * BinScale, xsize * BinScale).g = [0 0];
AvgFG = double(zeros(ysize * BinScale, xsize * BinScale));

E1E2Vec = [];
E2E3Vec = [];
fv = [];
gv = []; 
for y = 1:fgesize
    
    if (fge(1) == 0 && fge(2) == 0)
        continue;
    end
    
    E1E2 = ceil((log((fge(y,5) + 1)/(fge(y,6) + 1)) + 5) * BinScale);
    E2E3 = ceil((log((fge(y,6) + 1)/(fge(y,7) + 1)) + 5) * BinScale);
    
    f = fge(y,3);
    g = fge(y,4);
    
    curr = LookUpTable(E2E3, E1E2);
    if (isempty(curr.f))
        % Populate Container
        LookUpTable(E2E3, E1E2).f = f;
        LookUpTable(E2E3, E1E2).g = g;
        AvgFG(E2E3, E1E2) = double((f + g)/2);
        
        % build sample data
        E1E2Vec = [E1E2Vec; E1E2];
        E2E3Vec = [E2E3Vec; E2E3];
        fv = [fv;  double(f)];
        gv = [gv;  double(g)];
    else
        %% Check values and averge or new spot
        epsilon = 15;       
        tmpf = LookUpTable(E2E3, E1E2).f(1);
        tmpg = LookUpTable(E2E3, E1E2).g(1);
        deltaf = double(abs(tmpf) - abs(f));
        deltag = double(abs(tmpg) - abs(g));
        
        % check value for f
        if ((f >= 0 && tmpf >= 0 && deltaf < epsilon) || (f <= 0 && tmpf < 0 && deltaf < epsilon) || (abs(f) + abs(tmpf)) < epsilon)
           % Values are similar
            LookUpTable(E2E3, E1E2).f(1) = (tmpf + double(f))/2;
        else
           % store new value
           [t s] = size(LookUpTable(E2E3, E1E2).f);
           if (s < 2)
               LookUpTable(E2E3, E1E2).f(2) = 0;
           end
           
           tmpf = LookUpTable(E2E3, E1E2).f(2);
           LookUpTable(E2E3, E1E2).f(2) = (tmpf + double(f))/2;
        end
           
         % check value for g
        if ((g >= 0 && tmpg >= 0 && deltag < epsilon) || (g <= 0 && tmpg < 0 && deltag < epsilon) || (abs(g) + abs(tmpg)) < epsilon)
           % Values are similar
            LookUpTable(E2E3, E1E2).g(1) = (tmpg + double(g))/2;
        else
           % store new value
          [t s] = size(LookUpTable(E2E3, E1E2).g);
          if (s < 2)
               LookUpTable(E2E3, E1E2).g(2) = 0;
          end
          
          tmpg = LookUpTable(E2E3, E1E2).g(2);
          LookUpTable(E2E3, E1E2).g(2) = (tmpg + double(g))/2;
        end
        
    end
    
end


%% Create Grid and intrpolate sparse matrix. 
[h, w] = size(LookUpTable);
[gridx, gridy] = meshgrid(1:w, 1:h);
interpFV = griddata(E1E2Vec, E2E3Vec, fv, gridx, gridy, 'cubic');
interpGV = griddata(E1E2Vec, E2E3Vec, gv, gridx, gridy, 'cubic');
figure(2)
mesh(interpFV)

% Smooth interpolation with gaussian filters
interpFV = imgaussfilt3(interpFV, 6);
interpGV = imgaussfilt3(interpGV, 6);

%% Fill data into lookup table. 
for i = 1:w
   for j = 1:h
      curr = LookUpTable(j, i);
        if (isempty(curr.f))
         % Fill in with interpolated data
         LookUpTable(j, i).f = interpFV(j, i);
         LookUpTable(j, i).g = interpGV(j, i);
        end
   end
end



