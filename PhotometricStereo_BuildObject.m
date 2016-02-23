%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assingment 2 - Photostereo Imaging: Building 3D Objects. 

% Load TableLookUp
% load('LookUpTable.mat');
% load ('radius.mat');
% Load image data
cylinderdata = {'Photostereo_RealImages/cylinder-lamp1.tif'; 'Photostereo_RealImages/cylinder-lamp2.tif'; 'Photostereo_RealImages/cylinder-lamp3.tif'};
hexlightdata = {'Photostereo_RealImages/hex1-lamp1.tif'; 'Photostereo_RealImages/hex1-lamp2.tif'; 'Photostereo_RealImages/hex1-lamp3.tif'};
ellipsoiddata = {'Photostereo_RealImages/ellipsoid-lamp1.tif'; 'Photostereo_RealImages/ellipsoid-lamp2.tif'; 'Photostereo_RealImages/ellipsoid-lamp3.tif'};
spheredata = {'Photostereo_RealImages/sphere-lamp1.tif'; 'Photostereo_RealImages/sphere-lamp2.tif'; 'Photostereo_RealImages/sphere-lamp3.tif'};

img1 = rgb2gray(imread(spheredata{1}));
img2 = rgb2gray(imread(spheredata{2}));
img3 = rgb2gray(imread(spheredata{3}));

%% Build 3D mesh

% Iterate over all three images. 
[h, w] = size(img1);
th = 40;
BinScale = 30; % TODO: Create Global Static Vars to share. 

TDMap = [];
for i = 1:h
   for j = 1:w
       
      val = max([img1(i, j), img2(i, j), img2(i, j)]);
      if (val > th)
         % Process pixel. 
         E1 = img1(i, j);
         E2 = img2(i, j);
         E3 = img3(i, j); 
         
         % Create index for lookup table
         E1E2 = ceil(double(E1 + 1)/double(E2 + 1)) * BinScale;
         E2E3 = ceil(double(E2 + 1)/double(E3 + 1)) * BinScale;
         
         % Search lookup table
         f = LookUpTable(E2E3, E1E2).f(1);
         g = LookUpTable(E2E3, E1E2).g(1);
         
         x = ceil((((2 * f)/(1 + f^2 + g^2)) * radius) + ccx);
         y = ceil((((2 * g)/(1 + f^2 + g^2)) * radius) + ccy);
         z = ceil(((-1 + f^2 + g^2)/(1 + f^2 + g^2)) * radius);
         
         TDMap(i, j) = z;
         
      end
   end
end