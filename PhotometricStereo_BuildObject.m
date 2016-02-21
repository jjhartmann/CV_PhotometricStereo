%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assingment 2 - Photostereo Imaging: Building 3D Objects. 

% Load TableLookUp
load('LookUpTable.mat');

% Load image data
cylinderdata = {'Photostereo_RealImages/cylinder-lamp1.tif'; 'Photostereo_RealImages/cylinder-lamp2.tif'; 'Photostereo_RealImages/cylinder-lamp3.tif'};

img1 = rgb2gray(imread(cylinderdata{1}));
img2 = rgb2gray(imread(cylinderdata{2}));
img3 = rgb2gray(imread(cylinderdata{3}));


