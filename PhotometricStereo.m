%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assingment 2 - Photostereo Imaging

%% Init calabration Sphere. 
imgSphere = imread('Photostereo_RealImages/sphere-lamp1.tif');

%Convert to gray
imgSphereGray = rgb2gray(imgSphere);
imshow(imgSphereGray)
[x,y] = ginput(1)

