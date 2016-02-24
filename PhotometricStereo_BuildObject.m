%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assingment 2 - Photostereo Imaging: Building 3D Objects. 

% Load TableLookUp
% load('LookUpTable.mat');
% load ('radius.mat');
% Load image data
cylinderdata = {'Photostereo_RealImages/cylinder-lamp1.tif'; 'Photostereo_RealImages/cylinder-lamp2.tif'; 'Photostereo_RealImages/cylinder-lamp3.tif'};
hexlightdata = {'Photostereo_RealImages/hex1-lamp1.tif'; 'Photostereo_RealImages/hex1-lamp2.tif'; 'Photostereo_RealImages/hex1-lamp3.tif'};
hexlight2data = {'Photostereo_RealImages/hex2-lamp1.tif'; 'Photostereo_RealImages/hex2-lamp2.tif'; 'Photostereo_RealImages/hex2-lamp3.tif'};
ellipsoiddata = {'Photostereo_RealImages/ellipsoid-lamp1.tif'; 'Photostereo_RealImages/ellipsoid-lamp2.tif'; 'Photostereo_RealImages/ellipsoid-lamp3.tif'};
spheredata = {'Photostereo_RealImages/sphere-lamp1.tif'; 'Photostereo_RealImages/sphere-lamp2.tif'; 'Photostereo_RealImages/sphere-lamp3.tif'};
coneLightdata = {'Photostereo_RealImages/cone-lamp1.tif'; 'Photostereo_RealImages/cone-lamp2.tif'; 'Photostereo_RealImages/cone-lamp3.tif'};
conedarkdata = {'Photostereo_RealImages/cone2-lamp1.tif'; 'Photostereo_RealImages/cone2-lamp2.tif'; 'Photostereo_RealImages/cone2-lamp3.tif'};


img1 = rgb2gray(imread(hexlightdata{1}));
img2 = rgb2gray(imread(hexlightdata{2}));
img3 = rgb2gray(imread(hexlightdata{3}));

%% Build 3D mesh

% Iterate over all three images. 
[h, w] = size(img1);
[lw, lw] = size(LookUpTable);
th = 35;
BinScale = 100; % TODO: Create Global Static Vars to share. 

TDMap = [];
P = [];
Q = [];
Z = [];
for i = 1:h
   for j = 1:w
       
      val = max([img1(i, j), img2(i, j), img2(i, j)]);
      if (val > th)
         % Process pixel. 
         E1 = img1(i, j);
         E2 = img2(i, j);
         E3 = img3(i, j); 
         
         % Create index for lookup table
         E1E2 = ceil((log(double(E1 + 1)/double(E2 + 1)) + 5) * BinScale);
         E2E3 = ceil((log(double(E2 + 1)/double(E3 + 1)) + 5) * BinScale);
         
         % Get previous values
         previ = E1E2;
         prevj = E2E3;
         if (previ ~= 1)
             previ = i - 1;
         elseif (prevj ~= 1)
             prevj = j - 1;
         end
         
         prevP = P(previ, prevj);
         prevQ = Q(previ, prevj);
         
         %% Search lookup table
         f = LookUpTable(E2E3, E1E2).f;
         g = LookUpTable(E2E3, E1E2).g;
         
         % Find similar value
         [ft, fs] = size(f);
         % TODO: FIND APPROPRIATE VALUE
         
         x = ceil((((2 * f)/(1 + f^2 + g^2)) * radius));
         y = ceil((((2 * g)/(1 + f^2 + g^2)) * radius));
         z = ceil(((-1 + f^2 + g^2)/(1 + f^2 + g^2)) * radius);
         
         P(i, j) = double(x/z);
         Q(i, j) = double(y/z);
         Z(i, j) = 0;
         
         TDMap(i, j) = z;
         
      end
   end
end


% Integrate along multiple paths. 
[h, w] = size(P);
P2 = zeros(size(P));
Q2 = zeros(size(Q));
Z2 = zeros(size(P));
for i = 1:h
    
  if (mod(i, 2) == 1)
       for j = 1:w
           prevj = j;
           previ = i;
           if (j ~= 1)
               prevj = j - 1;
           end
           
           if (j == w && i ~= 1)
               previ = i - 1;
           end
           
           tmp = P(i, j);
           P2(i, j) =  ((P(previ, prevj) + P2(i, j))/2) + tmp;
           
           tmpQ = Q(i, j);
           Q2(i, j) =  ((Q(previ, prevj) + Q2(i, j))/2) + tmpQ;
           
           tmpZ = Z2(i, j);
           Z2(i, j) =  ((P2(i, j) + Q2(i, j))/2) + tmpZ;
       end
  else
       for j = w:1
           prevj = j;
           previ = i;
           if (j ~= w)
               prevj = j + 1;
           end
           
           if (j == 1 && i ~= 1)
               previ = i - 1;
           end
           
           tmp = P(i, j);
           P2(i, j) =  ((P(previ, prevj) + P2(i, j))/2) + tmp;
           
           tmpQ = Q(i, j);
           Q2(i, j) =  ((Q(previ, prevj) + Q2(i, j))/2) + tmpQ;
           
           tmpZ = Z2(i, j);
           Z2(i, j) =  ((P2(i, j) + Q2(i, j))/2) + tmpZ;
       end
  end
end


% Other way
for j = 1:w
    
  if (mod(j, 2) == 1)
       for i = 1:h
           prevj = j;
           previ = i;
           if (j ~= 1)
               prevj = j - 1;
           end
           
           if (j == w && i ~= 1)
               previ = i - 1;
           end
           
           tmp = P(i, j);
           P2(i, j) =  ((P(previ, prevj) + P2(i, j))/2) + tmp;
           
           tmpQ = Q(i, j);
           Q2(i, j) =  ((Q(previ, prevj) + Q2(i, j))/2) + tmpQ;
           
           tmpZ = Z2(i, j);
           Z2(i, j) =  ((P2(i, j) + Q2(i, j))/2) + tmpZ;
       end
  else
       for i = h:1
           prevj = j;
           previ = i;
           if (j ~= w)
               prevj = j + 1;
           end
           
           if (j == 1 && i ~= 1)
               previ = i - 1;
           end
           
           tmp = P(i, j);
           P2(i, j) =  ((P(previ, prevj) + P2(i, j))/2) + tmp;
           
           tmpQ = Q(i, j);
           Q2(i, j) =  ((Q(previ, prevj) + Q2(i, j))/2) + tmpQ;
           
           tmpZ = Z2(i, j);
           Z2(i, j) =  ((P2(i, j) + Q2(i, j))/2) + tmpZ;
       end
  end
end