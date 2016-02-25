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
th = 30;
BinScale = 100; % TODO: Create Global Static Vars to share. 

TDMap = [];
P = [];
Q = [];
Z = [];
prevE1E2 = 0;
prevE2E3 = 0;
prevf = 0;
prevg = 0;
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
%          if (prevf > 0 && prevg > 0)
%             prevf = LookUpTable(E2E3, E1E2).f;
%             prevg = LookUpTable(prevE2E3, prevE1E2).g;
%          end
         
         %% Search lookup table
         f = LookUpTable(E2E3, E1E2).f;
         g = LookUpTable(E2E3, E1E2).g;
         
         % Find similar value
         [ft, fs] = size(f);
         if (fs > 1)
            deltaf1 = abs(f(1) - prevf);
            deltaf2 = abs(f(2) - prevf);
            if (deltaf1 < deltaf2)
                f = f(1);
            else
                f = f(2);
            end
         end
         
         [gt, gs] = size(g);
         if (gs > 1)
            deltag1 = abs(g(1) - prevg);
            deltag2 = abs(g(2) - prevg);
            if (deltag1 < deltag2)
                g = g(1);
            else
                g = g(2);
            end
         end
         
         %% Build p and q
         x = ceil((((2 * f)/(1 + f^2 + g^2)) * radius));
         y = ceil((((2 * g)/(1 + f^2 + g^2)) * radius));
         z = ceil(((-1 + f^2 + g^2)/(1 + f^2 + g^2)) * radius);
         
         if (z < 1)
             z = abs(z);
         end
         
         P(i, j) = double(x/z);
         Q(i, j) = double(y/z);
         Z(i, j) = z;
         
         TDMap(i, j) = z;
         
         % set Previous
         prevE1E2 = E1E2;
         prevE2E3 = E2E3;
         prevf = f;
         prevg = g;
         
      end
   end
end


%% EXPERIMENTS: Integrate along multiple paths. 
[h, w] = size(P);
P2 = zeros(size(P));
Q2 = zeros(size(Q));
Z2 = zeros(size(P));
previ = 0;
prevj = 0;
for i = 1:h
    for j = 1:w
        tmp = P(i, j);
        tmpPrevP = 0;
        tmpPrevZ = 0;
        
        if(previ >= 1 && prevj > 1)
            tmpPrevP = P2(previ, prevj);
            tmpPrevZ = Z2(previ, prevj);
        else
            tmp = 0;
        end
        
        % Check for inf or nans
        if (tmp == -inf || isnan(tmp) || tmp == inf)
           tmp = 0; 
        end
        P2(i, j) =  tmpPrevP - tmp;
        Z2(i, j) =  double(tmpPrevZ + tmpPrevP)/2 - tmp;
        
        prevj= j;
        previ = i;
    end
end

previ = 0;
prevj = 0;
P3 = zeros(size(P));
for i = 1:h
    for j = w:-1:1
        tmp = P(i, j);
        tmpPrevP = 0;
        tmpPrevZ = 0;
        if(previ >= 1 && prevj > 1)
            tmpPrevP = P3(previ, prevj);
            tmpPrevZ = Z2(previ, prevj);
        else
            tmp = 0;
        end
        
        % Check for inf or nans
        if (tmp == -inf || isnan(tmp) || tmp == inf)
           tmp = 0; 
        end
        P3(i, j) =  tmpPrevP + tmp;
        Z2(i, j) =  double(tmpPrevZ + tmpPrevP)/2 + double(tmp + Z2(i, j))/2;
        
        prevj = j;
        previ = i;
    end
end


%%%%%%%%%% Q
previ = 0;
prevj = 0;
for j = 1:w
    for i = 1:h
        tmp = Q(i, j);
        tmpPrevQ = 0;
        if(previ > 1 && prevj > 1)
            tmpPrevQ = Q2(previ, prevj);
        end
        Q2(i, j) =  tmpPrevQ + tmp;
        
        prevj= j;
        previ = i;
    end
end


%% Interpolate

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