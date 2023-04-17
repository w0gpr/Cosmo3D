function [X1,Y1,Z1] = surfaceModelSingle(abrasionDepth,xP,zP,xPreset,zPreset)

%% Surface model builder
% This script builds the surface models that are then output to the main
% program to be read into the ray trace model. 

buf = 300;  % This is 3 meters of buffer on either side of the sample.
y1 = [-buf-1 -buf buf buf+1];

% convert from percentage to a usable range in the model
xP = xPreset(2)+(xPreset(3)-xPreset(2))*xP; 
zP = zPreset(2)+(zPreset(3)-zPreset(2))*zP;

x = linspace(xPreset(2)+1,xPreset(3)-1,25);
% z = linspace(zPreset(3),zPreset(2),zP);

% [xx,zz] = meshgrid(x,z);

% runs = numel(xx)
% x1 = [-1 0 481 xx(i) 821 1126 1127];
% z1 = [0 100 100 zz(i) 220 220 0];
x1 = [-1 xPreset(1:2) xP xPreset(3:end) xPreset(end)+1];
z1 = [0 zPreset(1:2) zP zPreset(3:end) 0];

xq = [x1(1:3) x x1(5:end)];
zq = pchip(x1,z1,xq);
% z = [z1(1:3) zq z1(5:end)];

[X1,Y1] = meshgrid(xq,y1);
Z1 = zeros(size(X1));
Z1(2:3,:) = [zq + abrasionDepth; zq + abrasionDepth];
