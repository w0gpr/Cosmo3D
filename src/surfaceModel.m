function [X1,Y1,Z1] = surfaceModel(i,abrasionDepth,xL,zL,xPreset,zPreset)

%% Surface model builder
% This script builds the surface models that are then output to the main
% program to be read into the ray trace model. There are three initial end
% members, and then eventually intermediate scenarios to solve for the
% minimized chi^2 reduction. 

buf = 300;
y1 = [-buf-1 -buf buf buf+1];

% switch i
%     case 1  % Full Block
%         x1 = [-1 0 481 482 1126 1127];
%         z1 = [0 100 100 220 220 0];
%     case 2  % Simple Triangular Block
%         x1 = [-1 0 481 821 1126 1127];
%         z1 = [0 100 100 220 220 0];
%     case 3  % No Block
%         x1 = [-1 0 820 821 1126 1127];
%         z1 = [0 100 100 220 220 0];
% end

x = linspace(xPreset(2)+1,xPreset(3)-1,xL);
z = linspace(zPreset(3),zPreset(2),zL);

[xx,zz] = meshgrid(x,z);

% runs = numel(xx)
% x1 = [-1 0 481 xx(i) 821 1126 1127];
% z1 = [0 100 100 zz(i) 220 220 0];
x1 = [-1 xPreset(1:2) xx(i) xPreset(3:end) xPreset(end)+1];
z1 = [0 zPreset(1:2) zz(i) zPreset(3:end) 0];

xq = [x1(1:3) x x1(5:end)];
zq = pchip(x1,z1,xq);
% z = [z1(1:3) zq z1(5:end)];

[X1,Y1] = meshgrid(xq,y1);
Z1 = zeros(size(X1));
Z1(2:3,:) = [zq + abrasionDepth; zq + abrasionDepth];
