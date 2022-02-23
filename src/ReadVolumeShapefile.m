%% JakCosmo
% This script reads the shapefiles output from QGIS to get the volume data.
clear

%% Read in the data and organize
filename = 'file:///home/brandon/GIS/Greenland/JakCosmo/data/vector/BlockVolume3.shp';

data = shaperead(filename);
[numID, ~] = size(data);

%% Convert data to matrix
polygons = zeros(numID,4);
for i = 1:numID
    polygons(i,1) = data(i).fid;
    polygons(i,2) = data(i).BlockID;
    polygons(i,3) = data(i).volume;
    polygons(i,4) = data(i).area;
end
polygons = [polygons (polygons(:,2)-1)/2 + 1];
minVolume = sum(polygons(:,3));
maxVolume = minVolume * 2;

%% Range of volumes

n = 1e5;
% Volume = zeros(n,1);
quarryVolume = polygons(:,3).* polygons(:,5)...
    + polygons(:,3) .* 0.3 .* polygons(:,5) .* randn(numID,n);
quarryVolumeSum = sum(quarryVolume);



%% Read in area shapefile
filename2 = 'file:///home/brandon/GIS/Greenland/JakCosmo/data/vector/AbrasionArea.shp';

abrasionArea = shaperead(filename2);
abrasionArea = abrasionArea.area;

aD = [4.1 1.9]; % abrasion depth calculated from field site in cm
% abrasionDepth = [aD(1)-aD(2)*2 aD(1) aD(1)+aD(2)*2];    % cm
abrasionDepth = aD(1) + aD(2)*randn(1,n);
abrasionVolume = abrasionDepth/100*abrasionArea;    % m^3
time = 200;

totalVolume = abrasionVolume + quarryVolumeSum;
erosionRate = totalVolume/time/abrasionArea*1e3;
% abrasionVolume = uint16(abrasionDepth/100*abrasionArea) 

%% results
rangeValues = [.025 .50 .975];
qVolume = quantile(quarryVolumeSum,rangeValues);
aVolume = quantile(abrasionVolume,rangeValues);
tVolume = quantile(totalVolume,rangeValues);
eRate = quantile(erosionRate,rangeValues);
figure(1)
clf
subplot(1,2,1)
h1 = histogram(quarryVolumeSum);
hold on
h2 = histogram(abrasionVolume);

% h1.Normalization = 'count';
% h1.BinWidth = 5;
% h2.Normalization = 'count';
% h2.BinWidth = 50;

% subplot(2,2,3)
% histogram(totalVolume)
subplot(1,2,2)
histogram(erosionRate)

Ratio = aVolume(2)/qVolume(2);