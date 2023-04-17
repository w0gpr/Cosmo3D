%% JakCosmo
% This script reads the shapefiles output from QGIS to get the volume data.
% clear

%% Read in the data and organize
filename = '~/data/vector/BlockVolume3.shp';

data = shaperead(filename);
[numID, ~] = size(data);

percentError = 0.3;

%% Convert data to matrix
polygons = zeros(numID,4);
for i = 1:numID
    polygons(i,1) = data(i).fid;
    polygons(i,2) = data(i).BlockID;
    polygons(i,3) = data(i).volume;
    polygons(i,4) = data(i).area;
end
polygons = [polygons (polygons(:,2)-1)/2 + 1];
minVolume = sum(polygons(:,3))
% maxVolume = minVolume * 2
% Area = sum(polygons(:,4))

%% Range of volumes

n = 1e5;
% Volume = zeros(n,1);
quarryVolume = polygons(:,3).* polygons(:,5)...
    + polygons(:,3) .* percentError .* polygons(:,5) .* randn(numID,n);
quarryVolumeSum = sum(quarryVolume);

% figure(1)
% clf
% hist(quarryVolumeSum)
% hold on


%% Buffered
filename = '~/data/vector/BlockVolume-05Clip.shp';

data = shaperead(filename);
[numID, ~] = size(data);

% percentError = 0.5;

%% Convert data to matrix
polygons_buf = zeros(numID,4);
for i = 1:numID
    polygons_buf(i,1) = data(i).fid;
    polygons_buf(i,2) = data(i).BlockID;
    polygons_buf(i,3) = data(i).volume;
    polygons_buf(i,4) = data(i).area;
end
polygons_buf = [polygons_buf (polygons_buf(:,2)-1)/2 + 1];
minVolume_buf = sum(polygons_buf(:,3))
% maxVolume_buf = minVolume_buf * 2
% Area_buf = sum(polygons_buf(:,4))



%% Range of volumes

% n = 1e5;
% Volume = zeros(n,1);
quarryVolume_buf = polygons_buf(:,3).* polygons_buf(:,5)...
    + polygons_buf(:,3) .* percentError .* polygons_buf(:,5) .* randn(numID,n);
quarryVolumeSum_buf = sum(quarryVolume_buf);

%% Full Range of volumes

n=5e5;

volRange = [floor(polygons(:,3)) ceil(polygons_buf(:,3))] * 1000;
volRange = volRange + 1;
volRange = sort(volRange')';
volRand = zeros(numID,n);
for i = 1:numID
    volRand(i,:) = randi([volRange(i,1),volRange(i,2)],1,n)./1000;
end
quarryVolume_all = volRand .* polygons_buf(:,5) + volRand .* percentError .* polygons_buf(:,5) .* randn(numID,n);
quarryVolumeSum_all = sum(quarryVolume_all);

%% 
mean(sum(volRand))
std(sum(volRand))
mean(quarryVolumeSum_all)
std(quarryVolumeSum_all)



figure(2)
clf
h1 = histogram(sum(volRand));
h1.BinWidth = 10;
hold on
h2 = histogram(quarryVolumeSum_all);
h2.BinWidth = 10;




%% 
figure(1)
clf
h1 = histogram(quarryVolumeSum);
h1.BinWidth = 10;
hold on
h2 = histogram(quarryVolumeSum_buf);
h2.BinWidth = 10;
h3 = histogram(quarryVolumeSum_all);
h3.BinWidth = 10;

mean(quarryVolumeSum);
mean(quarryVolumeSum_buf);
mean(quarryVolumeSum_all)
std(quarryVolumeSum_all)

% 
% %% Read in area shapefile
% filename2 = 'file:///home/brandon/GIS/Greenland/JakCosmo/data/vector/AbrasionArea.shp';
% 
% abrasionArea = shaperead(filename2);
% abrasionArea = abrasionArea.area;
% 
% aD = [4.1 1.9]; % abrasion depth calculated from field site in cm
% % abrasionDepth = [aD(1)-aD(2)*2 aD(1) aD(1)+aD(2)*2];    % cm
% abrasionDepth = aD(1) + aD(2)*randn(1,n);
% abrasionVolume = abrasionDepth/100*abrasionArea;    % m^3
% time = 200;
% 
% totalVolume = abrasionVolume + quarryVolumeSum;
% erosionRate = totalVolume/time/abrasionArea*1e3;
% % abrasionVolume = uint16(abrasionDepth/100*abrasionArea) 
% 
% %% results
% rangeValues = [.025 .50 .975];
% qVolume = quantile(quarryVolumeSum,rangeValues);
% aVolume = quantile(abrasionVolume,rangeValues);
% tVolume = quantile(totalVolume,rangeValues);
% eRate = quantile(erosionRate,rangeValues);
% figure(1)
% clf
% subplot(1,2,1)
% h1 = histogram(quarryVolumeSum);
% hold on
% h2 = histogram(abrasionVolume);
% 
% % h1.Normalization = 'count';
% % h1.BinWidth = 5;
% % h2.Normalization = 'count';
% % h2.BinWidth = 50;
% 
% % subplot(2,2,3)
% % histogram(totalVolume)
% subplot(1,2,2)
% histogram(erosionRate)
% 
% Ratio = aVolume(2)/qVolume(2);
