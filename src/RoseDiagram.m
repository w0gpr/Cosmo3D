%% Rose Diagram program
% This program is written to visualize geologic data as a rose diagram or
% polar coordinate histogram.  The specific application for originally
% writing this script was to plot data from an outcropping surface
% recently exposed in Greenland of glacial flow indicators such as
% striations, single chatter marks, chatter mark trains, and varying levels
% of "freshness" of mark, assumed created during the Little Ice Age
% readvance.
clear
% load data
data = importdata('JAK2018Camp3Flow.csv');
bins = 24; % 15 degree bins = 24 bins, 10 deg = 36,
rad1 = 20;
rad2 = 15;
% test = data;
edges = deg2rad(0:360/bins:360);
bins = edges;
data(data(:,1)>270,1) = data(data(:,1)>270,1) - 180; % This transforms data to be in southern orientation.

%% This sorts the data to the various ages determined by field "freshness" or certainty
age = sortrows(data,2);
age1 = age(age(:,2)==1,1);
age2 = age(age(:,2)==2,1);
age3 = age(age(:,2)==3,1);
age4 = age(age(:,2)==4,1);      % Indicates picked from SfM DEM

%% This sorts the data by the type of mark
mark = sortrows(data,3);
Chat = mark(mark(:,3)==1,1);
ChatTrain = mark(mark(:,3)==2,1);
Striae = mark(mark(:,3)==3,1);
Gouge = mark(mark(:,3)==4,1);
Pluck = mark(mark(:,3)==5,1);   % Picked from the SfM DEM
Questionable = mark(mark(:,3)==0,1);
Joint = mark(mark(:,3)==6,1);   % Picked from the SfM DEM
Data = data;
data = data(data(:,3)~=6,1);

%% Figures
figure (1)
clf
% baseline(Data(:,1),bins,0)
PHist(Data(:,1),bins,0)
% title 'All Data'



%% Part 2
figure (2)
clf

subplot(2,2,1)
baseline(data(:,1),bins,rad2)
PHist(Striae,bins,rad2)
hold off
% title 'A. Striations'

subplot(2,2,2)
baseline(data(:,1),bins,rad2)
PHist([Chat;ChatTrain],bins,rad2)
hold off
% title 'B. Chatter Marks'

subplot(2,2,3)
baseline(data(:,1),bins,rad2)
PHist(Gouge,bins,rad2)
hold off
% title 'C. Gouges'

subplot (2,2,4)
baseline(data(:,1),bins,rad2)
PHist(Pluck,bins,rad2)
hold off
% title 'D. Quarried Cliffs'


hold off

function PHist(data,bins,Radius)
polarhistogram(deg2rad(data),bins)
ax = gca;
ax.ThetaTick = 0:45:360;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaTickLabels = {'N','NE','E','SE','S','SW','W','NW'};
if Radius>0
    ax.RLimMode = 'manual';
    ax.RLim = [0 Radius];
end
end

function baseline(data,bins,Radius)
polarhistogram(deg2rad(data),bins,'DisplayStyle','stairs')
ax = gca;
ax.ThetaTick = 0:45:360;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaTickLabels = {'N','NE','E','SE','S','SW','W','NW'};
hold on
if Radius>0
    ax.RLimMode = 'manual';
    ax.RLim = [0 Radius];
end
end


