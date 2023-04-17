%% individualShieldingDepth
% This is a utility to calculate the estimated shielding factor versus the measured skyline shielding factor for a sample, based on the in field measurement of shielding via azimuth/elevation and the measured concentration based on a known expsoure history. This can be used to directly estimate the depth of erosion for a sample.


clear
dataTable = readtable('Camp3Samples.csv');
data.conc = dataTable{:,2};
data.sigma = dataTable{:,3};
data.P = dataTable{1,4} + dataTable{1,5};
conc = [data.conc-data.sigma*2 data.conc data.conc+data.sigma*2];

%% 
Lambda = 160;
lambda = 4.99e-7; % Be-10 decay rate
tExpose = 7400; % duration of exposure - years, 
tBurial = 220; % years
rho = 2.65;      % sample density

% SF = (data.conc*lambda)/(data.P*(1-exp(-lambda*(tExpose-tBurial))));
SFall = (conc*lambda)/(data.P*(1-exp(-lambda*(tExpose-tBurial)))); % The calculated shielding factor based on the measured concentration.
  
temp = -log(SFall)*Lambda/rho;

depth = [temp (temp(:,1)-temp(:,3))/2];

%%  
ER2B = [0 70
    99 70
    100 30
    120 15
    150 0
    330 0
    350 7
    355 30];

SFER2B = skyline(ER2B(:,1),ER2B(:,2));

temp = -log(1-(SFER2B - SFall(end,:)))*Lambda/rho;

depthSFER2B = [temp (temp(:,1)-temp(:,3))/2];
