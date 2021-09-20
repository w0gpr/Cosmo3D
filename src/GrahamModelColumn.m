%% GrahamModelColumn

% This is to determine the best fit of the vertical column of samples as if
% they were taken from a core. I am assuming the base sample is inline with
% the column and not offset by a few cm.

clear
 
% rng('shuffle')

data = readtable('Camp3Samples.csv');
sampleID = data{1:8,1};
D = data{1:8,2:end};


% This section breaks the data into more recognizable variable names, sets
% constants, and a priori knowledge of parameters
conc = D(1:5,1);
sigma = D(1:5,2);
Pspal = D(1,3);
Pmuon = D(1,4);
xSamp = D(1:5,5:6);
zSamp = D(1:5,7:8);
lambda = 4.99e-7; % decay rate
tExpose = 7000:100:8000; % duration of exposure - years, 
% tBurial = 500; % years
rho = 2.65;      % sample density
LambdaP = 160;  % g/cm^2 from Gosse and Phillips 2001
LambdaM = 1150; % double check this... this is from memory for testing
eD = 1;
z = 0:5:max(zSamp);
S = exp(-(z+eD)'*rho./[LambdaP LambdaM]);
concSolve = Pspal*S(:,1)/lambda*(1-exp(-lambda*tExpose))...
    + Pmuon*S(:,2)/lambda*(1-exp(-lambda*tExpose));

figure(1)
clf
% h1 = axes;
% set(h1,'Ydir','normal')
errorbar(conc,mean(zSamp,2),sigma,'horizontal')
% plot(conc,mean(zSamp,2))
set(gca, 'YDir', 'reverse')
hold on
plot(concSolve,z')