function [out] = JakCosmo3DFull(xP,zP,LambdaP,abrasionDepth)
%% This is a function to run simulations through different inputs
tic
rng('shuffle')

% Set the constants
numRays = 1000;
lambda = 4.99e-7; % Be-10 decay rate
tExpose = 7200; % duration of exposure - years, 
% tBurial = 200; % years
rho = 2.65;      % sample density
% LambdaP = 208;  % g/cm^2 from Gosse and Phillips 2001
% LambdaM = 1145; % double check this... this is from memory for testing

%% Read in sample data
% This section reads in the data from the .csv and then rotates the datum.
% It is done in a way that is dynamic and allows for buffers to be
% increased or decreased
data = readtable('Camp3Samples.csv');
numSamples = 8;
D = data{1:8,2:end};
abradePluckContact = 340;   % cm from base
wallHeight = 120;           % cm

% transform the datum from upper right to lower left
aPC = (abradePluckContact-max(max(D(:,5:6))))*-1;
D(:,5:6) = D(:,5:6)-max(max(D(:,5:6)));
D(:,7:8) = D(:,7:8)-wallHeight;
D(:,5:8) = D(:,5:8)*-1;

% This section breaks the data into more recognizable variable names
conc = D(:,1);
sigma = D(:,2);
Pspal = D(1,3);
Pmuon = D(1,4);
xSamp = D(:,5:6);
zSamp = D(:,7:8);

% This section is for expanding the sample domain and to add buffers around
% the sample location to allow for rays (high energy neutron) to be
% projected through. 300 cm (3 m) is initially used due to neutrons being
% effectively completely attenuated at these distances, however the code is
% such that the buffer size doesn't effect the processing time, so it could
% be expanded without increasing calculation time.
xBuf = 300;     % cm - buffer around data
zBuf = 100;     % cm - buffer around data
xSamp = xSamp + xBuf;
zSamp = zSamp + zBuf;
dT2 = zBuf+wallHeight;
% height = dT2+zBuf;       % cm
width = max(xSamp(:))+xBuf;   % cm

XSamp = [min(xSamp,[],2) max(xSamp,[],2)];
ZSamp = [min(zSamp,[],2) max(zSamp,[],2)];

% This describes the 'pinning' points used to create the surface model
xPreset = [0 xBuf+aPC xSamp(1,1) width];
zPreset = [zBuf zBuf dT2 dT2];

%% Build the surface model

rLmean = zeros(1,numSamples);

[X1,Y1,Z1] = surfaceModelSingle(abrasionDepth,xP,zP,xPreset,zPreset);
tri = delaunay(X1,Y1);
[numT,~] = size(tri);

%% cycle through the samples and ray traces
for ii = 1:numSamples
    % initialize the ray input
    if numRays == 1
        theta = 0; phi = 0;
    else
        [theta,phi] = cosmoDistribution(numRays);
    end
    rL = zeros(numRays,1);
    rad = 1.2e3;    % the radius or ray distance
    Q1 = [XSamp(ii,1)+(XSamp(ii,2)-XSamp(ii,1)).*rand(numRays,1),...
        rL,...
        ZSamp(ii,1)+(ZSamp(ii,2)-ZSamp(ii,1)).*rand(numRays,1)];
    Q2 = [rad.*sin(theta).*cos(phi)+Q1(:,1),...
            rad.*sin(theta).*sin(phi)+Q1(:,2),...
            rad.*cos(theta)+Q1(:,3)];

    for k = 1:numRays
        Q4 = Q2(k,:);
        Q3 = Q1(k,:);
        X = zeros(numT,1); Y = zeros(numT,1); Z = zeros(numT,1); 
        rays = zeros(numT,1);
        n=0;
        surfaceSample = false;
        for j = 1:numT
            P1=[X1(tri(j,1)) Y1(tri(j,1)) Z1(tri(j,1))];
            P2=[X1(tri(j,2)) Y1(tri(j,2)) Z1(tri(j,2))];
            P3=[X1(tri(j,3)) Y1(tri(j,3)) Z1(tri(j,3))];
            N = cross(P2-P1,P3-P1); % Normal to the plane of the triangle
            P0 = Q3 + dot(P1-Q3,N)/dot(Q4-Q3,N)*(Q4-Q3); % The point of intersection
            if dot(cross(P2-P1,P0-P1),N)>=0 && ...  % Is P0 inside triangle?
               dot(cross(P3-P2,P0-P2),N)>=0 && ...
               dot(cross(P1-P3,P0-P3),N)>=0 && ...
               P0(3)-Q3(3)>=-0.0001    % Does the ray project below the sample?
                n = n + 1;
                X(n) = P0(1); Y(n) = P0(2); Z(n) = P0(3); % If so, store P0
                temp = P0-Q3;
                temp = sqrt(temp*temp');
                rays(n) = temp;
            end
            if P0(3)-Q3(3) == 0
                surfaceSample = true;
            end
        end

        % This removes duplicates created due to intersections with lines
        rays = unique(rays(rays>0),'rows'); 
        if surfaceSample
            rays = [0;rays];
        end
        r = length(rays);
        rayLength = rays(1);
        if r>2
            for i = 2:2:r
                rayLength = rayLength + rays(i+1)-rays(i);
            end
        end
        rL(k) = rayLength;
    end
    rLmean(ii) = mean(rL);
end

%% Calculate the shielding factor, concentration, chi^2 reduction
S = exp(-rLmean*rho/LambdaP);
concSolve = Pspal*S/lambda*(1-exp(-lambda*tExpose))...
    + Pmuon/lambda*(1-exp(-lambda*tExpose));
M = sum(((conc-concSolve')./sigma).^2);

%% Generate outputs
% out.X1 = X1;
% out.Y1 = Y1;
% out.Z1 = Z1;
% out.concSolve = concSolve;
% out.M = M;

out = M;

toc
end
