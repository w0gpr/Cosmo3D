clear
 
rng('shuffle')
% rng(4096)

lambda = 4.99e-7; % decay rate
tExpose = 7200; % duration of exposure - years, 
% tBurial = 500; % years
rho = 2.65;      % sample density
LambdaP = 208;  % g/cm^2 from Gosse and Phillips 2001
LambdaM = 1145; % double check this... this is from memory for testing

% tic
% sm = 1;
makeFigures = true; % Make the 3D and 2D models
plotResults = true; % Show the chi^2 results
%% Read in sample data
% This section reads in the data from the .csv and then rotates the datum.
% It is done in a way that is dynamic and allows for buffers to be
% increased or decreased
data = readtable('Camp3Samples.csv');
numSamples = 8;
sampleID = data{1:numSamples,1};
D = data{1:8,2:end};
abradePluckContact = 340;   % cm from base
wallHeight = 120;           % cm
% transform the datum from upper right to lower left
aPC = (abradePluckContact-max(max(D(:,5:6))))*-1;
% aPC = abradePluckContact;
D(:,5:6) = D(:,5:6)-max(max(D(:,5:6)));
D(:,7:8) = D(:,7:8)-wallHeight;
D(:,5:8) = D(:,5:8)*-1;

% This section breaks the data into more recognizable variable names, sets
% constants, and a priori knowledge of parameters
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
height = dT2+zBuf;       % cm
width = max(xSamp(:))+xBuf;   % cm

% This describes the 'pinning' points used to create the surface model
blockXEdges = [xBuf+aPC max(xSamp(1))];
xPreset = [0 blockXEdges(1) blockXEdges(2) width];
zPreset = [zBuf zBuf dT2 dT2];

% This is the closed rectangular cross section of each sample. A y component
% is not needed in this model as the y component does not vary and is set
% to 0.
% rSx = [xSamp(:,1) xSamp(:,2) xSamp(:,2) xSamp(:,1) xSamp(:,1)];
% rSz = [zSamp(:,2) zSamp(:,2) zSamp(:,1) zSamp(:,1) zSamp(:,2)];

% This creates a random distribution of sample points to calculate the
% concentration, This makes a pseudo integration (MC) of the volume of
% rock sample collected.
% nP = 1; % number of points to generate inside each sample volume
% if nP>1
%     XSampLim = min(rSx')+(max(rSx')-min(rSx')).*rand(nP,1);
%     ZSampLim = min(rSz')+(max(rSz')-min(rSz')).*rand(nP,1);
% else    % This finds the centroid of the sample
%     XSampLim = (max(rSx')+min(rSx'))/2;
%     ZSampLim = (max(rSz')+min(rSz'))/2;
% end
XSamp = [min(xSamp,[],2) max(xSamp,[],2)];
ZSamp = [min(zSamp,[],2) max(zSamp,[],2)];

% numSamples = length(rSzPoints);


%% Solve for S parameter of measurements
SpMeasured = (conc*lambda)/(Pspal*(1-exp(-lambda*tExpose)))-(Pmuon/Pspal);
abradedDepthSolved = log(SpMeasured)*-160/rho;

%% Init abrasion depth guess and figure 5
abrasionDepth = 2.75;  % depth of abrasion in cm


%% Build the surface model
% x1 = [-1 0 481 821 1126 1127]; y1 = [-301 -300 300 301];
% z1 = [100 100 220 220];
% [X1,Y1] = meshgrid(x1,y1);
% [a,b] = size(X1);
% Z1 = zeros(size(X1));
% Z1(2:a-1,2:b-1) = [z1+1.5;z1+1.5];
xL = 5; zL = 5;
rLmean = zeros(xL*zL,numSamples);
for jj = 1:xL*zL
    %% Call the surface model
[X1,Y1,Z1] = surfaceModel(jj,abrasionDepth,xL,zL,xPreset,zPreset);
tri = delaunay(X1,Y1);
[numT,~] = size(tri);

if makeFigures == true
    figure(1)
    subplot(xL,zL,jj)
    trisurf(tri,X1,Y1,Z1)
    hold on
    plot(X1(2,2:end-1),Z1(2,2:end-1))
    hold off

    figure(2)
    hold on
    plot(X1(2,2:end-1),Z1(2,2:end-1),'k','LineWidth',0.5)
    hold off
end
%%
numRays = 500;
% numRays = numRays * nP;
%% Compute the intersection of the points to the surfaces
% Make the triangles for the vertices to be 'draped' over the x,y
% coordinates.

% input variables into this section:
% X1, Y1, Z1, Q1, numRays
% output variables from this section:
% rayLength
% Q2 will eventually be calculated by cosmoDistro

%% cycle through the samples and ray traces
% figure(1)
% clf



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
%     if nP == 1
%         Q1 = [XSampLim(ii) 0 ZSampLim(ii)]; 
%         Q3 = [rad.*sin(theta).*cos(phi)+Q1(1),...
%             rad.*sin(theta).*sin(phi)+Q1(2) rad.*cos(theta)+Q1(3)];
%     elseif nP > 1
%         error('n>1, This functionality of multiple points per sample not ready yet. Read comment in code')
% %         Q1 = [rSxPoints(:,ii) zeros(size(rSzPoints(:,ii))) rSzPoints(:,ii)]; 
% %         Q3 = [rad.*sin(theta).*cos(phi)+Q1(1),...
% %             rad.*sin(theta).*sin(phi)+Q1(2) rad.*cos(theta)+Q1(3)];
% %         Q3 = permute(Q3,[numel(Q3),1]);
% % To get the ability of multiple points per sample location running, I need
% % to abstract the ray trace for loop into a function, then be able to loop
% % through the different Q1 points for each ray randomly-evenly distributed.
%     end
    
%     rL = zeros(numRays,1);
    for k = 1:numRays
        Q4 = Q2(k,:);
        Q3 = Q1(k,:);
        X = zeros(numT,1); Y = zeros(numT,1); Z = zeros(numT,1); 
        rays = zeros(numT,1);
        n=0;
%         temp = 0;
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
        %         disp(temp)
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
        X = X(1:r); Y = Y(1:r); Z = Z(1:r);

        if r>2
            for i = 2:2:r
                rayLength = rayLength + rays(i+1)-rays(i);
            end
        end
        
%         figure(13)
%         trisurf(tri,X1,Y1,Z1)
%         hold on
%         if rem(k,10)==0
%             plot3([Q3(1);X;Q4(1)],[Q3(2);Y;Q4(2)],[Q3(3);Z;Q4(3)],'r-*')
%         end

        rL(k) = rayLength;
    end
    rLmean(jj,ii) = mean(rL);

end
% disp(rLmean)
end

%% This plots the incoming rays
% Xray = [Q1(:,1)'; Q2(:,1)'];
% Yray = [Q1(:,2)'; Q2(:,2)'];
% Zray = [Q1(:,3)'; Q2(:,3)'];
% figure(13)
% trisurf(tri,X1,Y1,Z1)
% hold on
% plot3(Xray,Yray,Zray,'r')
% hold off

%% Calculate the shielding factor, concentration, chi^2 reduction
% Pmuon = 0;
S = exp(-rLmean*rho/LambdaP);
concSolve = Pspal*S/lambda*(1-exp(-lambda*tExpose))...
    + Pmuon/lambda*(1-exp(-lambda*tExpose));
M = sum(((conc-concSolve')./sigma).^2);


if plotResults == true
    [~,b] = mink(M,10);
    figure11 = figure(12);
    clf
    axes1 = axes('Parent',figure11,'XTickLabel',sampleID,'XDir','reverse');
    box(axes1,'on');
    hold(axes1,'all');
    % Create plot
    plot([conc concSolve(b,:)'],'LineWidth',2);
    hold on
    errorbar(conc,sigma,'ko','LineWidth',1.5)
    xlabel('Sample name')
    ylabel('Concentration (Atoms/gram Qtz)')
%     legend({'Measured Concentration',[TRun{1},' - R chi^2 =', num2str(round(M(1)))],...
%         [TRun{2},' - R chi^2 =', num2str(round(M(2)))],[TRun{3},' - R chi^2 =', num2str(round(M(3)))]},...
%         'Location','north')
end
%%
% toc
%% 
% figure(4)
% subplot(1,2,1)
% histogram(theta*360/(2*pi))
% title('Theta')
% subplot(1,2,2)
% title('Phi')
% histogram(phi*360/(2*pi))

% 
% function [theta,phi] = cosmoDistribution(num)
% %% determine a random distribution of cosmic radiation products
% % alpha is from 0 to pi, where zero is horizontal and pi/2 is 90 degrees or
% % vertical
% rng('shuffle')
% phi = rand(num,1)*2*pi;
% 
% % This section takes twice as long computationally compared to below
% if num<100
%     theta = zeros(num,1);
%     n = 0;
%     while n<num
%         temp = pi/2*rand; y = rand;
%         Power = 2.3;  % 2.3
%         if cos(temp)^Power>=y
%             n = n +1;
%             theta(n) = temp;
%         end
%     end
% else
%     N = 4*num;
%     theta = pi/2*rand(N,1);
%     y = rand(N,1);
% 
%     i = cos(theta).^2.3>=y;
%     theta = theta(i);
% 
%     % This is a check to make sure there are enought points, 
%     % but is usually not needed for less than 1000
%     if length(theta)>num
%         theta = theta(1:num);
%     else
%         disp('Alert!')
%         N = 10*num;
%         theta = pi*rand(N,1);
%         y = rand(N,1);
% 
%         i = sin(theta).^2.3>=y;
%         theta = theta(i);
%         if length(theta)<num
%             error('x is too short')
%         end
%     end
% end
% end