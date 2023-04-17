function out = JakReadData()
%% Read in sample data
% This section reads in the data from the .csv and then rotates the datum.
% It is done in a way that is dynamic and allows for buffers to be
% increased or decreased. 

% This was an attempt to break the "ReadData" into a separate function, but was later abandoned to minimize complexity.
data = readtable('Camp3Samples.csv');
out.numSamples = 8;
out.sampleID = data{1:out.numSamples,1};
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
out.conc = D(:,1);
out.sigma = D(:,2);
out.Pspal = D(1,3);
out.Pmuon = D(1,4);
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
out.xPresent = [0 blockXEdges(1) blockXEdges(2) blockXEdges(2) width];
out.zPresent = [zBuf zBuf zBuf dT2 dT2];

% This is the closed rectangular cross section of each sample. A y component
% is not needed in this model as the y component does not vary and is set
% to 0.
out.rSx = [xSamp(:,1) xSamp(:,2) xSamp(:,2) xSamp(:,1) xSamp(:,1)];
out.rSz = [zSamp(:,2) zSamp(:,2) zSamp(:,1) zSamp(:,1) zSamp(:,2)];

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
out.XSamp = [min(xSamp,[],2) max(xSamp,[],2)];
out.ZSamp = [min(zSamp,[],2) max(zSamp,[],2)];

