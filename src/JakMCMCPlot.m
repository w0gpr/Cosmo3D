%% JakMCMCPlotting.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code reads the input file JakMCMCData.mat created from the output of
% JakMCMCCall.m. It is supposed to read in the scenarios that minimized the
% chi^2 reduction of the parameters xP, zP, LambdaP, erosionDepth, then
% feeds back into the JakMCMC.m function to create the outputs for plotting
% the data. A better approach would have been to do this during the MCMC
% call, but I could not determine how to do this easily using the
% MCMCStats toolbox. It probably could be done, but this should turn out
% easier. It also allows me to parallize it, which is always the dream :)
% Brandon Graham, University at Buffalo, Geology Department
% MIT License
% 20210930
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all

%% Initialize parallel pool
% pool = 1;
% if ~exist('pool','var')
%     pool = parpool;
% end

%% Initialize Data and setup
% tic
% info = JakReadData;
% numIter = 500;
% numRays = 1000;
% load '/home/brandon/Matlab/Cosmo3D/data/Output/JakMCMCData.mat';
% aC = sortrows([sschain chain]);   % allChain
% concSolve = zeros(numIter,info.numSamples);
% 
% 
% %% Run the parfor loop
% toc
% parfor ii = 1:numIter
%     out{ii} = JakCosmo3DFull(dataTable, aC(ii,2), aC(ii,3), aC(ii,4), aC(ii,5), numRays);
% end
% toc
% %% break apart the data
% X1 = out{1,1}{1,3};
% Y1 = out{1,1}{1,4};
% Z1 = out{1,1}{1,5};
% tri = out{1,1}{1,6};

loadData = false; % This is to load the data from the full MCMC
if loadData
%     load /home/brandon/Matlab/Cosmo3D/data/Output/JakMCMCPlotData.mat;
end
% 
% eM = [0 1 208 0.1    % end Members
%     0 1 182 0.1
%     0 1 208 4.1
%     0 1 182 4.1
%     1 0 208 0.1
%     1 0 182 0.1
%     1 0 208 4.1
%     1 0 182 4.1];

% for ii = 8:-1:1
%     out2{ii} = JakCosmo3DFull(dataTable, eM(ii,1), eM(ii,2), eM(ii,3), eM(ii,4), numRays);
% end

out2 = GrahamModel3DBuilderRayTrace(); % Load specific end members scenarios 
X3 = out2{1,3}(2,2:end-1);
% Y1 = out{1,1}{1,4};
Z3 = zeros(4,length(X3));
for jj = 1:4
Z3(jj,:) = out2{1,5}{1,jj}(2,2:end-1);
end
    
% figure(1)
% trisurf(tri,X1,Y1,Z1)
% chi2 = zeros(numIter,1);
% Z2 = zeros(numIter,29);
% for ii = 1:numIter
%     concSolve(ii,:) = out{1,ii}{1,2}(:)';
%     chi2(ii) = out{1,ii}{1,1};
%     Z2(ii,:) = out{1,ii}{1,5}(2,2:end-1);
% end
% X2 = out{1,1}{1,3}(2,2:end-1);
%% Plotting the data
figure(2); clf
plot(X2(1,:),Z2,'Color',[0.8 0.8 0.8 0.01],'LineWidth',1)
hold on
plot(X3,Z3,'LineWidth',3)
for i=1:info.numSamples
    plot(info.rSx(i,:),info.rSz(i,:),'k-','LineWidth',3)
end
plot(info.xPresent,info.zPresent,'k','LineWidth',3)
set(gcf, 'Position', [39 371 2316 323]);


figure11 = figure(3);
clf
axes1 = axes('Parent',figure11,'XTickLabel',info.sampleID,'XDir','reverse');
box(axes1,'on');
hold(axes1,'all');
% Create plot
plot(concSolve','Color',[0.8 0.8 0.8, 0.05],'LineWidth',1);
hold on
plot(out2{1,2}','LineWidth',3)
plot(data.conc,'k','LineWidth',3)
errorbar(data.conc,data.sigma,'ko','LineWidth',3)
% errorbar(data.conc,data.sigma*2,'ko','LineWidth',1)
xlabel('Sample name')
ylabel('Concentration (Atoms/gram Qtz)')
set(gcf, 'Position', [1 1 1920 961]);
    
    
% save /home/brandon/Matlab/Cosmo3D/data/Output/JakMCMCPlotData.mat    
