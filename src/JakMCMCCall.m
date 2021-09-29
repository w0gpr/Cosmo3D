%% JakMCMCCall.m
% This code calls the JakMCMC.m code to run through the MCMCStat model. The
% code is maintained at: https://mjlaine.github.io/mcmcstat/

data = readtable('Camp3Samples.csv');
numSamples = 8;
data = data{1:numSamples,2:3};

options.nsimu = 100;      % This is the number of simulations
params = {
    {'xPoint', 0.6, 0, 1}
    {'zPoint', 0.6, 0, 1}
    {'LambdaP', 190, 150, 240}
    {'ErosionDepth', 2.75, 0, 10}
    };

% modelfun = @JakMCMC;
% model.ssfun = @(params,data) sum(((data(1,:)-modelfun(params,data))./...
%     data(2,:)).^2);
model.ssfun = @JakMCMC;
model.N = numSamples;
% options.adaptint = 200;
% options.waitbar = 1;
options.method = 'mh';
options.updatesigma = 1;


[res,chain,s2chain,sschain] = mcmcrun(model,data,params,options);

%% Plotting

figure(1); clf
mcmcplot(chain,[],res,'pairs');

figure(2); clf
mcmcplot(chain,[],res, 'denspanel', 2);

figure(3); clf
mcmcplot(chain,[],res,'chainpanel');
%% 
figure(4); clf
for i = 1:length(params)
    subplot(2,2,i)
    plot(chain(:,i),sschain,'*')
    title(params{i}{1})
end

figure(5); clf
plot(1:options.nsimu,sschain)

chainstats(chain,res)