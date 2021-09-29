%% JakMCMCCall.m
% This code calls the JakMCMC.m code to run through the MCMCStat model. The
% code is maintained at: https://mjlaine.github.io/mcmcstat/

options.nsimu = 1000;      % This is the number of simulations
params = {
    {'xPoint', 0.6, 0, 1}
    {'zPoint', 0.6, 0, 1}
    {'LambdaP', 190, 150, 240}
    {'ErosionDepth', 2.75, 0, 10}
    };

model.ssfun = @JakMCMC;

[res,chain,s2chain,sschain] = mcmcrun(model,[],params,options);

%% Plotting

figure(1); clf
mcmcplot(chain,[],res,'pairs');

figure(2); clf
mcmcplot(chain,[],res, 'denspanel', 2);

figure(3); clf
mcmcplot(chain,[],res,'chainpanel');

chainstats(chain,res)