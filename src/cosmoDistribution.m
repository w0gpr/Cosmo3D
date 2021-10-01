function [theta,phi] = cosmoDistribution(num)
%% determine a random distribution of cosmic radiation products
% alpha is from 0 to pi, where zero is horizontal and pi/2 is 90 degrees or
% vertical
rng('shuffle')
phi = rand(num,1)*2*pi;

% This section takes twice as long computationally compared to below
if num<100
    theta = zeros(num,1);
    n = 0;
    while n<num
        temp = pi/2*rand; y = rand;
        Power = 2.3;  % 2.3
        if cos(temp)^Power>=y
            n = n +1;
            theta(n) = temp;
        end
    end
else
    N = 4*num;
    theta = pi/2*rand(N,1);
    y = rand(N,1);

    i = cos(theta).^2.3>=y;
    theta = theta(i);

    % This is a check to make sure there are enough points, 
    % but is usually not needed for less than 1000
    if length(theta)>num
        theta = theta(1:num);
    else
        disp('Alert!')
        N = 10*num;
        theta = pi*rand(N,1);
        y = rand(N,1);

        i = sin(theta).^2.3>=y;
        theta = theta(i);
        if length(theta)<num
            error('x is too short')
        end
    end
end
end