function [out,data,ver] = skyline(az,el,strike,dip);

% skyline.m
%
% Syntax: [correction_factor horizon] = skyline(azimuths,elevations,plunge,dip);
%
% This function calculates the shielding factor for cosmogenic isotope 
% production on a surface which is either dipping, shielded by high 
% topography, or both. 
%
% The first two arguments are vectors of similar 
% length containing the azimuth (0-360) and elevation (0-90) or points 
% on the horizon. This presumes that in the field the horizon was 
% approximated by a series of points with straight lines between them. 
% Note that this is incompatible with field procedures in which 
% the horizon is approximated by the average rise angle in a series 
% of equal sectors. This latter procedure is inappropriate anyway 
% because the relationship between rise angle and shielding is nonlinear.
%
% The second two arguments are optional and are scalars describing the 
% plunge and dip of the surface sampled. Use the convention that 
% the dip is down to the right of the strike direction, i.e. strike = 0,
% dip = 45 means that the surface is dipping 45 down to the E. 
%
% If you just want the shielding for a dipping surface enter
% empty vectors ([]) for the first two arguments.
%
% Note that this function first compares the "horizon" created by a
% dipping surface (second 2 arguments) with the far-field horizon 
% (first 2 arguments) and calculates a composite horizon before 
% determining the shielding factor. Thus the shielding due to a 
% dipping surface is never "double-counted" with a distant horizon.
%
% Uses highly simplistic numerical integration -- approximates
% horizon by a series of 1-degree sectors and uses an integration
% formula to calculate the shielding for each of these sectors,
% then sums them. 
% 
% The output argument correction_factor is the ratio of nuclide production 
% at the shielded site to that at an unshielded site at the same location.
%
% The output argument 'horizon' is a vector containing the horizon angle in degrees
% at 1-degree increments. Useful for plotting purposes. 
% 
% The output argument 'ver' is a string giving the version number of this function. 
%
% Written by Greg Balco -- UW Cosmogenic Nuclide Lab
% balcs@u.washington.edu
% First version, Feb. 2001
% Revised March, 2006
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math
%
% Copyright 2001-2007, University of Washington
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).
 

% what's the version

ver = '1.1';

% error checks

if max(max(az)) > 360;
	error('azimuths 0 to 360 please');
elseif max(max(el)) > 90;
	error('elevations 0 to 90 please');
elseif length(az) ~= length(el);
	error('vectors the same length please');
end;

% set up angular framework -- 1-degree increments...

angles = [0:pi/180:(2*pi)];

% step 1. Determine if dipping surface needed; create horizon
% from such surface.

if nargin == 4;

	% convert to radians
	
	strikeR = (strike/360) * (2*pi);
	dipR = (dip/360) * (2*pi);
	
	% calculate "horizon"
	% updip direction = strike direction - pi/2 by convention
	% a relative horiz angle to dip direction
	
	a = angles - (strikeR - (pi/2)); 
	
	horiz1 = atan(tan(dipR) .* cos(a));
	
	% remove entries less than 0...
	
	horiz1(find(horiz1 < 0)) = zeros(size(find(horiz1 < 0)));
	
elseif nargin == 2;

	horiz1 = zeros(size(angles));

else;

	error('Wrong number of arguments?');
	
end;

if ~(isempty(az) & isempty(el));
	
	% step 2. Interpolate linearly between supplied horizon points.

	% convert to radians

	azR = (az/360) * (2*pi);
	elR = (el/360) * (2*pi);

	% sort in ascending order

	[azR i] = sort(azR);
	elR = elR(i);

	% pad for interpolation;

	azR2(2:(length(azR)+1)) = azR;
	elR2(2:length(elR)+1) = elR;
	azR2(1) = azR(length(azR)) - (2*pi); 
	elR2(1) = elR(length(elR));
	azR2(length(azR)+2) = azR(1) + (2*pi); 
	elR2(length(azR)+2) = elR(1);

	% interpolate;

	horiz2 = interp1(azR2,elR2,angles);

	% flip 

	horiz2 = horiz2;
	
else;

	horiz2 = zeros(size(angles));
	
end;

% step 3. make composite horizon. 

horiz = max([horiz1;horiz2]);

% step 4. Integrate shielding function.

% Note: here we do the simplest possible thing, integrate in 
% 1-degree increments and add them up.

% sector angle = 1 degree;

B = (1/360)*(2*pi);

% integration formula

S = (B/(2*pi)) .* (sin(horiz).^3.3);

% above formula gives percent shielding; subtract from 1 to give corr. factor --

out = 1 - sum(S);

data = horiz.*180/pi;


	
	
	

