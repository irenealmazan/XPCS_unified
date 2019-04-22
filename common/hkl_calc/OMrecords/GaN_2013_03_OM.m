elempath.Henke='/science/matlab/calculations/atomicformfactors/elements_HENKE';
elempath.WNK = '/science/matlab/calculations/atomicformfactors/F0_WaasmaierKirfel';

% THIS IS BACKGROUND INFO THAT IS NEEDED . 
% FOR HKL to ANGLE - ALL ITEMS BELOW ARE REQUIRED TO BE EXACTLY AS SPEC USED FOR CALCULATIONS
% FOR ANGLE to HKL - (sigma and tau and even alphatarget can be left as dummy values)

	% This OM was used for the entire m-plane set of data during 2013-03 nitride run.
	cparam = [3.186 3.186 5.178 90 90 120];
	cparamdocu = 'GaN lattice parameters';

	if 0   % OM from spec
	h0 = [0 0 2];
	a0 = [9.721 -169.783 0 0];
	h1 = [1 -2 0];
	a1 = [15.796 -76.8284 0 0];
	lambda = 0.438105327;
	Energy = fhc/lambda;
	end

	if 1
	% For angles to hkl, can use these back calculated from the CTR's (we do not have
	% bragg peaks lined up so not able to get better OM directly from prpeaks
	h10m4 = [1 0 -4];	a10m4 = [19.3179   18.0611    0.2480    9.0246];
	h20m4 = [2 0 -4]; 	a20m4 = [19.3813   24.6919    0.2934   18.3120];
 	h30m4 = [3 0 -4]; 	a30m4 = [18.9021   36.3926    0.3661   28.0500];
	h020  = [0 2  0];	a020  = [15.8221  106.7081    0.4959    8.5966];
	h120  = [1 2  0]; 	a120  = [15.8787  114.8931    0.4697   17.9757];
 	h220  = [2 2  0];	a220  = [15.1391  129.6301    0.4035   27.9336];
	
	%    [1 0 -4] [0.9950    0.0403   -3.9436]   [19.3179   18.0611    0.2480    9.0246]
	%    [2 0 -4] [1.9900    0.0420   -3.9222]   [19.3813   24.6919    0.2934   18.3120]
	%    [3 0 -4] [2.9800    0.0436   -3.9010]   [18.9021   36.3926    0.3661   28.0500]
	%    [0 2  0] [     0    1.9918    0.0633]   [15.8221  106.7081    0.4959    8.5966]
	%    [1 2  0] [1.0000    1.9905    0.0810]   [15.8787  114.8931    0.4697   17.9757]
	%    [2 2  0] [2.0000    1.9892    0.0986]   [15.1391  129.6301    0.4035   27.9336]
	UB1 = calc_UB(cparam,h10m4,a10m4,h020,a020,'zaxis_mocvd');
	UB2 = calc_UB(cparam,h20m4,a20m4,h120,a120,'zaxis_mocvd');
	UB3 = calc_UB(cparam,h30m4,a30m4,h220,a220,'zaxis_mocvd');

	lambda = 0.438105327;
	Energy = fhc/lambda;
	end

	% sigma and tau changes every so often whenever we redo sigmatau, 
	%  and alphatarget changes a lot (as things wander)
	% need to put in the correct alpha for a particular scan if using hkl2angle 
	%  look at the header of the scan - (see #G4 entry 5)

	% if using calc_hkl2angle, use the appropriate sigma and tau and alphatarget for the scan.
	% if using calc_angle2hkl, just need a valid orientation matrix, Energy, and cparams
	%
	if 0
	% for the +0.5 0 -2 growths (mar1513_4)
	sigma = 0.48497; tau = -174.66867; alphatarget = 0.1;
	%(alphatarget = varied from 0.15 to 0.09 depending on scan)
	end

	if 0
	% for the -0.5 2 0 growths (mar1613_1)
	sigma = 0.48497; tau = -174.66867; alphatarget = 0.09; 
	%(alphatarget = varied from 0.09 to 0.012, 
	end

	if 1
	% for the +0.5 0 -4 growths (mar1713_2)
	sigma = 0.43445; tau = -175.799690; alphatarget = 0.09;
	% alphatarget varied from 0.09 to 0.08 or similar
	end

	if 0
	% if pretending there is no miscut so that H00 is also exactly normal to the surface
	sigma = 0;	tau = 90;  alphatarget = 0.0;
	end

	NormalLC = [sind(sigma).*cosd(tau);-sind(sigma).*sind(tau);cosd(sigma)];
	% these all give very similar, using 
	HKLnorm1 = (UB1.UBinv*NormalLC);
	HKLnorm2 = (UB2.UBinv*NormalLC);
	HKLnorm3 = (UB3.UBinv*NormalLC);
	% mean gives (using sigma for the -0.5 2 0) 
	% HKLnorm = [0.4390  0.0003  -0.0030]  or   scaled a bit 
		% HKLnorma = [1   0.00065  -0.00681];
	%   using the sigma for the 0.5 0 -4)
		% HKLnorm  = [1   0.00028  -0.00837];

filter1 = 2.295;  %28300  and 0.125 Sn filter box for Nitride expt calculated so far
                    % jul2905 around 2pm trying from set to see. Have trouble
% filtercorr parameters are used when there are thickness variations or materials variations for foil sets used for filters. Can keep all at 0 for 12id-d MOCVD run because it is same foils, 
% If correction is needed, it is used as a correction to the exponent of filter1^(filternumber+filtercorr)
filter2corr = 0;
filter4corr = -0.083;  % used for Sn foils 4th foil seems a little off
%filter4corr = 0;
filter8corr = 0;
filtercorr = [filter2corr,filter4corr,filter8corr];

% From marcalc_OMbetter (which uses the CTR's and back calculations to get angles, to figure out
% Bragg peak positions to put into a better OM to get orientation 
if 0
% RUNNING ALL THIS ABOVE GIVES
%for HKL				Angles (Delta, Theta, Mu, Gamma)
%    [0.9950    0.0403   -3.9436]   [19.3179   18.0611    0.2480    9.0246]
%    [1.9900    0.0420   -3.9222]   [19.3813   24.6919    0.2934   18.3120]
%    [2.9800    0.0436   -3.9010]   [18.9021   36.3926    0.3661   28.0500]
%    [     0    1.9918    0.0633]   [15.8221  106.7081    0.4959    8.5966]
%    [1.0000    1.9905    0.0810]   [15.8787  114.8931    0.4697   17.9757]
%    [2.0000    1.9892    0.0986]   [15.1391  129.6301    0.4035   27.9336]

end
