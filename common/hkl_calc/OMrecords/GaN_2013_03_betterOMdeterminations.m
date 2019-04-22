

% THIS IS BACKGROUND INFO THAT IS NEEDED . 
% FOR HKL to ANGLE - ALL ITEMS BELOW ARE REQUIRED TO BE EXACTLY AS SPEC USED FOR CALCULATIONS
% FOR ANGLE to HKL - (sigma and tau and even alphatarget can be left as dummy values)

	% This OM was used for the entire m-plane set of data during 2013-03 nitride run.
	cparam = [3.186 3.186 5.178 90 90 120];
	cparamdocu = 'GaN lattice parameters';

	if 0   % OM from spec
	% NOTE - THIS MUST BE USED IF GOING FROM SPEC HKL to ANGLES
		h0 = [0 0 2];
		a0 = [9.721 -169.783 0 0];
		h1 = [1 -2 0];
		a1 = [15.796 -76.8284 0 0];
		lambda = 0.438105327;
		Energy = fhc/lambda;
	end

	if 1
	% For angles to hkl, can use these in an OM
		% There were back calculated from the CTR's scans by CT 
		% (we do not have bragg peaks lined up and documented)
	h10m4 = [1 0 -4];		a10m4 = [19.3179   18.0611    0.2480    9.0246];
	h20m4 = [2 0 -4]; 	a20m4 = [19.3813   24.6919    0.2934   18.3120];
 	h30m4 = [3 0 -4]; 	a30m4 = [18.9021   36.3926    0.3661   28.0500];
	h020  = [0 2  0];		a020  = [15.8221  106.7081    0.4959    8.5966];
	h120  = [1 2  0]; 		a120  = [15.8787  114.8931    0.4697   17.9757];
 	h220  = [2 2  0];		a220  = [15.1391  129.6301    0.4035   27.9336];

	%   This was the data that was used as obtained from the CTR scans and indicated above.
	%	The hkl_PEAK, hkl of peak in CTR scan, the angles back calculated using spec alpha and OM	
	%    [1 0 -4] [0.9950    0.0403   -3.9436]   [19.3179   18.0611    0.2480    9.0246]
	%    [2 0 -4] [1.9900    0.0420   -3.9222]   [19.3813   24.6919    0.2934   18.3120]
	%    [3 0 -4] [2.9800    0.0436   -3.9010]   [18.9021   36.3926    0.3661   28.0500]
	%    [0 2  0] [     0    1.9918    0.0633]   [15.8221  106.7081    0.4959    8.5966]
	%    [1 2  0] [1.0000    1.9905    0.0810]   [15.8787  114.8931    0.4697   17.9757]
	%    [2 2  0] [2.0000    1.9892    0.0986]   [15.1391  129.6301    0.4035   27.9336]


	UB1 = calc_UB(cparam,h10m4,a10m4,h020,a020,'zaxis_mocvd');
	UB2 = calc_UB(cparam,h20m4,a20m4,h120,a120,'zaxis_mocvd');
	UB3 = calc_UB(cparam,h30m4,a30m4,h220,a220,'zaxis_mocvd');

	end