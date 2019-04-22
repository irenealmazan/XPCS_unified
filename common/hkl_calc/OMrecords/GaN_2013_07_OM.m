

% THIS IS BACKGROUND INFO THAT IS NEEDED . 
% FOR HKL to ANGLE - ALL ITEMS BELOW ARE REQUIRED TO BE EXACTLY AS SPEC USED FOR CALCULATIONS
% FOR ANGLE to HKL - (sigma and tau and even alphatarget can be left as dummy values)

	% This OM was used for the entire m-plane set of data during 2013-07 nitride run.
	cparam = [3.186 3.186 5.178 90 90 120];
	cparamdocu = 'GaN lattice parameters';

	if 0   % OM from spec
	% NOTE - THIS MUST BE USED IF GOING FROM SPEC HKL to ANGLES
		h0 = [2 0 -4];
		a0 = [19.5730 23.4474 0.2816 18.4114];
		h1 = [1 0 -4];
		a1 = [19.4142 16.8378 0.2784 8.9932];
		lambda = 0.438105327;
		Energy = fhc/lambda;
	end

	if 0
	% These were the positions of peaks, only found once and never tweaked or checked again
	% just to compare to last time 2013-03
%	h10m4 = [1 0 -4];		a10m4 = [19.3179   18.0611    0.2480     9.0246];
%	h20m4 = [2 0 -4]; 		a20m4 = [19.3813   24.6919    0.2934   18.3120];
%	h00m4	

	h10m4 = [1 0 -4];		a10m4 = [19.4142   16.8378   0.2784  8.9932];
	h20m4 = [2 0 -4]; 		a20m4 = [19.5730  23.4474   0.2816 18.4114];
 	h00m4 = [0 0 -4]; 		a00m4 = [10.3842  14.6890  0.2768  -0.0225];

%	UB1 = calc_UB(cparam,h10m4,a10m4,h020,a020,'zaxis_mocvd');
%	UB2 = calc_UB(cparam,h20m4,a20m4,h120,a120,'zaxis_mocvd');
%	UB3 = calc_UB(cparam,h30m4,a30m4,h220,a220,'zaxis_mocvd');

	end
