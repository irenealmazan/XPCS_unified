% PRETEND DATA
%
% For the inputs - this is form to use
% HKL =[...
%	h1 k1 l1
%	h2 k2 l2
%	h3 k3 l3
%	h4 k4 l4];
%%
% ANGLES = [...
%	del1 th1 mu1 gam1
%	del2 th2 mu2 gam2
%	del3 th3 mu3 gam3
%	del4 th4 mu4 gam4];
%  


% just some test data to use
DATAhkl = [...
	-1  2 0
	-0.5 2 0
	0 0 2
	0 0 -2
	0.5 0 -4
	0.6 0 -4
	0.7 0 -4];

DATAangle = [...
	9.721 -169.783 0 0
	15.796 -76.8284 0 0
	19.5186   15.7037    0.2641   4.8185
   	19.5321   15.9556    0.2661   5.7342
   	19.5475   16.2515    0.2684   6.6516];

% THIS IS BACKGROUND INFO THAT IS NEEDED . 
% FOR HKL to ANGLE - ALL ITEMS BELOW ARE REQUIRED TO BE EXACTLY AS SPEC USED FOR CALCULATIONS
% FOR ANGLE to HKL - (sigma and tau and even alphatarget can be left as dummy values)

	% This OM was used for the entire m-plane set of data during 2013-03 nitride run.
	cparam = [3.186 3.186 5.178 90 90 120];
	cparamdocu = 'GaN lattice parameters';

	h0 = [0 0 2];
	a0 = [9.721 -169.783 0 0];
	h1 = [1 -2 0];
	a1 = [15.796 -76.8284 0 0];
	lambda = 0.438105327;
	Energy = fhc/lambda;

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

	if 1
	% for the -0.5 2 0 growths (mar1613_1)
	sigma = 0.48497; tau = -174.66867; alphatarget = 0.09; 
	%(alphatarget = varied from 0.09 to 0.012, 
	end

	if 0
	% for the +0.5 0 -4 growths (mar1713_2)
	sigma = 0.43445; tau = -175.799690; alphatarget = 0.09;
	% alphatarget varied from 0.09 to 0.08 or similar
	end

	if 1
	% if pretending there is no miscut so that H00 is also exactly normal to the surface
	sigma = 0;	tau = 90;  alphatarget = 0.0;
	end


% First calculate the UB's - ca_UB only needs to be done ONCE per any 
%	particular orientation matrix used. 
%	(if either of the 2 bragg peaks or their angles, or the lattice parameters/angles changes)
%	These calculations do NOT depend on the sigma, tau, or Energy.
	[UBexpt] = calc_UB(cparam,h0,a0,h1,a1,'zaxis_mocvd');


% Angles to hkl, use this method. 
%	(For 2D detector, only the del and gam change
%	the mu and theta (sample orientation) are the same for each point  
%		for zaxis_mocvd angle order  - Delta Theta Mu Gamma
if 1
	HKLcalc =  calc_angles2hkl(DATAangle,UBexpt,Energy,'zaxis_mocvd');
end


% hkl to angles, use this method. (it also is tricker all around)
if 1
	%  Must put in the exact alpha that SPEC thought it was using. 
	%  Note - this may vary scan to scan 
	%  The alpha that spec thought it was at is in the header of the scan in line # G4, entry 5
	spec_params.alphatarget = alphatarget;
	spec_params.tau 	= tau;
	spec_params.sigma	= sigma;

	ANGLEScalc = calc_hkl2angles(DATAhkl,UBexpt,Energy,'zaxis_mocvd',spec_params) ; 
end

% THE ARE THE FUNCTIONS THAT ARE REQUIRED 
%	(maybe put them in a common directory 
% Base level functions 
% fhc  % gives constant fhc = 12398.4 i.e., hc in eV-Angstr
% ca_hphi_zaxis_mocvd
% ca_angle2hkl_zaxis_mocvd

% Upper level ones 
% calc_UB
% calc_hkl2angles
% calc_angles2hkl
%
% The ones to use in your programs are the calc_ functions. 
%	Their output is also self-documenting so that
% 	it is easy to check later if the wrong OM, alpha, etc was used or changed in the middle.



