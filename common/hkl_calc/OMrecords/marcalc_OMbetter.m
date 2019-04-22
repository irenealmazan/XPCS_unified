function OMchoices = marcalc_OMbetter
%Calculating better OM from the CTRs
%but the CTR's did not actually line up on BRAGG peaks, so I am 
%finding peaks, backing out angles, then using

OMspec = omstuff;  % read in the orientation matrix stuff that Spec was using

specparms.sigma = OMspec.sigma;
specparms.tau = OMspec.tau;

% thse are from 'findpoint' of nearest point to peak
% mar1713_3 #12
%mar1713_3 #12  at start of scan of (0 0-4) CTR
 Angles_start_004 = [19.3674 26.1442 0.30293 19.782679];
 HKL_start_004 = [2.14375 0.0422168 -3.91894];
 [HKL_104] = [0.9950    0.0403   -3.9436];
 [HKL_204] = [1.9900    0.0420   -3.9222];
 [HKL_304] = [2.9800    0.0436   -3.9010];
	specparms_004=specparms;
	specparms_004.alphatarget = asind(calc_sinalpha_zaxis_mocvd(Angles_start_004,[OMspec.sigma OMspec.tau]));
	% this equals 0.0834 : is what alpha_spec would have been at start of scan
	% inspecting data file, I think #G4 entry 5 is the alphatarget

%mar1613_3 #23 at start of scan of (-1 2 0) CTR
HKL_start_m120 = [-0.492864 1.99244 0.0545862];
Angles_start_m120= [15.7694 104.7244 0.50097 4.0817857];
 [HKL_020] = [0    	1.9918    0.0633];
 [HKL_120] = [1.0000    1.9905    0.0810];
 [HKL_220] = [2.0000    1.9892    0.0986];
	specparms_m120=specparms;
	specparms_m120.alphatarget= asind(calc_sinalpha_zaxis_mocvd(Angles_start_m120,[OMspec.sigma OMspec.tau]));
	% this equals 0.090 : is what alpha_spec would have been at start of scan
	% inspecting data file, I think #G4 entry 5 is the alphatarget

% Now on to calculating the angles that spec went to for these peaks
% UB won't change for any of this
UBspec = calc_UB(OMspec.cparam,OMspec.h0,OMspec.a0,OMspec.h1,OMspec.a1,'zaxis_mocvd');

[Angles_004_set]  = calc_hkl2angles([HKL_104;HKL_204;HKL_304],UBspec,OMspec.Energy,'zaxis_mocvd',specparms_004)
[Angles_m120_set] = calc_hkl2angles([HKL_020;HKL_120;HKL_220],UBspec,OMspec.Energy,'zaxis_mocvd',specparms_m120)

% following the above prescriptions to get out the alternate HKL
OMchoices.hkl = [Angles_004_set.hkl
		 Angles_m120_set.hkl]
OMchoices.angles = [Angles_004_set.angles
		   Angles_m120_set.angles];

% RUNNING ALL THIS ABOVE GIVES
%for HKL				Angles (Delta, Theta, Mu, Gamma)
%    [0.9950    0.0403   -3.9436]   [19.3179   18.0611    0.2480    9.0246]
%    [1.9900    0.0420   -3.9222]   [19.3813   24.6919    0.2934   18.3120]
%    [2.9800    0.0436   -3.9010]   [18.9021   36.3926    0.3661   28.0500]
%    [     0    1.9918    0.0633]   [15.8221  106.7081    0.4959    8.5966]
%    [1.0000    1.9905    0.0810]   [15.8787  114.8931    0.4697   17.9757]
%    [2.0000    1.9892    0.0986]   [15.1391  129.6301    0.4035   27.9336]

end

%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%
function  OM = omstuff
	% This OM was used for the entire m-plane set of data during 2013-03 nitride run.
	OM.cparam = [3.186 3.186 5.178 90 90 120];
	cparamdocu = 'GaN lattice parameters';

	OM.h0 = [0 0 2];
	OM.a0 = [9.721 -169.783 0 0];
	OM.h1 = [1 -2 0];
	OM.a1 = [15.796 -76.8284 0 0];
	OM.lambda = 0.438105327;   % what spec was using
	OM.Energy = fhc/OM.lambda

% double checked - sigma tau stay same
	OM.sigma = 0.43445; OM.tau = -175.79969;
% where I thought alphatarget is says  0.2378390714 for both scans, so it cant be right
end





