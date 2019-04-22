function [ANGLES,DOCU_anglenames,ANGLES2,DOCU_angle2names,ERROR] = ca_h2a_zaxis_mocvd(hkl,UB,Energy,params)
% [ANGLES,DOCU_anglenames,ANGLES2,DOCU_angle2names,ERROR] = ca_h2a_zaxis_mocvd(hkl,UB,Energy,params)
% version 20150517 C. Thompson - added more outputs, error codes, still only valid for mode 0 (fixed alpha)
% version 20130327 C. Thompson
% Still need to check and add error codes
% calculated thus far (hardwired) for fixed alpha (incident grazing angle w.r.t. surface)
% OUTPUT 
%	ANGLES = [Rho, Nu, Delta, Mu, Eta, Chi, Phi] in degrees, one row per hkl
%	DOCU_anglenames = documentation of angle names
%	ANGLES2 = [Alpha, Beta, Two_Theta] degrees
%	DOCU_angle2names = documentation of angles2 names 
%	ERROR = structure of error indicators
%
%   ERROR = structure of error indicators
%       ERROR.Alpha_min = lowest possible Alpha
%       ERROR.Alpha_max = highest possible Alpha
%       ERROR.Alpha_low = true if Alpha < Alpha_min
%       ERROR.Alpha_hi = true if Alpha > Alpha_max
%
%
% INPUT 
%	[h k l] in r.l.u, one row per hkl
%	UB is at least UB.UB from structure from calc_UB function, but can just pass whole structure
%	Energy is energy in eV
%	params is structure of parameteres necessary for spectrometer zaxis_mocvd  
%		params.alphatarget
%		params.sigma
%		params.tau
%		params.betatarget (not  used)
%		params.mode (not used)  << only mode 0 implemented
%		params.phi  (not used) 

%   if there were implemented	
%   params.mode: integer from 0 to 4
%   mode for:        
%   ----------------------------------
%   fixed alpha         0      <<< only one implemented 
%   fixed beta          1       
%   alpha=beta          2       
%   fixed phi           3       
%
%   params.sign_nu controls choice of Hlabx root   (params.sign_nu = +1 or -1)
%   giving sign of nu, modes 0-2 and 4-6
%
%   the following are all in degrees, not all are used in every mode
%   params.sigma and params.tau:	specify surface orientation 
%   params.alphatarget, params.betatarget: specify the target alpha (mode 0,4) or beta (mode 1,5)
%   params.Rho, params.Mu:		fixed values (used in all modes)
%   params.Phi:				fixed values (used for modes 3 and 7)
%   params.Chi:				fixed value (used for mores 4-7)

% Note this program does not find alpha iteratively but in a graphical sense. 
% There were times when spec zaxis code thus says it has found angles, but they are not correct
%	(there ARE no good solutions. This typically will be when Delta close to zero
%	This code will output that there are no solutions (ERROR.

% Equations for gamma, delta, theta (CT independently verified)
% are  equivalent to Equations 32a (gamma), 24 (delta), 33 ('omega') 
% in Lohmeier and Vlieg, J Appl Cryst 1993 vol 26, 706
% Note zaxis 'mu' is the same as their alpha, but we
%  do not 'align' the surface normal with diffractometer so alpha=mu CANNOT be used
%  as constraint to find all other angles in one calculation consistent with target alpha. 

% Note that there is a sign typo in  L+V(1993) equation 31
% The [cos(del)cos(gam)-cos(alpha)] coeff should be identical in x and y component
% but they have a + sign in the x component.  This affects calculation of h_phi,
% which is critical in calculating the UB matrix as well as converting angles to hkl. 
% I note that zaxis spec program  uses correct  h_phi  calculation, and that  
% L+V 1993 Equation 33 (for 'omega'), based on Equation 31, is correct and 
% does not propagate the error so it is just a typo in their paper.

% Order of output angles for zaxis_mocvd (best practice is to use order from 'pa' command)
	DOCU_anglenames = '[Delta  Theta  Mu  Gamma] : zaxis_mocvd';  
	DOCU_angle2names ='[Alpha  Beta  Two_Theta]';

% and for calculations hkl should be in columns
	HKL = hkl';
	UBHKL = UB.UB * [HKL] ./ (2*pi*Energy./fhc); % Matrix of Hphi values, one per column
	
% unpack the spectrometer parameters
	Sigma	= params.sigma;
	Tau		= params.tau;
	alphatarget = params.alphatarget;
	sigmatau = [Sigma  Tau];
 
% 2013 CT currently, hardwired zaxis_mocvd to be for fixed alpha mode (mode 0)
%	in future, would add another spectometer parameter with mode (spectrometer_params.mode,
%	and it would be flag to which set of calculations to use for hkl2angle.

% make a range of mu to use calculate and home in on which mu and theta gets to correct alpha
%	Note - unlike spec, this program will work as if alphatarget was frozen
	mu_range = ( [-1:.05:1].*(sigma + 1e-2) + alphatarget)';   

	LEN = length(hkl(:,1));
	ANGLES = zeros(LEN,4);
	
	for ii=1:LEN
		
		UBHKLii = UBHKL(:,ii);
		
		% from these mu's calc what  gam, delta, theta to reach Q
		gamangle = calc_gam(mu_range,UBHKLii);
		[deltaangle] = calc_delta(mu_range,gamangle,UBHKLii);
		[thetaangle] = calc_theta(deltaangle, mu_range,gamangle,UBHKLii);
 
		% what is actual incident angle to surface given these angles
		% (depends only on theta and mu actually)
		sinalpharange = ...
			calc_sinalpha_zaxis_mocvd([deltaangle thetaangle mu_range gamangle],sigmatau);

		ERROR2 = (sinalpharange - sind(alphatarget)).^2;
		% for this range - 
		% fit to A(x-x0)^2 and find x_0    
		% A = Acoef(1), -2Ax_0 = Acoef(2), A(x_0)^2 = Acoef(3)
		Acoef = polyfit(mu_range,ERROR2,2);

		mu_best = -Acoef(2)./(Acoef(1).*2);
		gam_best = calc_gam(mu_best,UBHKLii);
		[delta_best,ERRORdel] = calc_delta(mu_best, gam_best,UBHKLii);
		theta_best = calc_theta(delta_best, mu_best, gam_best,UBHKLii);

		ANGLES(ii,:) = [delta_best,theta_best,mu_best,gam_best];	

		% if curious on methodology used to find best mu - 
		% change the following from "if 0" to  "if 1" 
		% (it makes as many figures as there are hkl input
		% so do this only when there are just a few being calculated!)
		if 0;
			figure
			ANGLEall = [deltaangle,thetaangle,mu_range, gamangle];
			plot(mu_range, ERROR2,'o',...
				mu_range,polyval(Acoef,mu_range),'b-',...
						[mu_best mu_best],[-.2 .2].*max(ERROR2),'r-');
			xlabel(['mu trials,  best mu at ' num2str(mu_best), 'deg']);
			ylabel('Error [(sin(calc alpha) - sin(target alpha)]^2');
			title(['for hkl at [',num2str([hkl(ii,:)]),'] and sigma tau [',...
				 num2str([sigmatau(:)']),']']);
		end

	end
	[ANGLES2,docuangle2names] = calc_angles2_zaxis_mocvd(ANGLES,params);
end

%%%%%%%%%%%  HELPER FUNCTIONS BELOW %%%%%%%%%%%%%%%%%%%

function [sinalpha,sinbeta] = calc_sinalpha_zaxis_mocvd(angles,sigmatau)
	%[]sinalpha sinbeta] = calc_sinalpha_zaxis_mocvd([Delta Theta Mu Gamma], [Sigma Tau])
	Delta = angles(:,1);
	Theta = angles(:,2);
	Mu=angles(:,3);
	Gamma=angles(:,4);
	Sigma = sigmatau(1);Tau=sigmatau(2);
	%Thetamax = -Tau -90 if Sigma is Positive definite, Thetamax = -Tau + 90 if Sigma Negative definite
	
	sinalpha = sind(Mu).*cosd(Sigma)  + cosd(Mu).*sind(Sigma).*sind(Tau+Theta);	
	sinbeta =  sind(Gamma).*cosd(Sigma) - cosd(Gamma).*sind(Sigma).*sind(Theta + Tau - Delta);
end


function [gamangle] = calc_gam(mu,UBHKL)
	
	gamangle = asind(UBHKL(3,:) - sind(mu)) ;
	
end

function [deltaangle,ERROR] = calc_delta(mu,gam,UBHKL)
% ERROR is -1 if we really cannot get here from there
% but there is a fix up so that it will not cough up too much

	arg1 = 1 - sum(UBHKL.^2)./2 + sind(gam).*sind(mu);
	arg2 = cosd(gam).*cosd(mu);
	% HAD SOME PROBLEM WITH ARG1/ARG2>1, such as for [1 0 6] so acos is imag
	% does spec use real part in this case, or abs does it go to -delta and invert args
	% However, in any case, it can no longer possible to get to hkl with sinalpha but it does not cough up errors
	% It keeps sinalpha pretty close, 
		deltaangle = sign(arg2-arg1).*abs(acosd(arg1./arg2));
		%deltaangle = abs(acosd(arg1./arg2)); 
		%deltaangle = real(acosd(arg1./arg2));   
		%deltaangle = acosd(arg1./arg2);
		ERROR = sign(arg2-arg1);    % -1 means final calculation will not really get to hkl specified
end

function [thetaangle] = calc_theta(del,mu,gam,UBHKL)
	argy = UBHKL(2,:).*sind(del).*cosd(gam) - UBHKL(1,:).*(cosd(del).*cosd(gam) - cosd(mu));
	argx = UBHKL(1,:).*sind(del).*cosd(gam) + UBHKL(2,:).*(cosd(del).*cosd(gam) - cosd(mu));
	
	thetaangle = atan2(argy,argx).*180./pi;
	
end




