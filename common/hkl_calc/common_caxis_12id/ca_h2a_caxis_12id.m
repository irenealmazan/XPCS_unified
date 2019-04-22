function [ANGLES,DOCUanglenames,ANGLES2] = ca_h2angles_caxis_12id(hkl,UB,Energy,params)
% [ANGLES] = calc_hkl2angles_caxis_12id(hkl,UB,Energy,spectrometer_params)
% version 20130327 C. Thompson   !!! I forget whether this one is working yet sigma tau (caxis spec code is somewhat messed up)
% OUTPUT 
%	ANGLES [Delta Theta Mu Gamma] in degrees, one row per hkl
%	DOCU_anglenames = documentation of angle names
%	ANGLES2 = [Alpha, Beta, Two_Theta] in degrees, one row per hkl
%
% INPUT 
%	[h k l] in r.l.u, one row per hkl
%	UB is at least UB.UB from structure from calc_UB function, but can just pass whole structure
%	Energy is energy in eV
%	params is structure of parameteres necessary for spectrometer zaxis_mocvd 
%		params.alphatarget, params.sigma, params.tau  are all in degrees, 
%		(for now, it is hardwired to mode 0 (fixed alpha))


% Order of output angles for caxis_12id (best practice is to use order from 'pa' command)
	DOCUanglenames = '[Nu Chi Eta Delta]  : caxis_12id'; 

% and for calculations hkl should be in columns
	HKL = hkl';
	UBHKL = UB.UB * [HKL] ./ (2*pi*Energy./fhc); % Matrix of Hphi values, one per column
	
% unpack the spectrometer parameters
	sigma	= params.sigma;
	tau		= params.tau;
	alphatarget	= params.alphatarget;
	
	sigmatau	= [sigma  tau];
% 	hardwired caxis_12id for fixed alpha mode (mode 0)

	LEN = length(hkl(:,1));
	ANGLES = zeros(LEN,4);


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
		sinalpharange = calc_sinalpha_caxis_12id([deltaangle thetaangle mu_range gamangle],sigmatau);

		ERROR2 = (sinalpharange - sind(alphatarget)).^2;
		% for this range - fit to A(x-x0)^2 and find x_0
		Acoef = polyfit(mu_range,ERROR2,2);

		mu_best = sqrt(Acoef(3)./Acoef(1));
		gam_best = calc_gam(mu_best,UBHKLii);
		delta_best = calc_delta(mu_best, gam_best,UBHKLii);
		theta_best = calc_theta(delta_best, mu_best, gam_best,UBHKLii);

		ANGLES(ii,:) = [delta_best,theta_best,mu_best,gam_best];
		
		[ALPHA(ii),BETA(ii)] = calc_sinalpha_caxis_12id(ANGLES(ii,:),sigmatau);
		[TWOTHETA(ii)= 1; %%%% still need to finish

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

[ANGLES2,docuangle2names] = calc_angles2_caxis_12id(angles,params)
		
end

%%%%%%%%%%%%%%  HELPER FUNCTIONS

function [sinalpha,sinbeta] = calc_sinalpha_caxis_12id(angles,sigmatau)
	%[sinalpha sinbeta] = calc_sinalpha_caxis_12id([Delta Theta Mu Gamma], [Sigma Tau])
	% for caxis_12id spectrometer, [Nu Chi Eta Delta Phi(ignored) Mu(ignored)]
	% Note - Chi is CCW about +z so define Tau = +Chi_max - 90 (when using +|sigma|)
	%		(when using -|sigma| value, Tau = +Chi + 90)
	Nu	 	= angles(:,1);
	Chi	 	= angles(:,2);
	Eta		= angles(:,3);
	Delta	= angles(:,4);
	Sigma = sigmatau(1);Tau=sigmatau(2);
	%Thetamax = -Tau -90 if Sigma is Positive definite, Thetamax = -Tau + 90 if Sigma Negative definite
	
	sinalpha = sind(Eta).*cosd(Sigma)  - cosd(Eta).*sind(Sigma).*sind(Chi-Tau);
	
	sinbeta =  cosd(Sigma).*(  sind(Delta).*cosd(Eta) - cosd(Delta).*sind(Eta).*cosd(Nu)  ) + ...
		sind(sigma)*(   sind(Chi-Tau).* ( sind(Delta).*sind(Eta) +  cosd(Delta).*cosd(Eta).*cosd(Nu) ) ...
			- cosd(Delta).*cosd(Chi-Tau).*sind(Nu)    );
end

function [gamangle] = calc_gam(mu,UBHKL)
	
	gamangle = asin(UBHKL(3,:) - sind(mu)) .*180 ./pi;
	
end

function [deltaangle] = calc_delta(mu,gam,UBHKL)
	arg1 = 1 - sum(UBHKL.^2)./2 + sind(gam).*sind(mu);
	arg2 = cosd(gam).*cosd(mu);
	% HAD SOME PROBLEM WITH ARG1/ARG2>1, such as for [1 0 6] so acos is imag
	% does spec use real part in this case, or abs does it go to -delta and invert args
		%deltaangle = sign(arg2-arg1).*abs(acos(arg1./arg2)) .* 180 ./ pi;
		%deltaangle = abs(acos(arg1./arg2)) .* 180 ./ pi; 
		deltaangle = real(acos(arg1./arg2)) .* 180 ./ pi;   
	
end

function [thetaangle] = calc_theta(del,mu,gam,UBHKL)
	argy = UBHKL(2,:).*sind(del).*cosd(gam) - UBHKL(1,:).*(cosd(del).*cosd(gam) - cosd(mu));
	argx = UBHKL(1,:).*sind(del).*cosd(gam) + UBHKL(2,:).*(cosd(del).*cosd(gam) - cosd(mu));
	
	thetaangle = atan2(argy,argx).*180./pi;
	
end




