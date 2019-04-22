function [ANGLES2,docuangle2names] = calc_angles2_zaxis_mocvd(angles,params)
%% [ANGLES2,docuangle2] = calc_angles2_zaxis_mocvd(angles,params)
%% 
%%  OUTPUT
%%		ANGLES2  (degrees) Alpha Beta TwoTheta
%%		docuangles2 (string with names of the data or angles in angles2 with 2 spaces between each name)
%%  INPUT
%%		angles  (each row a new set of angles, in the order required for spectometer)
%%					[Delta  Theta  Mu  Gamma]     is order for zaxis_mocvd
%%		params   (must have params.sigma and params.tau for this spectrometer in order to calculate the angles2
%%
%% for sevchex 
%%     ANGLES2 are [Alpha  Beta  Two_Theta  PolarizationCorrection]
%%			in degrees (and Polarization correction (0-1) gives indication of how my exptl intensity dropped
%%

%%			note two spaced between angles names in order to accommodate or distinguish if 
%%				there are any stupid names with a space, (like "Two Theta" instead of "TwoTheta'
    docuangle2names =  '[Alpha  Beta  Two_Theta]';
    
if isstruct(angles); angles=angles.angles;end % just in case someone used the output of calc_hkl2angles

    Sigma 	= params.sigma;
    Tau		= params.tau;
    
	Delta = angles(:,1);
	Theta = angles(:,2);
	Mu=angles(:,3);
	Gamma=angles(:,4);
	
	%Note in finding Sigma and Tau; Thetamax = -Tau -90 if Sigma is Positive definite, Thetamax = -Tau + 90 if Sigma Negative definite
	
	sinalpha = sind(Mu).*cosd(Sigma)    + cosd(Mu).*sind(Sigma).*sind(Tau+Theta);
	sinbeta =  sind(Gamma).*cosd(Sigma) - cosd(Gamma).*sind(Sigma).*sind(Theta + Tau - Delta);
	
% Numerically, these methods of calculating TwoTheta should come out the same, but there are errors in precision cropping up?
% The cos form is more consistent with calculations with Hphi
     H = sqrt(2 .* sind(Gamma) .* sind(Mu) - 2 .* cosd(Gamma) .* cosd(Delta) .* cos(Mu) + 2);
     sinTwo_Theta_over2 = H./2;
     cosTwo_Theta = cosd(Gamma).*cosd(Delta).*cosd(Mu) - sind(Gamma).*sind(Mu);
	
	
	Alpha 	= asind(sinalpha);
	Beta 	= asind(sinbeta);
	Two_Thetac = acosd(cosTwo_Theta);
	Two_Thetas = 2 .*asind(sinTwo_Theta_over2);

% We are not outputting nhpi, but it might be useful to make sure to have it here
% However, we do not use them here and only use it
%   nphi = [sind(Sigma).*cosd(Tau); -sind(Sigma).*sind(Tau); cosd(Sigma)];   
%   nphi_hkl = [(UB.UBinv*nphi)*2*pi.*Energy./fhc]';
	    
	    
ANGLES2 = [Alpha, Beta, Two_Thetac];
end
