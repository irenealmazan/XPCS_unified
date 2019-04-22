function [ANGLES2,docuangle2names] = calc_angles2_sevchex(angles,params)
%% [ANGLES2,docuangle2] = calc_angles2_sevchex(angles,params)
%%  2015 - need to fix two_theta (typo) CT CT CT
%%  OUTPUT
%%		ANGLES2  (calculated of other angles or parameters per hkl for angles as calculated  - typically Alpha Beta TwoTheta
%%		docuangles2 (string with names of the data or angles in angles2)
%%  INPUT
%%		angles  (each row a new set of angles, in the order required for spectometer
%%		params   (must have params.sigma and params.tau	
%%
%% for sevchex 
%%     ANGLES2 are [Alpha  Beta  Two_Theta  PolarizationCorrection]
%%			in degrees (and Polarization correction (0-1) gives indication of how my exptl intensity dropped
%%


    docuangle2names =  '[Alpha  Beta  Two_Theta  PolarizationCorrection]';
    
if isstruct(angles); angles=angles.angles;end % just in case someone used the output of calc_hkl2angles

    Sigma 	= params.sigma;
    Tau		= params.tau;

% There are a variety of methods using H_phi and nphi that are used as necessary to 
% calculate Alpha, Beta and Two_Theta, however, if all the angles are known, and sigma tau, then
% it can be done as follows.
% However, we do not use them here and only use it
%   nphi = [sind(Sigma).*cosd(Tau); -sind(Sigma).*sind(Tau); cosd(Sigma)];   
%   nphi_hkl = [(UB.UBinv*nphi)*2*pi.*Energy./fhc]';

    Rho = angles(:,1);
    Nu = angles(:,2);
    Delta = angles(:,3);
    Mu = angles(:,4);
    Eta = angles(:,5);
    Chi = angles(:,6);
    Phi = angles(:,7);
    
% There are a variety of ways to find alpha and beta and two theta
% The following are just with the angles
    
% Numerically, these methods of calculating TwoTheta should come out the same, but there are errors in precision cropping up?
% The cos form is more consistent with calculations with Hphi
     H = sqrt(2 .* sind(Delta) .* sind(Rho) - 2 .* cosd(Delta) .* cosd(Nu) .* cos(Rho) + 2);
     sinTwo_Theta_over2 = H./2;
     cosTwo_Theta = cosd(Delta).*cosd(Nu).*cosd(Rho) - sind(Delta).*sind(Rho);
    
%% Polarization correction  (PolarizationCorr*100 is percent (Theoretical*Corr = Experimental)
    PolarizationCorr = (cosd(Delta).*cosd(Nu)).^2 + sind(Delta).^2;
    
    
    sinalpha = ...
cosd(Sigma) .* ( ...
	+cosd(Chi) ...
		 .* (cosd(Eta) .* sind(Rho) .* (1-cosd(Mu))+sind(Eta+Rho) .* cosd(Mu)) ...
	+sind(Chi) .*  cosd(Rho) .*  sind(Mu)) ...
	+ ...
sind(Sigma) .* ( ...
	+cosd(Chi) .*  cosd(Phi+Tau) .* cosd(Rho) .*  sind(Mu) ...
	-sind(Chi) .*  cosd(Phi+Tau) ...
		 .* ( cosd(Eta) .* sind(Rho) .* (1-cosd(Mu)) +sind(Eta+Rho) .* cosd(Mu)) ...
	+sind(Phi+Tau) ...
	 .* (+sind(Eta) .* sind(Rho) .* (cosd(Mu)-1) + cosd(Eta+Rho) .* cosd(Mu)));


sinbeta = ...
(((-sind(Chi) .* cosd(Delta) .* sind(Eta) .* sind(Mu)-cosd(Chi) .* cosd(Delta) .* cosd(Mu)) .* sind(Nu)+(cosd(Chi) .* cosd(Delta) .* sind(Mu)-sind(Chi) .* cosd(Delta) .* sind(Eta) .* cosd(Mu)) .* cosd(Nu)+sind(Chi) .* sind(Delta) .* cosd(Eta)) .* ... 
sind(Phi)+(-cosd(Delta) .* cosd(Eta) .* sind(Mu) .* sind(Nu)-cosd(Delta) .* cosd(Eta) .* cosd(Mu) .* cosd(Nu)-sind(Delta) .* sind(Eta)) .* cosd(Phi)) .* sind(Sigma) .* sind(Tau)+( ...
(-cosd(Delta) .* cosd(Eta) .* sind(Mu) .* sind(Nu)-cosd(Delta) .* cosd(Eta) .* cosd(Mu) .* cosd(Nu)-sind(Delta) .* sind(Eta)) .* sind(Phi)+ ...
((sind(Chi) .* cosd(Delta) .* sind(Eta) .* sind(Mu)+cosd(Chi) .* cosd(Delta) .* cosd(Mu)) .* sind(Nu)+(sind(Chi) .* cosd(Delta) .* sind(Eta) .* cosd(Mu)-cosd(Chi) .* cosd(Delta) .* sind(Mu)) .* cosd(Nu)-sind(Chi) .* sind(Delta) .* cosd(Eta)) .* ... 
cosd(Phi)) .* sind(Sigma) .* cosd(Tau)+ ...
((sind(Chi) .* cosd(Delta) .* cosd(Mu)-cosd(Chi) .* cosd(Delta) .* sind(Eta) .* sind(Mu)) .* sind(Nu)+(-sind(Chi) .* cosd(Delta) .* sind(Mu)-cosd(Chi) .* cosd(Delta) .* sind(Eta) .* cosd(Mu)) .* cosd(Nu)+cosd(Chi) .* sind(Delta) .* cosd(Eta)) .* ...
cosd(Sigma);

  Alpha = asind(sinalpha);
  Beta = asind(sinbeta);
  Two_Thetac = acosd(cosTwo_Theta);
  Two_Thetas = 2 .*asind(sinTwo_Theta_over2);

% We are not outputting nhpi, but it might be useful to make sure to have it here
% However, we do not use them here and only use it
%   nphi = [sind(Sigma).*cosd(Tau); -sind(Sigma).*sind(Tau); cosd(Sigma)];   
%   nphi_hkl = [(UB.UBinv*nphi)*2*pi.*Energy./fhc]';
	    
	    
ANGLES2 = [Alpha, Beta, Two_Thetac, PolarizationCorr];
