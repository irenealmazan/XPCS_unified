function [ANGLES2,docuangle2names] = calc_angles2_caxis_12id(angles,params)
%% [ANGLES2,docuangle2] = calc_angles2_caxis_12id(angles,params)
%% 
%%  OUTPUT
%%		ANGLES2  (degrees) Alpha Beta TwoTheta
%%		docuangles2 (string with names of the data or angles in angles2 with 2 spaces between each name)
%%  INPUT
%%		angles  (each row a new set of angles, in the order required for spectometer)
%%					[Nu  Chi  Eta Delta]     is order for caxis_12id
%%		params   (must have params.sigma and params.tau for this spectrometer in order to calculate the ANGLES2
%%
%% for caxis 
%%     ANGLES2 are [Alpha  Beta  Two_Theta]
%%			in degrees 
%%

%%			note two spaced between angles names in order to accommodate or distinguish if 
%%				there are any stupid names with a space, (like "Two Theta" instead of "TwoTheta'
    docuangle2names =  '[Alpha  Beta  Two_Theta]';
    
if isstruct(angles); angles=angles.angles;end % just in case someone used the output of calc_hkl2angles

    Sigma 	= params.sigma;
    Tau		= params.tau;
    % expected order 'Nu  Chi  Eta  Delta'; 
	Nu	= angles(:,1);
	Chi	= angles(:,2);
	Eta	=angles(:,3);
	Delta	=angles(:,4);
	
	%Note in finding Sigma and Tau; 
	%  Thetamax = -Tau -90 if Sigma is Positive definite,   (zaxis)
	%  Thetamax = -Tau + 90 if Sigma Negative definite   (caxis)

	
	sinalpha = sind(Eta).*cosd(Sigma)  + cosd(Eta).*sind(Sigma).*sind(Tau-Chi);
		
	sinbeta =  cosd(Sigma).*(  sind(Delta).*cosd(Eta) - cosd(Delta).*sind(Eta).*cosd(Nu)  ) + ...
		sind(Sigma)*(   sind(Chi-Tau).* ( sind(Delta).*sind(Eta) +  cosd(Delta).*cosd(Eta).*cosd(Nu) ) ...
			- cosd(Delta).*cosd(Chi-Tau).*sind(Nu)    );
			

	
% Numerically, these methods of calculating TwoTheta should come out the same, but there are errors in precision cropping up?
% The cos form is more consistent with calculations with Hphi
     H = sqrt(2 - 2 .* cosd(Delta) .* cosd(Nu) );
     %cosTwo_Theta = cosd(Delta).*cosd(Nu);
	 %sinTwo_Thetaover2 = H ./2;
	
	Alpha 	= asind(sinalpha);
	Beta 	= asind(sinbeta);
	Two_Thetac = acosd(cosd(Delta).*cosd(Nu));
	Two_Thetas = 2 .*asind(H ./ 2);

% We are not outputting nhpi, but it might be useful to make sure to have it here
% However, we do not use them here and only use it
%   nphi = [sind(Sigma).*cosd(Tau); -sind(Sigma).*sind(Tau); cosd(Sigma)];   
%   nphi_hkl = [(UB.UBinv*nphi)*2*pi.*Energy./fhc]';
	    
	    
ANGLES2 = [Alpha, Beta, Two_Thetac];
end
