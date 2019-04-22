function [h_phi,DOCU_anglenames] =ca_hp_zaxis_mocvd(angles)
% [h_phi,DOCU_anglenames] = ca_hp_zaxis_mocvd([angles(rowvector per point)])
% version 20130327 C. Thompson 
%
% OUTPUT 
% h_phi (as column vectors, one point per column) for use in calculating the U matrix based on Busing Levy
%		hkl = UBinv*h_phi * (2pi/lambda) using UB calculated with B based on
%		a*.a=2pi/lambda 
%		(original Busing Levy based on a*.a=1/lambda);  hkl = UBinv*h_phi  * (1/lambda)
% DOCU_anglesnames (string with  name and order as assumed to calculate h_phi ) varies for each spectrometer
%		e.g., '[Rho  Nu  Delta  Mu  Eta  Chi  Phi] : sevchex';  
%		e.g., '[Delta  Theta  Mu  Gamma] : zaxis_mocvd';
%		e.g., '[Nu  Chi  Eta  Delta] : caxis_12id';
% 
% INPUT
% angles - spectrometer angles input as *row*. << important, as multiple will be 
%      one set per row. Angles must be in order as required by spectrometername
%	   convention is to code all programs to use same order as output from 'pa'
%		 command
%

% If a single set of angles, but in column instead of row, rotate
if length(angles(1,:))==1;angles=angles';end 
% Expectation for multiple angles sets MUST be one set per ROW

% use order as in output in the pa command, 
% also include spectrometer name after : as a self check in documentation
DOCU_anglenames = '[Delta  Theta  Mu  Gamma] : zaxis_mocvd';  

Delta 	= angles(:,1)';
Theta 	= angles(:,2)';
Mu	= angles(:,3)';
Gamma	= angles(:,4)';

% outputs matrix of 3xN where N is number of points, each column a different h_phi(angles)
%     h_phi = SAMP*(DET' -eye(size(DET')))* [0;1;0]; with DET and SAMP for this spectrometer plugged in

	h_phi	=	[
		cosd(Theta).*sind(Delta).*cosd(Gamma) - sind(Theta).*(cosd(Delta).*cosd(Gamma) - cosd(Mu))
		sind(Theta).*sind(Delta).*cosd(Gamma) + cosd(Theta).*(cosd(Delta).*cosd(Gamma) - cosd(Mu))
					sind(Gamma) + sind(Mu)];
end
	
% for mocvd Zaxis, this is as seen in Lohmeier and Vlieg 1993
% but corrected for a typo they have in the last factor of 1st row
% they have +cos(mu) instead of -cos(mu) (their alpha is our mu)
% their later calculation do not propagate this error, thankfully
% zaxis code at 12id in spec is also correct (C. Thompson checking 2011)
% h_phi = [
%	cosd(theta).*sind(delta).*cosd(gam) - sind(theta).*(cosd(delta).*cosd(gam) - cosd(mu))
%	sind(theta).*sind(delta).*cosd(gam) + cosd(theta).*(cosd(delta).*cosd(gam) - cosd(mu))
%					sind(gam) + sind(mu)];

%from zaxis (spec code) is correct, although ordered differently). 
% h_phi_speczaxis = [
%	sind(theta) .* (cosd(mu) - cosd(gam).* cosd(delta)) + (cosd(gam).*sind(delta).*cosd(theta))
% -1.*cosd(theta).* (cosd(mu) - cosd(gam).* cosd(delta)) + (cosd(gam).*sind(delta).*sind(theta))
%	sind(gam) + sind(mu)]; 	

