function [h_phi,DOCUanglenames] =ca_hp_caxis_12id(angles)
% h_phi = ca_hp_caxis_12id([angles(rowvector)])
% version 20120504 C. Thompson (should work, but I need to recheck signs of rotations)
%
% Computes h_phi (column vectors) for use in calculating the U matrix based on Busing Levy
%  hkl = UBinv*h_phi   * (2pi/lambda) using UB calculated with B based on
%  a*.a=2pi/lambda 
%	(original Busing Levy based on a*.a=1/lambda), in that convention  hkl = UBinv*h_phi  * (1/lambda)
% However, angles must be in order in its vector as required by spectrometername
%	convention is to use same order as output from 'pa' command since this
% 	is also used typically for order to enter in the orientation matrix.
%		e.g., for 'zaxis_mocvd',  [Delta, Theta, Mu, Gamma]
%		e.g., for 'caxis_12id', [Nu Chi Eta Delta Phi(ignored) Mu(ignored)]
%

% for caxis_12id spectrometer, [Nu  Chi  Eta  Delta  Phi(ignored)  Mu(ignored)]
Nu	 	= angles(:,1)';
Chi	 	= angles(:,2)';
Eta		= angles(:,3)';
Delta		= angles(:,4)';
DOCUanglenames = '[Nu  Chi  Eta  Delta  Phi(ignored)  Mu(ignored)] : caxis_12id';

% outputs matrix of 3xN where N is number of points, each column a different h_phi(angles)
%     h_phi = SAMP*(DET' -eye(size(DET')))* [0;1;0]; with DET and SAMP for this spectrometer plugged in
%			see calc_ROTOPER_caxis_12id for separate DET and SAMP operators defined
% ct 12-2012 checkng my last paper calculations, change sign 1st link to opposite -cosd(Chi) ,,,, +sind(Chi)
if 0
	h_phi	= [
	-cosd(Chi).* cosd(Delta).* sind(Nu)  ...
		+ sind(Chi).*  ( cosd(Eta).* ( cosd(Delta).* cosd(Nu) - 1 ) + sind(Delta).* sind(Eta) );
	sind(Chi).* cosd(Delta).* sind(Nu)  ...
		+ cosd(Chi).* ( cosd(Eta).* ( cosd(Delta).* cosd(Nu) - 1 ) + sind(Delta).* sind(Eta) );
	sind(Delta).* cosd(Eta) - sind(Eta).* ( cosd(Delta).* cosd(Nu) - 1 )
			];
end
if 1
	h_phi	= [
	+cosd(Chi).* cosd(Delta).* sind(Nu)  ...
		- sind(Chi).*  ( cosd(Eta).* ( cosd(Delta).* cosd(Nu) - 1 ) + sind(Delta).* sind(Eta) );
	sind(Chi).* cosd(Delta).* sind(Nu)  ...
		+ cosd(Chi).* ( cosd(Eta).* ( cosd(Delta).* cosd(Nu) - 1 ) + sind(Delta).* sind(Eta) );
	sind(Delta).* cosd(Eta) - sind(Eta).* ( cosd(Delta).* cosd(Nu) - 1 )
			];
end
if 0
% from Fister paper - this is (after reduction) identical to mine except that his 1st component of h_phi is the negative of mine.
% This is not correct (in my opinion) when the Nu and Chi rotation operators are properly set up (they are actually RH coordinates, and
% all axes (Nu, Chi, Eta, and Delta) on the front huber is CCW about all axes x,y,z. 
% his caxis program might put a negative in front of the eta and delta (as if CW about x even though it is not?)
        a1 = cosd(Eta) - cosd(Nu) .* cosd(Delta) .* cosd(Eta) ...
            - sind(Delta) .* sind(Eta);
        a2 = cosd(Delta) .* sind(Nu);
        a3 = cosd(Eta) .* sind(Delta) + sind(Eta) ...
            - cosd(Nu) .* cosd(Delta) .* sind(Eta);
        h_phi = [
            a1 .* sind(Chi) + a2 .* cosd(Chi)
            a1 .* sind(Chi) - a1 .* cosd(Chi)
            a3];
end


end



