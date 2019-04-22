function [h_phi,DOCU_anglenames] =ca_hp_sevchex(angles)
% [h_phi,DOCU_anglenames = ca_hp_sevchex([angles(rowvector per point)])
% modified from ca_hp_zaxis_mocvd version 20130327 C. Thompson
% 28-SEP-14 GBS
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


% use order as output in the pa command, 
% use string names as they are output as headers in data file
% also include spectrometer name after : as a self check in documentation
% require 2 spaces between names
%DOCU_anglenames = '[Rho  Nu  Delta  Mu  Eta  Chi  Phi] : sevchex';
DOCU_anglenames = 'Rho  Nu  Delta  Mu  Eta  Chi  Phi';  

% If a single set of angles, but in column instead of row, rotate
if length(angles(1,:))==1;angles',end 
% Expectation for multiple angles sets MUST be one set per ROW

Rho 	= angles(:,1)';
Nu      = angles(:,2)';
Delta 	= angles(:,3)';
Mu	= angles(:,4)';
Eta 	= angles(:,5)';     
Chi	= angles(:,6)';
Phi     = angles(:,7)';

% outputs matrix of 3xN where N is number of points, each column a different h_phi(angles)
%     h_phi = SAMP*(DET' -eye(size(DET')))* [0;1;0]; with DET and SAMP for this spectrometer plugged in

h_phi = [ ...
(-sind(Eta).*sind(Phi) - sind(Chi).*cosd(Eta).*cosd(Phi)).*sind(Rho) ...
+ (cosd(Eta).*cosd(Mu).*sind(Phi) + (cosd(Chi).*sind(Mu) - sind(Chi).*sind(Eta).*cosd(Mu)).*cosd(Phi)).*cosd(Rho) ...
+ (-cosd(Delta).*cosd(Eta).*sind(Mu).*sind(Nu) - cosd(Delta).*cosd(Eta).*cosd(Mu).*cosd(Nu) - sind(Delta).*sind(Eta)).*sind(Phi) ...
+ ((sind(Chi).*cosd(Delta).*sind(Eta).*sind(Mu) + cosd(Chi).*cosd(Delta).*cosd(Mu)).*sind(Nu) ...
+ (sind(Chi).*cosd(Delta).*sind(Eta).*cosd(Mu) - cosd(Chi).*cosd(Delta).*sind(Mu)).*cosd(Nu) ...
- sind(Chi).*sind(Delta).*cosd(Eta)).*cosd(Phi); ...
(sind(Eta).*cosd(Phi) - sind(Chi).*cosd(Eta).*sind(Phi)).*sind(Rho) ...
+ (-cosd(Eta).*cosd(Mu).*cosd(Phi) + (cosd(Chi).*sind(Mu) - sind(Chi).*sind(Eta).*cosd(Mu)).*sind(Phi)).*cosd(Rho) ...
+ (cosd(Delta).*cosd(Eta).*sind(Mu).*sind(Nu) + cosd(Delta).*cosd(Eta).*cosd(Mu).*cosd(Nu) + sind(Delta).*sind(Eta)).*cosd(Phi) ...
+ ((sind(Chi).*cosd(Delta).*sind(Eta).*sind(Mu) + cosd(Chi).*cosd(Delta).*cosd(Mu)).*sind(Nu) ...
+ (sind(Chi).*cosd(Delta).*sind(Eta).*cosd(Mu) - cosd(Chi).*cosd(Delta).*sind(Mu)).*cosd(Nu) ...
- sind(Chi).*sind(Delta).*cosd(Eta)).*sind(Phi); ...
cosd(Chi).*cosd(Eta).*sind(Rho)+ (sind(Chi).*sind(Mu) + cosd(Chi).*sind(Eta).*cosd(Mu)).*cosd(Rho) ...
+ (sind(Chi).*cosd(Delta).*cosd(Mu) - cosd(Chi).*cosd(Delta).*sind(Eta).*sind(Mu)).*sind(Nu) ...
+ (-sind(Chi).*cosd(Delta).*sind(Mu) - cosd(Chi).*cosd(Delta).*sind(Eta).*cosd(Mu)).*cosd(Nu) ...
+ cosd(Chi).*sind(Delta).*cosd(Eta)];

end
	

