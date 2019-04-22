function [hkl] = calc_angles2hkl(angles,UB,Energy)
%[hkl] = calc_angles2hkl(angles, UB, Energy in eV)
% version 20150522 C. Thompson (made sure output angles.UB)
% version 20130327 C. Thompson, tweaked comments and defaults
%
% Current versions assume Energy of OM vectors and of angles is same.
% To change, one would need to carry Energy of each OM vector and Energy
%	of current angles, and scale the h_phi.
%
% OUTPUT
% hkl is structure	
%	hkl.hkl  (rows of hkl per each angle request) 
%	hkl.angles (angles used)
%	hkl.docu   
%	hkl.UB
%
% INPUT (plain except for UB)
%	angles  (may be structure with field .angles, or a plain matrix, 
%				each row a new set of angles, in the order required for spectometer
%		e.g.,	zaxis_mocvd [Delta  Theta  Mu  Gamma]
%				caxis_12id	[Chi  Nu  Eta  Delta]
%				sevchex		[Rho  Nu  Delta  Mu  Eta  Chi  Phi]
%				in associated low level files, standard practice is to use angle order
%				that is the same as in a 'pa' output in spec
% 	UB is structure as calculated calc_UB program : [UB] = calc_UB(cparam,h1,a1,h2,a2,spectrometername)
% 	Energy (eV). For calculating hkl from angles, this may be different energy than that used for finding the OM
% 	spectrometername is a string that indicates the spectometer  name, exact name is important- 
%			i.e., there should be function "ca_hp_{spectrometername}" defined specifically
%			e.g., if spectrometername is 'zaxis_mocvd', then ca_hp_zaxis_mocvd.m must be present and correct!

% Defaults: Energy is 28300 (eV), spectrometername= 'sevchex'
% 


if nargin<3,Energy=28300;end

spectrometername = UB.spectrometername;

% and in case someone put in the output of calc_hkl2angles to test things
if isstruct(angles); angles = angles.angles; end

	% use the h_phi calculations as function of angle specific for the spectrometer
	calc_hphi = str2func(['ca_hp_',spectrometername]);
	
	[h_phi,anglenames] = calc_hphi(angles);

	% we carry hkl in and out as point per rows [h k l;h2 k2 l2;... (but matrix calculations will always need them as columns
	%  [h h2 h3... ; k k2 k3 ...; l l2 l3 ...]
	hkl.hkl = [(UB.UBinv*h_phi)*2*pi.*Energy./fhc]';
	hkl.angles = angles;
	hkl.angles_label = anglenames;
	hkl.spectrometername = spectrometername;
	hkl.Energy = Energy;
	hkl.EnergyOM = UB.OM.a0E;
	hkl.docu = char(...
		['hkl were calculated from given angles on ',date],...
		['Using Energy as ',num2str(Energy),' eV ; UB documentation as follows '],...
		[UB.docu]);
	hkl.UB = UB;
	

