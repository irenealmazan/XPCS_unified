function [angles] = calc_hkl2angles(hkl,UB,Energy,spectrometer_params)
% [ANGLES] = calc_hkl2angles(hkl,UB,Energy,ParamsStructure)
% version 20150522 C. Thompson (made sure output angles.UB)
% version 20130327 C. Thompson (comments and defaults tweaked 2015-05)
%
% OUTPUT  - ANGLES is structure	
%	ANGLES.angles  	(calculated geometry angles for hkl, one set of angles per row corresponding to hkl
%	ANGLES.docuparams = struct2table(spectrometer_params); (easy to 'read' names and values, 2 spaces between names)
%	ANGLES.docu  (includes Energy, spectrometer name and spectrometer_params)
%	ANGLES.UB    (carries along the UB to self document)
%	ANGLES.hkl  (carries along the hkl used for each set of angles to self document)
%	ANGLES.angles2  (values of additional angles or calculated parameters specifically of interest with spectrometer
%	ANGLES.docuangles2 (names of the data or angles in angles2, usually specific to spectrometer
%	ANGLES.docuparams = struct2table(spectrometer_params); (easy to 'read' names and values)
%
% INPUT
%	hkl (matrix or single row - but acan handle multiple queries at once, (one hkl per row) [h1 k1 l1; h2 k2 l2; ....] 
%%		(note one can also use a structure input as long as the '.hkl' field  hkl values
%		(thus, one could use output from a call to calc_angles2hkl)
%	UB is structure from [UB] = calc_UB(cparam,h1,a1,h2,a2,spectrometername)
%			or	(UB = calc_UB(OM) with OM from the readspecscan helper files)
%	Energy (eV)
%	ParamsStructure is a structure that passes spectrometer specific information
%		See what is required by looking at specific ca_h2a_{spectrometername} files
%		Use the OMspecial output of the helper files of readspecscan for your diffractometer
%					 [OM,scninfo,OMspecial] = readspec_{helperfile for your diffractometer}

%
%	ATTENTION - TO recreate the diffractometer's path during an experiment, use exact parameters 
%	the diffractometer program used (Energy, orientation matrix and lattice parameters, sigma tau)

% 
% Defaults: Energy is 28300 (eV), spectrometername= 'sevchex';


if isstruct(hkl); hkl=hkl.hkl;end % just in case someone used the output of calc_angles2hkl

	% use appropriate lowerlevel ca_hkl2angles_spectrometername for particular spectrometer called
	spectrometername = UB.spectrometername;

	calc_h2a = str2func(['ca_h2a_',spectrometername]);
	
	[ANGLES,anglenames,ANGLES2,angle2names] = calc_h2a(hkl,UB,Energy,spectrometer_params);

	% we carry hkl  outside here as each point per rows [h k l;h2 k2 l2;... (but matrix calculations need
	% them  as columns [h h2 h3... ; k k2 k3 ...; l l2 l3 ...]
	angles.angles_label = anglenames; 
	angles.angles = ANGLES;  
	angles.specialOM = spectrometer_params;
	angles.docuparams = struct2table(spectrometer_params,'AsArray',1);
	angles.hkl = hkl;
	angles.UB = UB;
	angles.Energy = Energy;
	angles.EnergyOM = UB.OM.a0E;
	angles.docu = char(...
		['Angles were calculated from given hkl on ',date],...
		['Using Energy input as ',num2str(Energy),' eV; UB documentation as follows '],...
		[UB.docu]);
	angles.angles2_label = angle2names;
	angles.angles2 = ANGLES2;
