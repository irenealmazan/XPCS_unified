function [ANGLES,DOCU_anglenames,ANGLES2,DOCU_angle2names] = ca_h2a_caxis_12id(hkl,UB,Energy,params)
% [ANGLES] = ca_h2a_caxis_12id(hkl,UB,Energy,spectrometer_params)
% version 20130327 C. Thompson  
% tweak 2016-02 C Thompson to fit with other schemes better
%  Only works for mode 0 (fixed alpha)  
% OUTPUT 
%	ANGLES [Delta Theta Mu Gamma] in degrees, one row per hkl
%	DOCU_anglenames = documentation of angle names e.g., 'Nu  Chi  Eta  Delta'  
%	ANGLES2 = [Alpha, Beta, Two_Theta] in degrees, one row per hkl
%	DOCU_angle2names = documentation about the angle names e.g., 'Alpha  Beta  Two_Theta';
%
% INPUT 
%	[h k l] in r.l.u, one row per hkl
%	UB is at least UB.UB from structure from calc_UB function, but can just pass whole structure
%			(e.g., UB = calc_UB(OM)   where OM is the orietnation matrix structure
%			see calc_UB for other options to create UB
%	Energy is energy in eV
%	params is structure of parameters necessary to calculate angles from hkl
%		params.mode, 
%		params.alphatarget
%		params.sigma
%		params.tau  
%		(for now, it is hardwired to mode 0 (fixed alpha))
%	you may use just go ahead and use the specialOM for this input
%		(from the  [OM,scninfo,OMspecial] = readspec_caxis_v3_helper...
% 
% version 20160213 

% Order of output angles for caxis_12id (best practice is to use same order as from 'pa' command)
% Label with two or more spaces between each name
    DOCU_anglenames 	=	'Nu  Chi  Eta  Delta'; 
    DOCU_angle2names 	=	'Alpha  Beta  Two_Theta';

% and for calculations hkl should be in columns
	HKL = hkl';
%		|Q_xyz> = UB |hkl[r.l.u]>  / (2pi/Lambda) (gives Dimless Q_xyz 
%		with xyz fixed to phi (i.e., last circle)
	UBHKL = UB.UB * [HKL] ./ (2*pi*Energy./fhc) % Matrix of Hphi values, one per column
	
% unpack the spectrometer parameters
	mode 	= params.mode;
		if mode~=0;disp('this program currently only calculates for mode 0 (fixed alpha)'); 
		ANGLES=[];DOCU_anglesnames=[];ANGLES2=[];DOCU_angle2names=[];
		return
		end
	sigma	= params.sigma;
	tau		= params.tau;
	alphatarget	= params.alphatarget;
		
				
%	sigmatau	= [sigma  tau];
% 	hardwired caxis_12id for fixed alpha mode (mode 0)

	LEN = length(hkl(:,1));
	ANGLES = zeros(LEN,4);


% make a range of mu to use calculate and home in on which mu and theta gets to correct alpha
%	Note - unlike spec, this program will work as if alphatarget was frozen
	eta_range = ( [-2:.05:2].*(sigma + 1e-2) + alphatarget)';   
	
	LEN = length(hkl(:,1));
	ANGLES = zeros(LEN,4);  
	
	for ii=1:LEN
		
		UBHKLii = UBHKL(:,ii);
		
		% from these mu's calc what  gam, delta, theta to reach Q
		% For caxis, chi is a theta rotation, and eta has role of mu in old caxis
		%  'Nu  Chi  Eta  Delta'; 
		delangle = calc_del(eta_range,UBHKLii);
		[nuangle] = calc_nu(delangle,UBHKLii);
		[chiangle] = calc_chi(eta_range, delangle,nuangle,UBHKLii);
 
		% what is actual incident angle to surface given these angles
		% (depends only on theta and mu actually)
		%sinalpharange = calc_sinalpha_caxis_12id([deltaangle thetaangle mu_range gamangle],sigmatau);
		%  'Nu  Chi  Eta  Delta'; 
		[ANGLESrange] = calc_angles2_caxis_12id([nuangle chiangle eta_range delangle],params);
		sinalpharange = sind(ANGLESrange(:,1));

		ERROR2 = (sinalpharange - sind(alphatarget)).^2;
		% for this range - fit to A(x-x0)^2 and find x_0
		Acoef = polyfit(eta_range,ERROR2,2);

		eta_best = -Acoef(2)./(Acoef(1).*2);   
		del_best = calc_del(eta_best,UBHKLii);
		nu_best = calc_nu(del_best,UBHKLii);
		chi_best = calc_chi(eta_best, del_best, nu_best,UBHKLii);

		ANGLES(ii,:) = [nu_best,chi_best,eta_best,del_best];
		
		[ANGLES2] 	= calc_angles2_caxis_12id(ANGLES(ii,:),params);
		ALPHA(ii) 	= ANGLES2(:,1);
		BETA(ii) 	= ANGLES2(:,2);

		% if curious on methodology used to find best chi and eta - 
		% change the following from "if 0" to  "if 1" 
		% (it makes as many figures as there are hkl input
		% so do this only when there are just a few being calculated!)
		% If region does not include minimum, may need to make range
		% larger above. (if sigma is really large)
		if 1;
			figure
			plot(eta_range, ERROR2,'o',...
				eta_range,polyval(Acoef,eta_range),'b-',...
						[eta_best eta_best],[-.2 .2].*max(ERROR2),'r-');
			xlabel(['eta trials,  best eta at ' num2str(eta_best), 'deg']);
			ylabel('Error [(sin(calc alpha) - sin(target alpha)]^2');
			title(['for hkl at [',num2str([hkl(ii,:)]),'] and sigma tau [',...
				 num2str([sigma tau]),']']);
		end

		
	end

[ANGLES2,docuangle2names] = calc_angles2_caxis_12id(ANGLES,params);
		
end

%%%%%%%%%%%%%%  HELPER FUNCTIONS  %%%%%%%%%%%%

function [delangle] = calc_del(eta,UBHKL)

	H2o2 = sum(UBHKL.*UBHKL)./2;

	delangle = asind( (UBHKL(3,:) - sind(eta).*H2o2) ./ cosd(eta) );
	% implicit for mode 0 have to choose + or - del, choose +


	

end

function [nuangle] = calc_nu(del,UBHKL)

	H2o2 = sum(UBHKL.*UBHKL)./2;

	sinnu = UBHKL(1,:) ./ cosd(del);
	cosnu = (1 - H2o2) ./  cosd(del);
%	nuangle = atan2d(sinnu, cosnu); % May need error check for not real
%				also needs sign constraint 
%	nuangle = sign(nuangle).*nuangle; 
	% implicit for mode 0 have to choose + or - nu, should choose +
%	nuangle = acosd(cosnu);
	
	nuangle = asind(sinnu);
 


	
end

function [chiangle] = calc_chi(eta,del,nu,UBHKL)
	
	H2o2 = sum(UBHKL.*UBHKL)./2;
if 1
        argy1 = UBHKL(1,:) .* ( sind(del).*sind(eta) - H2o2.*cosd(eta) );
	argy2 = UBHKL(2,:) .* ( cosd(del).*sind(nu) );
	argx1 = UBHKL(2,:) .* ( sind(del).*sind(eta) - H2o2.*cosd(eta) );
	argx2 = UBHKL(1,:) .* ( cosd(del).*sind(nu) );	
	SS = +1;
end	

% Fister paper (probably not correct)
if 0
        a1 = cosd(eta) - cosd(nu) .* cosd(del) .* cosd(eta) ...
            - sind(del) .* sind(eta);
        a2 = cosd(eta) .* sind(nu);
        a3 = cosd(eta) .* sind(del) + sind(eta) ...
            - cosd(nu) .* cosd(del) .* sind(eta);
        argy1 = a1 .* UBHKL(1,:);
        argy2 = a2 .* UBHKL(2,:);
        argx1 = a2 .* UBHKL(1,:);
        argx2 = a1 .* UBHKL(2,:);
        SS = -1;
end

	chiangle = SS.*atan2d((argy1 + argy2),(argx1 - argx2));	
end
%		/* Calculate chi from diffraction equation */
%		t1 = cos(del) * sin(nu);
%		t2 = cos(nu) * cos(del) * cos(eta) + sin(del) * sin(eta) - cos(eta);
%		t3 = hphi[1] * t1 - hphi[0] * t2;
%		t4 = hphi[0] * t1 + hphi[1] * t2;
%		if ( mode == ZAXIS_ALBE )
%			chi = CHI*RAD; 	// Do not change chi if using specular mode.
%		else 
%			chi = atan2(t3, t4);



