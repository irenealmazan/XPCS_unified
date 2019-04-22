function [sigmatau] = find_sigmatau(eta,phi,rho,ALPHAguess)
%function [sigmatau] = find_sigmatau(eta,phi,rho,alpha_assumed)
%
%  use other function CW for zaxis, sevchex: if sigma = |sigma| then, tau is -PHImax-90
%		(tau definition OK for zaxis, and sevchex) 
%  use CCW for caxis, sevchub: if sigma = |sigma| then, tau is +PHImax-90
%		
%
% Function to fit sin of experimental ETA positions and PHI and alpha
% to the corresponding TAU and SIGMA (from the position of the maximum
% of the sin function and the magnitude.
% NOTE IT WILL GIVE BACK THE BEST FIT ALPHA FROM DATA
%
% Equivalents (eta is rotation about X, phi is rotation about Z, and Phi circle must be on top of Eta)
%  eta	phi	rho	sevc (hex or hub)
%  mu	th	0	zaxis (old mocvd)
%  eta	chi	0	caxis (used on front huber)
%
% ******* EVERYTHING IS IN DEGREES ******
% INPUT ARGUMENTS
%	eta is the vector of eta values
%	phi is the vector of phi values
%	rho is the setting of rho in sevc (use value of 0 if using caxis or zaxis)
%	alpha_assumed is the guessed alpha we think we are using in the alignment setup
%
% OUTPUT ARGUMENTS   is the structure sigmatau
%	sigmatau.sigma		sigma (in degrees)
%	sigmatau.tau		tau (in degrees)
%	sigmatau.alphacalc	alpha (incident angle) calculated from the curve
%	sigmatau.phimax		the rotation angle that had maximum tilt (used to calculate tau)
%	sigmatau.docu		various documentation (date, etc)
%
% 2015_0426  C. Thompson

if nargin<3;rho=0;rhoFLAG=0;rhoDOCU=[];else rhoFLAG=1;end;
if nargin<4;ALPHAguess = rho + mean(eta);end

[AMPLITUDE,PHIMAX,ALPHA]=ca_tausig(eta,phi,rho,ALPHAguess);

if rhoFLAG;
	rhoDOCU = ['used a value of rho = ',num2str(rho),' deg'];
end

sigmatau.sigma		= AMPLITUDE;
sigmatau.tau		= PHIMAX-90;  %CCW
	if sigmatau.tau < -180 | sigmatau.tau > 180;
	%then make more pretty 
		if sigmatau.tau < -180; sigmatau.tau = sigmatau.tau+360;
		else sigmatau.tau = sigmatau.tau-360;
		end
	end
sigmatau.alphacalc	= ALPHA;
sigmatau.phimax		= PHIMAX;
sigmatau.docu		= strvcat(date,...
	[rhoDOCU],...
	['alpha (assumed for starting the fit) was ',num2str(ALPHAguess),' deg'],...
	['alpha calculated from data given is ',num2str(sigmatau.alphacalc),' deg'],...
	['tau calculated from PHI_max-90  (CCW rotations), sigma is tilt ']); 

%%% the rest is just pretty plotting
phifull = [-170:10:170];

%subplot(2,1,1)
%
%HH=plot(phi,eta + rho - ALPHA,'b+',phifull,AMPLITUDE.*sind(1.*phifull - (PHIMAX-90)),'r-');
%HHxaxis = line([-sigmatau.tau -sigmatau.tau],[-AMPLITUDE AMPLITUDE]);
%
%	set(gca,'xlim',[min(phifull) max(phifull)]);
%	xlabel(['phi']),ylabel(['eta+rho-alpha_{calc}']);
%	TITLEstr = strvcat(['Tau  = ',num2str(sigmatau.tau),' (=-90-Max, Max at ',num2str(PHIMAX),' )'],...
%		[' TILT (sigma) = ',num2str(AMPLITUDE)],...
%		[' alpha_{calc} = ',num2str(ALPHA), ': alpha_{guess} = ',num2str(ALPHAguess)]);
%	title(TITLEstr);
%	set(HH(1),'markersize',10);
%	set(HH(2),'linewidth',1.5);
%	grid on
%subplot(2,1,2)
%set(gcf,'paperposition',[1.5 2.5 4 7]);

HH=plot(phi,eta,'b+',phifull,AMPLITUDE.*sind(phifull - (PHIMAX-90)) + ALPHA - rho,'r-');
HHxaxis = line([-sigmatau.tau -sigmatau.tau],[-AMPLITUDE AMPLITUDE]+ALPHA-rho);

	set(gca,'xlim',[min(phifull) max(phifull)]);
	xlabel(['phi']),ylabel(['eta [w/alpha@ ',num2str(ALPHA),']']);
	TITLEstr = strvcat(['Tau  = ',num2str(sigmatau.tau),' (=MAX-90, Max at ',num2str(PHIMAX),' )'],...
		[' TILT (sigma) = ',num2str(AMPLITUDE)],...
		[' alpha_{calc} = ',num2str(ALPHA), ': rho = ',num2str(rho)]);
	title(TITLEstr);
	set(HH(1),'markersize',10);
	set(HH(2),'linewidth',1.5);
	grid on

set(gcf,'paperposition',[1.5 2.5 4 3.5]);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AMPLITUDE,PHIMAX, ALPHA] = ca_tausig(eta,phi,rho,alphaguess)

% ensure angle input is in columns [phi eta rho]
rhoinput = rho.*ones(size(eta(:)));
angles = [phi(:) eta(:) rhoinput];

%Fitting program finds the phase shift, amplitude, and offset of the sin fit to the data
%  (sigma = amplitude)
%  (alpha = offset)
%  (phaseshift is eventually related to value of Tau
%		TAU = -90-phi_max  <<< this is the important one
%		phaseshift = phi_max-90;
%		phaseshift = -TAU-180;
%

% Approximate starting values
%		phaseshift = PHI_max - 90 
%		first guess amplitude is eta(value at phi_max) - alpha
PHIMAX0  =  max(phi(find(eta==max(eta))));

Xamplitude = max(eta)-min(eta);
Xoffset = alphaguess;
Xphase = (PHIMAX0 - 90);

X0 = [Xamplitude;Xphase;Xoffset];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = fminsearch(@(x) sindiff(x,angles),X0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AMPLITUDE	= X(1);
PHASESHIFT	= X(2);
ALPHA		= X(3);
PHIMAX		= PHASESHIFT + 90;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function OUTPUT = sindiff(Xguess,angles)
%function fitdiff = sindiff([AMPLITUDE;PHASE])
% Finction used with fmins to calculate the sin wave
% of eta with respect to phi (assumes that eta already has
% alpha subtracted

phi = angles(:,1);
eta = angles(:,2);
rho = angles(:,3);

eta_wo_rho = eta + rho;

Xamplitude	= Xguess(1);	% related to sigma
Xphase		= Xguess(2);    % eventually related to tau
Xoffset		= Xguess(3);	% alpha

EXPT = eta_wo_rho - Xoffset;
CALC = Xamplitude.*sind(phi - Xphase);
fitdiff= (EXPT-CALC).^2;
OUTPUT = sum(fitdiff);
end
