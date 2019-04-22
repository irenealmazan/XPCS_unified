function [Energycalc,deltacalc] = energycalib_2Ddet(THETA, HKL)
% Choosing which crystal 2016_02 (since crystalstr gets interpreted as numbers, still
% need to work on getting ot to print out the name choice that was used
% [Energy(keV),NUoffset(deg)] = energycalib([PHI(in degrees)], [HKL], crystalstring)
% Energy calibration from surface peaks using 2 inplane peaks with 2D detector
%
% THETA is THETA motor (or equivalent rotation about z) position at the peaks. 
%	on old zaxis (mocvd) it would be th
%	on caxis (huber) it would be chi
%	on sevchex (mocvd) and sevchub (huber) it would be phi
%
% for H K L, input [2 0 0;6 0 0] (not as [200;600]
% Example: energcal([13.9034; 27.162],[2 0 0; 6 0 0],'gan');
%
% NOTE THESE SHOULD BE IN PLANE PEAKS OTHERWISE ANALYSIS
% THAT WE USE IS NOT VERY MEANINGFUL.
% THETA values should be very carefully determined! Be Accurate!


% YES, THETA is motor position at the peaks. 
% Note, this is not really 'theta' as needed in
%		nL = 2d sin(2theta/2). 
% It is offset by an unknown but constant angle 
% (however, it is constant offset only when comparing inplane peaks) 
% 
% We know this, because when we look at the sample at THETA=0, 
% edges are not parallel with the beam direction.
% Also, even if they were, the substrates are not cut perfectly.
% so there is still an offset.
%
% Thus; theta = THETA - offset
%    sin(theta) = sin(THETA)cos(offset) - cos(THETA)sin(offset)
% 
% Therefore, what we are doing here is finding Lambda and delta
% that satisfy two equations for i=1,2 (two peaks)
% 	Lambda = 2d_i * (sin(THETA_i)cos(offset) - cos(THETA_i)sin(offset))
% Unfortunately, doesn't quite work up into nice little matrix equation
% where all we have to do is invert matrix....
%
%	let t=Lambda; x=sin(offset), y=cos(offset)=sqrt(1-x^2);
%	Let A_i = 2d_i*sin(THETA_i); B_i = 2d_i*cos(THETA_i)
%
%       Use matlab to figure out intersection of two lines to find x and t
% 	t = A_i*sqrt(1-x^2) - B_i.*x;	
%
%    L (lambda) = hc/Energy ==> 12398.4/Energy(eV) = L (angstroms)
%

% First see if the user can follow instructions
if length(HKL(:))<=2 
    disp('check the HKL input - why do you think that I could understand');
    disp('that, for example,  200 (two hundred) means 2 0 0');
    disp(' I need them separated and in a column, e.g., [2 0 0;4 0 0]');
    return
end 


%define crystal calculation type
CUBIC = 'CUB';HEXAGONAL = 'HEX';

%% here is where to add more crystals
CRYSTALS = strvcat(...
'gan1997 = [3.186 3.186 5.178]; Lattice = HEXAGONAL; %a=b; c (hexagonal)',...
'gan2016 = [3.188 3.188 5.186]; Lattice = HEXAGONAL; %a=b; c (hexagonal)',...
'al2o3 = [4.758 4.758 12.99]; Lattice = HEXAGONAL; %a=b; c (hexagonal)',...
'sto = [3.905 3.905 3.905]; %a=b=c; Lattice = CUBIC; %a=b=c (cubic)');


% Put in the choices for later
for ii=1:length(CRYSTALS(:,1));
    eval(CRYSTALS(ii,:));
end

% User gets to choose the crystal
disp('The following are the current crystal choices');
    disp(CRYSTALS);
    crystalstr=input('input choice (e.g., >> gan) >> ');
	if isempty(crystalstr);return;end
    a=crystalstr(1);b=crystalstr(2);c=crystalstr(3);

% start the business, finally

H = HKL(:,1);K=HKL(:,2);L=HKL(:,3);


if Lattice==CUBIC
%%%%%%% For cubic lattice  %%%%%
Dspacing = a./(sqrt(H.^2 + K.^2 + L.^2));  
Ddocu = strvcat(['Crystal ',crystalstr],'d_{hkl} = a/sqrt(h^2 + k^2 + l^2)', ['a= ',num2str(a),' Ang used for calculating  d from r.l.u.s ']);

%%%%%%% For OM with GaN (hexagonal) at room temperature %%%%
elseif Lattice==HEXAGONAL
Dspacing = 1./sqrt( 4./3 .* (H.^2 + H.*K + K.^2)./a^2 + L.^2./c.^2);
Ddocu = strvcat(['Crystal ',crystalstr],'inplane d_{hkl} = 1./sqrt( 4./3 .* (h^2 + hk + k^2)./a^2 + l^2./c.^2)', ['a and b= ',num2str(a),' Ang used for calculating  d from r.l.u.s '],['c = ',num2str(c)]);
end

DOCUmain = 'Using 2D detector - detector angles are not accurate. Get good phi angles';

deltaestimated  = 15; %degrees offset, this is estimated theta=THETA(motor)-delta
%		to get it, look at sample when THETA=0, how does it need to be
%		rotated to get the azimuth parallel to the incident beam?

% probably do not need to change anything below here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = HKL(:,1);K=HKL(:,2);L=HKL(:,3);
EVAL_doc = strvcat(DOCUmain,...
	Ddocu,[' '],...
	'peaks used were,',...
	['d = ',num2str(Dspacing(1)),' Angstrom: HKL= (',int2str(HKL(1,:)),') : phi = ',num2str(THETA(1)),' deg '],...
	['d = ',num2str(Dspacing(2)),' Angstrom: HKL= (',int2str(HKL(2,:)),') :  phi = ',num2str(THETA(2)),' deg ']);


% convert to radians and make sure column vector
THETAr = reshape(THETA .*pi./180,2,1);
hc=12398.4;

%	let t=Lambda; x=sin(delta), y=cos(delta)=sqrt(1-x^2);
%	Let A_i = 2d_i*sin(THETA_i); B_i = 2d_i*cos(THETA_i)
% 	t = A_i*sqrt(1-x^2) - B_i.*x;	
Ai = 2.*Dspacing.*sin(THETAr);
Bi = 2.*Dspacing.*cos(THETAr);

estdelr = deltaestimated.*pi./180;
sinestdelr = sin(estdelr);

A1 = Ai(1);A2=Ai(2);B1=Bi(1);B2=Bi(2);
% I just learned how to use anonymous functions
% This provides input into fzero that finds intersection
% (or when the difference between two lines is zero)
anonfunction = @(x) ((A1-A2).*sqrt(1-x^2) - (B1-B2).*x);
sindelcalc = fzero(anonfunction,sinestdelr);

Energycalc = hc./((A1.*sqrt(1-sindelcalc.^2) - B1.*sindelcalc).*1000); %keV
deltacalc = asin(sindelcalc).*180./pi;

% to see it graphically (but only works if initial estimate of delta was close
if 0
XXdeg = deltaestimated.*[.8:.001:1.2];
XX = sin(XXdeg.*pi./180);
YY = sqrt(1-XX.^2);
plot(XXdeg,hc./(A1.*YY-B1.*XX).*10,'r-',XXdeg,hc./(A2.*YY-B2.*XX).*1000,'g-');
xlabel(['delta[deg] (theta=PHI-delta)']);
ylabel(['Energy (keV)']);
end   

if 1;
EVAL_doc = strvcat(datestr(now),EVAL_doc,...
	['  '],...
	['theta(peak) = phi(motor)-OFFSET; OFFSET = ',num2str(deltacalc),' degree'],...
	[' (OFFSET can be useful for determining HKL directions on sample if edges squared at phi@zero) '],...
	['Energy = ',num2str(Energycalc),' keV: Wavelength (ang) = ',num2str(hc*1e-3/Energycalc)]);
delete information.txt
diary('information.txt')
disp(EVAL_doc);
diary off
disp(' ');
disp('****');
disp([' to print the file information.txt to document out the above']);
disp([' on linux box, the following probably works within matlab']);
disp([' ! lpr information.txt']);
end
