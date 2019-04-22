function [UB] = calc_UB(cparam,h0,a0,h1,a1,spectrometername)
% [UB]	= calc_UB(cparam,h0,a0(or0 angles),h1,a1(or1 angles),spectrometername)
% [UB]	= calc_UB(OM)	
%			with OM structure (fields .cparam, .h0, .a0, etc
% [UB]	= calc_UB(cparam)
%			to only calculate the B matrix and rparam (reciprocal space)
% version 20160213 C. Thompson
% 
%
% OUTPUT
	% [UB] is structure  of UB.B, UB.UBinv, UB.UB, UB.U, UB.docu, UB.cparam, UB.rparam
	% 	(if called with only cparam, then [UB] is structure of UB.B, UB.cparam, UB.rparam, UB.docu
%
% INPUT
	% cparam is a,b,c,alpha, beta, gamma of unit cell (in angstroms and degrees)
	% h0,a0 and h1,a1 are the [h k l] and [angles (degrees)...] for the two orientation matrix vectors
	% expectation, cparam, a0, h0, a1, h1 are ROW vectors 
	% 	spectrometer name is a string parameter for spectrometer, e.g., 'sevchex'; 
%
	% At present calc_UB assumes that both angles sets for the vectors were obtained at same Energy, and that the U (and UB) desired is at the same energy. 
%
%  Requires function ca_h2p_{spectrometername}.m for each spectrometer 
%
% Spectrometername - (*) not implemented yet
%	zaxis_mocvd 	[Delta  Theta  Mu  Gamma] 
%	(*)zaxis		classic zaxis; surface aligned to'omega'; 4 circles
%	caxis_12id		[Nu  Chi  Eta  Delta  Phi=0  Mu=0] 
%	(*)psic_12id	[Nu  Chi  Eta  Delta  Phi  Mu]  (Chi 0 is rotated 90 degrees in Eta compared to caxis_12id)				
%	(*)4circle		[TwoTheta  Theta  Chi  Phi]  (rught handed)
%	sevchex			[Rho  Nu  Delta  Mu  Eta  Chi  Phi] 
%	(*)sevchub		[Rho  Nu  Delta  Mu  Eta  Chi=0  Phi] 
%
% Compute the matrix U and UB based on 
% Busing and Levy, 1967 Acta Cryst. 22, p457-464
%	(but using 2pi/d instead of 1/d)

% NOTE UB is not unitary because B is not unitary
% to retrieve hkl[r.l.u.] from a set of angles
%		|u_phi(angles)> = ca_hp_{spectrometername}  (phi the 'last geo circle')
%		|hkl[r.l.u]> = (UBinv) |u_phi(angles)>*2pi/lambda
%		|Q_phixyz> = UB|hkl(r.l.u)>
%				note the B operator has the dimensions (units of 2pi/lambda)
%				thus Binv operator has dimensions lambda/2pi
%		|H_phixyz> = |Q_phi>/(2pi/lambda) (dimless HxHyHz on last geo circle)
%
% 		If spectrometer is at 0,0,0,0 ..... then xyzPhi and xyzlab are typically aligned
%		with xyz is fixed to lab (i.e., spectrometer at 0,0,0,0,0,...)
%		+z vertical up, +y(incident beam downstream), +x(horizontal, RH axes)
% Use rotations (SAMP and DET to determine with arbitrary angles of spectrometer)
%		SAMPinverse | H_phixyz>  = |H_lab xyz>
%		etc
%
% If using the B calc identical to Busing Levy, the the UB calculated is used as follows
%		|hkl[r.l.u]> = (UBinv) |u_phi(angles)>/lambda
%		|Q_xyz> = UB |hkl[r.l.u.]>*2*pi 
% 
% B |hkl[r.l.u.]> = Q_x'y'z' where x'y'z' is fixed to crystal
%		x' is along h, 
%		y' is in hk plane but perp to x'
%		z' is perpendicular to hk plane
%		If using our B, Q will be normal |Q| = 4pi sin(2theta/2)/lambda, 
%			if using Busing and Levy defined B, |Q| is 2 sin(2theta/2)/lambda units

% Can pass the OM structure from the readspec helper scans
% Note - While energy of each set of angles (or0 and or1) is carried via this OM structure - at prsent calc_UB assumes that both vectors were obtained at same Energy - it does not relatively scale.
if isstruct(cparam)  % use the OM structure if passed
	OM=cparam;
	cparam = OM.cparam;
	a0 = OM.a0;
	a1 = OM.a1;
	h0 = OM.h0;
	h1 = OM.h1;
	spectrometername = OM.spectrometername;
		if isfield(OM,'a0E');
			a0E = OM.a0E;
			a1E = OM.a1E;
		end	
elseif nargin==2 % using calc_UB(cparam,spectrometername)	
	spectrometername = h0;
elseif nargin==6
% create OM from input parameters
	OM.cparam = cparam;
	OM.h0	= h0;
	OM.h1	= h1;
	OM.a0	= a0;
	OM.a1	= a1;
	OM.spectrometername = spectrometername;
else  % just cpararm
	OM = []; % not a structure
end

UB.B = calc_B(cparam);

if (isstruct(OM)| nargin>2); 

	T_cryst 	= calc_Tcryst(cparam,h0,h1);
	[T_phi,anglenames] = calc_Tphi(cparam,a0,a1,spectrometername);
	
	if ~isfield(OM,'angledocu')
		OM.angledocu = anglenames;
	end

	UB.U	= T_phi*(T_cryst');
	UB.UB	= UB.U*UB.B;

	% B is not unitary, nor is UB so must take inverse! and not shortcut via transpose!!!
	UB.UBinv	= UB.UB\eye(3,3);   %UB.UB\eye(3,3) marginally faster than inv(UB.UB)
	UB.cparam	= cparam;
	UB.rparam	= calc_recip(cparam);
	UB.spectrometername = spectrometername;
	UB.OM		= OM;

		D0 = ['UB requested for ',spectrometername];
		D1 = [' '];
		D2 = ['UnitCell = [ ',num2str(cparam(:)'), ']'];
		D3 = ['or0(hkl) = ( ',num2str(h0(:)') ' )'];	
		D4 = ['or1(hkl) = ( ',num2str(h1(:)') ' )'];
		D5 = anglenames;	
		D6 = ['or0(angles) = ',num2str(a0(:)')];
		D7 = ['or2(angles) = ',num2str(a1(:)')];
	UB.docu = strvcat(D0,D1,D2,D3,D4,D5,D6,D7);

else   % if only needed the B matrix, only would put in a cpararm
	UB.U	= [];
	% UB.B already exists, calculated with cparam 
	UB.UB	= [];
	UB.UBinv	= [];
	UB.cparam	= cparam;
	UB.rparam	= calc_recip(cparam);
	UB.OM 		= OM;		%% empty
	UB.docu		= ...
		['B and rparam only  - calculated for UnitCell [',num2str(cparam(:)'),']'];
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rparam=calc_recip(cparam)
% calculates the reciprocal lattice vectors from the direct
% cparam = [a b c alpha beta gamma]   in angstrom and degrees
% rparam = [a* b* c* alpha* beta* gamma*]  in inverse angstrom and degrees
% 	note - using a dot a* = 2pi definitions (a.b* or a.c* = 0)
% OK works if cparam row or column vector

	cosalphap = ...
		(cosd(cparam(5))*cosd(cparam(6))-cosd(cparam(4)))/(sind(cparam(5))*sind(cparam(6)));
	cosbetap = ...
		(cosd(cparam(4))*cosd(cparam(6))-cosd(cparam(5)))/(sind(cparam(4))*sind(cparam(6)));
	cosgammap= ...
		(cosd(cparam(5))*cosd(cparam(4))-cosd(cparam(6)))/(sind(cparam(5))*sind(cparam(4)));

% 2015 error found CT these work for hexagonal and cubic, orthorhombic, tetragonal
% within our group (SRS) our papers using the calc_hkl have been on hexagonal, and cubic rlu
%	astar	=2*pi/(cparam(1)*sind(cparam(5))*sind(cparam(6)));
%	bstar	=2*pi/(cparam(2)*sind(cparam(6))*sind(cparam(4)));
%	cstar	=2*pi/(cparam(3)*sind(cparam(4))*sind(cparam(5)));
	V = sqrt(1 - cosd(cparam(4)).^2 - cosd(cparam(5)).^2 - cosd(cparam(6)).^2 ...
		+ 2 .* cosd(cparam(4)).*cosd(cparam(5)).*cosd(cparam(6)));
	
	astar	= 2*pi.*sind(cparam(4)) ./ (cparam(1) .* V);
	bstar	= 2*pi.*sind(cparam(5)) ./ (cparam(2) .* V);
	cstar	= 2*pi.*sind(cparam(6)) ./ (cparam(3) .* V);
	
	alphastar	= acosd(cosalphap);
	betastar	= acosd(cosbetap);
	gammastar	= acosd(cosgammap);

rparam = [astar bstar cstar alphastar betastar gammastar];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T_cryst = calc_Tcryst(cparam,h0,h1)
% function T_cryst = calc_Tcryst(cparam,or0(hkl),or1(hkl))
% OK works if vectors row or column vector
% Compute the matrix Tcryst from  2 orientation vectors according Busing
% and Levy
% 	vector 1 is parallel to or0
%	vector 2 is orthogonal to or0, but in plane of or0 and or1
%	vector 3 is orthogonal to or0 and vector 2

	B=calc_B(cparam);

% put orientation (hkl) vectors into the cartesian reciprocal lattice representation
	h0_cryst=(B*(reshape(h0,3,1)));
	h1_cryst=(B*(reshape(h1,3,1)));
% Then calculate Triplet of vectors and create matrix from it
	T_cryst=calc_tri(h0_cryst,h1_cryst);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Barray,Barrayinv]=calc_B(cparam)
% function [Barray,Barrayinv]=calc_B(cparam)
% OK works if cparam row or column vector
% Compute the matrix B according to Busing and Levy, which
% represents reciprocal lattice vectors of unit cell by cartesian coordinate system
% with the x same as astar, and the y perpendicular to x, but in plane of astar bstar,
% and the z axis perpendicular to astar and cstar (or x and y)
%  |HKL> = B|hkl> will supposedly convert an hkl in hkl representation to this cartesian representation
% cparam=[a b c alpha beta gamma];
%
% B is not a unitary matrix THUS B' is not equal to inv(B)
% However, B retains the magnitudes in reciprocal space
%	qvector = B x hvector, where hvector is the indices of the reflection (hkl)
%		and ha*+kb*+lc*= the vector in reciprocal space to that hkl point.
% so,  the length of  vector |ha*+ kb*+lc*| = |qx+qy+qz|
% Busing Levy paper uses reciprocal lattice vectors units of 1/a not 2pi/a (a*.a=1, etc)
% We have altered B to be consistent with 'modern' convention, (a*.a=2pi/a).
% This  alters the B matrix - but visibly only for the 3,3 entry --> 2pi/a(3) 
% note for the other entries of matrix, the factor of 2pi is already present 
% due to it being used in the calculation for a*, b*, c* used in the calculations.

% real space (a)
	a = cparam(1:3);
	alpha = cparam(4:6);
% reciprocal space  (b)  (calc_recip 2pi/d units, not 1/d units)
	rparam 	= calc_recip(cparam);
	b(:,1:3) 	= rparam(1:3);
	beta	= rparam(:,4:6);

	Barray=	[	b(1)	b(2)*cosd(beta(3))	b(3)*cosd(beta(2));
		      		0 	b(2)*sind(beta(3))	-b(3)*sind(beta(2))*cosd(alpha(1));
				0	0	2.*pi./a(3)];
				
% To use output of calc_recip to make Barray identical to Busing and Levy
% 	b(:,1:3) = rparam(1:3)./(2*pi);    % if wanted to calcu		 
%	Barray=[ 	b(1)	b(2)*cosd(beta(3))	b(3)*cosd(beta(2));
%		      		0 	b(2)*sind(beta(3))	-b(3)*sind(beta(2))*cosd(alpha(1));
%			 	0	0	1./a(3)];			 

Barrayinv = Barray\eye(3,3);	%UB.B\eye(3,3) marginally faster than inv(UB.B)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T_phi,anglenames] = calc_Tphi(cparam,a0,a1,spectrometername)
% function [Tphi,anglenames] = calc_Tphi(cparam,h0(angles),h1(angles),spectrometername)
% OK works if vectors row or column vector
% Compute the matrix T_phi according to Busing and Levy, which
% transform a or0 and or1 (expressed in the experimentally determined angles
% to the representation fixed to the diffractometer
% 	Calculate T_phi unitary matrix based on orientation vectors (in diffractometer angles)
% 	in orthonormalized basis set attached to diffractometer
%	anglenames is string with names of angles and order coded e.g., '[Delta Theta Mu Gamma]'

% only need to get the angle names and order once for particular spectrometer
% The ca_hp functions REQUIRE angle inputs as row vectors to work correctly
calc_hphi = str2func(['ca_hp_',spectrometername]);
	[h0_phi,anglenames]	= calc_hphi(a0(:)');
	h1_phi			= calc_hphi(a1(:)');	
% calculate triplet	
	T_phi = calc_tri(h0_phi,h1_phi);

end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tri_array=calc_tri(vector1,vector2)
% function array=calc_tri(vector1,vector2)
% OK works if vectors row or column vector
% procedure triple noted in Busing and Levy
% Generates a rotation matrix with columns equal to the orthogonal unit
% vectors t1,t2 and t3 which in a right-handed system chosen 
% Note, the input (vector1 and vector2) must be expressed in units of an 
%	orthonormal coordinate system, which is done through operation of B on the hkl.
% t1 is parallel to vector1
% t2 lies in the plane of vector1 and vector2
% t3 is perpendicular to that plane.

% make the three vectors (orthonormal, and t1, t2, t3, right handed)
	t1=vector1./norm(vector1); 
	t3=cross(vector1,vector2);
		t3=t3./norm(t3); 
	t2=cross(t3,t1);
		t2=t2./norm(t2);

% ensure  vectors are  columns no matter whether initial as column or row
	t1=reshape(t1,3,1);
	t2=reshape(t2,3,1);
	t3=reshape(t3,3,1);

% create unitary matrix with orthonormal set of vectors as each column
	tri_array = [t1 t2 t3];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

