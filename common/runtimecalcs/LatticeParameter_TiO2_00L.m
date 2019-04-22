function [c_angstroms, Delta_offset,MyTable] = LatticeParamTiO2(L,Phi,Lambda)
%[c, Delta] = LatticeParamTiO2(L,Phi,Lambda)   find c lattice parameter of TiO2
%
%	OUTPUT
%		c (in Angstroms) (is [c c] but both should be same assuming errors small
%		Delta (in Degrees - basically offset)
%
%	INPUT (measure 2 inplane peaks of 00L )
%		L [L1 L2]  e.g., (002) and (004) would be [2 4]
%		Phi [Phi1 Phi2;]	the phi positions for those in plane peaks
%							can have more than one set for different temperatures
%							[Phi1 Phi2; Phi1 Phi2... each row a diffent set)
%		Lambda (wavelength in Angstroms or Energy in keV
%				(I basically assume that if number is between 0 and 4 it is wavelength in Angsroms
%				and between 4 and 50 it is Energy in keV
%
%	For comparison, in the orientation matrix we said that c is 2.958, a is 4.594

if nargin<3;Lambda = 0.675663;end

if Lambda > 4;
	EkeV = Lambda;
	Lambda = 12.3984 ./ EkeV;
	disp(['Assume that you entered ' num2str(EkeV) ' in keV, so Lambda[Angstroms] is ' num2str(Lambda)])
else
	disp(char(['Assume that you entered ' num2str(Lambda) ' in Angstroms'],' '));
end

%	Thetaest = Nu./2;

	Num = L(2).*sind(Phi(:,1)) - L(1).*sind(Phi(:,2));
	Den = L(1).*cosd(Phi(:,2)) - L(2).*cosd(Phi(:,1));

	Delta_offset = atand(Num./Den);
	c_angstroms = mean(Lambda .* L./(2 .* sind (Phi + Delta_offset)),2);
	Theta_calc = Phi + Delta_offset;
	Phi_measured = Phi;

	
	STRING = char(['Calculated c[Ang]  Calculated Theta1[deg]  Calculated Theta2[deg] Calculated Offset[deg]']);
	
	MyTable = table(c_angstroms, Theta_calc, Delta_offset, Phi_measured);
	disp(MyTable)
%
