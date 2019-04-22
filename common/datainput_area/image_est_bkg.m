function [Muse] = findbkg(M,dim2fit,ndx_shoulder)
% [Muse] = findbkg(M,dim2fit,ndx_shoulder)
% dim2fit   every row is a profile to get shoulders and line 'X'
% dim2fit   every col is a profile to get shoulders and line 'Y'
% M  3d matrix (3rd dimension to mean over (only have the background)
% ndx_shoulder = [NL1 NL2; NR1 NR2]   left shoulder, right shoulder
% assume pixel in imageJ/spec notation

%if nargin<3;
	dim2mean=3;
%end
% Assume third dimension is the dimension to sum
% ndx are those to mean
% dim2mean typically is 3 usually and M3d already parsed

% convert to matlab convention
ndx_shoulder = ndx_shoulder+1;
ndx_full = [min(ndx_shoulder(:)):max(ndx_shoulder(:))];
ndx_s = [ndx_shoulder(1,1):ndx_shoulder(1,2) ndx_shoulder(2,1):ndx_shoulder(2,2)];

Muse = mean(M,dim2mean,'omitnan');

[NYcol,NXrow,steps] = size(Muse);

if strcmp(dim2fit,'X');  %fit line under peak through rows, X)
	for ii=1:NYcol
		pii = polyfit(ndx_s,Muse(ii,ndx_s),2);
		Muse(ii,ndx_full) = polyval(pii,ndx_full);
	end
else %fit line under peak through columns, Y)
	for ii=1:NXrow
		pii = polyfit(ndx_s',Muse(ndx_s,ii),2);
		Muse(ndx_full,ii) = polyval(pii,ndx_full');
	end	
end

% now, need to 'fit through' the peak
