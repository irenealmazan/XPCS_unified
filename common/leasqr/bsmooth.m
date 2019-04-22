function [spsdsum]= bsmooth(psdsum,npts);
% [spsdsum]= smooth(psdsum,npts)
% Smooth data in columns of psdsum
% by convolving with Gaussian of FW npts
% 05-jan-95 gbs

% Generate Gaussian;

gauss=[-2*npts:2*npts];
gauss=((2/npts)*(log(2)/pi)^(1/2))*exp(-((2*(log(2))^(1/2))*(gauss/npts)).^2);

% Loop through columns of psdsum;

[row,col] = size(psdsum);
spsdsum = zeros(row,col);

for i = 1:col

   rvector = psdsum(:,i)';

% Extend rvector by 2*npts on either end;

   evector = [rvector(1)*ones(1,2*npts) rvector rvector(row)*ones(1,2*npts)];

% Calculate convolution;

   svector = conv(evector,gauss);

% Remove extra points at ends;

   extra = 4*npts;
   spsdsum(:,i) = svector([extra+1:extra+row])';

end
