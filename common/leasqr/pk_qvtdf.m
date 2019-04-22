function prt=pk_qvtdf(x,f,lam,dp,func)
% function prt=pk_qvtdf(x,f,lam,dp,func)
% calculate partial derivatives of pseudo-voigt
% for use with leasqr.m
% 5-3-94 SMB (Bren@SLAC.stanford.edu)
leng=length(lam);
npts=length(x);
npeaks= floor(leng/4);
gu=zeros(npts,4);
fl=zeros(npts,4);
pl=zeros(npts,4);
g=zeros(npts,4);
f=zeros(npts,4);
pld=zeros(npts,1);
pgu=zeros(npts,1);
prt=zeros(npts,4);
for i=1:npeaks,
	base=(i-1)*4;
	gu(:,i)= (x-lam(base+1)) ./(0.600561*lam(base+3));
	fl(:,i)= (x-lam(base+1)) ./(0.5*lam(base+3));
	pl(:,i)= 1+fl(:,i).^2;
	g(:,i) = exp(-gu(:,i).^2);
	f(:,i)= lam(base+2) .*((lam(base+4)./pl(:,i))+((1-lam(base+4)).*g(:,i)));
	pld = (pl(:,i) - 2.*fl(:,i)*lam(base+4))./pl(:,i).^2;
	pgu = 2.*gu(:,i).*(1.-lam(base+4)).*g(:,i);
	prt(:,base+1) = lam(base+2).*(pgu./(0.600561*lam(base+3)) - 2.*pld./lam(base+3));
	prt(:,base+2) = f(:,i)./lam(base+2);
	prt(:,base+3) = lam(base+2).*(gu(:,i).*pgu - fl(:,i).*pld)./lam(base+3);
	prt(:,base+4) = lam(base+2).*(1./pl(:,i) - g(:,i));
end
prt(:,leng-1) = ones(npts,1);
prt(:,leng) = x;
