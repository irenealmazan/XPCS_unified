function f=pk_voigt(x,lam);
% f= pk_voigt(x,lam)
% actual function to calculate single voigt peak
% 5-2-94 SMB (Bren@SLAC.stanford.edu)
f=zeros(size(x));
delx=x-lam(1);
gau= exp(-(delx ./(0.600561*lam(3))).^2);
lor= 1+ (delx ./(0.5*lam(3))).^2;
f= f+ lam(2) .*((lam(4)./lor)+((1-lam(4)).*gau));
