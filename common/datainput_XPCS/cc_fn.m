function [ffit] = cc_fn(xfit,pp)
% function [ffit] = cc(xfit,pp)
% Used by 2-time correlation analysis to fit GaN m-plane diffuse data
% pp(1) constant in first factor
% pp(2) linear coef of i+j in first factor
% pp(3) exp decay constant of |i-j|
% pp(4) parameter to set matrix size
% 30-JUN-17 GBS


[ii,jj] = ind2sub([pp(4),pp(4)],xfit);

tau = max(0.01,pp(3));
ffit = (pp(1) + pp(2)*(ii+jj)).*exp(-abs(ii-jj)/tau);

return