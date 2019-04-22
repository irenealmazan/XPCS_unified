function [y] = NB_fn(x,p)
% function NB_fn(x,p)
%   Negative n=binomial distribution
%   p(1) = kbar
%   p(2) = M

kb = max(0,p(1));
M = max(1,p(2));
y = (gamma(x+M)./(gamma(M).*gamma(x+1))) .* (1 + M/kb).^(-x) .* (1 + kb/M)^(-M);

end

