function [f1,p1] = pk_voigtfit(x,p,data,wt1,varyp);
%function [FIT,params] = pk_voigtfit(x,params_init,DATA,wt1,varyp);
% varyp column vector like p, 0 if don't vary param,1 if do
% x, p are initial params, DATA is vector, wt1 is weigtht
% changed for piezo stuff 0
% does leastsqaures fit 
% params = (position peakheight  FWHM  percentage_of_lorenztian)
% (that is params(4)=1 all lorenz , if 0 all gaussian 
% I don't see a background

global verbose
verbose=1;

if nargin<4, wt1 = sqrt(data); end
if nargin<5, varyp = ones(size(p)); end

% I am having trouble with letting FWHM vary
options=[ 0.0001 .01
                  0.01      0.1
                  0.001    0.08
                  0.01      0.1];
                  
dp = 0.001.*varyp;

[f1,p1,kvg1,iter1,corp1,covp1,covr1,stdresid1,Z1,r21]=...
  leasqr(x,data,p,'pk_voigt',.0001,100,wt1,dp,...
  'pk_qvtdf',options);
%  'pk_qvtdf',options);


% information about leasqr routine
%function[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]=
%                   leasqr(x,y,pin,{func,stol,niter,wt,dp,dfdp,options})
%
% Version 3.beta
% Levenberg-Marquardt nonlinear regression of f(x,p) to y(x), where:
% x=vec or mat of indep variables, 1 row/observation: x=[x0 x1....xm]
% y=vec of obs values, same no. of rows as x.
% wt=vec(dim=1 or length(x)) of statistical weights.  These should be set
%   to be proportional to (sqrts of var(y))^-1; (That is, the covaraince
%   matrix of the data is assumed to be proportional to diagonal with diagonal
%   equal to (wt.^2)^-1.  The constant of proportionality will be estimated.),
%   default=1.
% pin=vector of initial parameters to be adjusted by leasqr.
% dp=fractional incr of p for numerical partials,default= .001*ones(size(pin))
%   dp(j)>0 means central differences.
%   dp(j)<0 means one-sided differences.
% Note: dp(j)=0 holds p(j) fixed i.e. leasqr wont change initial guess: pin(j)
% func=name of function in quotes,of the form y=f(x,p)
% dfdp=name of partials M-file in quotes default is prt=dfdp(x,f,p,dp,func)
% stol=scalar tolerances on fractional improvement in ss,default stol=.0001
% niter=scalar max no. of iterations, default = 20
% options=matrix of n rows (same number of rows as pin) containing 
%   column 1: desired fractional precision in parameter estimates.
%     Iterations are terminated if change in parameter vector (chg) on two
%     consecutive iterations is less than their corresponding elements
%     in options(:,1).  [ie. all(abs(chg*current parm est) < options(:,1))
%      on two consecutive iterations.], default = zeros().
%   column 2: maximum fractional step change in parameter vector.
%     Fractional change in elements of parameter vector is constrained to be 
%     at most options(:,2) between sucessive iterations.
%     [ie. abs(chg(i))=abs(min([chg(i) options(i,2)*current param estimate])).],
%     default = Inf*ones().
%
%          OUTPUT VARIABLES
% f=vec function values computed in function func.
% p=vec trial or final parameters. i.e, the solution.
% kvg=scalar: =1 if convergence, =0 otherwise.
% iter=scalar no. of interations used.
% corp= correlation matrix for parameters
% covp= covariance matrix of the parameters
% covr = diag(covariance matrix of the residuals)
% stdresid= standardized residuals
% Z= matrix that defines confidence region
% r2= coefficient of multiple determination
%
%  {}= optional parameters
% ss=scalar sum of squares=sum-over-i(wt(i)*(y(i)-f(i)))^2.