
% Plotting with resolution function

% Getting the number of points of resolution
% dtheta is angular resolution of spectrometer
% dL/L = dq/q = dtheta*cos(theta)/sin(theta)
% Therefore, dL = dtheta*L*cos(theta)/sin(theta) (NOTE L=1 for our scans)
%			also note that sin(theta)=lambda/2asub;

%%dtheta=.007*pi/180;	% dtheta in radians
%%Nres=round(dtheta*1*2*asub*sqrt(1-lambda*lambda/(4*asub*asub))/(lambda*del) );

% Then smoothing and plotting
%%Eaveres = bsmooth(abs(Eave),Nres);
%%figure(ip+3),cvdplot(Eaveres,QZZ,Fpos(II),N,Asub*1e8,Afilm*1e8);
%%text(.825,1e8,['lres = ',num2str(del*Nres),'  of Eave']);
%%text(.825,.2e8,['theta res = ',num2str(dtheta),' rads']);
