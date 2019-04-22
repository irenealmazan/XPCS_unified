function PERIOD = findperiod(A);
% PERIOD = findperiod (to get average period from successive clicks on x axis)
% Use graphics input to get period by user clicking on successive minima (or maxima)

if nargin<1; A = gcf; end

figure(A)

disp(['  Designed for any plot to find period between oscillations by graphics input'])
disp(['  Click on successive positions that define the period; when done, hit ENTER'])
disp([' '])
[X,Y] = ginput;

PERIOD = mean(diff(X));
STD = sqrt(  sum(diff(X)-PERIOD)^2)./(length(X)-1);
disp(['Average period in x axis units as calculated from these successive periods: ']),
disp(num2str(round(diff(X').*1000)./1000));
disp([' '])
disp([num2str(PERIOD),' (error ',num2str(STD),' xunits/oscill  (average)']);
disp([num2str(1./PERIOD),' oscill/ xunits (from average)']);



