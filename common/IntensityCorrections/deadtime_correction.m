function factor = deadtime_correction(prob)
% factor = deadtime_correction(prob)
% calculates dead-time correction factor 
% using approx. pulse pile-up formula
% prob = count rate * dead time

factor = 1 + prob + (3/2).*prob.^2 + (8/3).*prob.^3; 
