function gerr= inbounds(xpos,ypos);
% gerr= inbounds(xpos,ypos);
% Check to see whether graphics input is in axis limits
% gerr is true if cursor is within axis limits.

% Copyright 1996 Anneli Munkholm & Sean M. Brennan.
% Stanford Synchrotron Radiation Laboratory
% Stanford Linear Accelerator Center, Stanford CA 94309
% Bren@slac.stanford.edu

limits= axis;
gerr= 1;       % Assume position is in bounds
if xpos > limits(2) | xpos < limits(1), gerr= 0; end
if ypos > limits(4) | ypos < limits(3), gerr= 0; end
