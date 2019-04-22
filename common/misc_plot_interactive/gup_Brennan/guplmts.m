function guplmts
% limits for gup

% Copyright 1996 Anneli Munkholm & Sean M. Brennan.
% Stanford Synchrotron Radiation Laboratory
% Stanford Linear Accelerator Center, Stanford CA 94309
% Bren@slac.stanford.edu

limits=axis;
set(findobj(gcf,'tag','gup_xmin'),'string',num2str(limits(1)));
set(findobj(gcf,'tag','gup_xmax'),'string',num2str(limits(2)));
set(findobj(gcf,'tag','gup_ymin'),'string',num2str(limits(3)));
set(findobj(gcf,'tag','gup_ymax'),'string',num2str(limits(4)));
if length(limits)==6,  % 3D-figure
   set(findobj(gcf,'tag','gup_zmin'),'string',num2str(limits(5)));
   set(findobj(gcf,'tag','gup_zmax'),'string',num2str(limits(6)));
end
