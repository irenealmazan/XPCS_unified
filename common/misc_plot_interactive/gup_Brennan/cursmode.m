function cursmode(program)
% cursmode(program)
% Controls the mouse mode (cursor when within axes else pointer)

% Copyright 1996 Anneli Munkholm & Sean M. Brennan.
% Stanford Synchrotron Radiation Laboratory
% Stanford Linear Accelerator Center, Stanford CA 94309
% Bren@slac.stanford.edu

p=get(gca,'CurrentPoint');
if strcmp(get(gca,'visible'),'off')
	return
end
if strcmp(get(gcf,'pointer'),'arrow')   % mouse in pointer mode
	if inbounds(p(1,1),p(1,2)), 
		set(gcf,'pointer','crosshair')
		set(gcf,'WindowButtonDownFcn',[program '(''curs'')'])
		if strcmp(program,'peakfit')
		    global PEAK_XPOS
			PEAK_XPOS=[];
		end
	end
else								    % mouse in cursor mode
    if ~inbounds(p(1,1),p(1,2)), 
		set(gcf,'pointer','arrow')
		set(gcf,'WindowButtonDownFcn','')
	end
end
