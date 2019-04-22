function gup(command)
% Put up gui buttons to control plot features
% such as x or y scales linear or log, zoom feature
% or set limits using edit controls.  To find position
% of cursor, click mouse with cursor in plot window.
% to remove gui controls, type gup again.

% Copyright 1996 Anneli Munkholm & Sean M. Brennan.
% Stanford Synchrotron Radiation Laboratory
% Stanford Linear Accelerator Center, Stanford CA 94309
% Bren@slac.stanford.edu

global GUP_BUFF
if nargin == 0
	if isempty(findobj(gcf,'tag','gup_xlog'));
		command = 'init';
	else
		command = 'exit';
	end
end 
if strcmp(command,'init')
	if size(axis,2)==6
		gupinit3;         % 3D plot
	else
		gupinit;          % 2D plot
	end 
	GUP_BUFF=[];
elseif strcmp(command,'xlog') 
	current= get(findobj(gcf,'tag','gup_xlog'),'string');
	if current(4) == 'i',			% as in linear
		set(gca,'xscale','linear');
		set(findobj(gcf,'tag','gup_xlog'),'string','X Log');
	else
		set(gca,'xscale','log');
		set(findobj(gcf,'tag','gup_xlog'),'string','X Linear');
	end
	guplmts;
elseif strcmp(command,'ylog') 
	current= get(findobj(gcf,'tag','gup_ylog'),'string');
	if current(4) == 'i',			% as in linear
		set(gca,'yscale','linear');
		set(findobj(gcf,'tag','gup_ylog'),'string','Y Log');
	else
		set(gca,'yscale','log');
		set(findobj(gcf,'tag','gup_ylog'),'string','Y Linear');
	end
	guplmts;
elseif strcmp(command,'zlog') 
	current= get(findobj(gcf,'gup_tag','zlog'),'string');
	if current(4) == 'i',			% as in linear
		set(gca,'zscale','linear');
		set(findobj(gcf,'tag','gup_zlog'),'string','Z Log');
	else
		set(gca,'zscale','log');
		set(findobj(gcf,'tag','gup_zlog'),'string','Z Linear');
	end
	guplmts;
elseif strcmp(command,'zmax')
	limits= axis;
	limits(6)= str2num(get(findobj(gcf,'tag','gup_zmax'),'string'));
	axis(limits); 
elseif strcmp(command,'zmin') 
	limits= axis;
	limits(5)= str2num(get(findobj(gcf,'tag','gup_zmin'),'string'));
	axis(limits); 
elseif strcmp(command,'ymax')
	limits= axis;
	limits(4)= str2num(get(findobj(gcf,'tag','gup_ymax'),'string'));
	axis(limits); 
elseif strcmp(command,'ymin') 
	limits= axis;
	limits(3)= str2num(get(findobj(gcf,'tag','gup_ymin'),'string'));
	axis(limits); 
elseif strcmp(command,'xmax') 
	limits= axis;
	limits(2)= str2num(get(findobj(gcf,'tag','gup_xmax'),'string'));
	axis(limits); 
elseif strcmp(command,'xmin') 
	limits= axis;
	limits(1)= str2num(get(findobj(gcf,'tag','gup_xmin'),'string'));
	axis(limits); 
elseif strcmp(command,'curs') 
	cu_pos=get(gca,'CurrentPoint');
	s=sprintf('X= %-8.5g   Y= %-8.5g',cu_pos(1,1),cu_pos(1,2));
	disp(s)
	GUP_BUFF=[GUP_BUFF; cu_pos(1,1) cu_pos(1,2) ];
elseif strcmp(command,'zum')
	gupzoom;
elseif strcmp(command,'zout')
	axis('auto');
	guplmts;
elseif strcmp(command,'azm') 
	az=get(findobj(gcf,'tag','gup_azm'),'value');
	set(findobj(gcf,'tag','guptxt_azm'),'string',int2str(round(az)));
	el=get(findobj(gcf,'tag','gup_elv'),'value');
	view(az,el)
elseif strcmp(command,'elv')
	az=get(findobj(gcf,'tag','gup_azm'),'value');
	el=get(findobj(gcf,'tag','gup_elv'),'value');
	set(findobj(gcf,'tag','guptxt_elv'),'string',int2str(round(el)));
	view(az,el)
elseif strcmp(command,'exit')
	uicontrol=findobj(get(gcf,'children'),'Type','uicontrol');
	for j=1:length(uicontrol)
		tag=get(uicontrol(j),'tag');
		if strcmp(tag(1:3),'gup'), delete(uicontrol(j)); end
	end
	set(gcf,'pointer','arrow', ...
			'WindowButtonMotionFcn','',...
			'WindowButtonDownFcn','');
else 
% probably false initialize 
end 

