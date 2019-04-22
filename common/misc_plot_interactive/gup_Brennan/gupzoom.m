function gupzoom(action)
% rubberband zoom version for gup

% Copyright 1996 Anneli Munkholm & Sean M. Brennan.
% Stanford Synchrotron Radiation Laboratory
% Stanford Linear Accelerator Center, Stanford CA 94309
% Bren@slac.stanford.edu

if nargin<1, action='start'; end;

if strcmp(action,'start');
    ax=axis;   % get the extents of the current axis
    set(gcf,'pointer','cross', ...
			'WindowButtonMotionFcn','', ...
            'WindowButtonDownFcn','gupzoom(''down'');');
    % Now the figure will trap the buttondown, and zooming will happen
elseif strcmp(action,'down');
    % User has pressed down the mouse button.
		ax=axis; % get the extents of the current axis
        cpa=get(gca,'CurrentPoint');
        x1=cpa(1,1);y1=cpa(1,2);
        % We've got the first point, now the rbbox
        % We have to jump through hoops here, since the
        % RBBOX function returns nothing and updates
        % only the root object's PointerLocation.
        CurUnitsA=get(gca,'Units');     % store for later
        set(gca, 'Units', 'pixels');    % Want axis pos iin pixels
        fpa=get(gca, 'Position');       % Got it
        set(gca, 'Units', CurUnitsA);   % Reset axis units
        CurUnitsF=get(gcf,'Units');     % Same thing for figure.
        set(gcf, 'Units', 'pixels');
        fpf=get(gcf, 'Position');
        ppux=fpa(3)/(ax(2)-ax(1));      % Pixels per unit x-axis
        ppuy=fpa(4)/(ax(4)-ax(3));      % Pixels per unit y-axis
        cpf=get(gcf, 'CurrentPoint');   % Original zoom box
        set(gcf, 'Units', CurUnitsF);   % Restore figure units
        cp0=get(0, 'PointerLocation');  % Always in pixels
        rbbox([cpf 0 0], [cpf]);        % Call the RBBOX
        cp0fin=get(0, 'PointerLocation');   % get other box extent
		dxfin=cp0fin(1)-fpf(1)-fpa(1); 
		dyfin=cp0fin(2)-fpf(2)-fpa(2);
		gxs=get(gca,'Xscale');			% Get xscale (liniear or log)
		gys=get(gca,'Yscale');			% Get yscale (liniear or log)
		if gxs(1:3)=='lin' 
			x2=dxfin ./ ppux + ax(1); 
        else							% xaxis: Log scale
			delta=dxfin ./ (ppux*(ax(2)-ax(1)));
			x2=10^(delta*(log10(ax(2))-log10(ax(1)))+log10(ax(1)));
		end 							% yaxis: Log scale
		if gys(1:3)=='lin' 
			y2=dyfin ./ ppuy + ax(3); 
        else
			delta=dyfin ./ (ppuy*(ax(4)-ax(3)));
			y2=10^(delta*(log10(ax(4))-log10(ax(3)))+log10(ax(3)));
		end 
        x=[x1 x2];x1=min(x);x2=max(x);  % sort in min / max order
        y=[y1 y2];y1=min(y);y2=max(y);
		ax=[x1 x2 y1 y2];
		axis(ax);            % Set axis
        set(gcf,'pointer','arrow', ...
				'WindowButtonMotionFcn','cursmode(''gup'')', ...
			    'WindowButtonDownFcn','');
%	end;
	guplmts;
end; 
