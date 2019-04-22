function ROIS = makeroi(N,HF)
% MYROIS = myroichoice(Number,figuretouse)
% only does rectanular - still need to make to do trapezoid
% assumes round to integer 

if nargin<2;HF = gcf; end
if nargin<1;N=1;end

figure(HF)
magenta = [1 0 1];
cyan = [0 1 1];

ROIS=ones(N,4);
disp(['to input an ROI - use the cursor and click on the two corners']);
jj = 1;FLAG=0;
for ii = 1:N;
	if FLAG==0;
		[X,Y]=ginput(2);
		% Assume pixel counting so round to integer when drawing and creating
		if ~isempty(X);
			X=round(X);Y=round(Y);
			HL = line([X(1) X(2) X(2) X(1) X(1)],[Y(1) Y(1) Y(2) Y(2) Y(1)]);
			set(HL,'linewidth',1,'color',magenta);

			ROIS(jj,:) = [min(X) max(X)  min(Y) max(Y)];
			jj=jj+1;
		else
			FLAG=1;
		end
	end
end

N = length(ROIS(:,1));

if 0
for ii=1:N
	disp(['ROI (' int2str(ii) ')']) ;
	disp(['   X(col) (pixels): ',int2str(MYROIS(ii,1:2)),']']);
	disp(['   Y(row) (pixels): ',int2str(MYROIS(ii,3:4)),']']);
end
end
STRING = strvcat('ROIs chosen are the following',...
' minX maxX  minY maxY  ',...
' ',...
['MYROIS = [...'],...
int2str(ROIS),'  ];');

disp(STRING)

