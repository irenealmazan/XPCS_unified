% Analysis of Aug 2017 MPX flatfield data
% Read summed data created by make_flatfield.m  (small alterations 
% 29APR17 GBS

%% CT change here to use the flatfield from 2017-08
load FlatField_2017_0809_2_num64;
% Flatfield_2500 is the init (lo) thresholds, 10 sec, avg of 2500
% Flatfield_alt_10 is the Alt (hi) thresholds, 10 sec, avg of 500
% Flatfield cs 20 is the CS thresholds, 10 sec, avg of 500

%Fimage = Flatfield_alt_10.flatfield;
%disp('[ Using the flatfield with alt thresholds] ');

Fimage = FlatField.flatfield;
disp('[ Using the flatfield with thresholds as set for Aug 2017 run] ');

% plot some histograms to choose thresholds for bad pixels

ImageJ = 1; % use imageJ indexing [start at 0] when denoting ROIs and when plotting pixels
%	This applies to AXISdet and ROI values so must be consistent.
% Also in time, first image is 0

YROWpts = [1:size(Fimage,1)]-ImageJ;  
XCOLpts = [1:size(Fimage,2)]-ImageJ;
CLIM = [0 1000];

figure
HS = pcolor(XCOLpts,YROWpts,Fimage);
shading flat
xlabel('XCOL');ylabel('YROW');
title('Fimage');
prettyplot(gca);
colormap jet;
colorbar;
set(gca,'clim',CLIM);

% Flatfield ROI
ROIXCOLpts = [280 : 450]; % plotting X values
ROIYROWpts = [280 : 420]; % plotting Y values
colsndx = ROIXCOLpts + ImageJ; % matlab subscripts
rowsndx = ROIYROWpts + ImageJ;
	
ROIFimage = Fimage(rowsndx,colsndx);

figure
HS = pcolor(ROIXCOLpts,ROIYROWpts,ROIFimage);
shading flat
xlabel('XCOL');ylabel('YROW');
title('ROIFimage');
prettyplot(gca);
colormap jet;
colorbar;
set(gca,'clim',CLIM);

% Intensity histogram of ROI

%[Nim,Xim] = hist(ROIFimage(:),50);
figure;
hist(ROIFimage(:),50);
set(gca,'YScale','log');
ylim1 = [0.1 1e4];
set(gca,'Ylim',ylim1);
xlabel('photons');
ylabel('instances');
title('full time sum, ROI');

% Choose intensity limits for good pixels
BadLow = 180; % about 0.4 times average Fgmean 
% (Use 0.4 rather than 0.5 since top seems less to have lower illum)
BadHigh = 900; % about 2.0 times average Fgmean
hl = line(BadLow*[1 1],ylim1);
hl = line(BadHigh*[1 1],ylim1);

% Logical images of bad pixels

% Choose "Other" regions not to use
% e.g. boundaries of quadrants, and not illuminated corner
BadOimage = zeros(size(Fimage));
% "Fake" pixels with duplicate values at internal boundaries:
BadOimage(257:260,:) = 1;
BadOimage(:,257:260) = 1;
% not illuminated corner:
icorn = 70;
for ii = 1:icorn
    jj = 516-icorn+ii:516;
    BadOimage(ii,jj) = 1;
end

BadLimage = (Fimage < BadLow) & ~BadOimage;  
BadHimage = (Fimage > BadHigh) & ~BadOimage;
Badimage = BadLimage | BadHimage | BadOimage;

% Make version where bad pixels are NaN, good are 1

Badnan = Badimage + 1;
Badnan(Badnan>1) = NaN;

% Vector indices of bad pixels
[BadPixO] = find(BadOimage(:));
[BadPixL] = find(BadLimage(:));
[BadPixH] = find(BadHimage(:));

% Matrix subscripts of bad pixels
[BadPixYO,BadPixXO]=ind2sub(size(Fimage),BadPixO);
[BadPixYH,BadPixXH]=ind2sub(size(Fimage),BadPixH);
[BadPixYL,BadPixXL]=ind2sub(size(Fimage),BadPixL);

% X and Y plotting values of bad pixels
% If use ROI have to shift to be on original XY
%Yshift = min(ROI(3:4));
%Xshift = min(ROI(1:2));
Yshift = 0;
Xshift = 0;
BadPixYROWO = BadPixYO - ImageJ + Yshift;
BadPixXCOLO = BadPixXO - ImageJ + Xshift;
BadPixYROWH = BadPixYH - ImageJ + Yshift;
BadPixXCOLH = BadPixXH - ImageJ + Xshift;
BadPixYROWL = BadPixYL - ImageJ + Yshift;
BadPixXCOLL = BadPixXL - ImageJ + Xshift;

figure
HS = pcolor(XCOLpts,YROWpts,Fimage);
shading flat
xlabel('XCOL');ylabel('YROW');
title('Fimage w/ bad');
prettyplot(gca);
colormap jet;
colorbar;
% Note that in 'flat' shading, the pixel 0 is put between 0 and 1 (square) so to 
% get the circle on top of the pixel, need to shift by 1/2 
% If interpolated shading, would not need this shift in the cirles
HLo = line(BadPixXCOLO+.5, BadPixYROWO+.5);
HLh = line(BadPixXCOLH+.5, BadPixYROWH+.5);
HLl = line(BadPixXCOLL+.5, BadPixYROWL+.5);
set(HLo,'linestyle','none','marker','o','markeredgecolor','w');
set(HLh,'linestyle','none','marker','o','markeredgecolor','k');
set(HLl,'linestyle','none','marker','o','markersize',10,'markeredgecolor','w');
set(gca,'clim',CLIM);
set(gcf,'Paperposition',[2 2 5 4]);

figure
HS = pcolor(ROIXCOLpts,ROIYROWpts,ROIFimage);
shading flat
xlabel('XCOL');ylabel('YROW');
title('ROIFimage w/ bad');
prettyplot(gca);
colormap jet;
colorbar;
HLo = line(BadPixXCOLO+.5, BadPixYROWO+.5);
HLh = line(BadPixXCOLH+.5, BadPixYROWH+.5);
HLl = line(BadPixXCOLL+.5, BadPixYROWL+.5);
set(HLo,'linestyle','none','marker','o','markeredgecolor','w');
set(HLh,'linestyle','none','marker','o','markeredgecolor','k');
set(HLl,'linestyle','none','marker','o','markersize',10,'markeredgecolor','w');
set(gca,'clim',CLIM);
set(gcf,'Paperposition',[2 2 5 4]);


Fgood = Fimage.*~Badimage; % Bad pixels will be zero
Nbad = sum(Badimage(:));
Ngood = length(Fimage(:)) - Nbad;
Fgmean = sum(Fgood(:))/Ngood;

if 1
figure
HS = pcolor(XCOLpts,YROWpts,Fgood);
shading flat
xlabel('XCOL');ylabel('YROW');
title('Fimage, Not Bad Pixels');
prettyplot(gca);
colormap jet;
colorbar;
set(gca,'clim',CLIM);
end

% Normalize by average value of good pixels
imnorm = Fimage/Fgmean; % has 1101 zeros
imnormnan = imnorm.*Badnan; % has no zeros, all were bad pixels

FlatField.imnormnan = imnormnan;
FlatField.Badimage = Badimage;
FlatField.Docu = char(FlatField.Docu,...
[dateandtime,' : run through flatfield_analysis_MPX3.m']);

%save flatfield_MPX3.mat imnormnan Badimage
%save FlatField_2017_0809_2_num64 FlatField


