%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==8); % Calculate two-time correlations of pixels in ROIs over all frames

% ----------------- these pieces likely belong elsewhere - also look in datainput_AREA for ones to move
% Get X in terms of ML 
if ~isempty(lastframes);
% hardwire to use 'specpoint'
Xamount = tamount*(Xsteps - (lastframes(1) + delay))/(diff(lastframes) - delay);
else
Xamount = tamount*(Xsteps);
end

POINTSUMS=[1 Nt;
            40 70; % 1/2 ML
            71 100; % 1 ML
            lastframes(2)+delay Nt] - ImageJ;
POINTSUMS=[1 Nt;
            40 70; % 1/2 ML
            169 176; % 2.5 ML
            1 10; % 0 ML
            lastframes(2)+delay Nt] - ImageJ;
POINTSUMS=[1 Nt;
            1 10] - ImageJ; % 0 ML
% ------------------------------------------------------



disp('Calculating 2-time correlations, please wait.');
tic;
   
% First bin scans in time
if tbin > 1
    tbin = floor(tbin);
    Ntb = floor(Nt/tbin);
    IInormbe = reshape(IInormb(:,:,1:Ntb*tbin),Nr,Nc,tbin,Ntb);
    IInormbb = squeeze(sum(IInormbe,3));
    Xamounte = reshape(Xamount(1:Ntb*tbin),tbin,Ntb);
    Xamountb = squeeze(mean(Xamounte,1));
    POINTSUMSB(:,1) = floor((POINTSUMS(:,1)-1+ImageJ)/tbin)+1-ImageJ;
    POINTSUMSB(:,2) = floor((POINTSUMS(:,2)+ImageJ)/tbin) - ImageJ;
else
    Ntb = Nt;
    IInormbb = IInormb;
    Xamountb = Xamount;
    POINTSUMSB = POINTSUMS;
end

if isempty(POINTSUMSB); POINTSB=[1 Ntb]-ImageJ;
else POINTSB = POINTSUMSB;
end

% Correct for bad pixels and threshold variation, for input to smoothing
%Icorr = NaN*ones(size(IInormbb));
%for kk = 1:Ntb    
%    Icorr(:,:,kk) = IInormbb(:,:,kk)./imnormnanall; % normalize and set bad pixels = NaN
%end

Icorr = IInormbb; % already corrected

dIroiv = []; % vector to hold deviations from average, all good pixels in ROIs used

tic;
disp('Calculating mean for 2-time correlations, please wait.');

% Smoothing method for getting average

dj = -jjj:jjj;
di = -iii:iii;

% Calculate CC2 averages over integer and half-integer correlation regions
% Base on ML values in Xamountb

% For integer and half-integer during growth, use only values between 0 and 4.5 ML
Xn = Xamountb;
Xn(Xn<0) = NaN;
Xn(Xn>4.5) = NaN;
% Make into matrices, size of CC2
Xnc = ones(size(Xn'))*Xn;
Xnr = Xn'*ones(size(Xn));
% Set diagonal and lower triangle to NaN
an = 1 + tril(NaN*ones(size(Xnc)),0);
Xnc2 = an.*Xnc;
Xnr2 = an.*Xnr;
% Make sum and difference of rows and columns
Xns = Xnc2 + Xnr2;
Xnd = Xnc2 - Xnr2;
% Round to digitize
Xnsi = round(Xns);
Xndi = round(Xnd);

% Make matrices showing regions to check logic
CCIK = zeros(Ntb,Ntb);
CCHK = zeros(Ntb,Ntb);

%for ii = 1:5
%    for jj = ii:5
%        CCIK(Xnsi == ii+jj-2 & Xndi == jj-ii) = 5*(ii-1) + jj;
%        CCHK(Xnsi == ii+jj-1 & Xndi == jj-ii) = 5*(ii-1) + jj;
%    end
%end

for ii = 1:5
    for jj = ii:5
        CCIK(abs(Xns - (ii+jj-2)) < CWID & abs(Xnd - (jj-ii)) < CWID) = 5*(ii-1) + jj;
        CCHK(abs(Xns - (ii+jj-1)) < CWID & abs(Xnd - (jj-ii)) < CWID) = 5*(ii-1) + jj;
    end
end

% Plot logic checks
figure;
set(gcf,'Paperposition',[1 1 4 3]);
set(gcf,'Name','integer check');
HS = pcolor(Xamountb,Xamountb,CCIK);
shading flat;
%set(gca,'clim',[0 0.10]);
hx = xlabel('Growth Amount (ML)');
hy = ylabel('Growth Amount (ML)');
set(gca,'position',[.17 .15 .65 .8]);
axis equal
colorbar
colormap 'jet';

figure;
set(gcf,'Paperposition',[1 1 4 3]);
set(gcf,'Name','half integer check');
HS = pcolor(Xamountb,Xamountb,CCHK);
shading flat;
%set(gca,'clim',[0 0.10]);
hx = xlabel('Growth Amount (ML)');
hy = ylabel('Growth Amount (ML)');
set(gca,'position',[.17 .15 .65 .8]);
axis equal
colorbar
colormap 'jet';

% Loop over a sequence of ROIs for 2-time
% Use a different ROI list for 2-time


% Variables to save values for each step:
CC2M = NaN*ones(Ntb,Ntb,NRX);
CHM = NaN*ones(5,5,NRX);
CIM = NaN*ones(5,5,NRX);
CS0M = NaN*ones(1,NRX);
CSM = NaN*ones(5,NRX);
CE0M = NaN*ones(1,NRX);
CEM = NaN*ones(5,NRX);
CSEM = NaN*ones(1,NRX);
CCHM = NaN*ones(1,NRX);
CCIM = NaN*ones(1,NRX);

DXavg = NaN*ones(1,NRX);
    
for ix = 1:NRX
    TTROISI = TTROIS;
    TTROISI(:,1) = TTROISI(:,1) + (ix-1)*RXINC;
    TTROISI(:,2) = TTROISI(:,2) + (ix-1)*RXINC;
    
    DXavg(ix) = mean(TTROISI(1,1:2)) - XCEN;
    
    if abs(DXavg(ix)) < DXavgwin
        ttcroi = ttcroiin;
    else
        ttcroi = ttcroiout;
    end
    
    dIroiv = [];

% Use ROIs in ttcroi
    for iroi = ttcroi

% ROI indices
        ROIi = [TTROISI(iroi,3):TTROISI(iroi,4)]+ImageJ;
        ROIj = [TTROISI(iroi,1):TTROISI(iroi,2)]+ImageJ;

% expanded ROI subscripts to include window range around ROI
        ROIPi = [(TTROISI(iroi,3)-iii):(TTROISI(iroi,4)+iii)]+ImageJ;
        ROIPj = [(TTROISI(iroi,1)-jjj):(TTROISI(iroi,2)+jjj)]+ImageJ;

% first fix up NaN's, then use conv
% data in expanded ROI
        IcorrER = Icorr(ROIPi,ROIPj,:);

% bad pixels in expanded ROI
        ROIPnans = Badimage(ROIPi,ROIPj);

% Matrix subscripts of bad pixels in expanded ROI
        [nani,nanj] = find(ROIPnans);

% Correct to be subscripts in whole image
        naniw = nani + ROIPi(1) - 1;
        nanjw = nanj + ROIPj(1) - 1;

% Replace bad pixel values by local mean of good pixels
        IcorrER2 = IcorrER;
        for ii = 1:length(nani)   
            for kk = 1:Ntb    
                III = Icorr(naniw(ii)+di,nanjw(ii)+dj,kk);
                III(isnan(III))=[]; % remove bad pixels
                IcorrER2(nani(ii),nanj(ii),kk) = mean(III(:));  
            end
        end

% convolute - return only ROI:

%clear Iconv
%for kk = 1:Ntb
%    Iconv(:,:,kk) = conv2(IcorrER2(:,:,kk),ones(length(di),length(dj)),'valid')./(length(di)*length(dj));
%end

% Use 2-D Savitzky-Golay smoothing filter, cubic
        filt = sgsf_2d(di,dj,3,3,0);
        clear Iconv;
        for kk = 1:Ntb
            Iconv(:,:,kk) = conv2(IcorrER2(:,:,kk),filt,'valid');
        end

% Calculate difference from mean
        dIroi = Icorr(ROIi,ROIj,:) - Iconv;
% Loop for each time step
        BadimageROI = Badimage(ROIi,ROIj);
        clear dIroivj
        for jj = 1:Ntb
            dIroij = dIroi(:,:,jj); 
            dIroij(BadimageROI)=[]; % remove bad pixels
            dIroivj(:,jj) = dIroij(:);
        end

 % append to multi-roi list
        dIroiv = [dIroiv; dIroivj];
 
 
 % make some plots of smoothing results
 % ----------------------
        if plotsmooth
        for ii = 1:length(POINTSB(:,1))
            Ni = [POINTSB(ii,1): POINTSB(ii,2)]+ImageJ;
            NP = length(Ni);
	
            figure;clf;
            set(gcf,'Name','Iconv mean over scan pts');	
	
            HS = imagesc(sum(Iconv(:,:,Ni),3)./NP);
            %if ~isempty(AXISdet), axis(AXISdet);end
            DOCUInt = ['Iconv mean over SPEC scan pts [',int2str(POINTS(ii,:)),'])'];
            title(char(TITLEstr1,DOCUInt));
            %xlabel(XCOL);ylabel(YROW);
        	(gca);
            colorbar;
            colormap 'jet';
            set(gca,'ydir','Normal');
        end    

        % Plot intensity distributions in ROI, in POINTSUM

        if ~isempty(POINTSUMSB)   % sum over points sets in POINTSUMS
	
            figure; clf;
            set(gcf,'Name','SoX 1st ROI and Sum select points');
            for ii=1:length(POINTSUMSB(:,1));
                Ni = [POINTSUMSB(ii,1): POINTSUMSB(ii,2)]+ImageJ;
                HL(ii) = line(ROIi-ImageJ,sum(sum(IcorrER2(iii+1:end-iii,jjj+1:end-jjj,Ni),3),2)','color','b');
                line(ROIi-ImageJ,sum(sum(Iconv(:,:,Ni),3),2)','color','r');    
            end
            set(gca,'Yscale','lin')
            xlabel(YROW),ylabel('Int (sum)'),
            title(char(TITLEstr1,['summed over X in 1st ROI and over sel frames']));
            prettyplot(gca);
            legend(addnames2matrix('sum between [', int2str(POINTSUMS),']'))
            HLine.SoX_SoSelectPt = HL;clear HL;
            colormap 'jet';
            set(gca,'ydir','Normal');
    
            figure; clf;
            set(gcf,'Name','SoY 1st ROI and Sum selected points');
            for ii=1:length(POINTSUMSB(:,1));
                Ni = [POINTSUMSB(ii,1): POINTSUMSB(ii,2)]+ImageJ;
                HL(ii) = line(ROIj-ImageJ,sum(sum(IcorrER2(iii+1:end-iii,jjj+1:end-jjj,Ni),3),1)','color','b');
                line(ROIj-ImageJ,sum(sum(Iconv(:,:,Ni),3),1)','color','r');    
            end
            set(gca,'Yscale','lin')
            xlabel(XCOL),ylabel('Int (sum)'),
            title(char(TITLEstr1,['summed over Y in 1st ROI and between selected Spec Points']));
            prettyplot(gca);
            %legend(addnames2matrix('sum between [', int2str(POINTSUMS),']'));
            HLine.SoY_SoSelectPt = HL;clear HL;
            colormap 'jet';	
            set(gca,'ydir','Normal');
        end


        for ii = 1:length(POINTSB(:,1))
            Ni = [POINTSB(ii,1): POINTSB(ii,2)]+ImageJ;
            NP = length(Ni);
	
            figure;clf;
            set(gcf,'Name','dIroi mean over scan pts');	
	
            HS = imagesc(sum(dIroi(:,:,Ni),3)./NP);
            %if ~isempty(AXISdet), axis(AXISdet);end
            DOCUInt = ['dIroi mean over scan pts [',int2str(POINTS(ii,:)),'])'];
            title(char(TITLEstr1,DOCUInt));
            %xlabel(XCOL);ylabel(YROW);
        	(gca);
            colorbar;
            colormap 'jet';
            set(gca,'ydir','Normal');
	
        end  
        end


    toc;

    end
    


% Calc 2-time using this mean
% Mean is different for each pixel and time

    IIM2 = NaN*ones(Ntb,Ntb);
    for ii = 1:Ntb
        for jj = 1:ii
            IIM2(ii,jj) = mean(dIroiv(:,ii).*dIroiv(:,jj));
            IIM2(jj,ii) = IIM2(ii,jj);
        end
    end
    IID2 = diag(IIM2);
    CC2 = IIM2./sqrt(IID2*IID2'); % Normalized to make diagonal unity

    toc

%--------------------------------------------------------
    if 0 % Plot vs frame number
    figure;
    set(gcf,'Paperposition',[1 1 4 3]);
    set(gcf,'Name','2-time 2nd method');
    HS = pcolor(CC2);
    shading flat;
    %set(gca,'clim',[0 0.03]);	
    xlabel('Frame 1');
    ylabel('Frame 2');
    DOCUlast = ['2nd method 2-time correlation'];
    %title(char(TITLEstr1,[DOCUlast]));
    set(gca,'position',[.17 .11 .65 .8]);
    %prettyplot(gca);
    colorbar
    colormap 'jet';
    end
%--------------------------------------------------------
% Plot vs growth amount
    if plotall
    figure;
    set(gcf,'Paperposition',[1 1 4 3]);
    set(gcf,'Name','2-time 2nd method');
    HS = pcolor(Xamountb,Xamountb,CC2);
    shading flat;
    set(gca,'clim',[0 0.30]);
    hx = xlabel('Growth Amount (ML)');
    hy = ylabel('Growth Amount (ML)');
    DOCUlast = ['2-time correlation 2nd method'];
    title(['DXavg = ' num2str(DXavg(ix))]);
    %set(gca,'position',[.17 .12 .65 .68]);
    set(gca,'position',[.17 .15 .65 .8]);
    axis equal
    %set(gca,'FontSize',12);
    %set(hx,'FontSize',12);
    %set(hy,'FontSize',12);
    colorbar
    colormap 'jet';
    % Put some lines showing start and end of growth
    LW = 1.5;
    LC = 'w';
    pa = axis;
    hl = line(0*[1,1],pa(3:4));
    set (hl,'color',LC,'LineWidth',LW);
    hl = line(4.5*[1,1],pa(3:4));
    set (hl,'color',LC,'LineWidth',LW);
    hl = line(pa(1:2),0*[1,1]);
    set (hl,'color',LC,'LineWidth',LW);
    hl = line(pa(1:2),4.5*[1,1]);
    set (hl,'color',LC,'LineWidth',LW);
    set(gca,'Xtick',[0:4]);
    set(gca,'Ytick',[0:4]);
    axis([min(Xamountb) max(Xamountb) min(Xamountb) max(Xamountb)]);
    set(gca,'TickDir','out');
    if 0
    % CTR (0-1ML)
    ROIS_show  = [...                
    4.5  5.04  min(Xamountb)  max(Xamountb)];    % #49
    [HROI,COLORORDER] = showrois(ROIS_show,gcf);
    end
    end
%--------------------------------------------------------
% Calculate some metrics of the CC2 matrix

% Horizontal and vertical means of upper off-diagonal

    CX2 = sum(triu(CC2,1),1)./[0:Ntb-1];

    CY2 = sum(triu(CC2,1),2)./[Ntb-1:-1:0]';

% Radial and transverse averages

    for ii = 1:Ntb
        CR2(ii) = mean(diag(CC2,ii-1));
        CT2(ii) = mean(diag(flipud(CC2-diag(ones(Ntb,1))),ii-1))/2;
    end
%--------------------------------------------------------
    if 0
    P=[97:112];  %#49   <<< %%%% ATTENTION ATTENTION SEEMS TO BE PLACE WHERE CHANGE MANUALLY FOR DIFFERENT SCANS CT
    [row,column]=size(CC2);
    CCs = sum(CC2(:,P),2);
    for ii = P
        CCs(ii) = CCs(ii) - 1;
    end

    figure;
    hax=axes('Box','on');
    HL=plot(Xamountb, CCs);
    set(HL,'Marker','o','LineStyle','-');
    xlabel('Growth Amount (ML)'); ylabel('Intensity');
    title(char(TITLEstr1,['Summed over 2-time correlation']));
    xlim([min(Xamountb) max(Xamountb)]);
    set(gca,'position',[.17 .10 .65 .68]);
    prettyplot(gca)
    end

% --------------------------------------------------------


% Get integer and half-integer regions to average
    CI = NaN*ones(5,5);
    CH = NaN*ones(5,5);
    %for ii = 1:5
    %    for jj = ii:5
    %        CI(ii,jj) = mean(CC2(Xnsi == ii+jj-2 & Xndi == jj-ii));
    %        CH(ii,jj) = mean(CC2(Xnsi == ii+jj-1 & Xndi == jj-ii));
    %   end
    %end

    for ii = 1:5
        for jj = ii:5
            CI(ii,jj) = mean(CC2(abs(Xns - (ii+jj-2)) < CWID & abs(Xnd - (jj-ii)) < CWID));
            CH(ii,jj) = mean(CC2(abs(Xns - (ii+jj-1)) < CWID & abs(Xnd - (jj-ii)) < CWID));
        end
    end

% Also calculate correlations with regions before growth start and after end

    Xsi = Xamountb<0;
    % Upper triangle pre-growth
    CS0 = sum(sum(triu(CC2(Xsi,Xsi),1)))./(0.5*sum(Xsi)*(sum(Xsi)-1));
    % Integer correlations with pre-growth
    CS = NaN*ones(1,5);
    for ii = 1:5
        Xii = Xamountb>0 & round(Xamountb) == ii - 1;
        CS(ii) = mean(mean(CC2(Xsi,Xii)));
    end

    Xei = Xamountb>4.5;
    % Upper triangle post-growth
    CE0 = sum(sum(triu(CC2(Xei,Xei),1)))./(0.5*sum(Xei)*(sum(Xei)-1));
    % Half-integer correlations with post-growth
    CE = NaN*ones(1,5);
    for ii = 1:5
        Xii = Xamountb<4.5 & round(4.5 - Xamountb) == ii - 1;
        CE(ii) = mean(mean(CC2(Xii,Xei)));
    end

    % correlation of pre- and post-growth
    CSE = mean(mean(CC2(Xsi,Xei)));

% create a single mean value for half-integer, leave out first ML and diag
    CCH = 0;
    for ii = 2:4
        for jj = ii+1:5
            CCH = CCH + CH(ii,jj);
        end
    end
    CCH = CCH/6;

% create a single mean value for integer, leave out diag
    CCI = 0;
    for ii = 1:4
        for jj = ii+1:5
            CCI = CCI + CI(ii,jj);
        end
    end
    CCI = CCI/10;

    CC2M(:,:,ix) = CC2;
    CHM(:,:,ix) = CH;
    CIM(:,:,ix) = CI;
    CS0M(ix) = CS0;
    CSM(:,ix) = CS;
    CE0M(ix) = CE0;
    CEM(:,ix) = CE;
    CSEM(ix) = CSE;
    CCHM(ix) = CCH;
    CCIM(ix) = CCI;

end

% Summary 2-time plots over main Q regions

DXHmin = 32;
DXHmax = 56;

idxh = abs(DXavg) > DXHmin & abs(DXavg) < DXHmax;

CC2H = mean(CC2M(:,:,idxh),3);

DXImin = 3;
DXImax = 21;

idxi = abs(DXavg) > DXImin & abs(DXavg) < DXImax; 

CC2I = mean(CC2M(:,:,idxi),3);

figure;
set(gcf,'Paperposition',[1 1 4 3]);
set(gcf,'Name',['DX ' num2str(DXHmin) ' to ' num2str(DXHmax)]);
HS = pcolor(Xamountb,Xamountb,CC2H);
shading flat;
set(gca,'clim',[0 0.30]);
hx = xlabel('Growth Amount (ML)');
hy = ylabel('Growth Amount (ML)');
DOCUlast = ['2-time correlation'];
title(['DX ' num2str(DXHmin) ' to ' num2str(DXHmax)]);
%set(gca,'position',[.17 .12 .65 .68]);
set(gca,'position',[.17 .15 .65 .8]);
axis equal
%set(gca,'FontSize',12);
%set(hx,'FontSize',12);
%set(hy,'FontSize',12);
colorbar
colormap 'jet';
% Put some lines showing start and end of growth
LW = 1.5;
LC = 'w';
pa = axis;
hl = line(0*[1,1],pa(3:4));
set (hl,'color',LC,'LineWidth',LW);
hl = line(4.5*[1,1],pa(3:4));
set (hl,'color',LC,'LineWidth',LW);
hl = line(pa(1:2),0*[1,1]);
set (hl,'color',LC,'LineWidth',LW);
hl = line(pa(1:2),4.5*[1,1]);
set (hl,'color',LC,'LineWidth',LW);
set(gca,'Xtick',[0:4]);
set(gca,'Ytick',[0:4]);
axis([min(Xamountb) max(Xamountb) min(Xamountb) max(Xamountb)]);
set(gca,'TickDir','out');
    
figure;
set(gcf,'Paperposition',[1 1 4 3]);
set(gcf,'Name',['DX ' num2str(DXImin) ' to ' num2str(DXImax)]);
HS = pcolor(Xamountb,Xamountb,CC2I);
shading flat;
set(gca,'clim',[0 0.30]);
hx = xlabel('Growth Amount (ML)');
hy = ylabel('Growth Amount (ML)');
DOCUlast = ['2-time correlation'];
title(['DX ' num2str(DXImin) ' to ' num2str(DXImax)]);
%set(gca,'position',[.17 .12 .65 .68]);
set(gca,'position',[.17 .15 .65 .8]);
axis equal
%set(gca,'FontSize',12);
%set(hx,'FontSize',12);
%set(hy,'FontSize',12);
colorbar
colormap 'jet';
% Put some lines showing start and end of growth
LW = 1.5;
LC = 'w';
pa = axis;
hl = line(0*[1,1],pa(3:4));
set (hl,'color',LC,'LineWidth',LW);
hl = line(4.5*[1,1],pa(3:4));
set (hl,'color',LC,'LineWidth',LW);
hl = line(pa(1:2),0*[1,1]);
set (hl,'color',LC,'LineWidth',LW);
hl = line(pa(1:2),4.5*[1,1]);
set (hl,'color',LC,'LineWidth',LW);
set(gca,'Xtick',[0:4]);
set(gca,'Ytick',[0:4]);
axis([min(Xamountb) max(Xamountb) min(Xamountb) max(Xamountb)]);
set(gca,'TickDir','out');
  

figure;
hax = axes('Box','on');
HL = line(DXavg, CCHM);
set(HL,'Marker','o','LineStyle','-');
xlabel('DX (pixels)'); 
ylabel('Half Integer Correlation');
title(char(TITLEstr1,['Half integer correlation']));
set(gca,'position',[.17 .10 .65 .68]);
set(gca,'Xlim',[-100 100]);
prettyplot(gca)

figure;
hax = axes('Box','on');
HL = line(DXavg, CCIM);
set(HL,'Marker','o','LineStyle','-');
xlabel('DX (pixels)'); 
ylabel('Integer Correlation');
title(char(TITLEstr1,['Integer correlation']));
set(gca,'position',[.17 .10 .65 .68]);
set(gca,'Xlim',[-100 100]);
prettyplot(gca)

% Fit CIM, CHM at each Q

% Fitting routine variables 
global verbose;
verbose = [0 0];
stol = 0.000001;
niter = 100;
%       const slope decay size
pin = [1   0   2   5]'; % Zero fixes parameter
dp =   [1   1   1   0]'.*0.0001; % Fraction to vary for deriv

for ix = 1:NRX
    
    CH = CHM(:,:,ix);
    xfit = find(~isnan(CH));
    if isempty(xfit)
        pouthm(:,ix) = NaN;
        covphm(:,ix) = NaN;
        xisqh(ix) = NaN;
    else
        
        yfit = CH(xfit);

    % Choose weights
        wfit = ones(size(xfit));
    
        [ffit,pout,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(xfit,yfit,pin,'cc_fn',stol,niter,wfit,dp);
        xisqh(ix) = sum((wfit.*(ffit-yfit)).^2);
        pouthm(:,ix) = pout;
        sighm(:,ix) = sqrt(diag(covp));
    end
    
    CI = CIM(:,:,ix);
    xfit = find(~isnan(CI));
    if isempty(xfit)
        poutim(:,ix) = NaN;
        covpim(:,ix) = NaN;
        xisqi(ix) = NaN;
    else
        
        yfit = CI(xfit);

    % Choose weights
        wfit = ones(size(xfit));
    
        [ffit,pout,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(xfit,yfit,pin,'cc_fn',stol,niter,wfit,dp);
        xisqi(ix) = sum((wfit.*(ffit-yfit)).^2);
        poutim(:,ix) = pout;
        sigim(:,ix) = sqrt(diag(covp));
    end
end

%----------------------------------------------------------

 % -------------  end of checked  (goes with if WhattoDo==8
end
