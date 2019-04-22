% Script ptlfit.m
% Reads PbTiO3 normalized data
% Fits to model of CT and BS
% 29-aug-96 gbs

%close;
%close;
%clear all;
%clear global;

%scrnam = input('Enter data set name (default=ptl100)> ','s');
%if isempty(scrnam)
%     scrnam='ptl100';
%end

% control tweak function saves the
% scrname, dp, and params in controldata when the FIT button is
% hit. The params and dp will replace the vary and guess from the
% guess startup file.
% scrname is the name of the files to look at.

load ptlcontrolinfo

% Get initial guess fitting values and descriptor

eval([scrnam 'guess']);

% Get normalized data

eval(['load ' scrnam '.data']);

vary=dp;	%dp is output of control panel for things to vary
guess=params;	%params is inital guess from control panel

if exist('fighanf') == 1
   figure(fighanf);
   set(gcf,'Position',[300 100 630 500]);
else
   fighanf = figure('Position',[300 100 630 500]);
end
clf reset;
set(gcf,'PaperPosition',[0.25 1.5 8 8]);

% Set up data to fit: 

ptlxa = eval([scrnam '(:,1)']);
ptlya = eval([scrnam '(:,sigcol)']);
ptleba = eval([scrnam '(:,3)']);
ptlwa = 1 ./ ptleba;

ptlx = ptlxa;
ptly = ptlya;
ptlw = ptlwa;

if exist('exclmin') == 1
  nofitindx = ptlxa >= exclmin & ptlxa <= exclmax;
  ptlx(nofitindx) = [];
  ptly(nofitindx) = [];
  ptlw(nofitindx) = [];
end

Energy = 8051;
ptlinit(ptlx,Energy);

global verbose
verbose(1) = 0;
verbose(2) = 1;
stol = .00005;
niter = 30;

ptlg = guess;
ptlc = 0.001*vary;
%ptlc = [0.001 0.001 0.001 0.001 0.0001 0.001 0.001].*vary;
ptlo = [0 0.01; ...
        0 0.05; ...
	0 0.05; ...
        0 0.50; ...
        0 0.50; ...
        0 0.01; ...
        0 0.50; ...
        0 10.0];


[ptlf,ptlp,kvg,iter,corp,covp,covr,stdresid,Z,r2]= ...
leasqrl(ptlx,ptly,ptlg,'ptlfunc',stol,niter,ptlw,ptlc,'dfdp');
%leasqrl(ptlx,ptly,ptlg,'ptlfunc',stol,niter,ptlw,ptlc,'dfdp',ptlo);
ptls = zeros(size(ptlp));
ptls(ptlc~=0) = sqrt(diag(covp)); 
ptlchi2 = mean((ptlw.*(ptlf-ptly)).^2);

ptlinit(ptlxa,Energy);
ptlfa = ptlfunc(ptlxa,ptlp);
semilogy(ptlxa,ptlya,'go');
hold on;
%errorbar(ptlxa,ptlya,ptleba,ptleba,'go');
semilogy(ptlxa,ptlfa,'w-');
axis(plotaxis);
hold off;
  
xlabel('L (rlu)');
ylabel('Normalized Intensity (cts/mon)');
title([date,' : ',rdesc]);

disp(['fit parameters: ' sprintf('%10.3g',ptlp(:))]);
disp(['r2 = :          ' sprintf('%10.3g',r2)]);
disp(['red. chi2 = :   ' sprintf('%10.3g',ptlchi2)]);
disp(['kvg= :(1=good)  ' sprintf('%10.3g',kvg)]);
disp(['iter=:          ' sprintf('%10.3g',iter)]);
pax = axis;
%avary = setstr('V'.*vary + 'C'.*(~vary));
txtx = pax(1) + 0.05*(pax(2)-pax(1));
txty = pax(3) .* (pax(4)./pax(3)).^0.90;
dtxty = (pax(4)./pax(3)).^0.05;
for ii = [1:length(ptlp)]
  text(txtx,txty,[paramnames(ii,:) ...
              sprintf('%10.3g +-%10.3g',ptlp(ii),ptls(ii))]);
  txty = txty ./ dtxty;
end

txtx = pax(1) + 0.65*(pax(2)-pax(1));
txty = pax(3) .* (pax(4)./pax(3)).^0.90;
text(txtx,txty,['Reduced Chi^2: ' sprintf('%10.3g',ptlchi2)]);
%txty = pax(3) .* (pax(4)./pax(3)).^0.85;
%text(txtx,txty,['Corr. Coef. r^2: ' sprintf('%10.3g',r2)]);
if exist('exclmin') == 1
  txty = pax(3) .* (pax(4)./pax(3)).^0.80;
  text(txtx,txty,[sprintf('Excluded: %5.3f to %5.3f',exclmin,exclmax)]);
end
% eval('print -dps');
%eval(['print -dps transps/' rnam '.g0.ps']);
% eval(['! lpr -h -s -Pps transps/' rnam '.g0.ps']);



