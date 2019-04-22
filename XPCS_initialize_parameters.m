classdef XPCS_initialize_parameters
    
    properties(Constant)
    end
    
    methods(Static)
        
        function [plotsmooth,plotall,plotorig,plotnew,ixplotmin,ixplotmax,...
                pcolorCC2,plothalf,plotfull,plotcorrorig,plotcorrnew,...
                plotcontrast,pcolorextra,posflag,fitorig,fitnew, altnorm,plotdeltat]...
                = TTsput_flags()
            % Flags to control behavior of the master script
            
            
            % Plot flags
            plotsmooth = 0; % use to control whether to plot smoothing results if skip = 0
            
            % Note above will generate several plots per NRX, keep NRX small
            plotall = 1; % use to control whether to plot every 2-time in Q loop
            plotorig = 1; % plot original normalization
            plotnew = 1; % plot new normalization
            
            % Pixel range: Limit which ix to plot within 1:NRX
            ixplotmin = 25;
            ixplotmax = 50;
            ixplotmin = 1;
            ixplotmax = 5;
            
            % Color: max of clim for these plots
            %CMAX = 0.10;
            pcolorCC2 = 0; % use to control whether to plot 2-time mean over ttroi regions
            plothalf = 1;
            plotfull = 0;
            plotcorrorig = 0;
            plotcorrnew = 0;
            plotcontrast = 0;
            pcolorextra = 0; %use to control whether to pcolor Per pixel counts in ROIs, contrast, excess contrast
            
            % Fitting flags:
            posflag = 1; % average 2-time only on positive side
            fitorig = 0;
            fitnew = 0;
            
            altnorm = 0; % control calc/plot of alternative normalizations
            
            plotdeltat = 1; % control calc/plot of delta t analysis
            
        end
        
        function [ImageJ,Xstepcol,SINGLE,BKG,...
                scanflag,imname,p_image,ending,POINTSUMS] = TTsput_read_ini()
            
            
            ImageJ = 1; % use imageJ indexing [start at 0] when denoting ROIs and when plotting pixels
            %	This applies to AXISdet and ROI values so must be consistent.
            % Also in time, first image is 0

            Xstepcol = 0; % use pt number
            
            %%% SPEC POINT CONVENTION HERE  %%%%%%%%%%%
            SINGLE 		= 5;	% required for [0,1,2]
            
            BKG			= [];	% =[] no background subtracted,          
           
            % to feed the SCR function
            scanflag = 1;
            imname = [];%'pil_';
            p_image = 'p_image';
            ending = '.tif';
            
             POINTSUMS = [];%[145 2592]-ImageJ;
            
        end
        
        function [DOCU0,DOCU1,LOGFLAG,delay,XFRAC] = TTparameters_allruns()
            % Set parameters that stay same for each run
            
            DOCU0			= 'sput_180315a 110 TiO2';
            DOCU1 = [''];
            LOGFLAG = 1;   % 1 for log, 0 for lin need to change CLIM to match if using it
                        
            delay = 3; % delay in frames of growth start after valve switch
            XFRAC = 0.1; % use a range that is a fraction of xmax
        end
        
        function [TCV,nT,iiT,specfilenameM,SCNstrM,XCENV,YCENV,XWIDV,YWIDV,...
                ymax,tminv,tmaxv,clrs, DXImin, ...
                DXImax,NRY,PWV] = TTparameters_singlerun(flag_equil_or_growth)
            % Set parameters that varied for each scan
            
            
            switch flag_equil_or_growth
                
                case 'equilibrium'
                    % temperatures
                    TCV = [700;750;800;850;900]; % Equilibrium Thermocouple (new heater)
                    
                      % power in Watts
                    PWV = [0;0;0;0;0];
                    
                    % Filenames
                    specfilenameM = ['2018_0316_1';'2018_0317_1';'2018_0317_1';'2018_0318_1';'2018_0318_1'];
                    
                    % scans number (see lab book)
                    SCNstrM = ['56';'26'; '55' ; '37';'59']; % Equilibrium
                    
                    % Center pixels in between the CTRS
                    XCENV = [240;240; 238;232;230];% del
                    YCENV = [119;119;120;120;120]; % nu
                    
                    % Width of the larger ROI
                    XWIDV = [30;30;30;30;30];
                    % Pilatus detector on sevchex arm (X is del and Y is nu)
                    YWIDV = [30;30;30;30;30];
                    
                    % Positions of CTRs for ROIS
                    ymax = [8;8;8;8;8];
                    
                    % Min and max on time range for delta-time average
                    tminv = [200;1500;1;2160;720];
                    tmaxv = [3700;4800;4000;3700;1400];
                    
                    
                case 'growth'
                    TCV = [700;800;800;850;900;900;900]; % Growth Thermocouple (new heater)
                    
                     % power in Watts
                    PWV = [10; 10 ; 10 ; 10 ; 10 ; 20 ; 50];
                    
                    % Growth:
                    specfilenameM = ['2018_0316_1';'2018_0317_1';'2018_0317_1';'2018_0318_1';'2018_0318_1';'2018_0318_1';'2018_0318_1'];
                    
                    
                    SCNstrM = ['098';'045';'061';'041';'073';'102';'103']; % Growth
                    
                    % Center pixels in between the CTRS
                    %%{
                    XCENV = [240;238;238;253;236;240;240];% del
                    YCENV = [118;119;120;120;118;120;120]; % nu
                    %}
                    XWIDV = [30;30;30;30;30;30;30];
                    % Pilatus detector on sevchex arm (X is del and Y is nu)
                    YWIDV = [30;30;30;30;30;30;30];
                    
                    % Offset of two of the ROIS around the Positions of CTRs
                    ymax = [8;8;8;8;8;8;8];
                    
                    % Min and max on time range for delta-time average
                    tminv = [1000 ;150 ;200 ;1000; 700; 200;1050];
                    tmaxv = [5000;3600;3400;3800;3500;3600;3000];
                    
                case 'power_series'
                    
                    % temperatures in C
                    TCV = [800;800;800;800;800]; % Equilibrium Thermocouple (new heater)
                    
                    % power in Watts
                    PWV = [0; 10 ; 15 ; 20 ; 25];
                    
                    % Filenames
                    specfilenameM = ['2018_0813_1';'2018_0813_1';'2018_0813_1';'2018_0813_1';'2018_0813_1'];
                    
                    % scans number (see lab book)
                    SCNstrM = ['021';'018';'016';'032';'030']; % Equilibrium
                    
                    % Center pixels in between the CTRS
                    %XCENV = [230;236;240;240];% del
                    %YCENV = [120;118;120;120]; % nu
                    
                    XCENV = [93;93;93;93;93];% del
                    YCENV = [138;138;138;140;138]; % nu
                    
                    
                    % Width of the larger ROI
                    %XWIDV = [30;30;30;30];
                    % Pilatus detector on sevchex arm (X is del and Y is nu)
                    %YWIDV = [30;30;30;30];
                    
                    % Width of the larger ROI
                    XWIDV = [25;25;25;25;25];
                    % Pilatus detector on sevchex arm (X is del and Y is nu)
                    YWIDV = [45;45;45;45;45];
                    
                    
                    % Offset of two of the ROIS around the Positions of CTRs
                    ymax = [8;8;8;8;8];
                    
                    % Min and max on time range for delta-time average
                    tminv = [100; 100; 100;100;1000];
                    tmaxv = [4000;4000;4000;4000;3000];
                    %tminv = [720; 2500;1110;200];
                    %tmaxv = [1400;3500;2000;600];
                    
                    
                case 'only_data'
                    TCV = [800]; % Equilibrium Thermocouple (new heater)
                    
                    % power in Watts
                    PWV = [25];
                    
                    % Filenames
                    specfilenameM = ['2018_0813_1'];
                    
                    % scans number (see lab book)
                    SCNstrM = ['030'];
                    
                    % Center pixels in between the CTRS
                    %XCENV = [230;236;240;240];% del
                    %YCENV = [120;118;120;120]; % nu
                    
                    XCENV = [93];% del
                    YCENV = [135]; % nu
                    
                    
                    
                    % Width of the larger ROI
                    XWIDV = [25];
                    % Pilatus detector on sevchex arm (X is del and Y is nu)
                    YWIDV = [45];
                    
                    
                    % Offset of two of the ROIS around the Positions of CTRs
                    ymax = [8];
                    
                    % Min and max on time range for delta-time average
                    tminv = [1700];
                    tmaxv = [3200];
                    %tminv = [720; 2500;1110;200];
                    %tmaxv = [1400;3500;2000;600];


                    
                    
            end
            
          
            % index for the temperatures
            nT = length(TCV);% Do only first            
            
            % numbers of scans per temperature or number of total scans
            iiT = 1:nT; % to plot all
            %iiT = 1; % to plot only one
           
            %clrs = 'rgbcmkrgbcmkrgbcmkrgbcmk';
            clrs = 'orogobocomokxrxgxbxcxmxk*r*g*b*c*m*k.r.g.b.c.m.k';
            
            DXImin = 5;
            DXImax = 15;            
            NRY = 5; % Number of small TTROIS            
           
        end
        
        function [fitrange_time_iiT,pin_iiT,dp_iiT,qfitrange]= TTparameters_fit_singlerun(flag_equil_or_growth)
            % This function sets the fitrang, the initial guess and the parameters to 
            % maintain fix during the fit of the CC2NS functions in
            % SPCS_analysis.fit_CCN2S_with_leasqr
            
            switch flag_equil_or_growth
                
                case 'equilibrium'
                    
                    % set the fitrange of the CC2NS in time for each temperature:
                    fitrange_time_iiT = [2/3;2/5;2/5;2/5;2/5];
                    
                    pin_iiT = [0 1e-2 50 0;
                        0 1e-2 50 0;
                        0 1e-2 50 0;
                        0 1e-2 50 0;
                        0 1e-2 50 0];
                    
                    
                    dp_iiT =  [[1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001];
                    
                    % fit range for the tau time constant versus q
                    qfitrange.nu = [1.5e-3 1.5e-2;
                        1.5e-3 1.5e-2;
                        1.5e-3 1.5e-2;
                        1.5e-3 1.5e-2;
                        1.5e-3 1.5e-2];
                    
                case 'growth'
                    
                    % set the fitrange of the CC2NS in time for each temperature:
                    fitrange_time_iiT = [2/3;2/5;2/5;2/5;2/5;2/5;2/5];
                    
                    pin_iiT = [0 1e-2 50 0;
                        0 1e-2 50 0;
                        0 1e-2 50 0;
                        0 1e-2 50 0;
                        0 1e-2 50 0;
                        0 1e-2 50 0;
                        0 1e-2 50 0];
                    
                    
                    dp_iiT =  [[1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001];
                    
                    qfitrange.nu = [1.5e-3 1.5e-2;
                        1.5e-3 1.5e-2;
                        1.5e-3 1.5e-2;
                        1.5e-3 1.5e-2;
                        1.5e-3 1.5e-2;
                        1.5e-3 1.5e-2;
                        1.5e-3 1.5e-2];
                    
                case 'power_series'
                    
                    % set the fitrange of the CC2NS in time for each temperature:
                    fitrange_time_iiT = [1/2;1/2;1/2;1/2;1/2];
                    
                    
                    pin_iiT = [0 1e-3 500 0;
                        0 1e-3 500 0;
                        0 1e-3 500 0;
                        0 1e-3 500 0;
                        0 1e-3 500 0];
                    
                    dp_iiT =  [[1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001;
                        [1 1 1 0]*0.0001];
                    
                    qfitrange.nu = [2.5e-3 1.5e-2;
                        2.5e-3 1.5e-2;
                        2.5e-3 1.5e-2;
                        2.5e-3 1.5e-2;
                        2.5e-3 1.5e-2];
                    
                case 'only_data'
                    
                    % set the fitrange of the CC2NS in time for each temperature:
                    fitrange_time_iiT = [1/2.5];
                    
                    pin_iiT = [0 5e-3 50 0];
                    
                    
                    dp_iiT =  [[0 1 1 0]*0.0001];
                    
                    % fit range for the tau time constant versus q
                    qfitrange.nu = [1.5e-3 1.5e-2];
                    
            end
            
        end
        

        
        function [hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,tbin_allT...
                ,CWID_allT,N_degree_allT,mask] = TTparameters_2timecorr_calc(flag_equil_or_growth)
            % Set of parameters to calculate the area where the 2 times correlation
            % function is calculated:
            
            
             switch flag_equil_or_growth
                
                case 'equilibrium'
                    
                     
                     hwttr_allT = [1 1 1 1 1]; % row half width (pixels)=> box of 2*hwttr+1 pixels
                     hwttc_allT = [16 16 16 16 16]; % col half width (pixels)
                     
                     wrq_allT = [8 8 8 8 8];
                     wcq_allT = [0 0 0 0 0];
                     
                     % tuning the center of the reciprocal space for the CC2avg
                     % analysis (it can be larger than Ncs and Nrs respectively
                     offsetcc_allT = [0 0 0 0 0];
                     offsetrc_allT = [0 0 0 0 0];
                     
                     tbin_allT = [10 10 10 10 10]; % number of scans to bin together for 2-time calcs
                     % Note sumrois_wNaN doesn't work with single-pixel ROIs
                     
                     CWID_allT = [0.3 0.3 0.3 0.3 0.3]; % Parameter for integer/half-integer ML integration
                     
                     % degree of the polynomial when we fit the IInormbb reference
                     % in the new way of analysis
                     N_degree_allT = [1 1 1 1 1];
                    
                 case 'growth'
                     
                     hwttr_allT = [1 1 1 1 1 1 1]; % row half width (pixels)=> box of 2*hwttr+1 pixels
                     hwttc_allT = [16 16 16 16 16 16 16]; % col half width (pixels)
                     
                     wrq_allT = [8 8 8 8 8 8 8];
                     wcq_allT = [0 0 0 0 0 0 0];
                     
                     % tuning the center of the reciprocal space for the CC2avg
                     % analysis (it can be larger than Ncs and Nrs respectively
                     offsetcc_allT = [0 0 0 0 0 0 0];
                     offsetrc_allT = [0 0 0 0 0 0 0];
                     
                     tbin_allT = [10 10 10 10 10 10 10]; % number of scans to bin together for 2-time calcs
                     % Note sumrois_wNaN doesn't work with single-pixel ROIs
                     
                     CWID_allT = [0.3 0.3 0.3 0.3 0.3 0.3 0.3]; % Parameter for integer/half-integer ML integration
                     
                     % degree of the polynomial when we fit the IInormbb reference
                     % in the new way of analysis
                     N_degree_allT = [1 1 1 1 1 1 1];
                     
                 case 'power_series'
                     
                     hwttr_allT = [0 0 0 0 0]; % row half width (pixels)=> box of 2*hwttr+1 pixels
                     %hwttc_allT = [2 2 2 2 2]; % col half width (pixels)
                     hwttc_allT = [2 2 2 4 2]; % col half width (pixels)
                     
                     % number of boxes in the rows and columns directions
                     wrq_allT = [10 10 10 10 10];
                     wcq_allT = [4 4 4 2 4];
                     
                     % tuning the center of the reciprocal space for the CC2avg
                     % analysis (it can be larger than Ncs and Nrs respectively
                     offsetcc_allT = [0 0 0 0 0];
                     offsetrc_allT = [0 0 0 0 0];
                     
                     % number of scans to bin together for 2-time corr calcs
                    tbin_allT = [10 10 10 10 10];
                     % Note sumrois_wNaN doesn't work with single-pixel ROIs
                     
                     CWID_allT = [0.3 0.3 0.3 0.3 0.3]; % Parameter for integer/half-integer ML integration
                     
                     % degree of the polynomial when we fit the IInormbb reference
                    
                     N_degree_allT = [1 1 1 1 1];
                     
                 case 'only_data'
                     
                     
                     hwttr_allT = [0]; % row half width (pixels)=> box of 2*hwttr+1 pixels
                     hwttc_allT = [24]; % col half width (pixels)
                     
                     % number of boxes:
                     wrq_allT = [6];
                     wcq_allT = [0];
                     
                     % tuning the center of the reciprocal space for the CC2avg
                     % analysis (it can be larger than Ncs and Nrs respectively
                     offsetrc_allT = [-2];%[10];
                     offsetcc_allT = [0];
                     
                     % indexes for masking certain regions and excluding
                     % them from the two-time correlation function average:
                     mask.hwttr_allT = hwttr_allT; % row half width (pixels)=> box of 2*hwttr+1 pixels
                     mask.hwttc_allT = [5]; % col half width (pixels)
                     
                     % number of boxes:
                     mask.wrq_allT = wrq_allT;%[10];
                     mask.wcq_allT = wcq_allT;%[0];
                     
                     % tuning the center of the reciprocal space for the CC2avg
                     % analysis (it can be larger than Ncs and Nrs respectively
                     mask.offsetrc_allT = offsetrc_allT;%[-3];
                     mask.offsetcc_allT = offsetcc_allT;%[0];
                     
                     
                     
                     tbin_allT = [10]; % number of scans to bin together for 2-time calcs
                     % Note sumrois_wNaN doesn't work with single-pixel ROIs
                     
                     CWID_allT = [0.3]; % Parameter for integer/half-integer ML integration
                     
                     % degree of the polynomial when we fit the IInormbb reference
                     % in the new way of analysis
                     N_degree_allT = [1];
                     
                     
             end
            
           
            
        end
        
        function [POSITION,PAPERPOSITION,FONTSIZE,CMAX,CLIM,XCOLlabel,YROWlabel,...
                AXISdet,DOCUclim,INFOstr] = TTplot_parameters()
            % This script contains the parameters to plot the graphs.
            
            POSITION = [.17 .10 .65 .8];
            PAPERPOSITION = [1 1 5 4];
            FONTSIZE = 15;
            CMAX = 0.02;
            CLIM = [-1.5 2];
            
            XCOLlabel = 'pix(X)';
            YROWlabel = 'pix(Y)';
            AXISdet = [];
            DOCUclim=[];
            %INFOstr = char(INFOinterest,DUM);   % currently added onto 0 and 5
            INFOstr = [];    % turn this off quickly if desired
           

        end
        
        
        function [D_ds,kvector,lambda,pixel_size,th_Bragg,SplitQ_CTR,t_osc] = TT_experiment()
            
            % This function initializes the relevant parameters of the experiment:
            
            % Sample to front detector distance (p 117, book 208)
            %{ 
            %August 2018 
            D_ds = 2.54; % in [m]
            E_keV = 18; % in [keV]
            kvector = (2*pi/12.398)*E_keV*1e10; % in [1/m]
            lambda = 2*pi/kvector; % in [m]
            pixel_size = 172e-6; % in [m]
            th_Bragg = 4.7685/2; % in [deg], taken from lab book 208, page 117
            SplitQ_CTR = 3.1e-3; % taken from IInormb figure in 1/Angstrom
            t_osc = 700; % taken from integrated ROIS fluctuations in [s]
            %}
            
            % March 2018
             D_ds = 4.03; % in [m]
            kvector = 9.12e10; % in [1/m]
            lambda = 2*pi/kvector; % in [m]
            pixel_size = 172e-6; % in [m]
            th_Bragg = 4.8491/2; % in [deg], taken from lab book 206, page 55
             SplitQ_CTR = 3.1e-3; % taken from IInormb figure in 1/Angstrom
            t_osc = 700; % taken from integrated ROIS fluctuations in [s]
            
        end
        
    end
end