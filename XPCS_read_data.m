
classdef XPCS_read_data
    % This library contains all the functions which allow us to read the
    % data taken for an XPCS measurement and to prepare them for the
    % analysis
    properties(Constant)
    end
    
    
    methods(Static)
        
        
        function [ROIS_struct] = TTsput_prepare_ROIS(iT,XCENV,YCENV,XWIDV,YWIDV,ymax)
            
            XCEN = XCENV(iT);
            YCEN = YCENV(iT);
                        
            XWID = XWIDV(iT);
            YWID = YWIDV(iT);
            
            
             ROIS_struct.ROIS = [...
                (XCEN + [-XWID XWID]) (YCEN + [-YWID YWID])
                
                (XCEN + [-XWID XWID]) (YCEN-1 + [-1 1])
                
                (XCEN + [-XWID XWID]) (YCEN+2 - ymax(iT) + [-2 2])
                
                (XCEN + [-XWID XWID]) (YCEN + ymax(iT) + [-5 5])];
            
            
%             ROIS_struct.ROIS = [...
%                 (XCEN + [-XWID XWID]) (YCEN + [-YWID YWID])
%                 
%                 (XCEN + [-XWID XWID]) (YCEN + [-2 2])
%                 
%                 (XCEN + [-XWID XWID]) (YCEN - ymax(iT) + [-5 5])
%                 
%                 (XCEN + [-XWID XWID]) (YCEN + ymax(iT) + [-5 5])];
                                   
            ROIS_struct.AXISdet = [(XCEN + [-50 50]) (YCEN + YWID) ];
            
            ROIS_struct.TTROIS = [...
                (XCEN + [-XWID XWID]) (YCEN - mod(YCEN,4) + [-148 -145]) ];
            
            ROIS_struct.TTROIS = [...
                (XCEN + [-XWID XWID]) (YCEN + ymax(iT) + [-10 -7]) ];
            
            ROIS_struct.ttcroiout = [1]; % which TTROIs to use together to get 2-time outside X window
            ROIS_struct.DYavgwin = 25;
            ROIS_struct.ttcroiin = [1]; % which TTROIs to use together to get 2-time inside X window
            
            ROIS_struct.RYINC = 4; % Amount to increment ROI in X
            
            %DXHmin= round(ymax(iT)*0.75);
            %DXHmax= round(ymax(iT)*1.5);
        end
        
       
        function  [IIstruct] = TTsput_read(iT,TCV,specfilenameM,SCNstrM,DOCU0,DOCU1,ImageJ,Xstepcol,BKG,scanflag,imname,p_image,ending,POINTSUMS)
               
              %tamount = tamountv(iT);
            TC = TCV(iT);
            DOCUscan = [num2str(TC) 'C ' DOCU1];
                                    
            disp(['File index: ' num2str(iT)]);
            specfilename = specfilenameM(iT,:);
            SCNstr = SCNstrM(iT,:);
            
            STR = XPCS_read_data.read_paths_prepare_STR(scanflag,imname,p_image,ending);
            
            % read data and calculate the normalization
            index_SCN = 1; % if multiple SCNs, write array
            [II_orig,sdata,timestampX,TITLEstuct] = XPCS_read_data.read_data_MPX3(specfilename,STR,SCNstr,index_SCN,DOCU0,DOCUscan);
            
            for jjj = 1:size(II_orig,3)
               II_transp = squeeze(II_orig(:,:,jjj))';
               II(:,:,jjj) = II_transp;
            end
            
            [Norm] = XPCS_read_data.calc_Norm(sdata);
            
            % read time
            timestampX_flag = 0;   % timestamp flag from tif is not great, keep use spec
            lastframes_ini = [];
            [timeX,timestampX,lastframes, Xsteps,Xamount,SCNXLABEL] = XPCS_read_data.calc_TimeX(sdata,timestampX,timestampX_flag,lastframes_ini,Xstepcol,ImageJ);
            
            % correct data: normalization,background and flat field corrections
            BKG_FF_Flag = 0;
            imnormnan = [];
            [IInormb] = XPCS_read_data.from_II_to_IInorm(II,Norm,BKG,BKG_FF_Flag,ImageJ,imnormnan);
            %{
Flag_pixels = 1;
[ii_rows,jj_cols]= XPCS_read_data.find_high_counts_pixels(IInorm,Flag_pixels);
            %}
            
            
            % store data in a structure
            IIstruct.IInormb = IInormb;
            IIstruct.TITLEstuct = TITLEstuct;
            IIstruct.timeX = timeX;
            IIstruct.timestampX = timestampX;
            IIstruct.lastframes = lastframes;
            IIstruct.Xsteps = Xsteps;
            IIstruct.Xamount = Xamount;
            IIstruct.SCNXLABEL = SCNXLABEL;
            IIstruct.Nr = size(II,1); % detector size (rows), 195 Pixirad (nu)
            IIstruct.Nc = size(II,2); % detector size (rows), 487 Pixirad (del)
            IIstruct.Nt = size(II,3); % detector size (rows), 487 Pixirad (del)
            IIstruct.ROISfull = [1 IIstruct.Nc 1 IIstruct.Nr] - ImageJ;
            IIstruct.SPECpts = [1:IIstruct.Nt] - ImageJ;
            IIstruct.YROWpts = [1:IIstruct.Nr] - ImageJ;
            IIstruct.XCOLpts = [1:IIstruct.Nc] - ImageJ;
            IIstruct.DOCUInt = '[No FF]';
            
            if isempty(POINTSUMS)
                IIstruct.POINTSUMS=[1 IIstruct.Nt] - ImageJ;
            end
            
            
            
        end
            
        function  [IIstruct] = TTsimulation_read(iT,TCV,specfilenameM,SCNstrM,DOCU0,DOCU1,ImageJ,POINTSUMS)
               
              %tamount = tamountv(iT);
            TC = TCV(iT);
            DOCUscan = [num2str(TC) 'C ' DOCU1];
                                    
            disp(['File index: ' num2str(iT)]);
            specfilename = specfilenameM(iT,:);
            SCNstr = SCNstrM(iT,:);
            
            filepath = [specfilename SCNstr '/'];
            
            load([filepath specfilename SCNstr '_corr_dt.mat']);
            
          
            [Norm] = XPCS_read_data.calc_Norm_simulation(III,1e5);
            
             % correct data: normalization,background and flat field corrections
            BKG_FF_Flag = 0;
            BKG = [];
            imnormnan = [];
            [IInormb] = XPCS_read_data.from_II_to_IInorm(III,Norm,BKG,BKG_FF_Flag,ImageJ,imnormnan);
            
            
            % retain useful indexes for time (throw the 2 first ML)
             idt = damono > dt_minML;
            
            
             TITLEstuct.TITLEstr1 = char(...
                ['simulation #', SCNstr]);
            
            TITLEstuct.TITLEstrshort = char(...
                ['simulation #', SCNstr]);
            
            TITLEstuct.TITLEstr2 = char(...
                ['simulation #', SCNstr]);
            
            % calculate the ML/frame:
            MLdata = dtdata*gasconc;   % 
            
            
            % store data in a structure
            IIstruct.IInormb = IInormb;
            IIstruct.TITLEstuct = TITLEstuct;
            IIstruct.timeX = dtdata*[1:sum(idt)];%damono(idt);%./MLdata;
            IIstruct.timestampX =  dtdata*[1:sum(idt)];%damono(idt);
            IIstruct.lastframes = [1 sum(idt)- ImageJ];
            IIstruct.Xsteps =  [1:sum(idt)];%damono(idt);%./MLdata;
            IIstruct.Xamount =  dtdata*[1:sum(idt)];%damono(idt);%./MLdata;
            IIstruct.SCNXLABEL =  'spec point # ';
            IIstruct.Nr = size(III,1); % detector size (rows), 195 Pixirad (nu)
            IIstruct.Nc = size(III,2); % detector size (rows), 487 Pixirad (del)
            IIstruct.Nt = size(III,3); % detector size (rows), 487 Pixirad (del)
            IIstruct.ROISfull = [1 IIstruct.Nc 1 IIstruct.Nr] - ImageJ;
            IIstruct.SPECpts = [1:IIstruct.Nt] - ImageJ;
            IIstruct.YROWpts = [1:IIstruct.Nr] - ImageJ;
            IIstruct.XCOLpts = [1:IIstruct.Nc] - ImageJ;
            IIstruct.DOCUInt = '[No FF]';
            IIstruct.MLdata = MLdata;
            
            if isempty(POINTSUMS)
                IIstruct.POINTSUMS=[1 IIstruct.Nt] - ImageJ;
            end
            
            
            
        end
        
        
        function STR = read_paths_prepare_STR(scanflag,imname,p_image,ending)
            % this function prepares the STR structure which generates the
            % paths where the AREA detector images are stored.
            
            [SPECpath,AREApath,COMMONpath,HOMEpath] = pathdisplay;
            
            STR.scanflag=scanflag;
            STR.imname = imname;
            STR.SPECpath = SPECpath;
            STR.AREApath = AREApath;   % need to change pathdisplay if need pilatus
            STR.p_image = p_image;   %% [] if it uses the preferred points instead
            STR.ending = ending;

            
            
         end
        
        
        function [II,sdata,timestampX,TITLEstuct] = read_data_MPX3(specfilename,STR,SCNstr,index_SCN,DOCU0,DOCUscan)
            % reads the actual AREA detector files taken at each time step
            % in an XPCS mesurement and stores it in the matrix II.
            
            
            
            % FORCE only one scan at a time
            SCNs = eval(SCNstr);
            SCNs = SCNs(index_SCN);
            
            [NameMatrix,sdata] = make_imnames_2018_08(specfilename,SCNs,STR);
            FullNameArea = addnames2matrix([STR.AREApath,filesep],NameMatrix.fullfilenames);
            
            [II,timestampX] = load_MPX3(FullNameArea);
            
            % prepare the titles of the figures
            
            TITLEstuct.TITLEstr1 = char(...
                [pfilename(specfilename),' #', SCNstr,' : ', sdata.SCNDATE{1}],...
                [sdata.SCNTYPE{1}],...
                [DOCU0,' ',DOCUscan]);
            
            TITLEstuct.TITLEstrshort = char(...
                [pfilename(specfilename),' #',SCNstr, ' : ',sdata.SCNTYPE]);
            
            TITLEstuct.TITLEstr2 = char(...
                [pfilename(specfilename),' #',SCNstr]);
            
            
        end
        
        
        function [Norm] = calc_Norm(sdata)
            % This function reads information about sdata and calculates
            % the normalizing factor
            
            hubmon	= sdata.DATA(:,chan2col(sdata.collabels,'hexmon'));
            MONave 	= mean(hubmon);
            %valve_p	= sdata.DATA(:,chan2col(sdata.collabels,'valve_p'));
            %valve_v	= sdata.DATA(:,chan2col(sdata.collabels,'valve_v'));
            %secperpoint	= sdata.DATA(:,chan2col(sdata.collabels,'Seconds'));
            secperpoint	= sdata.DATA(:,chan2col(sdata.collabels,'Seconds'));
           
            
            % for norm use mean hexmon - since slit size changes a lot
            % we do not use the stated hex100mA from run params
            hex100mA = MONave./mean(secperpoint);
            
            if size(hubmon,2) > 1
                Norm	= hex100mA .* secperpoint ./ hubmon(:,1);   % sometimes extra columns appear with hexmon
            else
                Norm	= hex100mA .* secperpoint ./ hubmon;   % sometimes extra columns appear with hexmon                
            end
        end 
        
        function [Norm] = calc_Norm_simulation(III,cntrate)
            % This function reads information about sdata and calculates
            % the normalizing factor
            
            
            max_III = max(max(max(III)));
            
            Norm = (cntrate/max_III)*ones(size(III,3),1);
           
            
        end 
        
        
        function [timeX,timestampX,lastframes, Xsteps,Xamount,SCNXLABEL] = calc_TimeX(sdata,timestampX,timestampX_flag,lastframes,Xstepcol,ImageJ)
            % Note the timestamp for the medipix is only to 1 second, not to milliseconds
            % So need to use the spec information for time.
            %hubmon	= sdata.DATA(:,chan2col(sdata.collabels,'hubmon'));
            hubmon	= sdata.DATA(:,chan2col(sdata.collabels,'hexmon'));
            timestampSpec	= sdata.DATA(:,chan2col(sdata.collabels,'Time'));
            timestampEpoch  = sdata.DATA(:,chan2col(sdata.collabels,'Epoch'));
            
            if timestampX_flag==0
                timestampX = timestampSpec + timestampEpoch(1);
            end
            
            timeX = timestampX - timestampX(1); % from start of scan
            
            % Find valve switch points  (note that this output 1st point is 0th point (like spec)
            %disp(['Last frames (Spec point number convention) before valve switches: ' num2str(lastframes')]);
            NumOfLastFrames = numel(lastframes);
            
            switch NumOfLastFrames
                
                case 0
                    lastframes = [1 length(hubmon)-1];
                    
                case 1
                    lastframes = [lastframes length(hubmon)-1];
            end
            
             % use scan variable
            
            if isempty(Xstepcol)
                Xstepcol=1;
            end   
            
            if Xstepcol<1
                SCNXLABEL = 'spec point # ';
                Xsteps	= [0:length(timestampX)-ImageJ];
            else
                SCNXLABEL	= col2chan(sdata.collabels,Xstepcol);
                 Xsteps	= sdata.DATA(:,Xstepcol);
            end
            
            if strncmp(SCNXLABEL,'Time',4)
                SCNXLABEL	= [SCNXLABEL,' [\Delta sec] from SpecFile '];  % in pixirad,
                Xsteps	= timeX;
            elseif Xstepcol<1
                SCNXLABEL = 'spec point # ';
                Xsteps	= [0:length(timestampX)-ImageJ];               
            end
            
            % Express time axis in sec from start
            Xamount = timeX;
            
            
        end
        
        function [IInorm] = from_II_to_IInorm(II,Norm,BKG,BKG_FF_Flag,ImageJ,imnormnan)
            
            % This function normalizes, corrects for the bkg, and the flatfield
            
            Nt = size(II,3); % number of images in scan
            
            IInorm = NaN*ones(size(II));
            for ii=1:Nt % could be Nt
                IInorm(:,:,ii)	= II(:,:,ii) .* Norm(ii);
            end
            % Calculate mean background image
            if isempty(BKG)
                IIimbkg = 0;
            else
                IIimbkg	= mean(IInorm(:,:,[BKG(1):BKG(2)]+ImageJ),3);
            end
            
            % correct for bkg and flatfield
            if BKG_FF_Flag == 1
                
                %load flatfield_MPX3_2017_08

                % For MPX used flatfield / per-pixel normalization
                % Run flatfield_analysis_MPX3.m to make flatfield_MPX3.mat
                % contains 2 images: 
                % Badimage (bad pixels); 
                %imnormnan (pixel normalization, bad = NaN)
                
                IInormb = NaN*ones(size(II));
                for ii=1:Nt
                    IInormb(:,:,ii)	= (IInorm(:,:,ii)- IIimbkg)./imnormnan;
                end
                IInorm = IInormb;
                
            else 
                disp('[ Using no flatfield] ');
            end
            
            
        end
        
        function [ii_rows,jj_cols]= find_high_counts_pixels(IInormb,Flag_pixels)
           
            if Flag_pixels == 0
                % Pixels with high counts where i is y (row), j is x (col)
                ii_rows = [251 256 261 261 261 261 261 261 261 261 256 327 261 261 261 261 261 261 261 261 256 261 261 261 261 235 290 261 348 261 277 285 295 298 372 378 388 409 420 428 457 464 494 513 250 252 262 268 287 305 306 308 310 311 312 313 315 316 317 325 326 328 330 331 332 343 358 367 406 325 330 256 329 320 335 256 309 261 256 261 256 256 256 411   3 256 261 419 307 256 261   3   8 256 516 511];
                jj_cols = [  7  31  33  42  47  55  62  64  71  76  88  94 108 118 136 137 140 141 148 152 155 159 183 184 223 233 238 240 254 256 256 256 256 256 256 256 256 256 256 256 256 256 256 256 261 261 261 261 261 261 261 261 261 261 261 261 261 261 261 261 261 261 261 261 261 261 261 261 261 262 262 263 263 264 264 270 270 278 280 315 316 355 358 359 393 391 395 398 401 424 424 435 446 499 502 505];
            else
                % Can find these by e.g.
                bp = find(max(IInormb,[],3)>1500);
                [ii_rows,jj_cols] = ind2sub([516 516],bp);
            end
            
        end
        
        function Nfluct = find_pixel_high_fluctuations(IInormb,flimsd_flag)
            
            % Getting some pixels that fluctuate high for a single time point, e.g.
            % 414 232
            % Locate by looking for big changes in time
            disp('Fixing single-point fluctuations');
            dI = diff(IInormb,1,3);
            Nt = size(II,3); % number of images in scan
            
            
            if flimsd_flag
                % Could assume pattern fairly stable in time, normal fluctuations are counting
                % statistics
                mI = mean(IInormb,3);
                flim = flimsd*sqrt(mI);
            else
                % Alternately, just look for big ones
                flim = flimc;
            end
            
            Nfluct = 0;
            
            for ii = 1:Nt
                if ii == 1
                    fluct = dI(:,:,1)<-flim;
                    Ibase = IInormb(:,:,1);
                    Ifix = IInormb(:,:,2);
                    Ibase(fluct) = Ifix(fluct);
                    IInormb(:,:,1) = Ibase;
                elseif ii == Nt
                    fluct = dI(:,:,Nt-1)>flim;
                    Ibase = IInormb(:,:,Nt);
                    Ifix = IInormb(:,:,Nt-1);
                    Ibase(fluct) = Ifix(fluct);
                    IInormb(:,:,Nt) = Ibase;
                else
                    fluct = dI(:,:,ii)<-flim | dI(:,:,ii-1)>flim;
                    Ibase = IInormb(:,:,ii);
                    Ifix = (IInormb(:,:,ii-1)+IInormb(:,:,ii+1))/2;
                    Ibase(fluct) = Ifix(fluct);
                    IInormb(:,:,ii) = Ibase;
                end
                Nfluct = Nfluct + sum(fluct(:));
            end
            
            disp([num2str(Nfluct) ' found.']);
            
            
        end
        
       
            
            
    end
    
end