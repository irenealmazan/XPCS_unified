
classdef XPCS_analysis
    % This library contains all the functions which allow us to analyze the
    % time correlatin functions
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [ROIS_struct] = TTsput_prepare_ROIS(iT,XCENV,YCENV,XWIDV,YWIDV,y_or_xmax)
            
            XCEN = XCENV(iT);
            YCEN = YCENV(iT);
                        
            XWID = XWIDV(iT);
            YWID = YWIDV(iT);
            
            %{
            ROIS_struct.ROIS = [...
                (XCEN + [-XWID XWID]) (YCEN + [-YWID YWID])
                
                (XCEN + [-XWID XWID]) (YCEN + [-2 2])
                
                (XCEN + [-XWID XWID]) (YCEN - ymax(iT) + [-5 5])
                
                (XCEN + [-XWID XWID]) (YCEN + ymax(iT) + [-5 5])];
            %}
            
            %%{
            ROIS_struct.ROIS = [...
                (XCEN + [-XWID XWID]) (YCEN + [-YWID YWID])
                
                (XCEN + [-2 2]) (YCEN + [-YWID YWID]) %(YCEN + [-2 2])
                
                (XCEN - y_or_xmax(iT) + [-5 5]) (YCEN + [-YWID YWID]) %(YCEN - ymax(iT) + [-5 5])
                
                (XCEN  + y_or_xmax(iT) + [-5 5]) (YCEN + [-YWID YWID])]; %(YCEN + ymax(iT) + [-5 5])];
              %}         
            %{
            ROIS_struct.ROIS = [...
                (XCEN + [-XWID XWID]) (YCEN + [-YWID YWID])
                
                (XCEN + [-2 2]) (YCEN + [-YWID YWID]) %(YCEN + [-2 2])
                
                (86 + [-5 5]) (204+ [-5 5]) %(YCEN - ymax(iT) + [-5 5])
                
                (XCEN  + y_or_xmax(iT) + [-5 5]) (YCEN + [-YWID YWID])]; %(YCEN + ymax(iT) + [-5 5])];
            %}
            ROIS_struct.AXISdet = [(XCEN + [-50 50]) (YCEN + YWID) ];
            
            ROIS_struct.TTROIS = [...
                (XCEN + [-XWID XWID]) (YCEN - mod(YCEN,4) + [-148 -145]) ];
            
            ROIS_struct.TTROIS = [...
                (XCEN + [-XWID XWID]) (YCEN + y_or_xmax(iT) + [-10 -7]) ];
            
            ROIS_struct.ttcroiout = [1]; % which TTROIs to use together to get 2-time outside X window
            ROIS_struct.DYavgwin = 25;
            ROIS_struct.ttcroiin = [1]; % which TTROIs to use together to get 2-time inside X window
            
            ROIS_struct.RYINC = 4; % Amount to increment ROI in X
            
            %DXHmin= round(ymax(iT)*0.75);
            %DXHmax= round(ymax(iT)*1.5);
        end
          
        function Qval_struct = calculate_qval(XCEN_II,YCEN_II,Ncq_vect,Nrq_vect,sim_flag)
            
            if sim_flag
                Qval_struct.nu = (Ncq_vect-XCEN_II);% del direction  in 1/Angstroms
                Qval_struct.del = (Nrq_vect-YCEN_II ); %nu direction, in 1/Angstroms
                
            else
                [D_ds,kvector,lambda,pixel_size,th_Bragg,SplitQ_CTR,t_osc] = XPCS_initialize_parameters.TT_experiment();
                Qval_struct.nu = (kvector*(Ncq_vect-XCEN_II )*pixel_size/D_ds)*1e-10; % nu direction, in 1/Angstroms
                Qval_struct.del = (kvector*(Nrq_vect-YCEN_II)*pixel_size/D_ds)*sind(th_Bragg)*1e-10; % del direction  in 1/Angstroms
                
            end
            
        end        
        
        function [itt_range_struct] = prepare_subrange_for2corr(iT,timeX,XCENV,YCENV,XWIDV,YWIDV,tminv,tmaxv)
            % subrange for 2-time correlation function calculation
            
            itt_range_struct.xx = [-XWIDV(iT):XWIDV(iT)]; % Pixel locations relative to center 
            itt_range_struct.yy = [-YWIDV(iT):YWIDV(iT)]; % Pixel locations relative to center
            itt_range_struct.ixtt = XCENV(iT) + itt_range_struct.xx;
            itt_range_struct.iytt = YCENV(iT) + itt_range_struct.yy;
            itt_range_struct.Ncs = length(itt_range_struct.xx);
            itt_range_struct.Nrs = length(itt_range_struct.yy);
            itt_range_struct.ittt = find(timeX > tminv(iT) & timeX < tmaxv(iT));
            itt_range_struct.Nts = length(itt_range_struct.ittt);
            
            
        end       
      
        function [IIbin_struct] = bin_scans_in_time(Singlescan_struct,iT,tbinV,ImageJ,POINTSUMS)
            
            % Use sub-range of pixels and time
            
            iytt = Singlescan_struct.itt_range_struct.iytt;
            ixtt = Singlescan_struct.itt_range_struct.ixtt;
            ittt = Singlescan_struct.itt_range_struct.ittt;
            
          
            
            IInormbs = Singlescan_struct.IIstruct.II(iytt,ixtt,ittt);%Read_Singlescan_struct.IIstruct.IInormb(iytt,ixtt,ittt);
            
            Nrs = size(IInormbs,1);
            Ncs = size(IInormbs,2);
            Nts = size(IInormbs,3);
            
            %timestampXs = Read_Singlescan_struct.IIstruct.timestampX(ittt);
            timeX = Singlescan_struct.IIstruct.timeX(ittt);%Read_Singlescan_struct.IIstruct.timeX(ittt);
            %Xsteps = Read_Singlescan_struct.IIstruct.Xsteps(ittt);
            
            tbin = tbinV(iT);
            
            % First bin scans in time
            if tbin > 1
                tbin = floor(tbin);
                Ntb = floor(Nts/tbin);
                IIbin_struct.Ntb =  Ntb;
                IInormbe = reshape(IInormbs(:,:,1:Ntb*tbin),Nrs,Ncs,tbin,Ntb);
                IIbin_struct.II = squeeze(sum(IInormbe,3));
                timeX_bin = reshape(timeX(1:Ntb*tbin),tbin,Ntb);
                %Xstepse = reshape(Xsteps(1:Ntb*tbin),tbin,Ntb);
                IIbin_struct.timeX = squeeze(mean(timeX_bin,1));
                %IIbin_struct.Xstepsb = squeeze(mean(Xstepse,1));
                %timestampe = reshape(timestampXs(1:Ntb*tbin),tbin,Ntb);
                %IIbin_struct.timeb = squeeze(mean(timestampe,1));
                
                if ~isempty(POINTSUMS)
                    POINTSUMSB(:,1) = floor((POINTSUMS(:,1)-1+ImageJ)/tbin)+1-ImageJ;
                    POINTSUMSB(:,2) = floor((POINTSUMS(:,2)+ImageJ)/tbin) - ImageJ;
                else
                    POINTSUMSB = POINTSUMS;
                end
                
            else
                Ntb = Nts;
                IIbin_struct.Ntb = Nts;
                IIbin_struct.II = IInormbs;
                %IIbin_struct.Xamountb = Xamounts';
               % IIbin_struct.Xstepsb = Xsteps;
                %IIbin_struct.timeb = timestampXs;
                %POINTSUMSB = POINTSUMS;
            end
%             
%             if isempty(POINTSUMSB)
%                 IIbin_struct.POINTSB=[1 IIbin_struct.Ntb]-ImageJ;
%             else
%                 IIbin_struct.POINTSB = POINTSUMSB;
%             end
%                         
            IIbin_struct.YROWpts = [1:Nrs] - ImageJ;
            IIbin_struct.XCOLpts = [1:Ncs] - ImageJ;
            IIbin_struct.TITLEstuct = Singlescan_struct.IIstruct.TITLEstuct;%Read_Singlescan_struct.IIstruct.TITLEstuct;
        end
        
        function [IInormbb_ref_avg] = calc_IInormb_avg_inpixels(IIbin_struct,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,sim_flag,mask)
           % This function fits the temperal dependence of each pixel in
           % IInormb to a polynomial function whose degree is specified by
           % N_degree
           
            IInormbb = IIbin_struct.II;
            timeXb = IIbin_struct.timeX;
            
            Nrs = size(IInormbb,1);
            Ncs = size(IInormbb,2);
            Ntb = size(IInormbb,3);
            
            hwttr = hwttr_allT(iT);
            hwttc = hwttc_allT(iT);
            wrq = wrq_allT(iT);
            wcq = wcq_allT(iT);
            offsetcc = offsetcc_allT(iT);
            offsetrc = offsetrc_allT(iT);
            
            mask_struct.offsetr = mask.offsetrc_allT(iT);
            mask_struct.hwttr = ones(1,numel([-hwttr:hwttr]));
            if mask.hwttr_allT(iT) ~= 0
                midindex_hwttr = round(numel(mask_struct.hwttr)/2);%+mask_struct.offsetr;
                mask_struct.hwttr(midindex_hwttr-mask.hwttr_allT(iT):midindex_hwttr+mask.hwttr_allT(iT)) =  zeros(numel([-mask.hwttr_allT(iT):mask.hwttr_allT(iT)]),1);
            end
            
            mask_struct.offsetc = mask.offsetcc_allT(iT);
            mask_struct.hwttc = ones(1,numel([-hwttc:hwttc]));
            if mask.hwttc_allT(iT) ~= 0
                midindex_hwttc = round(numel(mask_struct.hwttc)/2);%+mask_struct.offsetc;
                mask_struct.hwttc(midindex_hwttc-mask.hwttc_allT(iT):midindex_hwttc+mask.hwttc_allT(iT)) = zeros(numel([-mask.hwttc_allT(iT):mask.hwttc_allT(iT)]),1);
            end
            
            
            
            
           
            
            ittccen =1 + (Ncs - 1)/2 ;% index of col center
            ittrcen = 1 + (Nrs - 1)/2; % index of row center        

            
            Ncq = 2*wcq + 1;
            Nrq = 2*wrq + 1;
         
            for icq = 1:Ncq
                 offttc = (icq-wcq-1)*(2*hwttc(icq)+1)+ offsetcc;  
                %indexc_non_masc = [-hwttc:hwttc] - mask_struct.hwttc;
                ittc_index = round(ittccen + offttc + [-hwttc:hwttc]) ;
                ittc = ittc_index(find(mask_struct.hwttc~=0));
               
                for irq = 1:Nrq
                    offttr = (irq - wrq - 1)*(2*hwttr+1)+ offsetrc;
                    ittr = round(ittrcen + offttr + [-hwttr:hwttr]);
                    
                    Qval_struct = XPCS_analysis.calculate_qval(ittccen,ittrcen,ittccen + offttc,ittrcen + offttr,sim_flag);
                    
                    IInormbb_ref_avg.scancq(icq).scanrq(irq).II = squeeze(mean(mean(IInormbb(ittr,ittc,:),1),2));
                    
                    IInormbb_ref_avg.qvector.nu(irq) = Qval_struct.nu;
                    
                    IInormbb_ref_avg.scancq(icq).scanrq(irq).timex = timeXb;
                end
                IInormbb_ref_avg.qvector.del(icq) = Qval_struct.del;
            end
            
            IInormbb_ref_avg.Ncq_Nrq = [Ncq Nrq];
            IInormbb_ref_avg.TITLEstruct = IIbin_struct.TITLEstuct;

       end
        
        function [CCN2_struct,IIbin_struct] = calc_2time_corr(IIbin_struct,iT,N_degree,flag_mean_or_poly)
                     
            IInormbb =  IIbin_struct.II;
            
            Nrs = size(IInormbb,1);
            Ncs = size(IInormbb,2);
            Ntb = size(IInormbb,3);
            
            switch flag_mean_or_poly
                case 'mean'
                    Ibar = mean(IInormbb,3);
                case 'poly'
                    Ibar = FittingFunctions.fit_IInormb_intime(IIbin_struct, N_degree(iT));
                    IIbin_struct.N_degree = N_degree(iT);
            end
            
            IIbin_struct.Ibar = Ibar;
            
            
            dI = IInormbb - Ibar;
            dlnI = IInormbb./Ibar - 1;
             
            % Calc 2-time using time average mean, no ensemble
            
            %CCN2_struct.IIM2 = NaN*ones(Nrs,Ncs,Ntb,Ntb);
            CCN2_struct.CCN2 = NaN*ones(Nrs,Ncs,Ntb,Ntb);
            
            for ii = 1:Ntb
                for jj = 1:ii
                    %CCN2_struct.IIM2(:,:,ii,jj) = dI(:,:,ii).*dI(:,:,jj);
                    %CCN2_struct.IIM2(:,:,jj,ii) = CCN2_struct.IIM2(:,:,ii,jj);
                   
                    CCN2_struct.CCN2(:,:,ii,jj) = dlnI(:,:,ii).*dlnI(:,:,jj);
                    CCN2_struct.CCN2(:,:,jj,ii) = CCN2_struct.CCN2(:,:,ii,jj);
                end
            end
            %IID2 = diag(IIM2);
            %CC2 = IIM2./sqrt(IID2*IID2'); % Normalized to make diagonal unity
            CCN2_struct.TITLEstuct = IIbin_struct.TITLEstuct;
        end
        
        function [IIbin_struct] = calc_Ibar_sgsf2D_space(IIbin_struct,maxpd,jjj,iii,do_mean)
                     
            III =  IIbin_struct.II;
            %maxpd = 2; % determines maximum degree of smoothing polynomial
            %iii = 1;%2;%3; %5;% half-width of pixel range in columns
            %jjj = 3;%2;%3; %5;% half-width of pixel range in rows
        
            if do_mean
                % Smooths the time average of III
                [Ibar,filt] = Functions_sos.smooth_sg_space_mean(III,maxpd,iii,jjj);%mean(III,3);
            else
                % Smooths each time frame of III
                [Ibar,filt] = Functions_sos.smooth_sg_space(III,maxpd,iii,jjj);%mean(III,3);
            end
            
            IIbin_struct.Ibar = Ibar;
            IIbin_struct.filt = filt;            
            %IIbin_struct.timeXb =  IIbin_struct.timeXb(1+jjj:Ntb-jjj);
            
            nrow = size(Ibar,1);
            ncol = size(Ibar,2);               
            
            IIbin_struct.II = III(jjj+1:jjj+nrow,iii+1:iii+ncol,:);
            
            
            
        end     
        
        function [IIbin_struct] = calc_Ibar_sgsf2D_space_and_time(IIbin_struct,maxpd,jjj,iii)
                     
            III =  IIbin_struct.II;
            Nr = size(III,1);           
            Nc = size(III,2);
            Ntb = size(III,3);
            %maxpd = 2; % determines maximum degree of smoothing polynomial
            %iii = 2; % half-width of pixel range in columns/time
            %jjj = 2; % half-width of pixel range in rows/space
        
            [Ibar,filt] = Functions_sos.smooth_sg_space_and_time(III,maxpd,iii,jjj);%mean(III,3);
        
            
            IIbin_struct.Ibar = Ibar;
            IIbin_struct.filt = filt;            
            IIbin_struct.timeX =  IIbin_struct.timeX(1+iii:Ntb-iii);
            
            
            
            IIbin_struct.II = III(jjj+1:Nr-jjj,1:Nc,1+iii:Ntb-iii);
            IIbin_struct.Ntb = size(IIbin_struct.II,3);
            
            
            
        end     
        
        
        function [CCN2_struct,IIbin_struct] = calc_2time_corr_sgsf2D(IIbin_struct)
                     
            III =  IIbin_struct.II;
            Ibar = IIbin_struct.Ibar;
            
            Nrs = size(III,1);
            Ncs = size(III,2);
            Ntb = size(III,3);
            
            % calculate speckels and use the total intensity to reescale
            % and remove Ibar:
            
            III_norm = squeeze(sum(sum(III,1),2)); 
            Ibar_norm = squeeze(sum(sum(Ibar,1),2));
        
            
            for tt = 1:Ntb
                dlnI(:,:,tt) = (Ibar_norm(tt).*III(:,:,tt))./(III_norm(tt).*Ibar(:,:,tt)) - 1;
            end
            
            
            %dI = IInormbbc - IInormbb_ref;
            %dlnI =  III./Ibar - 1;
             
            % Calc 2-time using time average mean, no ensemble ??
                        
            CCN2_struct.CCN2 = NaN*ones(Nrs,Ncs,Ntb,Ntb);
            
            for ii = 1:Ntb
                for jj = 1:ii
                   
                    CCN2_struct.CCN2(:,:,ii,jj) = dlnI(:,:,ii).*dlnI(:,:,jj);
                    CCN2_struct.CCN2(:,:,jj,ii) = CCN2_struct.CCN2(:,:,ii,jj);
                end
                
                % correction of the Poisson noise on the diagonal
                CCN2_struct.CCN2(:,:,ii,ii) = CCN2_struct.CCN2(:,:,ii,ii)- 1./Ibar(:,:,ii);
                
            end            
            CCN2_struct.TITLEstuct = IIbin_struct.TITLEstuct;
        end
        
       
        function [CCN2avg_struct,qvector_struct] = from_CCN2V_to_CCN2avg(Single_scan_struct,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,sim_flag)
            
            hwttr = hwttr_allT(iT);
            hwttc = hwttc_allT(iT);
            wrq = wrq_allT(iT);
            wcq = wcq_allT(iT);
            offsetcc = offsetcc_allT(iT);
            offsetrc = offsetrc_allT(iT);
            
            
            Nts = size(Single_scan_struct.CCN2_struct.CCN2,3);
            Ncs = size(Single_scan_struct.CCN2_struct.CCN2,2);
            Nrs = size(Single_scan_struct.CCN2_struct.CCN2,1);
            
            ittccen =1 + (Ncs - 1)/2 ;% index of col center
            ittrcen = 1 + (Nrs - 1)/2; % index of row center        

            CCN2V = Single_scan_struct.CCN2_struct.CCN2;
            III_avg = zeros(size(Single_scan_struct.IIbin_struct.II));
            
            
            Ncq = 2*wcq + 1;
            Nrq = 2*wrq + 1;
         
            for icq = 1:Ncq
                offttc = (icq-wcq-1)*(2*hwttc+1)+ offsetcc;               
                ittc = round(ittccen + offttc + [-hwttc:hwttc]) ;
               
                for irq = 1:Nrq
                    offttr = (irq - wrq - 1)*(2*hwttr+1)+ offsetrc;
                    ittr = round(ittrcen + offttr + [-hwttr:hwttr]);
                    
                    Qval_struct = XPCS_analysis.calculate_qval(ittccen,ittrcen,ittccen + offttc,ittrcen + offttr,sim_flag);
                    
                   % correction of the Poisson noise on the diagonal
                   III_avg = squeeze(mean(mean(Single_scan_struct.IIbin_struct.II(ittr,ittc,:),1),2));
                   CCN2avg = squeeze(mean(mean(CCN2V(ittr,ittc,:,:),1),2));
%                    for tt = 1:Nts
%                       CCN2avg(tt,tt) = CCN2avg(tt,tt) - 1./III_avg(tt); 
%                    end
%                    
                    CCN2avg_struct.scancq(icq).scanrq(irq).CCN2avg = CCN2avg ;%squeeze(mean(mean(CCN2V(ittr,ittc,:,:),1),2));
                    CCN2avg_struct.scancq(icq).scanrq(irq).timex = Single_scan_struct.IIbin_struct.timeX;
                    CCN2avg_struct.scancq(icq).scanrq(irq).TITLEstr2V = Single_scan_struct.IIbin_struct.TITLEstuct.TITLEstr2;
                    CCN2avg_struct.scancq(icq).scanrq(irq).nu = Qval_struct.nu;
                    CCN2avg_struct.scancq(icq).scanrq(irq).del = Qval_struct.del;
                    CCN2avg_struct.qvector.nu(irq) = Qval_struct.nu;
                    CCN2avg_struct.boxcenterrc.offttr(irq) = ittrcen + offttr;
                end
               CCN2avg_struct.qvector.del(icq) = Qval_struct.del;
               CCN2avg_struct.boxcenterrc.offttc(icq) = ittccen + offttc;
            end
            CCN2avg_struct.Nts_Ncs_Nrs = [Nts Ncs Nrs];
            CCN2avg_struct.Ncq_Nrq = [Ncq Nrq];
            CCN2avg_struct.TITLEstruct = Single_scan_struct.CCN2_struct.TITLEstuct;
            CCN2avg_struct.ittccen = ittccen ;
            CCN2avg_struct.ittrcen = ittrcen ;
            CCN2avg_struct.hwttr = hwttr;
            CCN2avg_struct.hwttc = hwttc;
            CCN2avg_struct.wrq = wrq;
            CCN2avg_struct.wcq = wcq;
            CCN2avg_struct.offsetcc = offsetcc;
            CCN2avg_struct.offsetrc = offsetrc;
            
            qvector_struct = CCN2avg_struct.qvector;
        end
       
        function [CCN2avg_struct,qvector_struct] = from_CCN2V_to_CCN2avg_several_regions(Single_scan_struct,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,D_ds,kvector,pixel_size,th_Bragg,mask)
            
            hwttr = hwttr_allT(iT);
            hwttc = hwttc_allT(iT);
            wrq = wrq_allT(iT);
            wcq = wcq_allT(iT);
            offsetcc = offsetcc_allT(iT);
            offsetrc = offsetrc_allT(iT);
            
            mask_struct.offsetr = mask.offsetrc_allT(iT);
            mask_struct.hwttr = ones(1,numel([-hwttr:hwttr]));
            if mask.hwttr_allT(iT) ~= 0
                midindex_hwttr = round(numel(mask_struct.hwttr)/2);%+mask_struct.offsetr;
                mask_struct.hwttr(midindex_hwttr-mask.hwttr_allT(iT):midindex_hwttr+mask.hwttr_allT(iT)) =  zeros(numel([-mask.hwttr_allT(iT):mask.hwttr_allT(iT)]),1);
            end
            
            mask_struct.offsetc = mask.offsetcc_allT(iT);
            mask_struct.hwttc = ones(1,numel([-hwttc:hwttc]));
            if mask.hwttc_allT(iT) ~= 0
                midindex_hwttc = round(numel(mask_struct.hwttc)/2);%+mask_struct.offsetc;
                mask_struct.hwttc(midindex_hwttc-mask.hwttc_allT(iT):midindex_hwttc+mask.hwttc_allT(iT)) = zeros(numel([-mask.hwttc_allT(iT):mask.hwttc_allT(iT)]),1);
            end
            
            
            
            
            
            Nts = size(Single_scan_struct.CCN2_struct.CCN2,3);
            Ncs = size(Single_scan_struct.CCN2_struct.CCN2,2);
            Nrs = size(Single_scan_struct.CCN2_struct.CCN2,1);
            
            ittccen =1 + (Ncs - 1)/2 ;% index of col center
            ittrcen = 1 + (Nrs - 1)/2; % index of row center        

            CCN2V = Single_scan_struct.CCN2_struct.CCN2;
            
            Ncq = 2*wcq + 1;
            Nrq = 2*wrq + 1;
      
            for icq = 1:Ncq
                offttc = (icq-wcq-1)*(2*hwttc(icq)+1)+ offsetcc;  
                ittc_index = round(ittccen + offttc + [-hwttc:hwttc]) ;
                ittc = ittc_index(find(mask_struct.hwttc~=0));
                
                for irq = 1:Nrq
                    offttr = (irq - wrq - 1)*(2*hwttr+1)+ offsetrc;
                    %indexr_non_masc = [-hwttr:hwttr] - mask_struct.hwttr;
                    %ittr_index = round(ittrcen + offttr + [-hwttr:hwttr]) ;
                    %ittr = ittr_index(find( mask_struct.hwttr~=0));
                    ittr = round(ittrcen + offttr + [-hwttr:hwttr]);
                    
                    Qval_struct = XPCS_analysis.calculate_qval(ittccen,ittrcen,ittccen + offttc,ittrcen + offttr,D_ds,kvector,pixel_size,th_Bragg);

                   
                    CCN2avg_struct.scancq(icq).scanrq(irq).CCN2avg = squeeze(mean(mean(CCN2V(ittr,ittc,:,:),1),2));
                    CCN2avg_struct.scancq(icq).scanrq(irq).timex = Single_scan_struct.IIbin_struct.Xamountb;
                    CCN2avg_struct.scancq(icq).scanrq(irq).TITLEstr2V = Single_scan_struct.IIbin_struct.TITLEstuct.TITLEstr2;
                    CCN2avg_struct.scancq(icq).scanrq(irq).nu = Qval_struct.nu;
                    CCN2avg_struct.scancq(icq).scanrq(irq).del = Qval_struct.del;
                    CCN2avg_struct.qvector.nu(irq) = Qval_struct.nu;
                    CCN2avg_struct.boxcenterrc.offttr(irq) = ittrcen + offttr;
                end
               CCN2avg_struct.qvector.del(icq) = Qval_struct.del;
               CCN2avg_struct.boxcenterrc.offttc(icq) = ittccen + offttc;
            end
            CCN2avg_struct.Nts_Ncs_Nrs = [Nts Ncs Nrs];
            CCN2avg_struct.Ncq_Nrq = [Ncq Nrq];
            CCN2avg_struct.TITLEstruct = Single_scan_struct.CCN2_struct.TITLEstuct;
            CCN2avg_struct.ittccen = ittccen ;
            CCN2avg_struct.ittrcen = ittrcen ;
            CCN2avg_struct.hwttr = hwttr;
            CCN2avg_struct.hwttc = hwttc;
            CCN2avg_struct.wrq = wrq;
            CCN2avg_struct.wcq = wcq;
            CCN2avg_struct.offsetcc = offsetcc;
            CCN2avg_struct.offsetrc = offsetrc;
            CCN2avg_struct.mask = mask_struct;
            
            qvector_struct = CCN2avg_struct.qvector;
        end       
        
        function [CCN2S_struct] = from_CCN2avg_to_CCN2S(CCN2avg_struct)
            % This function builds the CCN_structure where we store the
            % integrated 2-times correlation function. The inputs are:
            %   Nrq: the number of pixels/number of q's
            %   CCN_struct: the structure containing the integrated 2-times correlation function
            %   TITLEstr2V: the title string which specifies the name of
            %   the scans in iiT
            
            
            
            Ndt = CCN2avg_struct.Nts_Ncs_Nrs(1);
            idt = 1:Ndt;
            
            
            Ncq = CCN2avg_struct.Ncq_Nrq(1);
            Nrq = CCN2avg_struct.Ncq_Nrq(2);
            
            for icq = 1:Ncq
                
                CCNdt = NaN*ones(Nrq,Ndt);
                
                for irq = 1:Nrq
                    time_scan = CCN2avg_struct.scancq(icq).scanrq(irq).timex - CCN2avg_struct.scancq(icq).scanrq(irq).timex(1); % Delta in time from region start
                    CCN2S = CCN2avg_struct.scancq(icq).scanrq(irq).CCN2avg(idt,idt);
                    for ii = idt
                        CCNdt(irq,ii) = mean(diag(CCN2S,ii-1));
                    end
                    CCN2S_struct.scancq(icq).scanrq(irq).CCNdtV = CCNdt(irq,:);
                    CCN2S_struct.scancq(icq).scanrq(irq).time_1D = time_scan;
                    CCN2S_struct.scancq(icq).scanrq(irq).nu = CCN2avg_struct.scancq(icq).scanrq(irq).nu;
                    CCN2S_struct.scancq(icq).scanrq(irq).del = CCN2avg_struct.scancq(icq).scanrq(irq).del;
                end
            end
            CCN2S_struct.Ncq_Nrq = CCN2avg_struct.Ncq_Nrq;
            CCN2S_struct.Ndt = Ndt;
            CCN2S_struct.ittccen= CCN2avg_struct.ittccen;
            CCN2S_struct.ittrcen =  CCN2avg_struct.ittccen ;
            
            CCN2S_struct.TITLEstruct = CCN2avg_struct.TITLEstruct;
        end
        
        function [CCN2Sfit_struct] = fit_CCN2S(CCN2S_struct,fitfunc_param_str,fitfunc_str,param_legend,fit_range,fig_plot_fit_results)
            
            % initialize fitting results
            a1 = zeros(numel(CCN2S_struct),numel(CCN2S_struct(1).scanq),1);
            b1 = zeros(numel(CCN2S_struct),numel(CCN2S_struct(1).scanq),1);
            
          
            figh = figure(fig_plot_fit_results);
            clf;
            
            for iT = 1:numel(CCN2S_struct)
                for irq = 1:numel(CCN2S_struct(iT).scanq)
                    p_lower =[0 0 0]; 
                    p_upper =[7E-1 100 1e-2]; 
                    p_start = [6e-2 50 1e-3];
                    
                    [CCN2Sfit_struct(iT).scanq(irq).CCNdtV_fit] = FittingFunctions.fit_2time_corr(CCN2S_struct(iT).scanq(irq),fitfunc_param_str,param_legend,fitfunc_str,fit_range,p_lower,p_upper,p_start);                              

                    %DisplayFunctions_XPCS.display_fit_result(CCN_struct,iT,irq,figh);
                    
                end
               
             CCN2Sfit_struct(iT).TITLEstruct = CCN2S_struct.TITLEstruct;                        
            end
        end
        
        function [CCN2S_struct] = fit_CCN2S_with_leasqr(CCN2S_struct,iT,pin_iiT,dp_iiT,fitfunc_str,fit_range,indexq,flag_row_or_col,fig_plot_fit_results)
            
            figure(fig_plot_fit_results);
            clf;
            
            switch flag_row_or_col
                
                case 'row'
                    qtolook = [ 1:CCN2S_struct.Ncq_Nrq(2)];
                case 'col'
                    qtolook = [ 1:CCN2S_struct.Ncq_Nrq(1)];
                otherwise
                    disp('please select rows or cols')
                    return;
            end
            
            for iq = qtolook
                
                for ii_indexq = indexq
                    
                    switch flag_row_or_col
                        
                        case 'row'
                            CCN2S_struct_singlescan = CCN2S_struct.scancq(ii_indexq).scanrq(iq);
                            
                        case 'col'
                            CCN2S_struct_singlescan = CCN2S_struct.scancq(iq).scanrq(ii_indexq);
                            
                        otherwise
                            disp('please select rows or cols')
                            return;
                    end
                    
                    pin = pin_iiT(iT,:);
                    dp =  dp_iiT(iT,:);
                    w = ones(length(fit_range),1);
                    
                    [ fitres] = FittingFunctions.fit_2time_corr_with_leasqr(CCN2S_struct_singlescan,fitfunc_str,fit_range,pin,dp, w);
                    
                    %                     CCN2Sfit_struct(iT).TITLEstruct = CCN2S_struct(iT).TITLEstruct;
                    %
                    switch flag_row_or_col
                        case 'row'
                            CCN2S_struct.scancq(ii_indexq).scanrq(iq).CCNdtV_fit =  fitres;
                             CCN2S_struct.scancq(ii_indexq).scanrq(iq).CCNdtV_fit.fitfunc_str = fitfunc_str;
                        case 'col'
                            CCN2S_struct.scancq(iq).scanrq(ii_indexq).CCNdtV_fit =  fitres;
                            CCN2S_struct.scancq(iq).scanrq(ii_indexq).CCNdtV_fit.fitfunc_str = fitfunc_str;
                        otherwise
                            disp('please select rows or cols')
                            return;
                    end
                end
                
                
            end
            
            %DisplayFunctions_XPCS.display_fit_result(CCN_struct,iT,figh);
            
            
        end
       
        function [CCN2S_struct] = fit_CCN2S_with_leasqr_globalfit(CCN2S_struct,iT,pin_iiT,dp_iiT,fitfunc_str,fit_range,qtolook,indexq,flag_row_or_col,fig_plot_fit_results)
            
            figure(fig_plot_fit_results);
            clf;

            qvector = zeros(numel(qtolook),1); % the qvector to calculate the contrast as a function of q
            
            counter_q = 1;
            
            for iq = qtolook
                
                for ii_indexq = indexq
                    
                    switch flag_row_or_col
                        
                        case 'row'
                            CCN2S_struct_singlescan(counter_q) = CCN2S_struct.scancq(ii_indexq).scanrq(iq);
                             qvector(counter_q) =  CCN2S_struct_singlescan(counter_q).del;
                            
                        case 'col'
                            CCN2S_struct_singlescan(counter_q) = CCN2S_struct.scancq(iq).scanrq(ii_indexq);
                            qvector(counter_q) =  CCN2S_struct_singlescan(counter_q).nu;
                        otherwise
                            disp('please select rows or cols')
                            return;
                    end
                    
                end
                counter_q = counter_q + 1;
            end
            
            [pin_array] = FittingFunctions.prepare_coefs_forglobalfit(pin_iiT,iT,qvector);%pin_iiT(iT,:);            
            dp =  dp_iiT(iT,:);
            w = ones(length(fit_range),1);
            
            %[ fitres] = FittingFunctions.fit_2time_corr_with_leasqr(CCN2S_struct_singlescan,fitfunc_str,fit_range,pin,dp, w);
            [ fitres] = FittingFunctions.fit_2time_corr_with_leasqr_globalfit(CCN2S_struct_singlescan,fitfunc_str,fit_range,pin_array,dp, w);
            
            %                     CCN2Sfit_struct(iT).TITLEstruct = CCN2S_struct(iT).TITLEstruct;
            %
            counter_q = 1;
            for iq = qtolook
                
                for ii_indexq = indexq
                    
                    switch flag_row_or_col
                        case 'row'
                            CCN2S_struct.scancq(ii_indexq).scanrq(iq).CCNdtV_fit =  fitres(counter_q);
                            CCN2S_struct.scancq(ii_indexq).scanrq(iq).CCNdtV_fit.fitfunc_str = fitfunc_str;
                        case 'col'
                            CCN2S_struct.scancq(iq).scanrq(ii_indexq).CCNdtV_fit =  fitres(counter_q);
                            CCN2S_struct.scancq(iq).scanrq(ii_indexq).CCNdtV_fit.fitfunc_str = fitfunc_str;
                        otherwise
                            disp('please select rows or cols')
                            return;
                    end
                    counter_q = counter_q + 1;
                end
                
                
            end
            
            
            
        end
       
        
        function [pout_struct] = fit_tau(pout_struct,fitfunc_str)
                      
            
            for iT = 1:numel(pout_struct)  
                
                    %pout_struct_fitrange.fitrange = qrange;
                    %pout_struct_fitrange.qvector = pout_struct(iT).qvector(pout_struct_fitrange.fitrange);
                    %sigma = pout_struct(iT).sigma(pout_struct(iT).fitrange);
                    %pout_struct_fitrange.tau = pout_struct(iT).tau(pout_struct(iT).fitrange);
                    %[~,low_index] = find(pout_struct.qvector>qrange(1));
                    fit_range = pout_struct.fitrange; %low_index;
                
                    pin = [1e5];
                    dp =  [1]'*1e-4;
                    w = 1./pout_struct.sigma;
                    
                    [ pout_struct(iT).pout] = FittingFunctions.fit_tau_with_leasqr(pout_struct,fitfunc_str,fit_range,pin,dp, w);
                
                   
            end
        end
        
        function [pout_avg_struct] = calc_average_row_col(pcolor_struct,ppvector)
            
            for pp = ppvector
             
                pout = pcolor_struct(pp).pout;
                sigma = pcolor_struct(pp).sigma;
                
             
                 
                pout_avg_struct(pp).weight = 1./sigma.^2;
                pout_avg_struct(pp).weight_pout = pout.*pout_avg_struct(pp).weight;
                
                pout_avg_struct(pp).sum_row = sum(pout_avg_struct(pp).weight_pout,1,'omitnan');
                pout_avg_struct(pp).sum_col = sum(pout_avg_struct(pp).weight_pout,2,'omitnan');
                
                pout_avg_struct(pp).sum_row_weight = sum(pout_avg_struct(pp).weight,1,'omitnan');
                pout_avg_struct(pp).sum_col_weight = sum(pout_avg_struct(pp).weight,2,'omitnan');
                
                pout_avg_struct(pp).weight_mean_row = pout_avg_struct(pp).sum_row./pout_avg_struct(pp).sum_row_weight;
                pout_avg_struct(pp).weight_mean_col = pout_avg_struct(pp).sum_col./pout_avg_struct(pp).sum_col_weight;
                
                pout_avg_struct(pp).weight_mean = mean(pout_avg_struct(pp).weight_mean_row);
             
                
                pout_avg_struct(pp).del = pcolor_struct(pp).del;
                pout_avg_struct(pp).nu = pcolor_struct(pp).nu;
                pout_avg_struct(pp).ptitle = pcolor_struct(pp).ptitle;
                
                pout_avg_struct(pp).TITLEstruct = pcolor_struct(pp).TITLEstruct;
            end
           
        end
            
   
    end
    
end