classdef Functions_sos
    % This library contains all the functions to analyze the simulations
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [Cdt] = calc_Cdti(dt_minML,damono,ftm,nrow,ncol)
            
            % Calc delta-t correlations for each pixel
            
            % Do not use first dt_minML of growth for time correlations
            %dt_minML = 2;%if in ML
            idt = damono > dt_minML;
            
            III = abs(ftm(:,:,idt)).^2;
            ndt = size(III,3);
            ddam = damono(idt);
            ddam = ddam - ddam(1); % delta time coord (ML) (assumes evenly spaced)
            
            % Use time average as Ibar
            
            Ibar = mean(III,3);
            
            dlnI = III./Ibar - 1;
            
            Cdt = NaN*ones(size(dlnI));
            I1 = ones(ndt,1);
            Idt0 = conv(I1,I1);
            for ii = 1:nrow
                for jj = 1:ncol
                    It = squeeze(dlnI(ii,jj,:));
                    Idt = conv(It,flip(It))./Idt0;
                    Cdt(ii,jj,:) = Idt(end-ndt+1:end);
                end
            end
            
            
        end
        
        function [damHW] = calc_damHW(nrow,ncol,Cdt,ndt,dindex)
            
            %dindex can be any vector among dindex,ddam or dtime
            
            %load([runname '_corr_dt.mat']);

            
            
            damHW = dindex(end)*ones(nrow,ncol);
            fitfunc_str = 'FittingFunctions.CCN2single_fit';

            
            for ii = 1:nrow
                for jj = 1:ncol
                    Cdti = squeeze(Cdt(ii,jj,:));
                    Cdti = Cdti/Cdti(1);
                    err1 = 1;
                    
                    %{
            fit_range = [1:1:round(1/2.5*length(Cdti))];
            fitfunc_str = 'FittingFunctions.CCN2single_fit';
            pin = [0 1 1 0];%pin_iiT(iT,:);
            dp =  [[0 1 1 0]*0.0001];%dp_iiT(iT,:);
            w = ones(length(fit_range),1);
            
            CCfunc.time_1D = ddam';
            CCfunc.CCNdtV = Cdti;
            
            [ fitres] = FittingFunctions.fit_2time_corr_with_leasqr(CCfunc,fitfunc_str,fit_range,pin,dp, w);
            
            damHW(ii,jj) = fitres.pout(3);
            
            figure(100);
            clf;
            plot(CCfunc.time_1D(fit_range),CCfunc.CCNdtV(fit_range),'ob');
            hold on;
            plot(CCfunc.time_1D(fit_range),fitres.fitfunc,'r','LineWidth',3.0);
            title(['jj = ' num2str(jj) 'ii = ' num2str(ii)]);
            pause(.1);
                    %}
                    
                    for kk = 1:ndt
                        if (Cdti(kk) < 0.5);
                            err1 = 0; 
                            break; 
                        end
                    end
                    
                    if (err1==0)
                        isp = [(kk-1):kk];
                        %dindex = [1:numel(ddam)];
                        damHW(ii,jj) = interp1(Cdti(isp),dindex(isp),0.5)/log(2);                        
                    end
                    
                    if mod(ii,20) == 0
                        if mod(jj,20) == 0
                            figure(100);
                            clf;
                            plot(dindex,Cdti,'ob');
                            hold on;
                            plot(dindex,feval(fitfunc_str,dindex,[0 1.0 damHW(ii,jj) 0]),'r','LineWidth',3.0);
                            title(['jj = ' num2str(jj) 'ii = ' num2str(ii) ' tau = ' num2str(damHW(ii,jj),'%e')]);
                            pause(.1);
                        end
                    end
                end
            end
            
        end
        
        function [damHW] = calc_damHW_byfit(nrow,ncol,Cdt,dindex)
            %dindex can be any vector among dindex,ddam or dtime
            
            %load([runname '_corr_dt.mat']);
            
            
            
            damHW = dindex(end)*ones(nrow,ncol);
            
            
            for ii = 1:nrow
                for jj = 1:ncol
                    Cdti = squeeze(Cdt(ii,jj,:));
                    Cdti = Cdti/Cdti(1);
                    err1 = 1;
                    
                    %%{
                    fit_range = [1:1:round(length(Cdti)/4)];%[1:1:round(length(Cdti)/4)];
                    fitfunc_str = 'FittingFunctions.CCN2single_fit';
                    pin = [0 1 1e8 0];%pin_iiT(iT,:);
                    dp =  [[0 1 1 0]*0.0001];%dp_iiT(iT,:);
                    w = ones(length(fit_range),1);
                    
                    CCfunc.time_1D = dindex';
                    CCfunc.CCNdtV = Cdti;
                    try
                        [ fitres] = FittingFunctions.fit_2time_corr_with_leasqr(CCfunc,fitfunc_str,fit_range,pin,dp, w);
                        
                        
                        figure(100);
                        clf;
                        plot(CCfunc.time_1D(fit_range),CCfunc.CCNdtV(fit_range),'ob');
                        hold on;
                        plot(CCfunc.time_1D(fit_range),fitres.fitfunc,'r','LineWidth',3.0);
                        %set(gca,'Xscale','log','Yscale','log');
                        title(['jj = ' num2str(jj) 'ii = ' num2str(ii)]);
                        pause(.01);
                        %}
                        
                    catch
                        for kk = 1:numel(Cdti)
                            if (Cdti(kk) < 0.5);
                                err1 = 0;
                                break;
                            end
                        end
                        
                        if (err1==0)
                            isp = [(kk-1):kk];
                            %dindex = [1:numel(ddam)];
                            fitres.pout(3) = interp1(Cdti(isp),dindex(isp),0.5)/log(2);
                        end
                        
                         figure(100);
                        clf;
                        plot(dindex,Cdti,'ob');
                        hold on;
                        plot(dindex,feval(fitfunc_str,dindex,[0 1 fitres.pout(3) 0 ]),'r','LineWidth',3.0);
                        %set(gca,'Xscale','log','Yscale','log');
                        title(['jj = ' num2str(jj) 'ii = ' num2str(ii)]);
                        pause(.01);
                        
                    end
                    damHW(ii,jj) = fitres.pout(3);
                    
                    
                    
                end
            end
            
        end
        
            
        function [fit_res] = calc_damHW_byfit_double_exp(nrow,ncol,Cdt,damHW0,dindex)
            %dindex can be any vector among dindex,ddam or dtime
            
            %load([runname '_corr_dt.mat']);
            
            
            
            damHW = dindex(end)*ones(nrow,ncol);
            
            
            for ii = 1:nrow
                for jj = 1:ncol
                    Cdti = squeeze(Cdt(ii,jj,:));
                    Cdti = Cdti/Cdti(1);
                    err1 = 1;
                    
                    %%{
                    fit_range = [1:1:round(length(Cdti)/4)];%[1:1:round(length(Cdti)/4)];
                    fitfunc_str = 'FittingFunctions.CCN2single_fit_double_exp';
                    
                    pin = [0 0.8 damHW0.damHW(ii,jj) 0.1 5e8  0];%pin_iiT(iT,:);
                    dp =  [[0 1 0 1 1 0]*0.0001];%dp_iiT(iT,:);
                    
                    for tt = 1:length(fit_range)
                        w(tt) = 1./sqrt(length(fit_range)-(tt-1));%ones(length(fit_range),1);
                    end
                    
                    CCfunc.time_1D = dindex';
                    CCfunc.CCNdtV = Cdti;
                    try
                        [ fitres] = FittingFunctions.fit_2time_corr_with_leasqr(CCfunc,fitfunc_str,fit_range,pin,dp, w);
                        
                        %{
                        figure(100);
                        clf;
                        errorbar(CCfunc.time_1D(fit_range),CCfunc.CCNdtV(fit_range),w,'ob');%plot(CCfunc.time_1D,CCfunc.CCNdtV,'ob');%
                        hold on;
                        plot(CCfunc.time_1D(fit_range),fitres.fitfunc,'r','LineWidth',3.0);
                        %set(gca,'Xscale','log','Yscale','log');
                        title(['jj = ' num2str(jj) 'ii = ' num2str(ii)]);
                        pause(.01);
                        %}
                        
                    catch
                        for kk = 1:numel(Cdti)
                            if (Cdti(kk) < 0.5);
                                err1 = 0;
                                break;
                            end
                        end
                        
                        if (err1==0)
                            isp = [(kk-1):kk];
                            %dindex = [1:numel(ddam)];
                            fitres.pout(3) = interp1(Cdti(isp),dindex(isp),0.5)/log(2);
                        end
                        
                        %{
                        figure(100);
                        clf;
                        plot(dindex,Cdti,'ob');
                        hold on;
                        plot(dindex,feval(fitfunc_str,dindex,[0 1 fitres.pout(3) 0 ]),'r','LineWidth',3.0);
                        %set(gca,'Xscale','log','Yscale','log');
                        title(['jj = ' num2str(jj) 'ii = ' num2str(ii)]);
                        pause(.01);
                        %}
                    end
                    
                    fit_res.back(ii,jj) = fitres.pout(1);
                    fit_res.contrast_fast(ii,jj) = fitres.pout(2);
                    fit_res.damHW_fast(ii,jj) = fitres.pout(3);
                    fit_res.contrast_slow(ii,jj) = fitres.pout(4);
                    fit_res.damHW_slow(ii,jj) = fitres.pout(5);
                    fit_res.slope(ii,jj) = fitres.pout(6);
                    
                    
                    
                end
            end
            
        end
      
        function [Ibar,filt] = smooth_sg_space(III,maxpd,iii,jjj)
            % Smoothing method for getting average
            % Use 2-D Savitzky-Golay smoothing filter
            % For symmetric di or dj, use even pdi, pdj; next higher odd gives same
            % answer
            % these values should be moved to outer program
            % need to optimize by overplotting smoothed function and data
            %maxpd = 2; % determines maximum degree of smoothing polynomial
            %iii = 2; % half-width of pixel range in columns
            %jjj = 2; % half-width of pixel range in rows
            
            Nrs = size(III,1);
            Ncs = size(III,2);
            Ntb = size(III,3);
            
            
            
            di = [-iii:iii];
            dj = [-jjj:jjj];
            pdi = min(maxpd,2*iii);
            pdj = min(maxpd,2*jjj);
            filt = sgsf_2d(di,dj,pdi,pdj,0);
            
            Ncsc = Ncs - 2*iii;
            Nrsc = Nrs - 2*jjj;
            
            
            
            %III_mean = mean(III,3);
            
            
            for kk = 1:Ntb
                mat_2D = conv2(squeeze(III(:,:,kk)),filt,'valid');
                Ibar(:,:,kk) = mat_2D;
            end
            
            
            
        end
     
        
        function [Ibar,filt] = smooth_sg_space_mean(III,maxpd,iii,jjj)
             % Smoothing method for getting average
            % Use 2-D Savitzky-Golay smoothing filter
            % For symmetric di or dj, use even pdi, pdj; next higher odd gives same
            % answer
            % these values should be moved to outer program
            % need to optimize by overplotting smoothed function and data
            %maxpd = 2; % determines maximum degree of smoothing polynomial
            %iii = 2; % half-width of pixel range in columns
            %jjj = 2; % half-width of pixel range in rows
            
            Nrs = size(III,1);
            Ncs = size(III,2);
            Ntb = size(III,3);
            
           
            
            di = [-iii:iii];
            dj = [-jjj:jjj];
            pdi = min(maxpd,2*iii);
            pdj = min(maxpd,2*jjj);
            filt = sgsf_2d(di,dj,pdi,pdj,0);
            
            Ncsc = Ncs - 2*iii;
            Nrsc = Nrs - 2*jjj;
            
            
            
            III_mean = mean(III,3);
            
            mat_2D = conv2(III_mean,filt,'valid');
            
            Ibar = zeros(Nrsc,Ncsc,Ntb);
    
            for kk = 1:Ntb
                %mat_2D = conv2(squeeze(III(:,:,kk)),filt,'valid');
                Ibar(:,:,kk) = mat_2D;
            end
    
            
            
        end
        
        function [Ibar,filt] = smooth_sg_space_and_time(III,maxpd,iii,jjj)
             % Smoothing method for getting average
            % Use 2-D Savitzky-Golay smoothing filter
            % For symmetric di or dj, use even pdi, pdj; next higher odd gives same
            % answer
            % these values should be moved to outer program
            % need to optimize by overplotting smoothed function and data
            %maxpd = 2; % determines maximum degree of smoothing polynomial
            %iii = 2; % half-width of pixel range in columns/time
            %jjj = 2; % half-width of pixel range in rows/space
            
            Nrs = size(III,1);
            Ncs = size(III,2);
            Ntb = size(III,3);
            
           
            
            di = [-iii:iii];
            dj = [-jjj:jjj];
            pdi = min(maxpd,2*iii);
            pdj = min(maxpd,2*jjj);
            filt = sgsf_2d(di,dj,pdi,pdj,0);
            
            Ncsc = Ncs;
            Nrsc = Nrs - 2*jjj;
            Ntsc = Ntb - 2*iii;
            
            Ibar = zeros(Nrsc,Ncsc,Ntsc);
            
            for kk = 1:Ncs
               III_single_col = squeeze(III(:,kk,:));
               mat_2D = conv2(III_single_col,filt,'valid');
               Ibar(:,kk,:) = mat_2D; 
            end
            
            
            III_mean = mean(III,3);
          
            
            
        end
        
      
      
        
    end
    
end