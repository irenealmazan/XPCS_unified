% This scripts display the set of intensities taken during a time scan,
% their integrated intensity in the different ROIS



Qvalflag = 0;
[Allscans(iT).HS_struct] = DisplayTTM_plot.display_II(Read_Allscans(iT).IIstruct,log10(Read_Allscans(iT).IIstruct.II(:,:,1000)),'Log and frame = 2500',Qvalflag,13,ImageJ,SINGLE,XCOLlabel,YROWlabel,AXISdet,INFOstr,sim_flag,[-.5  2]);%,CLIM);


[Allscans(iT).HS_struct] = DisplayTTM_plot.display_II_with_rois(Read_Allscans(iT),Read_Allscans(iT).IIstruct.II(:,:,Allscans(iT).itt_range_struct.ittt),'Average in Log 10',iT+13,ImageJ,SINGLE,XCOLlabel,YROWlabel,AXISdet,INFOstr,CLIM);

[Allscans(iT).HS_struct] = DisplayTTM_plot.show_rois_only(Read_Allscans(iT),Allscans(iT),iT+14,XCOLlabel,YROWlabel,AXISdet);

SoXFLAG = 1;SoYFLAG = 0;LOGFLAG = 1;
[Allscans(iT).HS_sumimag_struct] = DisplayTTM_plot.make_summed_images(Read_Allscans(iT),Allscans(iT).ROIS_struct,[tminv(iT) tmaxv(iT)],SoXFLAG,SoYFLAG,LOGFLAG,iT+15,ImageJ,XCOLlabel,YROWlabel,AXISdet);
[Allscans(iT).HS_sumimag_struct] = DisplayTTM_plot.display_sumpoints_1D(Read_Allscans(iT),Allscans(iT),SoXFLAG,SoYFLAG,iT+16,ImageJ,XCOLlabel,YROWlabel,AXISdet);

SoXFLAG = 0;SoYFLAG = 1;LOGFLAG = 1;
[Allscans(iT).HS_sumimag_struct] = DisplayTTM_plot.make_summed_images(Read_Allscans(iT),Allscans(iT).ROIS_struct,[tminv(iT) tmaxv(iT)],SoXFLAG,SoYFLAG,LOGFLAG,iT+17,ImageJ,XCOLlabel,YROWlabel,AXISdet);
[Allscans(iT).HS_sumimag_struct] = DisplayTTM_plot.display_sumpoints_1D(Read_Allscans(iT),Allscans(iT),SoXFLAG,SoYFLAG,iT+18,ImageJ,XCOLlabel,YROWlabel,AXISdet);

DisplayTTM_plot.plot_summed_images(Read_Allscans(iT),Allscans(iT).ROIS_struct,log10(Read_Allscans(iT).IIstruct.II),iT+19,ImageJ,CLIM,XCOLlabel,YROWlabel,AXISdet,INFOstr,POINTSUMS);

ROIS_to_plot = [1 2 3 4];
[Allscans(iT).ROIS_struct.sum] = DisplayTTM_plot.plot_ROIs_as_counters(Read_Allscans(iT),Allscans(iT),ROIS_to_plot,iT+27) ;

% display the FT of the ROIs integrated intensity:
%DisplayTTM_plot.calc_fft_sumROIs(Allscans,Read_Allscans,[1:500],4,2) ;