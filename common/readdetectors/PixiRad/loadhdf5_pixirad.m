function [ data,imtime,thresh] = loadhdf5_CT(imname,thresh_order)
% [data,imtime,thresh] = loadhdf5(imname,thresh_order)   (assume 'color')
%	output
%		data (2D or 3D array): I(low threshold) - I(high threshold) 
%				may be 3D if multiple imnames 
%				
%		imtime (scalar or vector): timestamp (in seconds) of when image acquired
%				may be vector if multiple imnames
%		thresh (structure of arrays): raw 3D arrays with Ilow and Ihigh raw
%			thresh.low and thresh.high
%		
%	input
%		imname (string or matrix string) filename (needs paths)
%		thresh_order [LOW HIGH] e.g., [2 1] which image is set with 
%				LOW and HIGH thresholds in 'color' or 'window' collection mode
%			(assume I(low)-I(high) is good
%			
%		
%		
% Wonsuk Cha, SRS, MSD, ANL
% Dec. 18, 2015
% This script loads multiple three-dimensional hdf5 type files.
% imname: name of a file, if you leave blank, then this load all hdf5 files in the directory.
%
% edit 2016 July   Carol Thompson
% CT turned off hidden XY flipping, we prefer manipulation of xy
% 	or rowcolumns done higher up in analyzing stream
% 	Note - in 'color' collection mode of pilatus, this program assumes
%		threshold 2 is used for the high energy (image has good + bad)
%		threshold 1 is used the the low energy (image has bad)
%		It assumes that within file, (:,:,2)-(:,:,1) is I2-I2=Igood
%
%
% Timestamp (in seconds) 
% 		h5read(imname,'/entry/instrument/NDAttributes/NDArrayTimeStamp')
%

if nargin<2; thresh_order = [1 2]; end  % [low threshold     high threshold]

%
%% Load data
if nargin==1;
    fileinfo=dir('*.h5');% This creates information of folders.
    for n=1:length(fileinfo)
        imname=fileinfo(n).name;
        temp=h5read(imname,'/entry/data/data');
        imtimetemp = h5read(imname,'/entry/instrument/NDAttributes/NDArrayTimeStamp');
        thresh.low(:,:,n)=temp(:,:,thresh_order(1));
        thresh.high(:,:,n)=temp(:,:,thresh_order(2));
        imtime(n) = imtimetemp;
    end
elseif nargin==2
	for n=1:length(imname(:,1))
		imnameii = deblank(imname(n,:));
        temp=h5read(imnameii,'/entry/data/data');
        imtimetemp = h5read(imnameii,'/entry/instrument/NDAttributes/NDArrayTimeStamp');
        thresh.low(:,:,n)=temp(:,:,thresh_order(1));
        thresh.high(:,:,n)=temp(:,:,thresh_order(2));
        imtime(n) = imtimetemp;
    end
end

thresh.low=double(thresh.low);
thresh.high=double(thresh.high);

% Rows = Y, Cols = X, 
% CT prefers do any rearranging at a higher level since it changes per expt
%thresh.low=fliplr(thresh.low);
%thresh.high=fliplr(thresh.high);

data=thresh.low-thresh.high;
data=double(data);

disp(['The dataset is loaded.']);
disp(['The size is ' num2str(size(data,2)) ' X(cols) ' num2str(size(data,1)) ' Y(rows) ' num2str(size(data,3)) ' colors '])
end
