function [ImageStack,imtime,imdocu] = load_MPX3(fullimname,MASK,badpixfile)
% [ImageStack,imtime,imdocu] = load_MPX3(fullimname,MASK,badpixfile)   (assume 'color')
%	output
%		data (2D or 3D array): 
%				may be 3D if multiple imnames 
%				
%		imtime (scalar or vector): timestamp (in seconds) of when image acquired
%				may be vector if multiple imnames
%				it is present if the tif header has a time stamp, empty if it does not
%		
%		imdocu (string array or cell): information pulled out from this file to go with the ImageStack
%				(at present, only if applying badpix or not
%		
%	input
%		fullimname (string or matrix string) filename (needs paths)
%				assumption is that tif extension is '.tif' (and not '.TIF)
%		MASK (only use a masked region - not implemented yet, also would need X Y output change
%		badpixfile  (if present, use to replace bad pix, uses the EPICS formulation, also needs path if elsewhere
%			[X Y] [Xreplace Yreplace]
%		
%		


NPTS = length(fullimname(:,1));
[Yrow,Xcol] = size(imread(fullimname(1,:)));
imtime = [];


ImageStack = zeros([Yrow Xcol NPTS]);

% The fix below (with warning('off',...) will now only give the error message for the first iminfo call
% I preferred to have some error messages come out to check on the pilatus tifs since I will forget if I don't see them
	disp('If you are reading Pilatus tif files, their info is non-standard and Matlab will throw a message')
	disp('when  we did an iminfo command. Ignore error messages that look like')
	disp('... [Warning: zero count for tag] in matlab.io, etc]');

%%%%%%%  read in the images
ImInfo = imfinfo(fullimname(1,:));    %%% note recently the pilatus tiffs are not standard, this imfinfo will throw an error as long as we turn errors back on
for jj=1:NPTS;
	ImageStackii = imread(fullimname(jj,:),'tif');
	
	% handle that tifs from pilatus have issue as noted below; turn turn off for the multiple file reads
	%	note - using     warning('on','verbose');   helped in figuring out what to turn off
	
	warning('off','MATLAB:imagesci:tifftagsread:zeroComponentCount');

	ImInfo = imfinfo(fullimname(1,:));    
		%%% note recently the pilatus tiffs are not standard so this imfinfo works but also spits out error message
	imtime = [imtime;ImInfo.FileModDate];	% assume all same string length
	ImageStack(:,:,jj) = ImageStackii;
	
	warning('on','MATLAB:imagesci:tifftagsread:zeroComponentCount');  
end	


if nargin ==3 & ~isempty(badpixfile)
	[ImageStack,BadPixDocument] = apply_badpixreplace(ImageStack,badpixfile,[1 0]);   % no plot, will replace bad pix with NaN
	imdocu = BadPixDocument;
else
	imdocu = [];
end

	% Note apply_badpixreplace works for 3D stack, but it does not have mask option yet


% imfinfo of MPX3 tif file acquired via epics external trigger 2017-04 gives following informations

%                     Filename: '/home/cthompson/.gitrepose/DATA/2017/2017_04_mocvd_nitrideXPC…'
%                  FileModDate: '08-Apr-2017 19:19:23'       << this may change if file gets 'resaved' with different date
%                     FileSize: 532785
%                       Format: 'tif'
%                FormatVersion: []
%                        Width: 516
%                       Height: 516
%                     BitDepth: 16
%                    ColorType: 'grayscale'
%              FormatSignature: [73 73 42 0]
%                    ByteOrder: 'little-endian'
%               NewSubFileType: 0
%                BitsPerSample: 16
%                  Compression: 'Uncompressed'
%    PhotometricInterpretation: 'BlackIsZero'
%                 StripOffsets: 8
%              SamplesPerPixel: 1
%                 RowsPerStrip: 516
%              StripByteCounts: 532512
%                  XResolution: []
%                  YResolution: []
%               ResolutionUnit: 'Inch'
%                     Colormap: []
%          PlanarConfiguration: 'Chunky'
%                    TileWidth: []
%                   TileLength: []
%                  TileOffsets: []
%               TileByteCounts: []
%                  Orientation: 1
%                    FillOrder: 1
%             GrayResponseUnit: 0.0100
%               MaxSampleValue: 65535
%               MinSampleValue: 0
%                 Thresholding: 1
%                       Offset: 532520
%                         Make: 'Unknown'
%                        Model: 'Unknown'
%                     Software: 'EPICS areaDetector'
%                 SampleFormat: 'Two's complement signed integer'
%                  UnknownTags: [4x1 struct]

% UnknownTags: are fields ID, Offset, Value (showing way to see one element of many)
%		values shown are typical for the MPX3 as acquired without any special mode (1 det, saturate 4095)
%			UnknownTags(1).ID  (4 of them, of order 65000, 65001, 65002, 65003
%			UnknownTags(1).Offset (4 of them, of order 532742, 0, 0, 0)
%			UnkonwnTags(1).Value  (4 of them, 8.6053e+08, 9, 0, 0
