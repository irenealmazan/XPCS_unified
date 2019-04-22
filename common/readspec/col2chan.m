function [colname,numofcols,channames_as_column] = col2chan(chan,colchoose)
% [Columnlabel,Numofcols,AllLabels_matrix] = col2chan(Speclabels,Columnnum)
%  Given Spec label string, it will output the Label name assoicated with column
%
% OUTPUT 
% 	Columnlabel (string or string matrix of Labels)
%	Numofcols   (total number of labels in the string)
%	AllLabels_matrix  (reshaped the labels from spec file into a matrix)
%						(this might be useful to user to work with otherwise)
%		
% INPUT  
%	Speclabels  (the labels line row from the spec scan)  
% 	 			example, ['TwoTheta  h  k  l  epoch  time  detm']  
%				note 2 spaces between each label in row is **required**
%				this is done in spec to accomodates double names like 'Two Theta' in config
%	Columnnum   (vector of column numbers whose label are desired)
%				example, [6 7]
%				will output ['time';'detm']
%
% Complementary function is *chan2col* to output the speclabel of a particular column
% 		chan2col(speclabels,Chooselabel) 
%

chan_rephrased = [''''   (strrep(chan,'  ', ''','''))    ''''];
	channames_as_column = eval(['strvcat(',chan_rephrased,')']);
	numofcols = length(channames_as_column(:,1));   %just the total number of variables in the LABELS at start of scan

if any(colchoose> numofcols);
	colchoose = reshape(colchoose,1,length(colchoose));
	disp(['Requested column [',int2str(colchoose(find(colchoose>numofcols))), '] is not available  -  the total number of labels is ',int2str(numofcols)]);
	colname = [];
else
	colname = strtrim(channames_as_column(colchoose,:));
end
