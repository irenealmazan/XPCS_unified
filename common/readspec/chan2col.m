function [CHANCOLUMN,channelchoosematch] = chan2col(chan,channelchoose,quietflag)
% Columnnumber = chan2col(Speclabels,Chooselabel)
% Given spec LABELS string list from scan, ['Phi  H  K  L  mon  det']
% 	chan2col will output column numbers of chosen label
%   If there are two repeated with same label, it will return repeats also
%
% OUTPUT 
%	Columnnumber(s) associated with particular named string in spec file
% INPUT  
%	Speclabels  (the labels line from the spec scan)
% 	 			example, ['TwoTheta  h  k  l  epoch  time  detm']  
%				note 2 spaces between each label in row is **required** 
%				In spec this accomodates names like 'Two Theta' in config
%	Chooselabel (choice of label to match (can find several at once))
%				example, use chooselabel= ['time';'detm'] 
%				will output columns [6 7]
%				or e.g., chooselabel = ['Two Theta'] would be [1] in output
%				strvcat(strA,strB,strC,...) will be useful for 
%				constructing a string matrix with different length strings.
%
% Complementary function is *col2chan* to output the speclabel of a particular column
% 		col2chan(speclabels,[5 6]) would outputs ['time';'detm']
%
% Note that Speclabels can be long space delimited string, or string matrix array
% (however, note array strings are left justified)
%

% 2014-December-13 C. Thompson, fixed so accommodates the scan header better (which can have many more spaces
% 2016-Jan C. Thompson added quietflag, if argument exists, it will be 'quiet' 
%					added namematch to figure out which entries were matched in case 
% 2016-Feb C. Thompson, speclabels (or labels list) can be row with double space or a column 


% Convert chan into a cell array with each row a different string
% This requires that at least two spaces between each differentiated string

% if anything in 3rd argument, it will be 'quiet'
%   not give the message that the motors are not in the labels
if nargin<3;quietflag=0;end

 % replace double spaces with commas and  turn into column of strings
if length(chan(:,1))==1; 
	chan_rephrased = [''''   (strrep(chan,'  ', ''','''))    ''''];
	channames_as_column = strtrim(eval(['strvcat(',chan_rephrased,')']));
else  % or it already is column array
	channames_as_column = chan;
end
    % Then turn the column into a cell array
	numofcols = length(channames_as_column(:,1));   %just the total number of variables in the LABELS at start of scan	
	chan_rephrased2col = [ones(numofcols,1)*[''''], channames_as_column,  ones(numofcols,1)*[''';']];
	chan_rephrased2 = char(reshape(chan_rephrased2col',1,length(chan_rephrased2col(:))));
	channames_as_cellarray = strtrim(eval(['{',chan_rephrased2,'}']));
	
% in case input was not in column but in row, such as ['Test  Two Theta  Gamma']
% put it into a matrix and continue
% Undocumented feature (requires input with 2 spaces or more between  differeniated entries
%	and for only one space, will assume it is a ill-posed name like 'Two Theta'
if length(channelchoose(:,1))==1;
	chanchoose_rephrased = [''''   (strrep(channelchoose,'  ', ''','''))    ''''];
	chanchoosenames_as_column = eval(['strvcat(',chanchoose_rephrased,')']);
	channelchoose = chanchoosenames_as_column;
end

CHANCOLUMN = [];
channelchoosematch = [];jj=1;
for ii = 1:length(channelchoose(:,1))
	THE_MATCH 	= find( strcmp(strtrim(channelchoose(ii,:)),channames_as_cellarray));
		if ~isempty(THE_MATCH)
			%THE_MATCH = THE_MATCH(1);  %%% here is where it kicks out multiples if desired
			channelchoosematch(jj) = ii;
			jj=jj+1;
		elseif quietflag ==0;  % if match is empty tell us unless quiet is on (not zero)
			disp([''''  strtrim(channelchoose(ii,:))  ''''  ' choice is not parsed within chan: ',chan]);
		end
	CHANCOLUMN = [CHANCOLUMN;THE_MATCH];
end

