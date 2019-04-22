function longernamematrix = addnames2matrix(prestring,StringMatrix,poststring)
% longmatrix = addnames2matrix(prestring,StringMatrix,poststring)
% 	Given a string matrix, takes the prename and the postname  and adds to ends
%	Useful with filename matrix, add path or extenstion. 
%	
%	OUTPUT (plain matrix)

%		longermatrix   with the added paths and extensions
%	INPUT (plain or struct)
%		prestring	(plain) the string to put in front of all the rows
%		poststring  (plain) the string to put at end of all the rows
%					be sure to include any fileseparators etc
%		StingMatrix (plain or structure)
%					if plain, it is a string matrix, e.g., rows of filenames
%					if struct, use input as follows (as 2 element cell) 
%						{MyStruct,'filenames'} will tell MyStruct.filenames
%						or, MyStruct.filename will also just work


if nargin<1;help addnames2matrix;return;end
if nargin<3;poststring=[];end

PREm=[];POSTm=[];

% if just a plain matrix
if ~isstruct(StringMatrix) & ~iscell(StringMatrix)

	thematrix = StringMatrix;
	
else
	
	thematrix = getfield(StringMatrix{1},StringMatrix{2});
	
end

	Npts = length(thematrix(:,1));
	
	if ~isempty(prestring);
		PREm 	= char(ones(Npts,1)*[prestring]);
	end
	if ~isempty(poststring);
		POSTm	= char(ones(Npts,1)*[poststring]);
	end

longernamematrix=[PREm thematrix POSTm];
