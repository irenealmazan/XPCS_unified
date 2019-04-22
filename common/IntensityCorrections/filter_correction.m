function ABSCORR = filter_correction(FILTER,filter1,filtercorr)
% ABSCORR = filter_correction(filterNUMBER,value_of_filter1,[filtercorrections])
% XIA filter number	1	2	4	8; absorbers (nominal) 1xfoil	2xfoils	4xfoils	8xfoils
%
% OUTPUT (plain)
%	ABSCORR (vector of same length as input filterNUMBER). It will be > 1
%			I(no filters) = I(with filterNUMBER in) .* ABSCORR
%
% INPUT (plain)
%	filterNUMBER (may be vector) Integers (values allowed 0-15)
%	value_of_filter1 (will be > 1) absorption strength of slot 1 filter; I(without)/I(with)
%	filtercorrections (optional)  1x3 vector, default is [0 0 0]; 
%						apply to effective thickness of filters 2, 4, 8 due to thickness
%						Look at file for how to determine the filter corrections if needed.


% XIA slot	 has 4 slots, enumerate as 		"1"	 "2" "4" "8"
%		combinations then are 0 to 15
% filtercorr is vector to correct the absorption calculations of slots 2, 3, 4 with respect to ideal
% It corrects for thickness errors or inhomogeneities in the filters with respect to 'ideal'
%
% Note filtercorr = [filter2corr; filter4corr; filter8corr]
%		filtercorr are corrections added to the exponents due to problems with uneven absorbers.
% 			Given I0, and I2 (measured with filter 2 in beam)
%		e.g., filter2corr :   log(I2/I0) = [2 + filter2corr]*log(filter1)
%			Given I2 (with filter 2 in beam), and I4 (with filter 4 in beam)
%		e.g., filter4corr:	log(I4/I2) = [4 + filter2corr +filter4corr]*log(filter1)
%		e.g., filter8corr:	log(I8/I4) = [8 + filter2corr + filter4corr]*log(filter1)
%			it goes without saying to do the absorber measurments carefully!! 

if nargin<2, 
	filter1 = 2.239; 
end
if nargin<3;
	filtercorr = [0 0 0]

end

filter2corr = filtercorr(1);
filter4corr = filtercorr(2);
filter8corr = filtercorr(3);

% these are used
% [1 0 0 0], [0 1 0 0]  (#1 and #2) can use filter1^n;
% for [0 0 1 0] in the filter mix, instead of multiplying by another 
% filter, use filter1^(1+filter4corr)
% Likewise, for [0 0 0 1] (#8)
% use filter1^(1+filter8corr);


% I'm sure this could be done more elegantly

FILTERNUM = FILTER+1;
pow = ...
[0
 1 
 2  + filter2corr
 3  + filter2corr
 4  + filter4corr
 5  + filter4corr
 6  + filter4corr + filter2corr
 7  + filter4corr + filter2corr
 8  + filter8corr
 9  + filter8corr
 10 + filter8corr + filter2corr
 11 + filter8corr + filter2corr
 12 + filter8corr + filter4corr
 13 + filter8corr + filter4corr
 14 + filter8corr + filter4corr + filter2corr
 15 + filter8corr + filter4corr + filter2corr];

ABSCORR =  filter1.^pow(FILTERNUM);
