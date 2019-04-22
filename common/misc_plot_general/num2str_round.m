function Xoutstr = num2str_round(Xin,DIG)
% Xoutstr = num2str_round(Xin,number_of_digits) :
%		Create string from number with decimal precision to number digits specified. 
%		Default decimal precision is to 4 digits; example num2str(round(5.0101234) outputs the string '5.0101';
%			Compare
%				num2str(25.0101234,4) 		outputs '25.01' (specifies total number of digits)
%				num2str_round(25.0101234,4)	outputs '25.0101'

	if nargin<2;DIG = 4;end
	
	% create the sprintf code for number of decimals shown
	Dcode = ['%0.' int2str(DIG) 'f'];
	
	Xoutstr = num2str(Xin,Dcode);	
	
end
