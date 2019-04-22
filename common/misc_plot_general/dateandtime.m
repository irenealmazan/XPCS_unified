function DateandTimeStr = dateandtime
% function DateandTimestr = dateandtime
% outputs string with spaces and [], e.g., ' [29-May-2018  14:20] ' suitable for titles

	D = clock;
	DateandTimeStr = ...
	[' [',datestr(D,1),'  ',datestr(D,15),'] '];
