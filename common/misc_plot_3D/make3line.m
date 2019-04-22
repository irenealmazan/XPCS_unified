function [HL] = make3line(matrixofpoints)
% [HL] = make3line(matrixofpoints) adds a line in 3D that connects points in the currently active set of axes
% OUTPUT  (structure) handle of the line
% INPUT  (plain) 3xN matrix where each column is a point [x;y;z]
		
		X	=	matrixofpoints(1,:);
		Y	=	matrixofpoints(2,:);
		Z	=	matrixofpoints(3,:);
		
		HL.H = line(X,Y,Z);
end

