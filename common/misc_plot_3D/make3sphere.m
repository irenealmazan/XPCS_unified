function [HL] = make3sphere(matrixofpoints,radius,COLOR,TRANSPARENCY)
% [HL] = make3sphere(XYZmatrix,radius,COLOR,TRANSPARENCY)
%  In the currently active axes:
%  Given a point xyz, draws a little sphere centered there
%		May control the radius, color and transparency
%		Try >> Hlight = light   to get some lighting going
%		To change all shinyness, try >> material [metal, shiny, dull]
%
% OUTPUT (structure)
%	[HL]  structure, HL.H is the handle for each sphere 
%			thus HL.H(1) gives properties of 1st one plotted, etc 
%
% INPUT  (plain)
%	XYZmatrix is 3xN matrix of positions of sphere centers, 
%			where N is number of spheres (each column is an [x;y;z] pos
%	radius   radius in xyz units, default is 0.001
%	COLOR	 is [R G B] value or 'r', 'g', etc
%	TRANSPARENCY  is alpha, default is opaque (1) (total transparent is 0)
%
% Note - the default also gives these spheres reflective surfaces so
%	if lighting is on, they will shine!  >> Hlight = light;
%	to indiscriminately change all surfaces shinyness, 
%			try >> material [metal, shiny, dull]
%	
        if nargin<2;radius = 0.1;end
        if nargin<3;COLOR = [0.8 .8 .8];end
	if nargin<4;TRANSPARENCY= [1];end
        
        %[Xs, Ys, Zs] = sphere(20);  %makes sphere, radius 1 at origin)
        %R0 radius
        %CEN = [X0 Y0 Z0];
        %X = R0.*Xs + X0; 
        %Y = R0.*Ys + Y0;
        %Z = R0.*Zs + Z0;
		
		X0	=	matrixofpoints(1,:);
		Y0	=	matrixofpoints(2,:);
		Z0	=	matrixofpoints(3,:);
        
        Ns = length(X0);
        N=15+1;
        
        [Xs,Ys,Zs]=sphere(N-1); %(default is already N=15) sphere(20) OK 
        Xs = Xs.*radius;
        Ys  = Ys.*radius;
        Zs = Zs.*radius;
       
        
        Xtot = reshape( ((ones(N.*N,1) * X0) + (Xs(:)*ones(1,Ns))),N,N,Ns);
        Ytot = reshape( ((ones(N.*N,1) * Y0) + (Ys(:)*ones(1,Ns))),N,N,Ns);
        Ztot = reshape( ((ones(N.*N,1) * Z0) + (Zs(:)*ones(1,Ns))),N,N,Ns);
        
        props = sphereproperties(ones(size(Xtot(:,:,1))),COLOR,TRANSPARENCY);
		clear HL;
       
        for ii=1:Ns
            
            HL.H(ii) = surface(Xtot(:,:,ii),Ytot(:,:,ii),Ztot(:,:,ii),props);
        end;
end

function [props]=sphereproperties(CDATA,COLOR,TRANSPARENCY)
	%for metallic somewhat transparent ball of CDATA colors that is lit up
	props.AmbientStrength = 0.3;
	props.DiffuseStrength = 0.6;
	props.SpecularStrength = 1;
	props.SpecularExponent = 25;
	props.SpecularColorReflectance = 0.8;

	props.FaceColor= 'flat';%'texture';%'texture';%'flat';  %'texture'
	props.EdgeColor = 'none';%COLOR;%'interp';
%    	props.Cdatamapping = 'direct';
        props.FaceLighting = 'phong';
	props.Cdata = ones(size(CDATA)).*.2; %ones(size(CDATA)).*0.01;	

	props.FaceAlpha = TRANSPARENCY;%0.2;% 0.3; %  'texture';
	props.EdgeAlpha = TRANSPARENCY;%'flat';%
	%props.Alphadatamapping = 'none';% default is 'scaled'
	%props.Alphadata = ones(size(CDATA)).*TRANSPARENCY;
	%FaceAlpha: {'flat'  'interp'  'texturemap'}
	%EdgeAlpha: {'flat'  'interp'  'texturemap'}


	%for metallic somewhat transparent ball of CDATA colors that is lit up
	%props.AmbientStrength = 0.3;
	%props.DiffuseStrength = 0.3;
	%props.SpecularStrength = 1;
	%props.SpecularExponent = 25;
	%props.SpecularColorReflectance = 0.8;
	
	%props.FaceColor= 'texture';  %'flat'
	%props.EdgeColor = 'none';%'interp';
	%props.FaceLighting = 'phong';
	%props.Cdata = ones(size(CDATA)).*0.5;
	
	%props.FaceAlpha = 0.05; %  'texture';
	%props.EdgeAlpha = 0.05;%'none';%'flat' need alphadata same size as zdata
	%%props.Alphadatamapping = 'none';% default is 'scaled'
	%%props.Alphadata = ones(size(DETECTOR)).*.8;

end

%	AlphaData
%	AlphaDataMapping: [ none | direct | {scaled} ]
%	EdgeAlpha: [ flat | interp ] -or- {an Alpha}.
%	FaceAlpha: [ flat | interp | texturemap ] -or- {an Alpha}.
%	CData
%	CDataMapping: [ direct | {scaled} ]
%	EdgeColor: [ none | flat | interp ] -or- {a ColorSpec}.
%       FaceColor: [ none | {flat} | interp | texturemap ] -or- a ColorSpec.
%	DisplayName
%	LineStyle: [ {-} | -- | : | -. | none ]
%	LineWidth
%	Marker: [ + | o | * | . | x | square | diamond | v | ^ | > | < | pentagram | hexagram | {none} ]
%	MarkerEdgeColor: [ none | {auto} | flat ] -or- a ColorSpec.
%	MarkerFaceColor: [ {none} | auto | flat ] -or- a ColorSpec.
%	MarkerSize
%	MeshStyle: [ {both} | row | column ]
%	XData
%	YData
%	ZData
%	FaceLighting: [ none | {flat} | gouraud | phong ]
%	EdgeLighting: [ {none} | flat | gouraud | phong ]
%	BackFaceLighting: [ unlit | lit | {reverselit} ]
%	AmbientStrength
%	DiffuseStrength
%	SpecularStrength
%	SpecularExponent
%	SpecularColorReflectance
%	VertexNormals
%	NormalMode: [ {auto} | manual ]
% 'AmbientStrength'    'DiffuseStrength   'SpecularStrength'   'SpecularExponent'   'SpecularColorReflectance'
%    'Shiny',	0.3, 	0.6, 	0.9,	20,		1.0
 %   'Dull',		0.3,	0.8,	0.0,	10,		1.0
  %  'Metal',	0.3,	0.3,	1.0,	25,		.5
