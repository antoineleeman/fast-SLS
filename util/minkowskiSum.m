function Psum = minkowskiSum(PA,PB)
% computes the minkowski sum of two general polygons
% usage: Psum = minkowskiSum(PA,PB)
%
% arguments: (input)
% PA - a 2-d polyshape, or an nx2 array of x,y coordinates traversing a polygon
% PB - a 2-d polyshape, or an mx2 array of x,y coordinates traversing a polygon
%
% arguments: (output)
% Psum - a 2-d polyshape that represents the resulting Minkowski sum.
%
% Example of use
% Pa = polyshape([-1 2;0 2;-0.75 2.25;-1 3]);
% Pb = polyshape([1 0;2 0;3 1;4 3;2 .3]);
% MSumShape = minkowskiSum(Pa,Pb)
% MSumShape = 
%   polyshape with properties:
% 
%       Vertices: [8×2 double]
%     NumRegions: 1
%       NumHoles: 0
%
% Author: John D'Errico
% Date: 12/1/2022


% check arguments
if nargin ~= 2
  error('Exactly two arguments are required')
end

if ~isa(PA,'polyshape')
  PA = polyshape(PA);
end

if ~isa(PB,'polyshape')
  PB = polyshape(PB);
end

% triangulate each polyshape
Atri = triangulation(PA);
Btri = triangulation(PB);

% how many triangles in each?
NA = size(Atri.ConnectivityList,1);
NB = size(Btri.ConnectivityList,1);

% loop over each pair of triangles in the triangulations
ind12 = dec2base(0:8,3) - '0' + 1;
for iA = 1:NA
  for iB = 1:NB
    % form the Minkowski sum of the two triangles iA and iB from each
    % triangulated polyshape
    tAxy = Atri.Points(Atri.ConnectivityList(iA,:),:);
    tBxy = Btri.Points(Btri.ConnectivityList(iB,:),:);
    
    % combinations of all possible vertices from the two triangles.
    Hxy = tAxy(ind12(:,1),:) + tBxy(ind12(:,2),:);
    
    % form the convex hull to resolve the resulting hexagon
    Hexhull = convhull(Hxy);
    Hexhull(end) = [];
    nhull = numel(Hexhull);
    
    if nhull > 3
      % but there may have been collinear points. This will happen if two
      % edges of the original triangles were parallel. I want to remove those
      % collinear points to get rid of a warning from polyshape.
      tolxy = norm(eps(max(Hxy,[],1) - min(Hxy,[],1))*100);

      % Euclidean distance between pairs
      D = @(xy1,xy2) sqrt(sum((xy1 - xy2).^2,2));

      % test for colllinearity using the triangle inequality
      collflag = -D(Hxy(circshift(Hexhull,-1),:),Hxy(circshift(Hexhull,1),:)) + ...
                     D(Hxy(Hexhull,:),Hxy(circshift(Hexhull,1),:)) + ...
                     D(Hxy(Hexhull,:),Hxy(circshift(Hexhull,-1),:));

      Hexhull(collflag <= tolxy) = []; % No need to test abs value. negative numbers are also bad
    end

    HexShape = polyshape(Hxy(Hexhull,:));
    
    if (iA == 1) && (iB == 1)
      Psum = HexShape;
    else
      Psum = union(Psum,HexShape);
    end
  end
end


