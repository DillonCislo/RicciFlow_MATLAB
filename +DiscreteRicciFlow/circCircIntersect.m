function [xout, yout] = circCircIntersect(x1, y1, r1, x2, y2, r2)
%CIRCCIRCINTERSECT This function calculates the point of intersection (if
%any) between two sets of circles.  Each circle is defined by its center in
%Cartesian (x,y)-coordinates and its radius.  Usually, two points per pair
%are returned.  NaNs are returned instead when the circles do not 
%intersect or are identical.  Two identical points are returned when the
%two circles are tangent
%
%   INPUT ARGUMENTS:
%
%       - x1:       Nx1 list of the x-coordinates of the centers of the
%                   first set of circles
%
%       - y1:       Nx1 list of the y-coordinates of the centers of the
%                   first set of circles
%
%       - r1:       Nx1 list of the radii of the first set of circles
%
%       - x2:       Nx1 list of the x-coordinates of the centers of the
%                   second set of circles
%
%       - y2:       Nx1 list of the y-coordinates of the centers of the
%                   second set of circles
%
%       - r2:       Nx1 list of the radii of the second set of circles
%
%   OUTPUT ARGUMENTS:
%
%       - xout:     Nx1 list of the x-coordinates of the intersection
%                   points between the two circle sets
%
%       - yout:     Nx1 list of the y-coordinates of the intersection
%                   points between the two circle sets
%
% by Dillon Cislo 01/22/2020

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------

% Validate first x-coordinate list
validateattributes(x1, {'numeric'}, ...
    {'vector', 'finite', 'nonnan', 'real'});

% The number of sets of circles for which to calculate the intersections
numCirc = numel(x1);

% Validate remaining coordinate lists
validateattributes(y1, {'numeric'}, ...
    {'vector', 'finite', 'nonnan', 'real', 'numel', numCirc});
validateattributes(x2, {'numeric'}, ...
    {'vector', 'finite', 'nonnan', 'real', 'numel', numCirc});
validateattributes(y2, {'numeric'}, ...
    {'vector', 'finite', 'nonnan', 'real', 'numel', numCirc});

% Validate radii lists
validateattributes(r1, {'numeric'}, ...
    {'vector', 'finite', 'nonnan', 'real', 'positive', 'numel', numCirc});
validateattributes(r2, {'numeric'}, ...
    {'vector', 'finite', 'nonnan', 'real', 'positive', 'numel', numCirc});

% Make sure that inputs are column vectors
if (size(x1,2) ~= 1), x1 = x1.'; end
if (size(y1,2) ~= 1), y1 = y1.'; end
if (size(r1,2) ~= 1), r1 = r1.'; end
if (size(x2,2) ~= 1), x2 = x2.'; end
if (size(y2,2) ~= 1), y2 = y2.'; end
if (size(r2,2) ~= 1), r2 = r2.'; end

%--------------------------------------------------------------------------
% Calculate Circle Intersection Points
%--------------------------------------------------------------------------

% The Euclidean distance between circle centers
d = sqrt( (x2-x1).^2 + (y2-y1).^2 );

% Determine if any bad intersection types exist ---------------------------

indx1 = d > (r1+r2); % Circles too far apart to intersect
indx2 = r2 > (d+r1); % Circle one is completely inside circle two
indx3 = r1 > (d+r2); % Circle two is completely inside circle one

% Circles are identical
indx4 = (d < (10*eps)) & (abs(r1-r2) < (10*eps));

% Combine indices - 
indx = indx1 | indx2 | indx3 | indx4;

% Calculate basic intersection points -------------------------------------

% The angle that the separation vector between the two circles makes with
% the x-axis
a0 = atan2( (y2-y1), (x2-x1) );

% The angle opposite r2 in the triangle formed by the radii vectors and the
% separation vector
a1 = acos( (d.^2 + r1.^2 - r2.^2) ./ (2 .* r1 .* d ) );

% Linear combinations of the previous angles
alpha1 = a0 + a1;
alpha2 = a0 - a1;

% The coordinates of the circle intersection points
xout = [x1 x1] + [r1 r1] .* cos([alpha1 alpha2]);
yout = [y1 y1] + [r1 r1] .* sin([alpha1 alpha2]);

% Replace complex results (bad intersections) with NaNs -------------------
xout(indx, :) = NaN;
yout(indx, :) = NaN;   

end

