function h = computePowerDists(L, r)
%COMPUTEPOWERDISTS Given a set of edge lengths of a triangulation and a
%set of circle radii for each vertex of that triangulation, this function
%performs a completely vectorized calculation of the perpendicular distance
%from each triangulation edge to the power center of the unique circle
%orthogonal to the three circles at each vertex of a face. Intended for
%internal use with the 'DiscreteRicciFlow' package.
%
%   INPUT PARAMETERS:
%
%       - L:    #Fx3 list of edge lengths.  L(f,i) is the length of the
%               edge opposite vertex i in face f.
%       - r:    #Fx3 list of circle packing radii.  r(f,i) is the radius of
%               the circle centered on vertex i in face f.
%
%   OUTPUT PARAMETERS:
%
%       - h:    #Fx3 list of power distances.  h(f,i) is the perpendicular
%               distance of the edge opposite vertex i in face f to the
%               power center of face f
%
% by Dillon Cislo 11/17/2019

%==========================================================================
% THE GEOMETRY OF THE CALCULATION
%
% To simplify the calculation we endow each face with its own set of
% abstract coordinates.  For a face defined as f = [ vi vj vk ]
%
%               vk=(xk,yk)
%                   /\
%                  /  \
%                 /    \
%             Lj /      \ Li
%               /        \
%              /          \
%   vi=(0,0)  o------------o  vj=(Lk,0)
%                   Lk
%
%
% The fact that the power circle is simultaneously orthogonal to all three
% vertex centered circles enforces the following equidistance condition:
%
%           di^2 - ri^2 = dj^2 - rj^2 = dk^2 -rk^2 = R^2
%
% where R is the radius of the power circle and di is the distance between
% the center of the power circle and the vertex i.  The formulae for the
% center of the power circle (x,y) can be derived by directly substituting
% the abstract coordinates above into the equidistance condition.
%==========================================================================

% Cache variables to eliminate redundant calculations
Li = L(:,1); Lj = L(:,2); Lk = L(:,3);
Li2 = Li.^2; Lj2 = Lj.^2; Lk2 = Lk.^2;

ri = r(:,1); rj = r(:,2); rk = r(:,3);
ri2 = ri.^2; rj2 = rj.^2; rk2 = rk.^2;

% Calculate (xk,yk) for each face
xk = ( Lj2 + Lk2 - Li2 ) ./ ( 2 .* Lk ); xk2 = xk.^2;
yk2 = Lj2 - xk2; yk = sqrt( yk2 );

% Compute the center of the power circle (x,y)
x = ( Lk2 + ri2 - rj2 ) ./ ( 2 .* Lk );
y = ( xk2 + yk2 - 2 .* xk .* x + ri2 - rk2 ) ./ ( 2 .* yk );

% Compute the distances bewteen the vertices and the power center
di = sqrt( x.^2 + y.^2 );
dj = sqrt( (Lk-x).^2 + y.^2 );
dk = sqrt( (xk-x).^2 + (yk-y).^2 );

% Calculate the areas of each of the triangle subdivisions
si = ( Li + dj + dk ) ./ 2; % Semi-perimieters for Heron's formula
sj = ( Lj + dk + di ) ./ 2;
sk = ( Lk + di + dj ) ./ 2;

Ai = sqrt( si .* (si-Li) .* (si-dj) .* (si-dk) );
Aj = sqrt( sj .* (sj-Lj) .* (sj-dk) .* (sj-di) );
Ak = sqrt( sk .* (sk-Lk) .* (sk-di) .* (sk-dj) );

% The power distances can be calculated as the heights from the triangle
% area formula
hi = 2 .* Ai ./ Li;
hj = 2 .* Aj ./ Lj;
hk = 2 .* Ak ./ Lk;

% If the power center lies on one of the edges, clip output to the real
% part of the calculation
h = real([ hi hj hk ]);

end

