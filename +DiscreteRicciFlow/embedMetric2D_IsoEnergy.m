function V2D = embedMetric2D_IsoEnergy(L, F, fID, E, feIDx)
%EMBEDMETRICIN2D_ISOENERGY
%Given a flat discrete Riemannian metric on a triangulation
%(a set of positive edge lengths obeying the triangule inequality in each
%face) this function provides a parameterization of that matric as a set of
%2D vertex coordinates. Embedding is done by minimizing the norm of a
%matrix operator derived from the face-based gradient operator acting on
%scalar vertex signals that must vanish for any valid 2D embedding.
%
%WARNING: Face connectivity list should be consistently oriented
%
%Mainly intended for internal use with the 'DiscreteRicciFlow' package
% 
%   INPUT PARAMETERS:
%
%       - L:            #Ex1 list of triangulation edge lengths
%
%       - F:            #Fx3 face connectivity list
%
%       - fID:          The face ID of the seed face whos position is
%                       constrained in the embedding
%
%       - E:            #Ex2 list of vertex IDs defining edges
%
%       - feIDx:        #Fx3 face-edge correspondence tool.  feIDx(f,i) is
%                       the ID of the edge opposite vertex i in face f
%
%   OUTPUT PARAMETERS:
%
%       - V2D:          #Vx2 2D vertex coordinates
%
% by Dillon Cislo 11/20/2019

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------

% Ensure mandatory input parameters are supplied
if (nargin < 1), error('Please supply target edge lengths'); end
if (nargin < 2), error('Please supply face connectivity list'); end

% Verify the properties of the input triangulation
validateattributes( F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive'} );
validateattributes( L, {'numeric'}, ...
    {'2d', 'ncols', 1, 'finite', 'nonnan', 'positive'} );

numFaces = size(F,1); % The number of faces in the triangulation
numVertex = max(F(:)); % The number of vertices in the triangulation

% Construct a MATLAB-representation of the input triangulation
TR = triangulation( F, zeros( max(F(:)), 2 ) );

% Check for optional input parameters
if (nargin < 3), fID = 1; end
if (nargin < 4), E = sort( TR.edges, 2 ); end

if (nargin < 5)
    
    e1IDx = sort( [ F(:,3), F(:,2) ], 2 );
    e2IDx = sort( [ F(:,1), F(:,3) ], 2 );
    e3IDx = sort( [ F(:,2), F(:,1) ], 2 );
    
    [~, e1IDx] = ismember( e1IDx, E, 'rows' );
    [~, e2IDx] = ismember( e2IDx, E, 'rows' );
    [~, e3IDx] = ismember( e3IDx, E, 'rows' );
    
    feIDx = [ e1IDx e2IDx e3IDx ];
    
end

% Verify the properties of the optional input parameters
validateattributes( E, {'numeric'}, ...
    {'2d', 'ncols', 2, 'integer', 'positive'} );
validateattributes( feIDx, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive'} );

% Map edge lengths to faces
L_F = L(feIDx);

%--------------------------------------------------------------------------
% Calculate Face Geometric Quantities
%--------------------------------------------------------------------------

% Calculate face areas ----------------------------------------------------

% Some convenience variables for clarity
% Li(f) is the length of the edge opposite vertex i in face f
Li = L_F(:,1); Lj = L_F(:,2); Lk = L_F(:,3);

% The semi-perimeter of each face
S = sum(L_F, 2) ./ 2;

A = sqrt( S .* (S-Li) .* (S-Lj) .* (S-Lk) );

% Calculate the position of the third vertex in the intrinsic geometry ----
Xk = ( Lj.^2 + Lk.^2 - Li.^2 ) ./ ( 2 .* Lk );
Yk = sqrt( Lj.^2 - Xk.^2 );

% Calculate area-normalized edge-vectors in the plane ---------------------
eXJ = Yk ./ ( 2 .* A );
eYJ = -Xk ./ ( 2 .* A );

eXK = zeros(numFaces,1);
eYK = Lk ./ ( 2 .* A );

%--------------------------------------------------------------------------
% Construct the Gradient Operator
%--------------------------------------------------------------------------

% The row indices of non-zero values
I = [ 0*numFaces + repmat(1:numFaces, 1, 4), ...
    1*numFaces + repmat(1:numFaces, 1, 4) ].';

% The column indices of non-zero values
J = repmat([ F(:,2); F(:,1); F(:,3); F(:,1) ], 2, 1);

% The actual non-zero values
V = [ eXJ; -eXJ; eXK; -eXK; eYJ; -eYJ; eYK; -eYK ];

% The sparse gradient operator
G = sparse(I, J, V, 2*numFaces, numVertex);

%--------------------------------------------------------------------------
% Construct and Solve the Linear Embedding Problem
%--------------------------------------------------------------------------

% Construct the basic problem ---------------------------------------------

% A matrix that performs a pi/2 CCW rotation about the normal
Irn = (1:(2*numFaces))';
Jrn = [(numFaces+1):(2*numFaces) 1:numFaces]';
Vrn = [ones(numFaces,1); -ones(numFaces,1)];
RN = sparse( Irn, Jrn, Vrn, 2*numFaces, 2*numFaces );

% The linear embedding operator
C = [ sparse(2*numFaces, numVertex), G ] - ...
    RN * [ G sparse(2*numFaces, numVertex) ];

% The (vanishing) RHS
d = zeros( 2*numFaces, 1 );

% Construct problem constraints -------------------------------------------
Ieq = (1:6);
Jeq = [ F(fID,:), numVertex + F(fID,:) ];
Aeq = sparse( Ieq, Jeq, 1, 6, 2*numVertex );

beq = [ 0; -Lk(fID); -Xk(fID); 0; 0; Yk(fID) ];

% Construct problem options -----------------------------------------------
options = optimoptions('lsqlin', 'Algorithm', ...
    'interior-point', 'Display', 'off');

% Solve the problem -------------------------------------------------------
V2D = lsqlin(C, d, [], [], Aeq, beq, [], [], [], options);

%--------------------------------------------------------------------------
% Format Output
%--------------------------------------------------------------------------

V2D = reshape(V2D, numVertex, 2);

% Re-orient normals upward if necessary -----------------------------------
TR = triangulation( F, [V2D, zeros(numVertex, 1)] );

% The face unit normals of the embedding
fN = TR.faceNormal;

% The sign of each face normal
pm = unique(sign(fN(:,3)));
if numel(pm) > 1
    warning(['Inconsistent face ordering. ' ...
        'Re-orient faces and attempt embeddding again']);
elseif pm < 0
    V2D = [ -V2D(:,1), V2D(:,2) ];
end

% Check edge length fidelity ----------------------------------------------
L2D = V2D(E(:,2),:) - V2D(E(:,1),:);
L2D = sqrt(sum(L2D.^2, 2));

edgeErr = max( abs(L2D-L) ./ L );
if edgeErr > 0.01
    warning(['Relative embedding edge length discrepancy of ', ...
        num2str(edgeErr)]);
end

end

