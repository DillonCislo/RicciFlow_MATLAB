function V2D = embedMetric2D_CircleQueue(L, F, fID, E, feIDx)
%EMBEDMETRICIN2D_CIRCLEQUEUE
%Given a flat discrete Riemannian metric on a triangulation
%(a set of positive edge lengths obeying the triangule inequality in each
%face) this function provides a parameterization of that matric as a set of
%2D vertex coordinates. Embedding is done by constraining the location of a
%seed face and then populating a queue of neighboring faces.  Faces in the
%queue are embedded by finding the (usually two) intersections of the
%circles defined by the edge lengths connecting the two already embedded
%vertices with the yet-to-be embedded vertex in the plane.  The position of
%the new vertex is chosen to be the intersection that maintains the normal
%direction of the planar triangulation.
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
% by Dillon Cislo 11/18/2019

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
% Create Embedding Queue and Embed the First Face
%--------------------------------------------------------------------------

% The edge lengths of the seed face
Li = L_F(fID,1); Lj = L_F(fID,2); Lk = L_F(fID,3);

% The internal angle of the first vertex
theta = acos( (Lj^2 + Lk^2 - Li^2) ./ (2 * Lj * Lk) );

% Embed the seed face in 2D
Vi = [0 0]; Vj = [Lk 0]; Vk = Lj .* [ cos(theta) sin(theta) ];

V2D = zeros( max(F(:)), 2 );
V2D( F(fID,:), : ) = [ Vi; Vj; Vk ];

% Create flags marking vertices that have already been embedded
embedFlags = F;
embedFlags(ismember(F, F(fID,:))) = 0;

% Create the embedding queue that holds neighboring faces to previously
% embedded faces
Que = [];

% A convenience variable to help populate the queue
fID_Counter = zeros(1, size(F,1));
fID_Counter(fID) = 1;

%--------------------------------------------------------------------------
% Embed All Remaining Faces
%--------------------------------------------------------------------------

while nnz(embedFlags) > 0
    
    % Update queue --------------------------------------------------------
    
    % Find the faces attached to the current face
    nbrs = TR.neighbors(fID);
    nbrs =  nbrs(~isnan(nbrs)); % Removing NaN's for boundary faces
    
    % Append neighbors to the end of the queue
    Que = [ Que nbrs ];
    
    % Pop off the next face to embed
    fID = Que(1); Que(1) = [];
    
    while nnz(embedFlags(fID,:)) == 0
        
        if fID_Counter(fID) == 0
            
            nbrs = TR.neighbors(fID);
            nbrs = nbrs(~isnan(nbrs));
            
            % Find unique face that have not yet been added to the queue
            nbrs = nonzeros(nbrs .* (~fID_Counter(nbrs)))';
            
            % Add unique new faces the end of the queue
            Que = [Que nbrs];
            
            % Mark the face as visited in the queue variable
            fID_Counter(fID) = 1;
            
        end
        
        fID = Que(1); Que(1) = [];
        
    end
    
    % Find the intersection of the circles centered at the previously
    % embedded vertices of the current face with a radius equal to those
    % edge lengths.  The intersection point that maintains the orientation
    % of the outward pointing face normal will be the location of the next
    % embedded vertex -----------------------------------------------------
    
    % Find the IDs of the previously embedded vertices
    embeddedInFace = find(embedFlags(fID,:) == 0);
    v1_IDx = F( fID, embeddedInFace(1) );
    v2_IDx = F( fID, embeddedInFace(2) );
    
    % The ID of the non-embedded vertex
    v3_IDx = nonzeros(embedFlags(fID,:));
    
    % Find the 2D coordinates of the previously embedded vertices
    x1 = V2D(v1_IDx, 1); y1 = V2D(v1_IDx, 2);
    x2 = V2D(v2_IDx, 1); y2 = V2D(v2_IDx, 2);
    
    % Find the radii of the corresponding circles (note indices used)
    r1 = L_F( fID, embeddedInFace(2) );
    r2 = L_F( fID, embeddedInFace(1) );
    
    % Find the intersection(s) of the circles
    [xout, yout] = DiscreteRicciFlow.circCircIntersect( x1, y1, r1, ...
        x2, y2, r2 );
    
    % Optional debugging for node positions
    if any(isnan([xout; yout]))
        warning('NaN in circle intersection calculations');
    end
    
    % Choose the which intersection to use as the new embedding location --
    
    % The the other face attached to the edge defined by the already
    % embedded vertices
    adjFace = cell2mat( TR.edgeAttachments( v1_IDx, v2_IDx ) );
    adjFace(adjFace == fID) = [];
    
    % The vertex positions of the alternative face
    pf = V2D(F(adjFace,:), :)';
    
    % Some convenience variables for the circle intersection points
    k1 = [ xout(1); yout(1) ];
    k2 = [ xout(2); yout(2) ];
    
    % The sum of the distances between the vertices and the intersection
    % points
    dk1 = sum(sqrt(sum((bsxfun(@minus, pf, k1)).^2)));
    dk2 = sum(sqrt(sum((bsxfun(@minus, pf, k2)).^2)));
    
    % Choose the point with the greatest sum of distances
    if dk1>dk2
        next_vertex = k1;
    else
        next_vertex = k2;
    end
    
    % Embedding new vertex coordinates into the list
    V2D(v3_IDx, :) = next_vertex';
    
    % Flagging new vertex as 'embedded'
    embedFlags(embedFlags == v3_IDx) = 0;
    fID_Counter(fID) = 1;
    
end

%--------------------------------------------------------------------------
% Format Output
%--------------------------------------------------------------------------

% Check normal orientation ------------------------------------------------
TR = triangulation( F, [V2D, zeros(size(V2D,1), 1)] );

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
    

