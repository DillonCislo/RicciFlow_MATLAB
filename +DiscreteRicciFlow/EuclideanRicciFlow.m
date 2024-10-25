function [ L, V2D, allL ] = EuclideanRicciFlow( F, V3D, varargin )
%EUCLIDEANRICCIFLOW An implementation of the discrete Euclidean Ricci flow
%for shape deformation and parameterization.  Given a 3D mesh triangulation
%this pipeline produces a new set of mesh edge lengths such that the
%discrete Gaussian curvature is entirely localized on mesh boundary
%vertices (the precise distribution of the curvature is determined by user
%specified options).  The user can also optionally embed this discrete
%metric in 2D to produce a conformal parameterization of the 3D mesh.
%
%   INPUT PARAMETERS:
%
%       - F:            #Fx3 face connectivity list
%
%       - V3D:          #Vx3 3D vertex coordinate list
%       
%       OPTIONAL INPUTS:
%
%       - 'Display':        The level of detail provided in user feedback
%                           - {'true'}  Displays details about Newton
%                           iterations and other processes
%                           - 'false' Displays nothing
%
%       - 'EdgeLengths':    #Ex1 list of initial edge lengths. This option
%                           is used when we wish to supply an intrinsic
%                           triangulation instead of an extrinsic one (i.e.
%                           with V3D = [])
%
%       - 'BoundaryType':   - {'Free'} Boundary is free to take any shape
%                           - 'Fixed'  Boundary is set to take one of a
%                           particular set of pre-defined shapes
%
%       - 'BoundaryShape':  The shape of a fixed boundary
%                           - {'Polygon'} Target curvature is distributed
%                           evenly over a supplied set of boundary
%                           singularities
%                           - 'Circles' All boundaries are set to be
%                           circles
%
%       - 'Cones':          The set of boundary vertices that will become
%                           the vertices of a 'Polygon' type fixed boundary
%
%       - 'Distortion':     Determines handling of area distortion
%                           - {'Basic'} No extra effort is made to reduce
%                           area distortion
%                           - 'Optimal' Area distortion is minimized as
%                           much as possible while maintaining the
%                           conformal property of the parameterization
%
%       - 'Tolerance':      Newton method minimization of the Ricci energy
%                           will terminate when the maximum absolute
%                           difference between the current curvature and
%                           the target curvature is less than tolerance
%
%       - 'CircTolerance':  Curvature flow towards circle shaped boundaries
%                           will terminate when the maximum absolute change
%                           in target curvatures between iterations is less
%                           than the circle tolerance
%
%       - 'PCGTolerance':   Tolerance for the PCG method used to calculate
%                           the update step for the Ricci energy
%                           minimization
%
%       - 'MaxIter':        The maximum number of iterations allowed for
%                           the Newton method minimization of the Ricci
%                           energy
%
%       - 'MaxCircIter':    The maximum number of iterations allowed for
%                           optimization of circle shape boundary
%                           curvatures
%
%       - 'MaxPCGIter':     The maximum number of iterations allowed for
%                           each calculation of the update step in the
%                           minimization of the Ricci energy
%
%       - 'Embedding':      The method used to embed the discrete metric in
%                           the plane
%                           - {'IsoEnergy'} Find the embedding that
%                           minimizes an energy that must vanish for a
%                           valid isometric embedding.  Faces must be
%                           consistently oriented to yield a valid result
%                           - 'CircIntersect' Iterative embed faces by
%                           checking for intersections of circles defined
%                           by edge lengths. Faces do not need to be
%                           consistently oriented, however this method is
%                           much slower and prone to self-intersections
%
%       - 'ScaleEmbedding': Whether or not to re-scale the embedding to the
%                           unit disk
%                           - {'true'} Vertex coordinates of the embedding
%                           are re-scaled to lie within the unit disk
%                           - 'false' vertex coordinates are not re-scaled
%
%       - 'ScaleMetric':    Whether or not to scale the output lengths
%                           produced by Ricci flow after embedding
%                           - {'true'} Output lengths are replaced by the
%                           edge lengths of the 2D embedding
%                           - 'false' Output lengths are returned as is
%
%       - 'ZeroID':         The ID of the point to conformally map to
%                           (0, 0). This is only used if the output
%                           embedding is mapped to the unit disk and if the
%                           input mesh is a simple topological disk
%
%       - 'OneID'           The ID of the point to conformally map to
%                           (1, 0). This is only used if the output
%                           embedding is mapped to the unit disk and if the
%                           input mesh is a simple topological disk
%                           
%
%   OUTPUT PARAMETERS:
%       
%       - L:            #Ex1 list of final edge lengths.  This is the
%                       representation of the discrete Riemannian metric
%                       with the target curvature
%
%       - V2D:          #Vx2 2D vertex coordinate list
%
%       - allL:         The edge lengths of each intermediate iteration
%                       along the flow.  Each entry in this cell array is
%                       an #Ex1 list of edge lengths
%
% by Dillon Cislo 11/16/2019

%==========================================================================
% INPUT PROCESSING
%==========================================================================

% Ensure mandatory input parameters are supplied
if (nargin < 1), error('Please supply face connectivity list'); end
if (nargin < 2), V3D = []; end

% Verify the properties of the input triangulation topology
validateattributes( F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive'} );

%--------------------------------------------------------------------------
% Determine mesh topology
%--------------------------------------------------------------------------

% Construct MATLAB-style representation of the input triangulation
if isempty(V3D)
    
    numVertex = max(F(:));
    TR = triangulation(F, zeros(numVertex,3));
    
else
    
    validateattributes( V3D, {'numeric'}, ...
        {'2d', 'ncols', 3, 'finite', 'nonnan'} );
    numVertex = size(V3D,1);
    TR = triangulation( F, V3D ); 
    
end

% Vertex IDs defining edges
E = sort( TR.edges, 2 );

% numFaces = size(F,1); % The number of faces
numEdges = size(E,1); % The number of edges

% The Euler characteristic of the mesh
% EulerChi = numVertex - numEdges + numFaces;

% Compute the boundaries of the mesh --------------------------------------

% bdy{i} is an ordered list of vertex IDs defining the ith boundary
bdy = DiscreteRicciFlow.compute_boundaries( F );

% A list of all boundary vertices for convenience
allBdyIDx = horzcat(bdy{:});

% A list of all vertices not on the boundary for convenience
interiorIDx = setdiff( 1:numVertex, allBdyIDx );

numBdy = numel(bdy); % The number of mesh boundaries

% The genus of the mesh
% genus = ( 2 - numBdy - EulerChi ) / 2;

%**************************************************************************
% TODO: Add subroutine to cut closed surfaces (medium priority)
assert( numBdy > 0, ...
    'Closed surfaces are not currently supported' );
%**************************************************************************

%--------------------------------------------------------------------------
% Process Optional Inputs
%--------------------------------------------------------------------------

% Parameter name sets - used to validate user input
allBdyTypes = { 'free', 'fixed' };
allBdyShapes = { 'polygon', 'circles' };
allDistTypes = { 'basic', 'optimal' };
allEmbedTypes = { 'circintersect', 'isoenergy' };

% Default parameters (see function documentation)
iterDisp = true;
bdyType = 'free';
bdyShape = 'polygon';
cones = [];
distHandling = 'basic';
tol = 10e-7;
circTol = 10e-8;
pcgTol = 10e-7;
maxIter = 1000;
maxCircIter = 1000;
maxPCGIter = numVertex;
embedMethod = 'isoenergy';
scaleEmbedding = true;
scaleMetric = true;
L0 = [];
zeroID = [];
oneID = [];

for i = 1:length(varargin)
    if isa(varargin{i},'double') 
        continue;
    end
    if isa(varargin{i},'logical')
        continue;
    end
    if ~isempty(regexp(varargin{i},'^[Dd]isplay','match'))
        iterDisp = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Ee]dge[Ll]engths','match'))
        L0 = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Bb]oundary[Tt]ype','match'))
        bdyType = lower(varargin{i+1});
    end
    if ~isempty(regexp(varargin{i},'^[Bb]oundary[Ss]hape','match'))
        bdyShape = lower(varargin{i+1});
    end
    if ~isempty(regexp(varargin{i},'^[Cc]ones','match'))
        cones = varargin{i+1};
        if ( numel(cones) <= 3 )
            warning(['Insufficient number of cones supplied. ', ...
                'Using maximal polygon']);
            cones = [];
        end 
    end
    if ~isempty(regexp(varargin{i}, '^[Dd]istortion', 'match'))
        distHandling = lower(varargin{i+1});
    end
    if ~isempty(regexp(varargin{i}, '^[Tt]olerance', 'match'))
        tol = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Cc]irc[Tt]olerance', 'match'))
        circTol = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Pp][Cc][Gg][Tt]olerance', 'match'))
        pcgTol = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Mm]ax[Ii]ter', 'match'))
        maxIter = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Mm]ax[Cc]irc[Ii]ter', 'match'))
        maxCircIter = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Mm]ax[Pp][Cc][Gg][Ii]ter', 'match'))
        maxPCGIter = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Ee]mbedding', 'match'))
        embedMethod = lower(varargin{i+1});
    end
    if ~isempty(regexp(varargin{i}, '^[Ss]cale[Ee]mbedding', 'match'))
        scaleEmbedding = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Ss]cale[Mm]etric', 'match'))
       scaleMetric = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'ZeroID')
        zeroID = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'OneID')
        oneID = varargin{i+1};
    end
end

% Process initial edge length input
if isempty(L0)
    assert(~isempty(V3D), ['You must supply either an embedded ' ...
        'triangulation or an initial discrete metric']);
else
    validateattributes(L0, {'numeric'}, {'vector', 'positive', ...
        'finite', 'real', 'numel', numEdges});
    if (size(L0,2) ~= 1), L0 = L0.'; end
end

% Check that boundary type is valid
if ~any(strcmp(bdyType, allBdyTypes))
    error('Invalid boundary type');
end

% Check that boundary shape is valid
if ~any(strcmp(bdyShape, allBdyShapes))
    error('Invalid boundary shape');
end

%**************************************************************************
% TODO: Free boundary conditions don't seem to work for multiply connected
% meshes. I have no idea why - figure out if your code is crap or if this
% is a meaningful restriction (medium priority)
if (numBdy > 1) && strcmp(bdyType, 'free')
    warning(['Free boundary conditions not currently supported ', ...
        'for multiply connected meshes. Using the maximal polygon ', ...
        'mapping instead.']);
    bdyType = 'fixed';
end
%**************************************************************************

% Check that distortion handling is valid
if ~any(strcmp(distHandling, allDistTypes))
    error('Invalid distortion handling type');
end

% Check that the embedding method is valid
if ~any(strcmp(embedMethod, allEmbedTypes))
    error('Invalid embedding method');
end

% Check that cone singularities are valid
numCones = numel(cones);
if numCones ~= 0
    
    if any( ~ismember( cones, allBdyIDx ) )
        error('Cone singularities must lie on mesh boundaries');
    end
    
    % Override any input cones if desired boundary shape is set to circles
    if strcmp(bdyShape, 'circles')
        warning(['Desired boundary shape is circles.' ...
            'Ignoring user supplied cones']);
        cones = [];
        numCones = 0;
    end
    
    %**********************************************************************
    % TODO: Add handling for sub-maximal polygonal boundaries for multiply
    % connected surfaces (low priority)
    if strcmp(bdyShape, 'polygon') && (numBdy > 1)
        warning(['Circles and maximal polygons are currently the ' ...
            'only supported boundary shapes for multiply connected ' ...
            'surfaces.  Ignoring user supplied cones']);
        cones = [];
        numCones = 0;
    end
    %**********************************************************************
    
end

% Validate minimization parameters
validateattributes( tol, {'numeric'}, {'scalar', 'real', 'positive'} );
validateattributes( circTol, {'numeric'}, {'scalar', 'real', 'positive'} );
validateattributes( pcgTol, {'numeric'}, {'scalar', 'real', 'positive'} );

validateattributes( maxIter, {'numeric'}, ...
    {'scalar', 'integer', 'positive'} );
validateattributes( maxCircIter, {'numeric'}, ...
    {'scalar', 'integer', 'positive'} );
validateattributes( maxPCGIter, {'numeric'}, ...
    {'scalar', 'integer', 'positive'} );

% Determine whether or not to record the intermediate edge lengths
if nargout > 2
    recordAllL = true;
else
    recordAllL = false;
end

% Process the zero point
if ~isempty(zeroID)
    if (numBdy ~= 1)
        zeroID = []; % Only use for topological disks
    else
        validateattributes(zeroID, {'numeric'}, {'scalar', 'integer', ...
            'finite', 'positive', 'real'});
        assert(~ismember(zeroID, allBdyIDx), ...
            'The zero point cannot lie on the mesh boundary');
    end
end

% Process the one point
if ~isempty(oneID)
    if (numBdy ~= 1)
        oneID = []; % Only use for topological disks
    else
        validateattributes(oneID, {'numeric'}, {'scalar', 'integer', ...
            'finite', 'positive', 'real'});
        assert(ismember(oneID, allBdyIDx), ...
            'The one point must lie on the mesh boundary');
    end
end


%==========================================================================
% CONSTRUCT CONNECTIVITY STRUCTURE TOOLS
%==========================================================================

% Display output
if iterDisp, fprintf('Calculating intial parameters... '); end

% Construct face-edge correspondence tool ---------------------------------
% Given a list of scalar edge quantities, 'EQ', the output of
% 'EQ(feIDx(f,i))' is that quantity corresponding to the edge opposite the
% ith vertex in face f

e1IDx = sort( [ F(:,3), F(:,2) ], 2 );
e2IDx = sort( [ F(:,1), F(:,3) ], 2 );
e3IDx = sort( [ F(:,2), F(:,1) ], 2 );

[~, e1IDx] = ismember( e1IDx, E, 'rows' );
[~, e2IDx] = ismember( e2IDx, E, 'rows' );
[~, e3IDx] = ismember( e3IDx, E, 'rows' );

feIDx = [ e1IDx e2IDx e3IDx ];

% Construct vertex face attachment list -----------------------------------
% vA{i} is the list of face IDs attached to the ith vertex
vA = vertexAttachments( TR );

% Construct vertex edge attachment list -----------------------------------
% A tool used to help sum edge based quantities over all of the edges
% attached to a particular vertex. V2E(:,1) is a list of (repeating) vertex
% IDs.  A vertex will be repeated as many times as there are edges attached
% to it. V2E(i,2) will return an edge ID attached to the vertex stored in
% V2E(i,1) (in no particular order)
V2E = [ E(:), repmat((1:numEdges)', 2, 1) ];

%==========================================================================
% CALCULATE INITIAL CIRCLE PACKING METRIC AND INVERSIVE DISTANCES
%==========================================================================

% Calculate initial mesh edge lengths
if isempty(L0)
    L0 = V3D( E(:,2), : ) - V3D( E(:,1), : );
    L0 = sqrt( sum(L0.^2, 2) );
end
L0_F = L0(feIDx);

assert( all(all( sum(L0_F, 2) - 2 .* L0_F > 0 )), ...
    'Input edge lengths violate the triangle inequality');

% Handle output
if recordAllL
    allL = {L0};
else
    allL = {};
end

% Compute circle packing radii (gamma) for each vertex in each face
gamma_F0 = ( repmat( sum(L0_F,2), [1 3] ) - 2 .* L0_F ) ./ 2;

% The circle packing radius for each vertex is chosen be the minimum gamma
% value from among all of the faces attached to that vertex
gamma_V0 = zeros( numVertex, 1 );
for v = 1:numVertex
    gamma_V0(v) = min( gamma_F0( vA{v}) );
end

% Update the radii list for each face and edge
gamma_E0 = gamma_V0(E);

% Calculate the inversive distance for each edge
invD_E = ( L0.^2 - sum(gamma_E0.^2, 2) ) ./ ...
    ( 2 .* gamma_E0(:,1) .* gamma_E0(:,2) );

%==========================================================================
% SET TARGET CURVATURE
%==========================================================================

% Free boundaries ---------------------------------------------------------
if strcmp( bdyType, 'free' )
    
    freeBoundary = true;
    
    if strcmp(distHandling, 'basic')
        
        % For the basic free boundary condition we set the target curvature
        % of all interior points to zero and only update the radii of
        % interior vertices during minimization.  In this way the curvature
        % is automatically distributed along the boundary
        Ktar = zeros(numVertex, 1);
        
    elseif strcmp(distHandling, 'optimal')
        
        %******************************************************************
        % TODO: Add optimal distortion handling for free boundaries (high
        % priority)
        %
        % I don't really understand how this works yet.  Apparently there
        % is some method of curvature diffusion, but maybe this is just the
        % basic free boundary condition in disguise.
        error(['Optimal area distortion for free boundaries' ...
            'is not currently supported']);
        %******************************************************************
        
    end
    
% Fixed boundaries --------------------------------------------------------
else
    
    freeBoundary = false;
    
    if strcmp(distHandling, 'basic')
    
        % Sub-Maximal Polygons --------------------------------------------
        % Currently only supported for topological disks.  In this case
        % the total curvature (2*pi) is evenly distributed among the user
        % supplied cones.
        if numCones > 0
            
            Ktar = zeros(numVertex, 1);
            Ktar(cones) = 2 * pi / numCones;
            
        % Maximal Polygons/Circles ----------------------------------------
        % Here a total curvature of 2*pi is assigned to the outermost
        % boundary and a total curvature of -2*pi is assigned to each inner
        % boundary.  These curvatures are distributed evenly among all
        % vertices belonging to that boundary component.  The outermost
        % boundary is set to be the one with the longest length in the 3D
        % configuration.
        else
            
            Ktar = zeros(numVertex, 1);
            
            % Determine the outermost boundary ----------------------------
            if numBdy > 1
                
                % Initial guess for outer boundary
                outerBdy = 1;
                
                % Length of the first boundary component
                bdyE = sort([bdy{1}' circshift(bdy{1}, -1)'], 2);
                bdyE = ismember( E, bdyE, 'rows' );
                maxBdyL = sum(L0(bdyE));
                
                for i = 2:numBdy
                    
                    % The length of the current boundary component
                    bdyE = sort([bdy{i}' circshift(bdy{i}, -1)'], 2);
                    bdyE = ismember( E, bdyE, 'rows' );
                    curBdyL = sum(L0(bdyE));
                    
                    % Outer boundary is set to longest boundary
                    if curBdyL > maxBdyL
                        maxBdyL = curBdyL;
                        outerBdy = i;
                    end
                    
                end
                
                % Shift outermost boundary to first position in boundary
                % component list
                bdy = circshift( bdy, -(outerBdy-1));
                
            end
            
            % Assign target curvatures ------------------------------------
            for i=1:numBdy
                
                % Maximal Polygons ----------------------------------------
                % The total curvature for each boundary is distributed
                % evenly among all vertices belonging to that boundary
                if strcmp(bdyShape, 'polygon')
                    if i == 1
                        Ktar(bdy{i}) = 2 * pi / numel(bdy{i});
                    else
                        Ktar(bdy{i}) = -2 * pi / numel(bdy{i});
                    end
                
                % Circles -------------------------------------------------
                % The total curvature is distributed as a ratio of the sum
                % of the lengths of the boundary edges attached to the
                % vertex to twice the total length of that boundary
                % component
                elseif strcmp(bdyShape, 'circles')
                    
                    % Compute the length of each boundary edge ------------
                    
                    % The IDs of the edges in the current boundary
                    bdyE = sort([bdy{i}' circshift(bdy{i}, -1)'], 2);
                    bdyE = find(ismember(E, bdyE, 'rows'));
                    
                    % The lengths of the edges in the current boundary
                    bdyL = L0(bdyE);
                    
                    % The total length of the current boundary
                    totalBdyL = sum(bdyL);
                    
                    % Calculate the sum of the two boundary edges attached
                    % to each boundary vertex -----------------------------
                    
                    % The vertex IDs defining each boundary edge
                    bdyEIDx = E(bdyE, :);
                    
                    newKV = zeros(numVertex, 1);
                    for v = bdy{i}
                        newKV(v) = sum(bdyL(any(bdyEIDx == v, 2)));
                    end
                    
                    % Update the target curvatures ------------------------
                    newKV = pi .* newKV ./ totalBdyL;
                    if (i > 1), newKV = -newKV; end
                    
                    Ktar = Ktar + newKV;
                    
                end
                
            end
            
        end
        
    elseif strcmp(distHandling, 'optimal')
        
        %******************************************************************
        % TODO: Add optimal distortion handling for fixed boundaries (high
        % priority)
        %
        % I definitely don't understand how this works. Apparently there is
        % some projected gradient method that can interatively flow towards
        % a superior Ktar that make area distortion uniform, but I'm not
        % sure if you have to run this during Ricci energy minimization or
        % if it is separate.  Could be very cool if I can get this working
        error(['Optimal area distortion for fixed boundaries' ...
            'is not currently supported']);
        %******************************************************************
        
    end
    
end

%==========================================================================
% CALCULATE INITIAL VERTEX CURVATURE
%==========================================================================

% Calculate internal angles from triangulation edge lengths ---------------

% Some convenience variables to vectorize the cosine law calculation
Gi = L0_F; Gj = circshift(L0_F, [0 -1]); Gk = circshift(L0_F, [0 -2]);

% The internal angles
intAng = ( Gj.^2 + Gk.^2 - Gi.^2 ) ./ ( 2 .* Gj .* Gk );
intAng = acos(intAng);

% Combine sum internal angles around vertex 1-ring to calculate Gaussian
% curvatures --------------------------------------------------------------
K0 = sparse( F(:), 1, intAng(:), numVertex, 1 );
K0(interiorIDx) = 2 * pi - K0(interiorIDx);
K0(allBdyIDx) = pi - K0(allBdyIDx);

if iterDisp, fprintf('Done\n'); end

%==========================================================================
% MINIMIZE RICCI ENERGY VIA NEWTON METHOD
%==========================================================================

iterNum = 1;
if iterDisp
    fprintf('Beginning Ricci Energy Minimization Round %d\n', iterNum);
end

[ K, L, gamma_V, allL ] = DiscreteRicciFlow.minimizeRicciEnergy( F, E, ...
    V3D, feIDx, V2E, allBdyIDx, ...
    K0, Ktar, L0, gamma_V0, invD_E, ...
    tol, pcgTol, maxIter, maxPCGIter, iterDisp, freeBoundary, ...
    allL, recordAllL );

% Update iteration number
iterNum = iterNum+1;

%--------------------------------------------------------------------------
% Flow Mesh Curvatures to Generate Circle Shaped Boundary
%--------------------------------------------------------------------------

if strcmp( bdyShape, 'circles' )
    
    while true % A shitty approximation of a do-while loop
        
        %------------------------------------------------------------------
        % Update target curvatures
        %------------------------------------------------------------------
        
        % The previous steps target curvatures
        prevKtar = Ktar;
        Ktar = zeros(numVertex, 1);
        
        for i = 1:numBdy
            
            % Compute the length of each boundary edge under the current
            % metric ------------------------------------------------------
            
            % The IDs of the edges in the current boundary
            bdyE = sort([bdy{i}' circshift(bdy{i}, -1)'], 2);
            bdyE = find(ismember( E, bdyE, 'rows' ));
            
            % The lengths of the edges in the current boundary
            bdyL = L(bdyE);
            
            % The total length of the current boundary
            totalBdyL = sum(bdyL);
            
            % Calculate the sum of the two boundary edges attached to each
            % boundary vertex ---------------------------------------------
            
            % The vertex IDs defining each boundary edge
            bdyEIDx = E(bdyE, :);
            
            newKV = zeros(numVertex, 1);
            for v = bdy{i}
                newKV(v) = sum(bdyL(any(bdyEIDx == v, 2)));
            end
            
            % Update the target curvatures --------------------------------
            newKV = pi * newKV ./ totalBdyL;
            if (i > 1), newKV = -newKV; end
            
            Ktar = Ktar + newKV;

        end
        
        %------------------------------------------------------------------
        % Check for convergence
        %------------------------------------------------------------------
        
        % The circle flow convergence parameter
        circleErr = max(abs(Ktar-prevKtar));
        
        if circleErr < circTol
            break;
        elseif iterNum > maxCircIter
            if iterDisp
                fprintf(['Maximum iteration number exceeded. ' ...
                    'Terminating circle flow\n']);
            end
            break;
        end
        
        %------------------------------------------------------------------
        % Minimize Ricci energy for new target curvatures
        %------------------------------------------------------------------
        
        if iterDisp
            fprintf(...
                [ 'Beginning Ricci Energy Minimization Round %d\n', ...
                'Maximum curvature difference = %f\n' ], ...
                [ iterNum, circleErr ] );
        end
        
        [ K, L, gamma_V, allL ] = ...
            DiscreteRicciFlow.minimizeRicciEnergy( F, E, V3D, ...
            feIDx, V2E, allBdyIDx, ...
            K, Ktar, L, gamma_V, invD_E, ...
            tol, pcgTol, maxIter, maxPCGIter, iterDisp, freeBoundary, ...
            allL, recordAllL );
        
        % Increase iteration count
        iterNum = iterNum+1;
        
    end
    
end

%--------------------------------------------------------------------------
% Check that Final Output is a Valid Metric
%--------------------------------------------------------------------------

if any( L <= 0 )
    warning('Final metric has non-positive edge lengths!');
end

L_F = L(feIDx);
if any( (sum(L_F, 2) - 2 .* L_F) <= 0 )
    warning('Final metric does NOT satisfy the triangle inequality!');
end
            
%==========================================================================
% EMBED OUTPUT METRIC IN R^2
%==========================================================================
% NOTE: Both supported embedding procedures require a seed face whose
% vertex positions will be constrained.  Internally this algorithm will
% alays use the first face in the connectivity list.  Alternative seed
% faces can be tried externally by calling the embedding functions on their
% own.

if nargout > 1
    
    if iterDisp, fprintf('Embedding metric in 2D '); end
    
    if strcmp(embedMethod, 'isoenergy')
        if iterDisp
            fprintf('using the isometric energy method... ');
        end
        V2D = DiscreteRicciFlow.embedMetric2D_IsoEnergy( L, F, 1, ...
            E, feIDx );
    elseif strcmp(embedMethod, 'circintersect')
        if iterDisp
            fprintf('using the circle intersection method... ');
        end
        V2D = DiscreteRicciFlow.embedMetric2D_CircleQueue( L, F, 1, ...
            E, feIDx );
    end
    
    if iterDisp, fprintf('Done\n'); end
    
    % Re-scale circle output to unit disk ---------------------------------
    if strcmp(bdyType, 'fixed') && (numCones ==  0) && scaleEmbedding
        
        if iterDisp, fprintf('Re-scaling embedding to unit disk... '); end
        
        V2D = ( V2D - min(V2D, [], 1) );
        V2D = V2D ./ max(V2D);
        V2D = 2 .* ( V2D - 0.5 );

        if (numBdy == 1)

            % Strongly pin outer boundary to the unit disk
            V2D = complex(V2D(:,1), V2D(:,2));
            V2D(allBdyIDx) = exp(1i .* angle(V2D(allBdyIDx)));
            
            % Map zero point to (0, 0)
            if ~isempty(zeroID)
                V2D = (V2D - V2D(zeroID)) ./ (1 - conj(V2D(zeroID)) .* V2D);
            end

            % Map one point to (1, 1)
            if ~isempty(oneID)
                V2D = exp(-1i .* angle(V2D(oneID))) .* V2D;
            end

            V2D = [real(V2D), imag(V2D)];

        end
        
        if iterDisp, fprintf('Done\n'); end
        
    end
    
    % Re-scale output metric to match embedding edge lengths --------------
    if scaleMetric
        
        if iterDisp, fprintf('Re-scaling output metric... '); end
        
        L = V2D(E(:,2), :)-V2D(E(:,1), :);
        L = sqrt(sum(L.^2, 2));
        
        if iterDisp, fprintf('Done\n'); end
        
    end
    
end

end
