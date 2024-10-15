function [ K, L_E, gamma_V, allL ] = minimizeRicciEnergy( F, E, V3D, ...
    feIDx, V2E, bdyIDx, ...
    K0, Ktar, L0, gamma_V0, invD_E, ...
    tol, pcgTol, maxIter, maxPCGIter, iterDisp, freeBoundary, ...
    allL, recordAllL )
%MINIMIZERICCIENERGY Minimizes the Ricci Energy of a mesh triangulation
%given a set of vertex target curvatures and an initial circle packing
%metric using the inversive distance scheme.  Intended for internal use
%with the 'DiscreteRicciFlow' package - WARNING: NO INPUT VALIDATION IS
%PERFORMED
%
%   INPUT PARAMETERS:
%
%       - F:            #Fx3 face connectivity list
%       
%       - E:            #Ex2 edge connectivity list
%
%       - V3D:          #Vx3 3D vertex coordinate list
%
%       - feIDx:        #Fx3 face-edge correspondence tool.  feIDx(f,i) is
%                       the ID of the edge opposite vertex i in face f
%
%       - V2E:          #VEx2 vertex-edge correspondence tool.  Used to sum
%                       edge-based quantites over all edges attached to a
%                       particular vertex.
%
%       - bdyIDx:       #BVx1 list of boundary vertex IDs (includes all
%                       boundary components)
%
%       - K0:            #Vx1 list of initial vertex curvatures
%
%       - Ktar:         #Vx1 list of vertex target curvatures
%
%       - L0:            #Ex1 list of initial edge lengths
%
%       - gamma_V0:      #Vx1 list of circle packing radii
%       
%
%       - invD_E:       #Ex1 list of edge inversive distances
%
%       - tol:          The tolerance of the minimization procedure.  The
%                       method will terminate when the maximum absolute
%                       difference between the current curvature and the
%                       target curvature is less than tol
%
%       - pcgTol:       The tolerance for the pre-conditioned conjugate
%                       gradient method used to solve for the Newton method
%                       update step
%
%       - maxIter:      The maximum number of allowed minimization
%                       iterations
%
%       - maxPCGIter:   The maximum number of allowed iterations for the
%                       pre-conditioned conjugate gradient method used to
%                       solve for the Newton method update step
%
%       - iterDisp:     The level of output detail supplied about the
%                       minimization
%
%       - freeBoundary: Boolean specifying free vs fixed boundary shape
%
%       - allL:         A cell array containing the intermediate edge
%                       lengths from each iteration along the flow
%
%       - recordAllL:   Boolean indicating whether or not to record
%                       intermediate edge lengths
%
%
%   OUTPUT PARAMETERS:
%
%       - K:            #Vx1 list of output vertex target curvatures
%
%       - L:            #Ex1 list of output edge lengths
%
%       - gamma_V:      #Vx1 list of output circle packing radii
%
%       - allL:         A cell array containing the intermediate edge
%                       lengths from each iteration along the flow
%
% by Dillon Cislo 11/17/2019

%--------------------------------------------------------------------------
% Minimization Pre-Processing
%--------------------------------------------------------------------------

% Set current mesh properties ---------------------------------------------

% The Gaussian curvature
K = K0;

% The edge lengths
L_E = L0;
L_F = L_E(feIDx);

% The circle packing radii
gamma_V = gamma_V0;
gamma_F = gamma_V0(F);

% The conformal factors
u_V = log(gamma_V);

% Set mesh connectivity structure tools -----------------------------------
% (Some potential to optimize at the expense of clarity)
numVertex = size(K0,1);
numEdges = size(E,1);

% Vertex IDs of interior vertices
intIDx = setdiff(1:numVertex, bdyIDx, 'stable');

% The number of interior vertices
numIntVertex = numel(intIDx);

% A version of the V2E array without any boundary vertex contributions
intV2E = V2E;
intV2E( ismember(V2E(:,1), bdyIDx), : ) = [];
[~, intV2E(:,1)] = ismember( intV2E(:,1), intIDx );

% Edge IDs of interior edges
intEIDx = find( ~any( ismember( E, bdyIDx ), 2 ) );

intE = E(intEIDx, :);
[~, intE] = ismember( intE, intIDx );

% Build Hessian construction tools ----------------------------------------
% The row and column indicies of the non-zero elements of the sparse
% Hessian matrix.  These will not change and should only be computed once

if false %freeBoundary
    
    % For free boundary conditions we ignore boundary vertex contributions
    hessRow = [ intV2E(:,1); intE(:,1); intE(:,1) ];
    hessCol = [ intV2E(:,1); intE(:,2); intE(:,1) ];
    
else
    
    hessRow = [ V2E(:,1); E(:,1); E(:,2) ];
    hessCol = [ V2E(:,1); E(:,2); E(:,1) ];
    
end

% Set minimization variables ----------------------------------------------

% Calculate the initial error
if freeBoundary
    err = max(abs(K(intIDx)-Ktar(intIDx)));
else
    err = max(abs(K-Ktar));
end

% The current interation number
iterNum = 0;

%--------------------------------------------------------------------------
% Run minimization
%--------------------------------------------------------------------------

% Display progress
if iterDisp
    fprintf('Iteration   Error\n');
    fprintf('%4d %15f\n', [ iterNum, err ] );
end

while err > tol
    
    % For each face calculate the distance from each edge to the power
    % circle center -------------------------------------------------------
    h_F = DiscreteRicciFlow.computePowerDists( L_F, gamma_F );
    
    % Calculate edge weights ----------------------------------------------
    h_F = h_F ./ L_F;
    w_E = sparse( feIDx(:), 1, h_F(:), numEdges, 1 );
    
    % Construct the Hessian of the Ricci energy ---------------------------
    % Hij = { -wij          [vi, vj] \in M
    %       { sum(wik,k)    i = j
    %       { 0             else
    
    if false %freeBoundary
        hessVal = [ w_E(intV2E(:,2)); -w_E(intEIDx); -w_E(intEIDx) ];
        H = sparse( hessRow, hessCol, hessVal, ...
            numIntVertex, numIntVertex );
    else
        hessVal = [ w_E(V2E(:,2)); -w_E; -w_E ];
        H = sparse( hessRow, hessCol, hessVal, numVertex, numVertex );
    end
    
    %**********************************************************************
    % Solve the linear system for the update step direction ---------------
    % TODO: Add a legitimate line search step size method based on the
    % strong Wolfe/Armijo conditions (medium priority)
    
    % Free boundary conditions - ignore boundary node entries
    if freeBoundary
        
        H(bdyIDx, :) = [];
        H(:, bdyIDx) = [];
        
        % Pre-conditioned conjugate gradient method
        HL = ichol(H);
        [du_int, flag] = pcg( H, (Ktar(intIDx)-K(intIDx)), ...
            pcgTol, maxPCGIter, HL, HL' );
        
        % Update conformal factors along step
        u_V(intIDx) = u_V(intIDx) + du_int;
        
    else
        
        % Pre-conditioned conjugate gradient method
        HL = ichol(H);
        [du, flag] = pcg( H, (Ktar-K), pcgTol, maxPCGIter, HL, HL' );
        
        u_V = u_V + du;
        u_V = u_V - mean(u_V);
        
    end
    %**********************************************************************
    
    % Update circle packing radii -----------------------------------------
    gamma_V = exp(u_V);
    gamma_F = gamma_V(F);
    
    % Update edge lengths from the new circle packing according to the
    % inversive distance scheme -------------------------------------------
    gammaI = gamma_V(E(:,1)); gammaJ = gamma_V(E(:,2));
    L_E = sqrt( gammaI.^2 + gammaJ.^2 + 2 .* gammaI .* gammaJ.* invD_E );
    L_F = L_E(feIDx);
    
    % Record intermediate edge lengths
    if recordAllL, allL = [ allL, L_E ]; end
    
    % Calculate internal angles from triangulation edge lengths -----------
    
    % Some convenience variables to vectorize the cosine law calculation
    Gi = L_F; Gj = circshift(L_F, [0 -1]); Gk = circshift(L_F, [0 -2]);
    
    % The internal angles
    intAng = ( Gj.^2 + Gk.^2 - Gi.^2 ) ./ ( 2 .* Gj .* Gk );
    
    % Some optional debugging
    if any(abs(intAng)  > 1)
        
        if ~isempty(V3D)
            
            badFace = any(abs(intAng) > 1, 2 );
            faceColors = 0.8 .* ones(size(F));
            faceColors(badFace, : ) = repmat([1 0 0], sum(badFace), 1);
            
            trisurf( triangulation(F, V3D), ...
                'FaceVertexCData', faceColors, ...
                'FaceColor', 'flat', 'EdgeColor', 'k' );
            axis equal
            
        end

        error('acos out of bounds: probably bad triangle quality');
        
    end
    
    intAng = acos(intAng);
    
    % Calculate the Gauss curvature on each vertex ------------------------
    K = full(sparse( F(:), 1, intAng(:), numVertex, 1 ));
    K(intIDx) = 2 * pi - K(intIDx);
    K(bdyIDx) = pi - K(bdyIDx);
    
    % Update minimization parameters --------------------------------------
    
    if freeBoundary
        err = max(abs(K(intIDx)-Ktar(intIDx)));
    else
        err = max(abs(K-Ktar));
    end
    
    iterNum = iterNum + 1;
    
    % Display progress
    if iterDisp
        fprintf('%4d %15f\n', [ iterNum, err ] );
    end
    
    % Terminate if maximum iteration number is exceeded
    if iterNum > maxIter
        if iterDisp
            disp(['Maximum number of iterations exceeded. ', ...
                'Terminating minimization of Ricci energy']);
        end
        break;
    end
    
end


end

