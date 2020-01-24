function boundaries = compute_boundaries(F)
%COMPUTE_BOUNDARIES Finds the boundary vertices of a potientially multiply
%connected mesh. Re-creates the functionality of 'compute_boundaries.m' in
%the ImSAnE package.
%
%   INPUT PARAMETERS:
%
%       - F:            #FxP polygonal face connectivity list (P=3 for a
%                       triangulation).  Faces do not need to have
%                       consistent orientation (even though they should and
%                       everyone's life would be better if they did)
%
%   OUTPUT PARAMETERS:
%
%       - boundaries:   1x#B cell array of boundary components.
%                       boundaries{i} is a 1x#BV row vector of the vertex
%                       IDs in the ith boundary componenent
%
% by Dillon Cislo 11/18/2019

% Build the unordered vertex adjacency list
A = sparse( F, F(:, [2:end 1]), 1 );
A = A+A';

% Find vertices in the edges that occur only once
[ bdyIDx, ~ ] = find(A == 1);
bdyIDx = unique(bdyIDx);

numBdy = 0; % The number of boundary components
boundaries = {}; % The output parameter
while ~isempty(bdyIDx)
    
    numBdy = numBdy + 1;
    
    % Choose the starting point of the current boundary component
    i = bdyIDx(1);
    
    % If a boundary vertex is connected to more than 2 other boundary
    % vertices, something is wrong (i.e. a triangle sharing only one
    % vertex with the rest of the mesh)
    u = find(A(i,:) == 1);
    if (numel(u) ~= 2), warning('Problem in boundary'); end
    
    boundary = [i u(1)];
    
    % Iterate over the rest of the boundary
    s = boundary(2);
    i = 2;
    while i <= size(A,2)
        
        u = find(A(s,:) == 1);
        if (numel(u) ~= 2), warning('Problem in boundary'); end
        
        if u(1) == boundary(i-1)
            s = u(2);
        else
            s = u(1);
        end
        
        if s ~= boundary(1)
            boundary = [ boundary s ];
        else
            break;
        end
        
        i = i+1;
       
    end
    
    bdyIDx = setdiff( bdyIDx, boundary );
    boundaries{numBdy} = boundary;
    
end

