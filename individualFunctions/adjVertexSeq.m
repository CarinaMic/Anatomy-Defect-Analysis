% Identify sorted adjacent faces around centre vertice (circle) 
% For open and closed meshes

% Input:    boundVerticesNr (for open meshes): row number of boundary vertices (parameter data from function boundSTL)
%           vertexMap (cell format): centre vertex (row number of cell corresponds to centre vertex) and 
%                                   its adjacent vertices in the cell (row number); 
%                                   without boundary vertices (for open meshes)
%                                   (parameter data from function adjVertex)
%           vertexFaceMap (cell format): centre vertex (row number of cell corresponds to centre vertex) and 
%                                       the row numbers of the faces containing adjacent vertices (in cell);
%                                       without boundary vertices (for open meshes)
%                                       (parameter data from function adjVertex)
%           verticesNrROI: numbering of vertices in ROI (region of interest)

% Output:   vertexAdjMap(cell format): Sorted adjacent faces around centre vertice (circle); without centre vertice
%                                       Column 1-2: row numbers of the vertices of the faces; without centre vertice
%                                       Column 3: row number of the faces = row number of the normals
%                                       without boundary vertices (for open meshes)
%           vertexComPoint(cell format): second common point next to the central vertice of the adjacent faces 
%                                       sorted as in vertexAdjMap (circle)
%                                       without boundary vertices (for open meshes)

% Developed by C.Micheler, 
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [vertexAdjMap,vertexComPoint] = adjVertexSeq(boundVerticesNr,vertexMap,vertexFaceMap,verticesNrROI)

    verticesCount = size(verticesNrROI,1);
    
    % Preallocation cell array (cells with different size)
    vertexAdjMap = cell(verticesCount,1); % Sorted faces structure of the centre vertex (circle of adjacent faces)
    vertexComPoint = cell(verticesCount,1); % Second common point (sorted)
    
    for i = 1:verticesCount
    
        % Skip if vertex is a boundary vertex
        if ismember(verticesNrROI(i,1), boundVerticesNr) % for open mesh
            continue;
        end
    
        % Add row number of the faces / row number of the normals
        % Column 1-2: row numbers of the adjacent vertices;
        % Column 3: row number of the faces = row number of the normals
        currVertexFaces = vertexMap{i,1};
        currVertexFaces(:,3) = vertexFaceMap{i,1};
        % Initialization
        sortedFaces = zeros(size(currVertexFaces));
        comPoints = zeros(size(currVertexFaces,1),1);
    
        % Start with the first face
        sortedFaces(1,:) = currVertexFaces(1,:);
        comPoints(1) = currVertexFaces(1,1);
        nextVertex = currVertexFaces(1,1);

        % Use a logical array to keep track of remaining faces
        remaining = true(size(currVertexFaces,1),1);
        remaining(1) = false;  % First face is already used

        % Reorder remaining faces based on adjacency
        for j = 2:size(currVertexFaces,1)
            nextRow = find((currVertexFaces(:,1) == nextVertex | currVertexFaces(:,2) == nextVertex) & remaining, 1);
            if isempty(nextRow)
                break; % Break the loop if no next face is found
            end
            sortedFaces(j,:) = currVertexFaces(nextRow,:);
            nextVertex = setdiff(currVertexFaces(nextRow, 1:2), nextVertex);
            comPoints(j,1) = nextVertex;
            remaining(nextRow) = false;  % Mark as used
        end

        % Store
        vertexAdjMap{i,1} = sortedFaces;
        vertexComPoint{i,1} = comPoints;
    end
    
    disp('Vertex face map sorted (sequence)');

end