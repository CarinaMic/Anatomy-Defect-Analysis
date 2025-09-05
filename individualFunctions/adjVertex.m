% Identify adjacent vertices of a centre vertex (vertex face map); 
% For open and closed meshes

% Input:    boundVerticesNr (for open meshes): row number of boundary vertices 
%                                               (parameter data from function boundSTL)
%           verticesNrROI: numbering of vertices in ROI (region of interest)
%           facesROI: faces structure of the vertices in ROI

% Output:   allVertexFaceMap (cell format): centre vertex (row number of cell corresponds to centre vertex) and 
%                                           the row numbers of the faces containing adjacent vertices (in cell)
%           vertexFaceMap (cell format): centre vertex (row number of cell corresponds to centre vertex) and 
%                                       the row numbers of the faces containing adjacent vertices (in cell);
%                                       without boundary vertices (for open meshes)
%           vertexMap (cell format): centre vertex (row number of cell corresponds to centre vertex) and 
%                                   its adjacent vertices in the cell (row number); 
%                                   without boundary vertices (for open meshes)

% Developed by C.Micheler, 
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function [allVertexFaceMap,vertexFaceMap,vertexMap] = adjVertex(boundVerticesNr,verticesNrROI,facesROI)

    % Precompute mapping from vertices to faces
    allVertexFaceMap = cell(max(facesROI(:)),1);
    for f = 1:size(facesROI,1)
        for v = facesROI(f,:)
            allVertexFaceMap{v} = [allVertexFaceMap{v},f]; 
        end
    end
    
    % Preallocation cell array (cells with different size)
    verticesCount = size(verticesNrROI,1);
    vertexFaceMap = cell(verticesCount,1);
    vertexMap = cell(verticesCount,1);
    
    % Map faces to vertices
    for i = 1:verticesCount
        % Face numbers (row) with corresponding vertice in it
        % Cell row = row number in verticesNrROI -> corresponding vertice
        vertex = verticesNrROI(i,1);
        if ismember(vertex, boundVerticesNr) % for open mesh 
            % Skip if vertex is a boundary vertex
            continue;
        end
        vertexFaceMap{i} = allVertexFaceMap{vertex};
        % Faces with the vertices number
        vertexFaces = facesROI(vertexFaceMap{i},:)';
        % Adjacent Vertices (Remove the central vertice)
        vertexFaces(vertexFaces == vertex) = [];
        vertexMap{i} = reshape(vertexFaces,2,[]).';
    end
    
    disp('Vertex face map determined');

end