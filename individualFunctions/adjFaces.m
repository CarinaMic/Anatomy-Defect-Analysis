% Adjacent faces: Find faces that share two vertices with the current face
% Function for (1) normal and (2) mean curvature
% For open and closed meshes

% Input:    boundFacesNr (for open meshes): faces structure of boundary faces
%                                           (parameter data from function boundSTL)
%           allVertexFaceMap (cell format): centre vertex (row number of cell corresponds to centre vertex) and 
%                                           the row numbers of the faces
%                                           containing adjacent vertices (in cell)
%                                           with boundary vertices
%                                           (parameter data from function adjVertexSeq)
%           facesROI: faces structure of the vertices in ROI
%           verticesROI: coordinates of the vertices in ROI 

% Output:   faceAdjMap: row number of faces which are adjacent to the certain face
%                       row number of faceAdjMap is the centre face
%                       3 adjacent faces to every face
%           edgeAdjMap: norm of the common edges of the adjacent faces
%                       common edges in one row
%           sumAreaFaces: sum of the area of the surrounding faces around the centre face
%           areaFaceMap: area which is one third of the area of the triangles around the centre face

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function [faceAdjMap,edgeAdjMap,sumAreaFaces,areaFaceMap] = adjFaces(boundFacesNr,allVertexFaceMap,facesROI,verticesROI)
    
    faceAdjMap = zeros(size(facesROI,1),3);
    edgeAdjMap = zeros(size(facesROI,1),3);
    areaFace = zeros(size(facesROI,1),1);
    crossProduct = zeros(size(facesROI,1),3);
    facesNrROI = (1:size(facesROI,1))';
    
    for i = 1:size(facesROI,1)
    
        % Skip if face is a boundary face
        if ismember(facesNrROI(i,1), boundFacesNr) % for open mesh
            continue;
        end
    
        currentRow = facesROI(i,:);
        % 3 cells for the 3 vertices of currentRow (one cell per vertex -> cell: faceNr which contains the vertex)
        vertexCells = allVertexFaceMap(currentRow); 
        % Concatenate cell contents
        allNumbers = [vertexCells{:}];
        % Count occurrences
        uniqueNumbers = unique(allNumbers);
        counts = histc(allNumbers, uniqueNumbers);
        % Adjacent faces: two common vertices -> two occurrences of the FaceNr
        faceNrTwice = uniqueNumbers(counts == 2);
        % Store adjacent faces
        % row N -> 3 facesNrs of adjacent faces to faceNr N
        faceAdjMap(i,:) = faceNrTwice;
    
        % Map of common edges (norm) of adjacent faces
        tempEdges = zeros(3,3);
        for j = 1:3
            allNumbersE = [currentRow facesROI(faceAdjMap(i,j),:)];
            uniqueNumbersE = unique(allNumbersE);
            countsE = histc(allNumbersE, uniqueNumbersE);
            % Common edge: two common vertices -> two occurrences of the verticeNr
            verticeNrTwice = uniqueNumbersE(countsE == 2);
            tempEdges(j,:) = verticesROI(verticeNrTwice(2),:) - verticesROI(verticeNrTwice(1),:);
            edgeAdjMap(i,j) = vecnorm(tempEdges(j,:));
        end
        crossProduct(i,:) = cross(tempEdges(1,:),tempEdges(2,:));
        areaFace(i,1) = 0.5*vecnorm(crossProduct(i,:),2,2);
    end
    
    % Sum of the adjacent faces
    for i = 1:size(facesROI, 1)
        if ismember(facesNrROI(i, 1), boundFacesNr)
            sumAreaFaces(i, 1) = 0; % boundary faces (for open mesh)
        else
            sumAreaFaces(i, 1) = sum(areaFace(faceAdjMap(i, :)));
        end
    end
    
    areaFaceMap = sumAreaFaces / 3;

    disp('Adjacent faces calculated');

end

