% Calculate vector of every edge in ROI and the area of the faces
% Function for mean (2) and gauss (3) curvature
% For open and closed meshes

% Input:    boundVerticesNr (for open meshes): row number of boudary vertices 
%                                               (parameter data from function boundSTL)
%           vertexAdjMap (cell format): Sorted adjacent faces around centre vertice (circle); without centre vertice
%                                        Column 1-2: row numbers of the vertices of the faces; without centre vertice
%                                        Column 3: row number of the faces = row number of the normals
%                                       (parameter data from function adjVertexSeq)
%           vertexComPoint (cell format): second common point next to the central vertice of the adjacent faces 
%                                           sorted as in vertexAdjMap (circle)
%                                           (parameter data from function adjVertexSeq)
%           verticesNrROI: numbering of the vertices in ROI (region of interest)
%           verticesROI: coordinates of the vertices in ROI 

% Output:   comEdge (cell format): coordinates of the common edges of the adjacent faces (coordinates)
%                                   sorted as in vertexAdjMap (circle)
%                                   without boundary vertices (for open meshes)
%           comEdgeNorm (cell format): norms of the common edges of the adjacent faces (coordinates)
%                                       sorted as in vertexAdjMap (circle)
%                                       without boundary vertices (for open meshes)
%           sumAreaFaces: sum of the area of the surrounding faces around the centre vertice
%                           without boundary vertices (for open meshes)
%           baryAreaFaces: barycentric area which is one third of the area of the triangles around the centre vertice
%                           without boundary vertices (for open meshes)

% Developed by C.Micheler, 
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [comEdge,comEdgeNorm,sumAreaFaces,baryAreaFaces]= edgesArea(boundVerticesNr,vertexAdjMap,vertexComPoint,verticesNrROI,verticesROI)

    verticesCount = size(verticesNrROI, 1);
    
    % Coordinates of the centre vertice
    vertices = verticesROI(verticesNrROI(:, 1), :);
    
    % Preallocate arrays
    % Coordinates of the adjacent vertices
    adjVertices = cell(verticesCount,1);
    % Vector of the edges
    edges = cell(verticesCount,1);
    % Common edge
    comEdge = cell(verticesCount,1);
    comEdgeNorm = cell(verticesCount,1);
    % Calculate the area of the faces
    crossProduct = cell(verticesCount,1);
    areaFaces = cell(verticesCount,1);
    sumAreaFaces = zeros(verticesCount,1);
    
    for i = 1:verticesCount
    
        % Skip if vertex is a boundary vertex
        if ismember(verticesNrROI(i,1), boundVerticesNr) % for open mesh
            continue;
        end
    
        % Coordinates of the adjacent vertices
        % Vertice for edge1 of the face
        adjVertices{i,1}(:,:,1) = verticesROI(vertexAdjMap{i,1}(:,1),:);
        % Vertice for edge2 of the face
        adjVertices{i,1}(:,:,2) = verticesROI(vertexAdjMap{i,1}(:,2),:);
        % Calculate vector of the adjacent edges
        % edge{:,:}(:,:,1) -> edge1 of the face
        % edge{:,:}(:,:,2) -> edge2 of the face
        edges{i,1} = adjVertices{i,1} - vertices(i,:);
    
        % Common edge and its norm
        % Circle back to the first element
        comEdge{i,1} = verticesROI(vertexComPoint{i,1}(:,1),:) - vertices(i,:);
        comEdgeNorm{i,1} = vecnorm(comEdge{i,1},2,2);
    
        % Calculate the area of the faces
        % Norm of the cross product of two vectors correspond to the area of the parallelogram
        % -> The area of the triangle is the half
        % |a x b| = |a||b|sin(phi)
        % normal vector = cross product
        % normals of the stl-file: length of the normals normalized to length 1
        % -> new calculation with the original edges
        crossProduct{i,1} = cross(edges{i,1}(:,:,1),edges{i,1}(:,:,2));
        areaFaces{i,1} = 0.5*vecnorm(crossProduct{i,1},2,2);
        sumAreaFaces(i,1) = sum(areaFaces{i,1});
    end
    
    baryAreaFaces = sumAreaFaces / 3;
    
    disp('Edges and face areas calculated');

end