% Identify boundary vertices / faces of a mesh (stl)

% Input:    importData with ConnectivityList and Points
%           Points = vertices
%           ConnectivityList = faces structure of the stl

% Output:   boundFacesNr: Row numbers of boundary face
%           boundFaces: Faces structure of boundary faces
%           boundVerticesNr: Row numbers of boundary vertices
%           boundVertices: Coordinates of boundary vertices 

% Developed by C.Micheler, 
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function [boundFaces,boundFacesNr,boundVertices,boundVerticesNr] = boundSTL(importData) % derived from function edging in Volume

    % Identify boundary vertices / faces
    % Boundary edge: triangular side (combination of two vertices) only once in row of faces list
    
    faces = importData.ConnectivityList;
    vertices = importData.Points;
    
    % Extract all edges from the faces
    edges = [faces(:, [1, 2]); faces(:, [2, 3]); faces(:, [3, 1])];
    sortedEdges = sort(edges, 2);
    
    % Identify unique edges and their counts
    [uniqueEdges, ~, edgeIndices] = unique(sortedEdges, 'rows');
    edgeCounts = accumarray(edgeIndices, 1);
    
    % Boundary edges appear only once
    isBoundaryEdge = edgeCounts == 1;
    boundaryEdges = uniqueEdges(isBoundaryEdge, :);
    
    % Precompute face edges for all faces
    faceEdges = [sort(faces(:, [1, 2]), 2), (1:size(faces, 1))'; ...
        sort(faces(:, [2, 3]), 2), (1:size(faces, 1))'; ...
        sort(faces(:, [3, 1]), 2), (1:size(faces, 1))'];

    % Find faces that include boundary edges
    boundaryFaces = ismember(faceEdges(:, 1:2), boundaryEdges, 'rows');
    boundaryFaceIndices = faceEdges(boundaryFaces, 3);

    % Extract boundary faces
    boundFacesNr = unique(boundaryFaceIndices);
    boundFaces = faces(boundFacesNr,:);

    % Store boundary vertices
    boundVerticesNr = unique(boundaryEdges(:)); % Flatten the array and then find unique vertices
    boundVertices = vertices(boundVerticesNr,:);
    
    disp('Boundary vertices and faces identified')

end