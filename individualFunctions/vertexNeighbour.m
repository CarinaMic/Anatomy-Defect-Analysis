%% Vertex-to-Nearest-Neighbour Distance

% Input:    verticesRef: 'neighbours' (reference vertices)
%           vertices: vertices to search for neighbours in verticesRef
%           faces: connectivity list of vertices

% Output:   nearVertex: vertex-to-nearest-neughbour distance 
%           nearVertexMax: maximum of vertex-to-nearest-neughbour distance 
%           nearVertexMean: mean of vertex-to-nearest-neughbour distance 
%           nearVertexStd: standard deviation of vertex-to-nearest-neughbour distance 
%           nearVertexFace: vertex-to-nearest-neighbour distance averaged on one face (values averaged from three vertices)
%           nearVertexFaceNorm: normalised nearVertexFace for colorised pelvis (colormap viridis)

% Developed by C.Micheler, 
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [nearVertex,nearVertexMax,nearVertexMean,nearVertexStd,nearVertexFace,nearVertexFaceNorm] = ...
    vertexNeighbour(verticesRef, vertices, faces)

% Vertex-to-Nearest-Neighbour calculates the distance from each vertex to its
% nearest neighbor in the vertices list.
% Using a KD-tree searcher for improved performance
% Minimum distances from each vertex to the nearest vertex in reference

% Create a KD-tree searcher for verticesB
MdlB = KDTreeSearcher(verticesRef);

% Find the nearest neighbor in verticesB for each vertex in verticesA
idx = knnsearch(MdlB,vertices);

% Calculate the distances based on the nearest neighbor indices
minDistances = sqrt(sum((vertices - verticesRef(idx,:)).^2, 2));

% Pelvis with vertex-to-nearest-neighbour color
% Retrieve nearVertex for each vertex of each face
nearVertexMatrix = minDistances(faces);
% Calculate mean for each face
nearVertexFace = mean(nearVertexMatrix,2);
% Normalisation
lowerBound = min(nearVertexFace);
upperBound = max(nearVertexFace);
% Normalize data within this range
normalData = (nearVertexFace - lowerBound) / (upperBound - lowerBound);
normalData(normalData < 0) = 0;   % Clamp values below the range to 0
normalData(normalData > 1) = 1;   % Clamp values above the range to 1
% Apply colormap (recommended: viridis)
% 16-bit colour: alpha: 1bit R: 5bit G: 5bit B: 5bit -> 32768 colors
colourNum = 32768;
colourMap = viridis(colourNum);
colourIdx = round(normalData * (colourNum-1)) + 1;
rgbColour = colourMap(colourIdx,:);

nearVertexFaceNorm = rgbColour;
nearVertex = minDistances;
nearVertexMax = max(minDistances);
nearVertexMean = mean(minDistances);
nearVertexStd = std(minDistances);

disp('Vertex-to-Nearest-Neighbour Distance');

end
