%% Hausdorff Distance

% Input:    verticesRef: 'neighbours' (reference vertices)
%           vertices: vertices to search for distances in verticesRef

% Output:   hausdorff: Hausdorff distances

% Developed by C.Micheler, 
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [hausdorffDist] = hausdorffDistance(verticesRef, vertices)

% Hausdorff-Distance calculates the Hausdorff distance between two point clouds.
% This is done by calculating the maximum of the minimum distances from each vertex in
% vertices to the nearest vertex in verticesRef.

% Create KD-tree searchers for both vertex sets
MdlRef = KDTreeSearcher(verticesRef);
MdlTest = KDTreeSearcher(vertices);

% First, calculate distance from each vertex in 'vertices' to its nearest neighbor in 'verticesRef'
idxRef = knnsearch(MdlRef, vertices);
distancesRef = sqrt(sum((vertices - verticesRef(idxRef,:)).^2, 2));

% Second, calculate distance from each vertex in 'verticesRef' to its nearest neighbor in 'vertices'
idxTest = knnsearch(MdlTest, verticesRef);
distancesTest = sqrt(sum((verticesRef - vertices(idxTest,:)).^2, 2));

% Calculate the Hausdorff Distance
% Maximum of the minimal distances in both directions
hausdorffDist = max(max(distancesRef), max(distancesTest));

% Store results in obj structure
hausdorffDist = hausdorffDist;

disp('Hausdorff Distance calculated');

end