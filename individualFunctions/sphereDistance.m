%% Vertex-to-Nearest-Neighbour Distance

% Input:    typeName: 'alpha' | 'refinedAlpha' | 'inside' | 'grid'
%           approx: struct to be updated (output container)
%           i: Numeric identifier used only for logging 
%           defIdx: Row index under which results are stored for this defect
%           verticesDefect: defect mesh vertices
%           facesDefect: defect mesh triangular connectivity 
%           bestModelType: suffix used to name fields (e.g., 'P50', 'P75').

% Output:   approx: updated struct

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function approx = sphereDistance(approx, i, defIdx, verticesDefect, facesDefect, typeName, bestModelType)

% Spheres
spherePts = approx.(typeName).bestModelDist.(['sphere' bestModelType 'Mesh']);
bestModelNr = approx.(typeName).bestModel.(['score' bestModelType 'ModelNr']);
bestModelPct = approx.(typeName).bestModel.(['score' bestModelType 'Pct']);
scaledSpheres = approx.(typeName).(bestModelPct).modelScale(bestModelNr).scaledSpheres; % scaled spheres

% Nearest Neighbour
Mdl = KDTreeSearcher(spherePts);
[~, dSurface] = knnsearch(Mdl, verticesDefect);

% Points in sphere model
centers = scaledSpheres(1:3,:)';
radii   = scaledSpheres(4,:)';
D = pdist2(verticesDefect, centers);
insideMask = any(D <= radii', 2);      % Points inside sphere model?
%Inside = 0, Outside = vertex-to-nearest neighbour distance
minDistances = dSurface;
minDistances(insideMask) = 0;

% Colourised nearest neighbour
nearVertexMatrix = minDistances(facesDefect);
faceMean = mean(nearVertexMatrix, 2);
minVal = min(faceMean); maxVal = max(faceMean);
normVals = (faceMean - minVal) / (maxVal - minVal);
normVals = max(0, min(1, normVals));
colourNum = 32768;
colourMap = viridis(colourNum);
colourIdx = round(normVals * (colourNum-1)) + 1;
rgbColour = colourMap(colourIdx,:);

% Save
approx.(typeName).bestModelDist.(['nearDist' bestModelType]){defIdx,1} = minDistances;
approx.(typeName).bestModelDist.(['nearDist' bestModelType 'Face']){defIdx,1} = faceMean;
approx.(typeName).bestModelDist.(['nearDist' bestModelType 'Norm']){defIdx,1} = normVals;
approx.(typeName).bestModelDist.(['nearDist' bestModelType 'FaceColour']){defIdx,1} = rgbColour;
approx.(typeName).bestModelDist.(['nearDist' bestModelType 'Mean'])(defIdx,1) = mean(minDistances);
approx.(typeName).bestModelDist.(['nearDist' bestModelType 'Max'])(defIdx,1) = max(minDistances);
approx.(typeName).bestModelDist.(['nearDist' bestModelType 'Std'])(defIdx,1) = std(minDistances);

disp(['Vertex-to-Nearest-Neighbour Distance: i=' num2str(i) ', typeName=' typeName ', defectNr=' num2str(defIdx)]);

end