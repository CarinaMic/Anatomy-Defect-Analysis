% Transformation (Kabsch): Translation and Rotation / Orientation (incl. landmark error)

% Input:    refLandmarks: reference landmarks for transformation
%           weight: weighting of the landmarks for transformation
%           importData with ConnectivityList and Points

% Output:   trafoR: rotation matrix
%           trafoT: translation vector
%           trafoLRMS: least root mean square distance between two paired sets
%           trafoVertices: transformed vertices
%           trafoFaces: transformed faces
%           trafoNormals: transformed normals
%           trafoCentreFaces: centre of transformed faces
%           trafoLandmarks: transformed landmarks (incl. landmark errror)
%           trafoEuler: euler angle or rotation matrix

% Developed by C.Micheler, 
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [trafoR,trafoT,trafoLRMS,trafoVertices,trafoFaces,trafoNormals,trafoCentreFaces,trafoLandmarks,trafoEuler] = ...
    trafoKabsch(refLandmarks, weight, importData)

allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', ...
    'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};

% Collect landmarks
Landmarks = zeros(length(weight), 3); % weight defines number of landmarks
for idx = 1:length(weight)
    if isfield(importData, allLandmarks{idx})
        Landmarks(idx, :) = importData.(allLandmarks{idx});
    else
        Landmarks(idx, :) = NaN;
    end
end
% Filter out rows with NaN values
validRows = ~any(isnan(Landmarks), 2);
validIndices = find(validRows); % Get the valid indices (those not NaN)
Landmarks = Landmarks(validRows, :);
refLandmarks = refLandmarks(validRows, :);
weight = weight(validRows);

% Kabsch
[trafoR, trafoT, trafoLRMS] = ...
    Kabsch(Landmarks', refLandmarks', weight);
rotated_Kabsch = trafoR * Landmarks' + trafoT;

% Euler angles (z,y,x)
trafoEuler = rad2deg(rotm2eul(trafoR));

% Assign to obj.rotation using the valid landmarks list
for k = 1:length(validIndices)
    idx = validIndices(k);  % Retrieve the actual landmark index
    trafoLandmarks.(allLandmarks{idx}) = rotated_Kabsch(:, k)';  % Assign the rotated value 
end
% For the landmarks that weren't valid, assign NaN
invalidIndices = setdiff(1:length(allLandmarks), validIndices);
for idx = invalidIndices
    trafoLandmarks.(allLandmarks{idx}) = NaN(1, 3);
end

% Transform vertices
rotatedVertices = trafoR * importData.vertices' + trafoT;
trafoVertices = rotatedVertices';
% Faces / Normals
trafoFaces = importData.faces;
trafoNormals = importData.normals;
% Calculate the centroid of the triangular faces in the mesh
V1 = trafoVertices(importData.faces(:,1),:);
V2 = trafoVertices(importData.faces(:,2),:);
V3 = trafoVertices(importData.faces(:,3),:);
faceVertices = cat(3, V1, V2, V3);
trafoCentreFaces = mean(faceVertices, 3);

% Landmark errors
for idx = 1:length(allLandmarks)
    landmark = allLandmarks{idx};
    if any(validIndices == idx)
        % Calculate Euclidean distance between transformed and reference landmarks
        transformedLandmark = trafoLandmarks.(landmark); 
        referenceLandmark = refLandmarks(validIndices == idx, :);
        trafoLandmarks.(['error_' landmark]) = norm(transformedLandmark - referenceLandmark); 
    else
        % Landmark was filtered out, assign NaN and note it
        trafoLandmarks.(['error_' landmark]) = NaN; 
    end
end

disp('model transformed');
end