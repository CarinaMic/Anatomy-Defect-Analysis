% Translation to reference point

% Input:    importData with ConnectivityList and Points
%           refPoint: reference Point for translation

% Output:   transVector: translation vector
%           transVertices: shifted vertices
%           transFaces: connectivity list of shifted vertices
%           transNormals: normals of the shifted faces
%           transCentreFaces: centre of the shifted faces
%           transLandmarks: shifted landmarks

% Developed by C.Micheler, 
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [transVector,transVertices,transFaces,transNormals,transCentreFaces,transLandmarks] = transSTL(refPoint,importData)

% Translation
transVector = refPoint - importData.acentre;
transVertices = importData.vertices + transVector;

% List of all landmarks
allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre',...
    'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};

% Translate available landmarks
transLandmarks = struct();
for idx = 1:length(allLandmarks)
    landmark = allLandmarks{idx};
    if isfield(importData, landmark)
        if strcmp(landmark, 'acentre')
            transLandmarks.acentre = refPoint; 
        else
            fieldName = ['trans', landmark];
            transLandmarks.(fieldName) = importData.(landmark) + transVector;
        end
    end
end

% Faces / Normals
transFaces = importData.faces;
transNormals = importData.normals;
% Calculate the centroid of the triangular faces in the mesh
V1 = transVertices(importData.faces(:,1),:);
V2 = transVertices(importData.faces(:,2),:);
V3 = transVertices(importData.faces(:,3),:);
faceVertices = cat(3, V1, V2, V3);
transCentreFaces = mean(faceVertices, 3);

disp('vertices shifted')
end