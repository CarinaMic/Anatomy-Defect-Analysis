% Scaling mesh (stl)

% Input:    importData with ConnectivityList and Points
%           scaleFactor

% Output:   scaledVertices: scaled vertices
%           scaledFaces: connectivity list of scaled vertices
%           scaledNormals: normals of the scaled faces
%           scaledCentreFaces: centre of the scaled faces
%           scaledLandmarks: scaled landmarks

% Developed by C.Micheler, 
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function [scaledVertices,scaledFaces,scaledNormals,scaledCentreFaces,scaledLandmarks] = scaleSTL(import,scaleFactor)

% Scale Mesh: https://de.mathworks.com/help/driving/ref/extendedobjectmesh.scale.html
% Attention: Scaling leads to new position of the geometry
mesh = extendedObjectMesh(import.vertices,import.faces);
scaledMesh = scale(mesh,scaleFactor);
scaledVertices = scaledMesh.Vertices;
scaledFaces = scaledMesh.Faces;

% List of all landmarks
allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', ...
    'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};

% Scale available landmarks (for transformation)
scaledLandmarks = struct();
for idx = 1:length(allLandmarks)
    landmark = allLandmarks{idx};
    if isfield(import, landmark)
        fieldName = ['scaled', landmark];
        scaledLandmarks.(fieldName) = import.(landmark) * scaleFactor;
    end
end

% Normals
scaledNormals = import.normals;

% Calculate the centroid of the triangular faces in the mesh
V1 = scaledMesh.Vertices(scaledMesh.Faces(:,1),:);
V2 = scaledMesh.Vertices(scaledMesh.Faces(:,2),:);
V3 = scaledMesh.Vertices(scaledMesh.Faces(:,3),:);
faceVertices = cat(3, V1, V2, V3);
scaledCentreFaces = mean(faceVertices, 3);

disp('model scaled')

end