%% Fill pelvis defect with points (Pre-filter with boundary box)

% Input:    pelvisNum: Numeric identifier used only for logging
%           inputPoints: Candidate points to be pre-filtered and tested
%           box:  Oriented bounding box with fields edgeVector,cornerpoints
%           meshFaces: Triangular faces of the target mesh
%           meshVertices: Vertex coordinates of the target mesh

% Output:   gridPoints: mask/indices/points inside mesh

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [gridPoints] = fillDefect(pelvisNum,inputPoints,box,meshFaces,meshVertices)

% Bounding box (not aligned to axis)
% Bounding box's axis directions
dir1 = box.edgeVector(1,:)/norm(box.edgeVector(1,:)); % x % One box edge; right-hand-rule
dir2 = box.edgeVector(2,:)/norm(box.edgeVector(2,:)); % y
dir3 = box.edgeVector(3,:)/norm(box.edgeVector(3,:)); % z

% Create the rotation matrix local -> global (box)
rotMatrix = [dir1; dir2; dir3];

% Rotate the input points and the box vertices to the aligned space
rotatedPoints = inputPoints * rotMatrix';
rotatedBox = box.cornerpoints * rotMatrix';
% Check if the point is inside the aligned box
minX = min(rotatedBox(:, 1));
maxX = max(rotatedBox(:, 1));
minY = min(rotatedBox(:, 2));
maxY = max(rotatedBox(:, 2));
minZ = min(rotatedBox(:, 3));
maxZ = max(rotatedBox(:, 3));
% Check if each coordinate of each point is within the bounds
isInsideX = rotatedPoints(:, 1) >= minX & rotatedPoints(:, 1) <= maxX;
isInsideY = rotatedPoints(:, 2) >= minY & rotatedPoints(:, 2) <= maxY;
isInsideZ = rotatedPoints(:, 3) >= minZ & rotatedPoints(:, 3) <= maxZ;
% Combine the checks for all coordinates
isInsideAll = isInsideX & isInsideY & isInsideZ; % logical
gridPoints.boxInsideMask = isInsideAll;
% Extract the points that are inside the box
rotatedPointsInside = rotatedPoints(isInsideAll,:);
gridPoints.boxInsideIdx = find(isInsideAll); % idx

% Transform points to global cosy
pointsInside = rotatedPointsInside * rotMatrix;
gridPoints.boxInside = pointsInside;

% Fill pelvis defects with points
% Check which points are inside the mesh
% Vertices/faces define the refined alpha shape mesh
insideDefect = inpolyhedron(meshFaces, meshVertices, ... % Mesh
    gridPoints.boxInside);
gridPoints.insideMask = false(size(gridPoints.boxInsideMask,1), 1); % Logical mask of point base (gridPoints)
gridPoints.insideMaskIdx = gridPoints.boxInsideIdx(insideDefect,:);
gridPoints.insideMask(gridPoints.insideMaskIdx) = true;
gridPoints.inside = gridPoints.boxInside(insideDefect,:);

disp(['pelvis defect filled with points (pre-filter bounding box): pelvis defect ', num2str(pelvisNum)]);

end