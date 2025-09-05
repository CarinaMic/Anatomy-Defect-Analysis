%% Intersection of boundary box and defect (Pre-filter)

% Input:    typeName: 'alpha' | 'refinedAlpha' | 'inside' | 'grid'
%           pointsNum: Index of point set (logging)
%           volNum: Index of volume (logging)
%           box.*: bounding box
%               - edgeVector  – edge vectors of the box (local COSY)
%               - cornerpoints  – sorted corner points of the box
%           inputPoints: all points to be checked (depending on typeName)
%           inputPointNums: number in inputPoints (for 'inside'/'grid')
%           inputMasksIdx: global indices of the respective points in inputPoints (for 'inside'/'grid')
%           inputMaskNums: Total number [Defect, Vertices, Points] (only for 'inside'/'grid')

% Output:   inter.(typeName).pairs.*: points inside box (pairwise comparison)

% Developed by C.Micheler, 
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [intersect] = intersectBox(intersect,typeName,pointsNum,volNum,box,inputPoints,inputPointNums,inputMasksIdx,inputMaskNums)

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
minBox = min(rotatedBox, [], 1);
maxBox = max(rotatedBox, [], 1);
% Combine the checks for all coordinates
isInsideAll = all(rotatedPoints >= minBox & rotatedPoints <= maxBox, 2); % logical
% Extract the points that are inside the box
switch typeName
    case {'alpha', 'refinedAlpha'}
        % Logical mask base:
        % alpha: defect Mesh
        % refinedAlpha: defect Mesh + refined patches
        inter.(typeName).pairs.boxDefectInsideMask = isInsideAll; % logical mask 
        inter.(typeName).pairs.boxDefectInsideMaskIdx = find(isInsideAll); % idx
    case {'inside', 'grid'}
        % Logical mask base:
        % inside: defect Mesh + refined patches - reference pelvis vertices - grid points
        isDefect = isInsideAll(1:inputPointNums(1),1);
        isVertices = isInsideAll((inputPointNums(1)+1):(inputPointNums(1)+inputPointNums(2)), 1);
        isPoints = isInsideAll((inputPointNums(1)+inputPointNums(2)+1):end, 1);
        inputMaskDefect = inputMasksIdx(1:inputPointNums(1),1);
        inputMaskVertices = inputMasksIdx((inputPointNums(1)+1):(inputPointNums(1)+inputPointNums(2)), 1);
        inputMaskPoints = inputMasksIdx((inputPointNums(1)+inputPointNums(2)+1):end, 1);
        inter.(typeName).pairs.boxDefectInsideMaskIdx = inputMaskDefect(isDefect,1); % idx
        inter.(typeName).pairs.boxVerticesInsideMaskIdx = inputMaskVertices(isVertices,1);
        inter.(typeName).pairs.boxPointsInsideMaskIdx = inputMaskPoints(isPoints,1);
        % Logical mask                                                             
        inter.(typeName).pairs.boxDefectInsideMask = false(inputMaskNums(1),1);
        inter.(typeName).pairs.boxVerticesInsideMask = false(inputMaskNums(2),1);
        inter.(typeName).pairs.boxPointsInsideMask = false(inputMaskNums(3),1);
        inter.(typeName).pairs.boxDefectInsideMask(inter.(typeName).pairs.boxDefectInsideMaskIdx) = true;
        inter.(typeName).pairs.boxVerticesInsideMask(inter.(typeName).pairs.boxVerticesInsideMaskIdx) = true;
        inter.(typeName).pairs.boxPointsInsideMask(inter.(typeName).pairs.boxPointsInsideMaskIdx) = true;
end

% Points in global cosy                                                          
pointsInside = inputPoints(isInsideAll,:);
inter.(typeName).pairs.boxInside = pointsInside;

disp(['intersection "', (typeName), '" defect - box: pelvis pair ', num2str(pointsNum), '-', num2str(volNum)]);

end