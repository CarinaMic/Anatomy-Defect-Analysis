%% Intersection of two anatomical voluem models

% Input:    typeName: 'alpha' | 'refinedAlpha' | 'inside' | 'grid'
%           pointsNum: index of point set (logging)
%           volNum: index of volume (logging)
%           inputIdx: global indices corresponding to rows of `inputPointBase`
%           inputNums: counts [#defect, #vertices, #points] within `inputIdx`
%           inputPointBase: concatenated base array [defect; vertices; points]
%           inputBaseNums: counts [#defect, #vertices, #points] within `inputPointBase`
%           inputVolume.*: struct with fields:
%               - faces    triangulated surface connectivity
%               - vertices vertex coordinates
%           inputMaskNums: Total number [Defect, Vertices, Points] (only for 'inside'/'grid')

% Output:   inter.(typeName).pairs.*: inside/outside masks; indices into global space, subset of inputIdx

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [intersect] = intersectPairs(intersect,typeName,pointsNum,volNum,inputIdx,inputNums,inputPointBase,inputBaseNums,inputVolume,inputMaskNums)

switch typeName
    case {'alpha', 'refinedAlpha'}
        inputAll = inputPointBase(inputIdx,:); % All
    case {'inside', 'grid'}
        defectBase = inputPointBase(1:inputBaseNums(1),:);
        verticesBase = inputPointBase((inputBaseNums(1)+1):(inputBaseNums(1)+inputBaseNums(2)), :);
        pointsBase = inputPointBase((inputBaseNums(1)+inputBaseNums(2)+1):end, :);

        defectIdx = inputIdx(1:inputNums(1),:);
        verticesIdx = inputIdx((inputNums(1)+1):(inputNums(1)+inputNums(2)), :);
        pointsIdx = inputIdx((inputNums(1)+inputNums(2)+1):end, :);

        inputDefect = defectBase(defectIdx,:);
        inputVertices = verticesBase(verticesIdx,:);
        inputPoints = pointsBase(pointsIdx,:);

        inputAll = [inputDefect; inputVertices; inputPoints];
end

% Inpolyhedron: points inside mesh
% https://de.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume
insideVol = inpolyhedron(inputVolume.faces,inputVolume.vertices,inputAll);

% Plot for control
figure;
hold on;
grid on;
axis equal;
title('Inside and Outside Points Test');
xlabel('X');
ylabel('Y');
zlabel('Z');
% Plot the volume (triangulated surface)
trisurf(inputVolume.faces, inputVolume.vertices(:,1), inputVolume.vertices(:,2), inputVolume.vertices(:,3), ...
    'FaceColor', 'cyan', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% Plot inside points (in blue)
scatter3(inputAll(insideVol,1), inputAll(insideVol,2), inputAll(insideVol,3), ...
    50, 'b', 'filled', 'DisplayName', 'Inside Points');
% Plot outside points (in red)
scatter3(inputAll(~insideVol,1), inputAll(~insideVol,2), inputAll(~insideVol,3), ...
    50, 'r', 'filled', 'DisplayName', 'Outside Points');
view(3)

% Extract the points that are inside
switch typeName
    case {'alpha', 'refinedAlpha'}
        % Logical mask base:
        % alpha: defect Mesh
        % refinedAlpha: defect Mesh + refined patches
        inter.(typeName).pairs.defectInsideMaskIdx = inputIdx(insideVol,1); % idx
        % Logical mask                                                             
        inter.(typeName).pairs.defectInsideMask = false(inputMaskNums,1);
        inter.(typeName).pairs.defectInsideMask(inter.(typeName).pairs.defectInsideMaskIdx) = true;
    case {'inside', 'grid'}
        % Logical mask base:
        % inside: defect Mesh + refined patches - reference pelvis vertices - grid points
        isDefect = insideVol(1:inputNums(1),1);
        isVertices = insideVol((inputNums(1)+1):(inputNums(1)+inputNums(2)), 1);
        isPoints = insideVol((inputNums(1)+inputNums(2)+1):end, 1);
        inputMaskDefect = inputIdx(1:inputNums(1),1);
        inputMaskVertices = inputIdx((inputNums(1)+1):(inputNums(1)+inputNums(2)), 1);
        inputMaskPoints = inputIdx((inputNums(1)+inputNums(2)+1):end, 1);
        inter.(typeName).pairs.defectInsideMaskIdx = inputMaskDefect(isDefect,1); % idx
        inter.(typeName).pairs.verticesInsideMaskIdx = inputMaskVertices(isVertices,1);
        inter.(typeName).pairs.pointsInsideMaskIdx = inputMaskPoints(isPoints,1);
        % Logical mask                                                           
        inter.(typeName).pairs.defectInsideMask = false(inputMaskNums(1),1);
        inter.(typeName).pairs.verticesInsideMask = false(inputMaskNums(2),1);
        inter.(typeName).pairs.pointsInsideMask = false(inputMaskNums(3),1);
        inter.(typeName).pairs.defectInsideMask(inter.(typeName).pairs.defectInsideMaskIdx) = true;
        inter.(typeName).pairs.verticesInsideMask(inter.(typeName).pairs.verticesInsideMaskIdx) = true;
        inter.(typeName).pairs.pointsInsideMask(inter.(typeName).pairs.pointsInsideMaskIdx) = true;
end


% Points
pointsInside = inputAll(insideVol,:);                                            
inter.(typeName).pairs.inside = pointsInside;

disp(['intersection "', (typeName), '" defect: pelvis pair ', num2str(pointsNum), '-', num2str(volNum)]);

end
