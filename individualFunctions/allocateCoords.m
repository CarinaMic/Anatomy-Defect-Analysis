%% Function to allocate coordinates to defects, vertices, and points

% Input:    intersect: struct to be updated (output container)
%           baseIntersect: label used for logging
%           typeName: 'alpha' | 'refinedAlpha' | 'inside' | 'grid'
%           inputIDs: IDs for the selected categorie
%           pelvisDefect: struct array providing required data, fields used:
%                       - .transform.trafo.(weightKabsch).vertices
%                       - .transform.trafo.(weightKabsch).faces
%                       - .volume.patches.all.comVertices / comFaces
%                       - .volume.shrink.inside.mainClusterVertices / mainClusterPoints
%                       - .volume.gridPoints.inside                (for 'grid')
%           refData: struct with reference pelvis mesh:
%                       - .vertices, .faces 
%           pointBaseGrid: grid points (for 'grid')
%           weightKabsch: fieldname selecting the transform in pelvisDefect(u).transform.trafo.*
%           currentPct: struct carrying precomputed filters from intersectCount, fields used:
%                       - .filteredFreqDefect{i}.PointBaseIndex
%                       - .filteredFreqVertices{i}.PointBaseIndex   (for 'inside'/'grid')
%                       - .filteredFreqPoints{i}.PointBaseIndex     (for 'inside'/'grid')
%           currentPctName: target fieldname under inter.(typeName).*

% Output:   intersect : updated struct with inter.(typeName).(currentPctName) 

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [intersect] = allocateCoords(intersect, baseIntersect, typeName, inputIDs, pelvisDefect, refData, ...
    pointBaseGrid, weightKabsch, currentPct, currentPctName)

numIDs = length(inputIDs); % Number of point bases

% Preallocate storage
currentPct.insideDefect = cell(numIDs, 1);
currentPct.usedDefect = zeros(numIDs, 1);
currentPct.usedOrigDefect = zeros(numIDs, 1);
currentPct.lengthInsideDefect = zeros(numIDs, 1);
currentPct.facesDefectInside = cell(numIDs, 1);
if ismember(typeName, {'inside', 'grid'})
    currentPct.insideVertices = cell(numIDs, 1);
    currentPct.usedVertices = zeros(numIDs, 1);
    currentPct.lengthInsideVertices = zeros(numIDs, 1);
    currentPct.facesVerticesInside = cell(numIDs, 1);
    currentPct.insidePoints = cell(numIDs, 1);
    currentPct.usedPoints = zeros(numIDs, 1);
    currentPct.lengthInsidePoints = zeros(numIDs, 1);
end

% Iterate through each point base
for i = 1:numIDs
    u = inputIDs(i); % Pelvis defect ID
    % Extract base elements
    switch typeName
        case 'alpha'
            pointBaseDefect = pelvisDefect(u).transform.trafo.(weightKabsch).vertices;
            pointBaseOrigDefect = pelvisDefect(u).transform.trafo.(weightKabsch).vertices;
            facesBaseDefect = pelvisDefect(u).transform.trafo.(weightKabsch).faces;
        case 'refinedAlpha'
            pointBaseDefect = pelvisDefect(u).volume.patches.all.comVertices;
            pointBaseOrigDefect = pelvisDefect(u).transform.trafo.(weightKabsch).vertices;
            facesBaseDefect = pelvisDefect(u).volume.patches.all.comFaces;
        case {'inside', 'grid'}
            pointBaseDefect = pelvisDefect(u).volume.patches.all.comVertices;
            pointBaseOrigDefect = pelvisDefect(u).transform.trafo.(weightKabsch).vertices;
            facesBaseDefect = pelvisDefect(u).volume.patches.all.comFaces;
            pointBaseVertices = refData.vertices;
            facesBaseVertices = refData.faces;
            pointBasePoints = pointBaseGrid;

            pointInputVertices = pelvisDefect(u).volume.shrink.inside.mainClusterVertices;
            if ismember(typeName, {'inside'})
                pointInputPoints = pelvisDefect(u).volume.shrink.inside.mainClusterPoints;
            else % grid
                pointInputPoints = pelvisDefect(u).volume.gridPoints.inside;
            end
    end

    % Defect allocation
    pointBaseDefectIdx = currentPct.filteredFreqDefect{i}.PointBaseIndex;
    insideDefect = pointBaseDefect(pointBaseDefectIdx, :);
    currentPct.insideDefect{i} = insideDefect;
    currentPct.usedDefect(i) = length(pointBaseDefectIdx) / length(pointBaseDefect);
    currentPct.lengthInsideDefect(i) = length(pointBaseDefectIdx);
    validOrigIdx = pointBaseDefectIdx(pointBaseDefectIdx <= length(pointBaseOrigDefect)); % without patches
    currentPct.usedOrigDefect(i,1) = length(validOrigIdx) / length(pointBaseOrigDefect);

    % Find defect faces
    verticesFacesDefect = false(length(pointBaseDefect), 1);
    verticesFacesDefect(pointBaseDefectIdx) = true;
    facesDefectInside = all(verticesFacesDefect(facesBaseDefect), 2); % Faces with all vertex indices inside
    currentPct.facesDefectInside{i} = facesBaseDefect(facesDefectInside, :); % Faces structure for points base structure (vertices)

    % Vertices and points allocation (for 'inside' and 'grid')
    if ismember(typeName, {'inside', 'grid'})
        % Vertices
        pointBaseVerticesIdx = currentPct.filteredFreqVertices{i}.PointBaseIndex;
        insideVertices = pointBaseVertices(pointBaseVerticesIdx, :);
        currentPct.insideVertices{i} = insideVertices;
        currentPct.usedVertices(i) = length(pointBaseVerticesIdx) / length(pointInputVertices);
        currentPct.lengthInsideVertices(i) = length(pointBaseVerticesIdx);

        % Find vertices faces
        verticesFacesVer = false(length(pointBaseVertices), 1);
        verticesFacesVer(pointBaseVerticesIdx) = true;
        facesVerticesInside = all(verticesFacesVer(facesBaseVertices), 2);
        currentPct.facesVerticesInside{i} = facesBaseVertices(facesVerticesInside, :); % Faces structure for points base structure (vertices)

        % Points
        pointBasePointsIdx = currentPct.filteredFreqPoints{i}.PointBaseIndex;
        insidePoints = pointBasePoints(pointBasePointsIdx, :);
        currentPct.insidePoints{i} = insidePoints;
        currentPct.usedPoints(i) = length(pointBasePointsIdx) / length(pointInputPoints);
        currentPct.lengthInsidePoints(i) = length(pointBasePointsIdx);
    end
end

% Aggregate statistics
currentPct.sumInsideDefect = sum(currentPct.lengthInsideDefect);
currentPct.usedDefectMean = mean(currentPct.usedDefect);
currentPct.usedDefectStd = std(currentPct.usedDefect);
currentPct.usedOrigDefectMean = mean(currentPct.usedOrigDefect);
currentPct.usedOrigDefectStd = std(currentPct.usedOrigDefect);

if ismember(typeName, {'inside', 'grid'})
    currentPct.sumInsideVertices = sum(currentPct.lengthInsideVertices);
    currentPct.usedVerticesMean = mean(currentPct.usedVertices);
    currentPct.usedVerticesStd = std(currentPct.usedVertices);
    currentPct.sumInsidePoints = sum(currentPct.lengthInsidePoints);
    currentPct.usedPointsMean = mean(currentPct.usedPoints);
    currentPct.usedPointsStd = std(currentPct.usedPoints);
end

% Return updated object
inter.(typeName).(currentPctName) = currentPct;
% Data storage of intersection data
sizePct = whos('currentPct');
inter.(typeName).(currentPctName).sizeIntersectAllo = sizePct.bytes;

disp(['defect intersection "', (typeName), '": coordinates for "', (baseIntersect), '" and "', (currentPctName), '"']);
end
