%% Comparison of intersection types

% Input:    inter: struct to be updated (output container)
%           type1: first intersection type name ('alpha' | 'refinedAlpha' | 'inside' | 'grid')
%           type2: second intersection type name ('alpha' | 'refinedAlpha' | 'inside' | 'grid')
%           baseName: base/category label
%           pctName: percentage label
%           numIDs: number of pelvis IDs / point bases to compare
%           type1Data: struct holding filtered frequency tables for type1
%           type2Data: struct holding filtered frequency tables for type2 

% Output:   inter: updated struct with intersect.comp with comparison of two types 

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [inter] = compareIntersect(inter, type1, type2, baseName, pctName, numIDs, type1Data, type2Data)

combi = [type1 '_' type2];

% Initialize storage
inter.comp.(combi).(baseName).(pctName).uniqueDefectType1 = cell(numIDs, 1);
inter.comp.(combi).(baseName).(pctName).uniqueDefectType2 = cell(numIDs, 1);
inter.comp.(combi).(baseName).(pctName).sharedDefect = cell(numIDs, 1);
inter.comp.(combi).(baseName).(pctName).lengthUniqueDefectType1 = zeros(numIDs, 1);
inter.comp.(combi).(baseName).(pctName).lengthUniqueDefectType2 = zeros(numIDs, 1);
inter.comp.(combi).(baseName).(pctName).lengthSharedDefect = zeros(numIDs, 1);
if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
    % Vertices
    inter.comp.(combi).(baseName).(pctName).uniqueVerticesType1 = cell(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).uniqueVerticesType2 = cell(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).sharedVertices = cell(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType1 = zeros(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType2 = zeros(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).lengthSharedVertices = zeros(numIDs, 1);
    % Points
    inter.comp.(combi).(baseName).(pctName).uniquePointsType1 = cell(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).uniquePointsType2 = cell(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).sharedPoints = cell(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).lengthUniquePointsType1 = zeros(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).lengthUniquePointsType2 = zeros(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).lengthSharedPoints = zeros(numIDs, 1);
end

% Iterate over PointBases/IDs
for b = 1:numIDs
    % Extract Defect indices
    defect1 = [];
    defect2 = [];
    if isfield(type1Data, 'filteredFreqDefect') && numel(type1Data.filteredFreqDefect) >= b
        defect1 = type1Data.filteredFreqDefect{b}.PointBaseIndex;
    end
    if isfield(type2Data, 'filteredFreqDefect') && numel(type2Data.filteredFreqDefect) >= b
        defect2 = type2Data.filteredFreqDefect{b}.PointBaseIndex;
    end
    defect1 = unique(defect1);
    defect2 = unique(defect2);

    % Compare Defects
    inter.comp.(combi).(baseName).(pctName).uniqueDefectType1{b} = setdiff(defect1, defect2);
    inter.comp.(combi).(baseName).(pctName).uniqueDefectType2{b} = setdiff(defect2, defect1);
    inter.comp.(combi).(baseName).(pctName).sharedDefect{b} = intersect(defect1, defect2);
    inter.comp.(combi).(baseName).(pctName).lengthUniqueDefectType1(b) = numel(inter.comp.(combi).(baseName).(pctName).uniqueDefectType1{b});
    inter.comp.(combi).(baseName).(pctName).lengthUniqueDefectType2(b) = numel(inter.comp.(combi).(baseName).(pctName).uniqueDefectType2{b});
    inter.comp.(combi).(baseName).(pctName).lengthSharedDefect(b) = numel(inter.comp.(combi).(baseName).(pctName).sharedDefect{b});

    % For inside and grid
    if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
        % Extract Vertices and Points indices
        vertices1 = [];
        vertices2 = [];
        points1 = [];
        points2 = [];
        if isfield(type1Data, 'filteredFreqVertices') && numel(type1Data.filteredFreqVertices) >= b
            vertices1 = type1Data.filteredFreqVertices{b}.PointBaseIndex;
        end
        if isfield(type2Data, 'filteredFreqVertices') && numel(type2Data.filteredFreqVertices) >= b
            vertices2 = type2Data.filteredFreqVertices{b}.PointBaseIndex;
        end
        if isfield(type1Data, 'filteredFreqPoints') && numel(type1Data.filteredFreqPoints) >= b
            points1 = type1Data.filteredFreqPoints{b}.PointBaseIndex;
        end
        if isfield(type2Data, 'filteredFreqPoints') && numel(type2Data.filteredFreqPoints) >= b
            points2 = type2Data.filteredFreqPoints{b}.PointBaseIndex;
        end

        vertices1 = unique(vertices1);
        vertices2 = unique(vertices2);
        points1 = unique(points1);
        points2 = unique(points2);

        if isempty(vertices1) || isempty(vertices2) % if alpha/refinedAlpha vs inside/grid
            inter.comp.(combi).(baseName).(pctName).sharedVertices{b} = [];
            inter.comp.(combi).(baseName).(pctName).sharedPoints{b} = [];
            inter.comp.(combi).(baseName).(pctName).lengthSharedVertices(b) = 0;
            inter.comp.(combi).(baseName).(pctName).lengthSharedPoints(b) = 0;

            inter.comp.(combi).(baseName).(pctName).uniqueVerticesType1{b} = vertices1;
            inter.comp.(combi).(baseName).(pctName).uniqueVerticesType2{b} = vertices2;
            inter.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType1(b) = numel(vertices1);
            inter.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType2(b) = numel(vertices2);

            inter.comp.(combi).(baseName).(pctName).uniquePointsType1{b} = points1;
            inter.comp.(combi).(baseName).(pctName).uniquePointsType2{b} = points2;
            inter.comp.(combi).(baseName).(pctName).lengthUniquePointsType1(b) = numel(points1);
            inter.comp.(combi).(baseName).(pctName).lengthUniquePointsType2(b) = numel(points2);
        else
            % Compare Vertices
            inter.comp.(combi).(baseName).(pctName).uniqueVerticesType1{b} = setdiff(vertices1, vertices2);
            inter.comp.(combi).(baseName).(pctName).uniqueVerticesType2{b} = setdiff(vertices2, vertices1);
            inter.comp.(combi).(baseName).(pctName).sharedVertices{b} = intersect(vertices1, vertices2);
            inter.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType1(b) = numel(inter.comp.(combi).(baseName).(pctName).uniqueVerticesType1{b});
            inter.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType2(b) = numel(inter.comp.(combi).(baseName).(pctName).uniqueVerticesType2{b});
            inter.comp.(combi).(baseName).(pctName).lengthSharedVertices(b) = numel(inter.comp.(combi).(baseName).(pctName).sharedVertices{b});

            % Compare Points
            inter.comp.(combi).(baseName).(pctName).uniquePointsType1{b} = setdiff(points1, points2);
            inter.comp.(combi).(baseName).(pctName).uniquePointsType2{b} = setdiff(points2, points1);
            inter.comp.(combi).(baseName).(pctName).sharedPoints{b} = intersect(points1, points2);
            inter.comp.(combi).(baseName).(pctName).lengthUniquePointsType1(b) = numel(inter.comp.(combi).(baseName).(pctName).uniquePointsType1{b});
            inter.comp.(combi).(baseName).(pctName).lengthUniquePointsType2(b) = numel(inter.comp.(combi).(baseName).(pctName).uniquePointsType2{b});
            inter.comp.(combi).(baseName).(pctName).lengthSharedPoints(b) = numel(inter.comp.(combi).(baseName).(pctName).sharedPoints{b});
        end
    end
end

% Summarize results
inter.comp.(combi).(baseName).(pctName).sumUniqueDefectType1 = sum(inter.comp.(combi).(baseName).(pctName).lengthUniqueDefectType1);
inter.comp.(combi).(baseName).(pctName).sumUniqueDefectType2 = sum(inter.comp.(combi).(baseName).(pctName).lengthUniqueDefectType2);
inter.comp.(combi).(baseName).(pctName).sumSharedDefect = sum(inter.comp.(combi).(baseName).(pctName).lengthSharedDefect);
if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
    inter.comp.(combi).(baseName).(pctName).sumUniqueVerticesType1 = sum(inter.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType1);
    inter.comp.(combi).(baseName).(pctName).sumUniqueVerticesType2 = sum(inter.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType2);
    inter.comp.(combi).(baseName).(pctName).sumSharedVertices = sum(inter.comp.(combi).(baseName).(pctName).lengthSharedVertices);
    inter.comp.(combi).(baseName).(pctName).sumUniquePointsType1 = sum(inter.comp.(combi).(baseName).(pctName).lengthUniquePointsType1);
    inter.comp.(combi).(baseName).(pctName).sumUniquePointsType2 = sum(inter.comp.(combi).(baseName).(pctName).lengthUniquePointsType2);
    inter.comp.(combi).(baseName).(pctName).sumSharedPoints = sum(inter.comp.(combi).(baseName).(pctName).lengthSharedPoints);
end

disp(['comparison intersect types "', (type1),'" - "',(type2), '": pelvis defect categorie ', ...
    (baseName), ' and ' ,(pctName)])
end