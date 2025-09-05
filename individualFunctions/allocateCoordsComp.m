%% Comparison of intersection types: allocate coordinates

% Input:    inter: struct to be updated (output container)
%           combi: combination name '<type1>_<type2>' used under inter.comp.*
%           baseName: base/category label
%           pctName: percentage label
%           inputIDs: pelvis IDs / point-base indices to process
%           pelvisDefect: struct array providing defect meshes per ID
%           pelvis: struct array with reference pelvis data
%                     - pelvis(1).import.processed.vertices
%           allPelvis: struct with global reference points
%                     - .refPoints.gridPoints.pointsBox
%           type1: first intersection type name ('alpha' | 'refinedAlpha' | 'inside' | 'grid')
%           type2: second intersection type name ('alpha' | 'refinedAlpha' | 'inside' | 'grid')
%           weightKabsch: fieldname selecting the transform in pelvisDefect(u).transform.trafo.*

% Output:   inter: updated struct with intersect.comp with comparison of two types 

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [inter] = allocateCoordsComp(inter, combi, baseName, pctName, inputIDs, pelvisDefect, pelvis, allPelvis, type1, type2, weightKabsch)
numIDs = length(inputIDs);

% Initialize storage in allDefect
inter.comp.(combi).(baseName).(pctName).insideDefectType1 = cell(numIDs, 1);
inter.comp.(combi).(baseName).(pctName).insideDefectType2 = cell(numIDs, 1);
inter.comp.(combi).(baseName).(pctName).insideDefectShared = cell(numIDs, 1);
if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
    inter.comp.(combi).(baseName).(pctName).insideVerticesType1 = cell(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).insideVerticesType2 = cell(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).insideVerticesShared = cell(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).insidePointsType1 = cell(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).insidePointsType2 = cell(numIDs, 1);
    inter.comp.(combi).(baseName).(pctName).insidePointsShared = cell(numIDs, 1);
end

% Process each PointBase
for l = 1:numIDs
    u = inputIDs(l); % pelvisDefect ID

    % Select data for type1
    switch type1
        case 'alpha'
            pointBaseDefectType1 = pelvisDefect(u).transform.trafo.(weightKabsch).vertices;
        case 'refinedAlpha'
            pointBaseDefectType1 = pelvisDefect(u).volume.patches.all.comVertices;
        case {'inside', 'grid'}
            pointBaseDefectType1 = pelvisDefect(u).volume.patches.all.comVertices;
            pointBaseVerticesType1 = pelvis(1).import.processed.vertices;
            pointBasePointsType1 = allPelvis.refPoints.gridPoints.pointsBox;
    end
    % Select data for type2
    switch type2
        case 'alpha'
            pointBaseDefectType2 = pelvisDefect(u).transform.trafo.(weightKabsch).vertices;
        case 'refinedAlpha'
            pointBaseDefectType2 = pelvisDefect(u).volume.patches.all.comVertices;
        case {'inside', 'grid'}
            pointBaseDefectType2 = pelvisDefect(u).volume.patches.all.comVertices;
            pointBaseVerticesType2 = pelvis(1).import.processed.vertices;
            pointBasePointsType2 = allPelvis.refPoints.gridPoints.pointsBox;
    end

    % Allocate defect indices
    defectIdx1 = inter.comp.(combi).(baseName).(pctName).uniqueDefectType1{l};
    defectIdx2 = inter.comp.(combi).(baseName).(pctName).uniqueDefectType2{l};
    defectShared = inter.comp.(combi).(baseName).(pctName).sharedDefect{l};
    if ~isempty(defectIdx1)
        inter.comp.(combi).(baseName).(pctName).insideDefectType1{l} = pointBaseDefectType1(defectIdx1, :);
    else
        inter.comp.(combi).(baseName).(pctName).insideDefectType1{l} = [];
    end
    if ~isempty(defectIdx2)
        inter.comp.(combi).(baseName).(pctName).insideDefectType2{l} = pointBaseDefectType2(defectIdx2, :);
    else
        inter.comp.(combi).(baseName).(pctName).insideDefectType2{l} = [];
    end
    if ~isempty(defectShared) % Assuming shared defects are same in both (points of refined patches at the end)
        inter.comp.(combi).(baseName).(pctName).insideDefectShared{l} = pointBaseDefectType1(defectShared, :); % Assuming shared indices are valid for type1
    else
        inter.comp.(combi).(baseName).(pctName).insideDefectShared{l} = [];
    end

    % Allocate vertices and points for inside/grid
    if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
        verticesIdx1 = inter.comp.(combi).(baseName).(pctName).uniqueVerticesType1{l};
        verticesIdx2 = inter.comp.(combi).(baseName).(pctName).uniqueVerticesType2{l};
        verticesShared = inter.comp.(combi).(baseName).(pctName).sharedVertices{l};
        pointsIdx1 = inter.comp.(combi).(baseName).(pctName).uniquePointsType1{l};
        pointsIdx2 = inter.comp.(combi).(baseName).(pctName).uniquePointsType2{l};
        pointsShared = inter.comp.(combi).(baseName).(pctName).sharedPoints{l};

        if ~isempty(verticesIdx1)
            inter.comp.(combi).(baseName).(pctName).insideVerticesType1{l} = pointBaseVerticesType1(verticesIdx1, :);
        else
            inter.comp.(combi).(baseName).(pctName).insideVerticesType1{l} = [];
        end
        if ~isempty(verticesIdx2)
            inter.comp.(combi).(baseName).(pctName).insideVerticesType2{l} = pointBaseVerticesType2(verticesIdx2, :);
        else
            inter.comp.(combi).(baseName).(pctName).insideVerticesType2{l} = [];
        end
        if ~isempty(verticesShared)
            inter.comp.(combi).(baseName).(pctName).insideVerticesShared{l} = pointBaseVerticesType1(verticesShared, :);
        else
            inter.comp.(combi).(baseName).(pctName).insideVerticesShared{l} = [];
        end
        if ~isempty(pointsIdx1)
            inter.comp.(combi).(baseName).(pctName).insidePointsType1{l} = pointBasePointsType1(pointsIdx1, :);
        else
            inter.comp.(combi).(baseName).(pctName).insidePointsType1{l} = [];
        end
        if ~isempty(pointsIdx2)
            inter.comp.(combi).(baseName).(pctName).insidePointsType2{l} = pointBasePointsType2(pointsIdx2, :);
        else
            inter.comp.(combi).(baseName).(pctName).insidePointsType2{l} = [];
        end
        if ~isempty(pointsShared)
            inter.comp.(combi).(baseName).(pctName).insidePointsShared{l} = pointBasePointsType1(pointsShared, :);
        else
            inter.comp.(combi).(baseName).(pctName).insidePointsShared{l} = [];
        end

    end
end
end