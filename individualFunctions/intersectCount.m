%% Intersection combinations (counts)

% Input:    baseIntersect: category label used for logging %%%
%           typeName: 'alpha' | 'refinedAlpha' | 'inside' | 'grid'
%           intersectPctNames: list of field name with intersection threshold, e.g., {'pct50','pct75',...};
%           referenceDefectIdx: per-pair cell array of indices (defect)
%           referenceVerticesIdx: per-pair cell array of indices (vertices)   [required for 'inside'/'grid']
%           referencePointsIdx: per-pair cell array of indices (points)     [required for 'inside'/'grid']
%           inputIDs: IDs for the selected categorie
%           machingIdx: row number of the k-th pair from pairs.IDs in the global 'all'-pair table

% Output:   inter: updated struct with inter.(typeName)

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function inter = intersectCount(inter, baseIntersect, typeName, intersectPctNames,...
    referenceDefectIdx, referenceVerticesIdx, referencePointsIdx, inputIDs, matchingIdx)

% Preallocation
numPointBases = size(inputIDs, 1); % num of point bases = num IDs
numMatching = numPointBases - 1;
inter.(typeName).pairs.frequencyDefect = cell(numPointBases, 1);
if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
    inter.(typeName).pairs.frequencyVertices = cell(numPointBases, 1);
    inter.(typeName).pairs.frequencyPoints = cell(numPointBases, 1);
end

% Calculate frequencies for defects, vertices, and points
for i = 1:numPointBases
    % Start and end indices of defect point input (first column in pairs)
    startIdx = (i - 1) * numMatching + 1;
    endIdx = i * numMatching;
    % Indices for this defect point input (point base)
    currentPairIdx = matchingIdx(startIdx:endIdx); % matchingIdx = idx of pairs 'all'

    % Collect indices
    currentPointBaseDefectIdx = vertcat(referenceDefectIdx{currentPairIdx});
    if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
        currentPointBaseVerticesIdx = vertcat(referenceVerticesIdx{currentPairIdx});
        currentPointBasePointsIdx = vertcat(referencePointsIdx{currentPairIdx});
    end

    % Count frequencies
    % Defects
    [uniqueIndicesDefect, ~, idxCountsDefect] = unique(currentPointBaseDefectIdx);
    pointBaseFreqDefect= accumarray(idxCountsDefect, 1);
    inter.(typeName).pairs.frequencyDefect{i,1} = table(uniqueIndicesDefect, pointBaseFreqDefect, ...
        'VariableNames', {'PointBaseIndex', 'Frequency'});
    if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
        % Vertices
        [uniqueIndicesVertices, ~, idxCountsVertices] = unique(currentPointBaseVerticesIdx);
        pointBaseFreqVertices = accumarray(idxCountsVertices, 1);
        inter.(typeName).pairs.frequencyVertices{i,1} = table(uniqueIndicesVertices, pointBaseFreqVertices, ...
            'VariableNames', {'PointBaseIndex', 'Frequency'});
        % Points
        [uniqueIndicesPoints, ~, idxCountsPoints] = unique(currentPointBasePointsIdx);
        pointBaseFreqPoints = accumarray(idxCountsPoints, 1);
        inter.(typeName).pairs.frequencyPoints{i,1} = table(uniqueIndicesPoints, pointBaseFreqPoints, ...
            'VariableNames', {'PointBaseIndex', 'Frequency'});
    end
end


% Filter frequencies based on intersection percentages
timeIntersectFilter = NaN(length(intersectPctNames), 1);
for k = 1:length(intersectPctNames)
    pctCount = inter.(typeName).(intersectPctNames{k}).count;
    if pctCount > 1
        t2 = tic;
        % Frequency filter
        filterThreshold = pctCount - 1; % minus 1, as point is also in own defect
        % Preallocation
        filteredFreqDefect = cell(numPointBases, 1);
        if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
            filteredFreqVertices = cell(numPointBases, 1);
            filteredFreqPoints = cell(numPointBases, 1);
        end
        for i = 1:numPointBases
            % Filter frequencyDefect
            filteredFreqDefect{i} = inter.(typeName).pairs.frequencyDefect{i}( ...
                inter.(typeName).pairs.frequencyDefect{i}.Frequency >= filterThreshold, :);
            % Filter frequencyVertices and frequencyPoints
            if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
                filteredFreqVertices{i} = inter.(typeName).pairs.frequencyVertices{i}( ...
                    inter.(typeName).pairs.frequencyVertices{i}.Frequency >= filterThreshold, :);
                filteredFreqPoints{i} = inter.(typeName).pairs.frequencyPoints{i}( ...
                    inter.(typeName).pairs.frequencyPoints{i}.Frequency >= filterThreshold, :);
            end
        end
        % Store filtered tables
        inter.(typeName).(intersectPctNames{k}).filteredFreqDefect = filteredFreqDefect;
        if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
            inter.(typeName).(intersectPctNames{k}).filteredFreqVertices = filteredFreqVertices;
            inter.(typeName).(intersectPctNames{k}).filteredFreqPoints = filteredFreqPoints;
        end
        timeIntersectFilter(k,1) = toc(t2);
    end
end
inter.(typeName).timeIntersectFilter = timeIntersectFilter;

disp(['defect intersection "', (typeName), '": frequency tables for pelvis defect categorie "', (baseIntersect), '"']);

end