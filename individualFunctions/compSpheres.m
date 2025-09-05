%% Comparison approxi models (point frequency, voxel-based)

% Input:    typeName: 'alpha' | 'refinedAlpha' | 'inside' | 'grid'
%           approx: struct to be updated (output container)
%           i: Numeric identifier used only for logging
%           pctNames: subfield labels to process (e.g., {'p50','p75'}).

% Output:   approx: updated struct

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [approx] = compSpheres(approx, i, typeName, pctNames)

% Global grid points
gridPoints = approx.(typeName).gridSpheres;
nPoints = size(gridPoints,1);

% Point frequency
pointFreq = zeros(nPoints,1);
modelList  = struct('pctName',{},'modelNr',{},'pointFreq',{},'scoreSum',{},'scoreMean',{},'scoreGlobal',{});

dupMap  = containers.Map('KeyType','char','ValueType','logical'); % for duplicates
for k = 1:numel(pctNames)
    pctName = pctNames{k};
    if isfield(approx.(typeName), pctName) && isfield(approx.(typeName).(pctName), 'modelScale')

        modelScale = approx.(typeName).(pctName).modelScale;
        for modelNr = 1:numel(modelScale)
            idx  = modelScale(modelNr).intersectIdxGlobal(:);
            if isempty(idx)
                continue
            end

            % Test for duplicates
            key = sprintf('%d,', sort(idx));
            if isKey(dupMap, key)
                continue
            end
            dupMap(key) = true;

            % Grid points inside sphere model and reference
            mask = false(nPoints,1);
            mask(idx) = true;
            % Frequency
            pointFreq  = pointFreq + mask;

            modelList(end+1).pctName   = pctName;
            modelList(end).modelNr   = modelNr;
        end
    end
end

%  Scoring (3 methods)
for m = 1:numel(modelList)
    pctName = modelList(m).pctName;
    modelNr = modelList(m).modelNr;
    idx = approx.(typeName).(pctName).modelScale(modelNr).intersectIdxGlobal(:);
    freqVals = pointFreq(idx);          % point frequency
    scoreSum = sum(freqVals);           % absolut
    scoreMean = mean(freqVals);         % mean frequency
    pInside = sum(pointFreq > 0);       % all points inside sphere models and reference
    scoreGlobal = scoreSum / pInside;   % global; reference: all points inside sphere models and reference
    modelList(m).pointFreq = freqVals ;
    modelList(m).scoreSum  = scoreSum;
    modelList(m).scoreMean = scoreMean;
    modelList(m).scoreGlobal = scoreGlobal;
end

% Save results
approx.(typeName).modelFit = modelList;

disp(['Compared sphere models: i=' num2str(i) ', typeName=' typeName]);

end