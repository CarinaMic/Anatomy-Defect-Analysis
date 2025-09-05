%% Volume approximation: intersection of sphere model and reference

% Input:    typeName: 'alpha' | 'refinedAlpha' | 'inside' | 'grid'
%           approx: struct to be updated (output container)
%           i: Numeric ID used only for logging
%           reference: Struct with reference surface
%           pctNames: List of subfield labels to process (e.g., {'p50','p75'}).
%           voxelSize: Cartesian grid spacing used for voxelization

% Output:   approx: updated struct
%           approx.(typeName).(pctName).spheres    % clumped-sphere assembly
%           approx.(typeName).(pctName).volSpheres % volume estimated from spheres

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [approx] = volSphereIntersect(approx, i, reference, typeName, pctNames, voxelSize)

% Reference
refF = reference.faces;
refV = reference.vertices;

% Global bounding box
globalMin = [ Inf  Inf  Inf];
globalMax = [-Inf -Inf -Inf];
for k = 1:numel(pctNames)
    pctName = pctNames{k};
    if isfield(approx.(typeName).(pctName), 'models')
        model = approx.(typeName).(pctName).models;
        for modelNr = 1:numel(model)
            spheres     = model(modelNr).spheres;  % [xyzr]
            minCorner = min(spheres(1:3,:) - spheres(4,:), [], 2)';
            maxCorner = max(spheres(1:3,:) + spheres(4,:), [], 2)';
            globalMin = min(globalMin, minCorner);
            globalMax = max(globalMax, maxCorner);
        end
    end
end
globalMin = floor(globalMin/voxelSize)*voxelSize;
globalMax = ceil( globalMax/voxelSize)*voxelSize;
% Voxel grid
[xg,yg,zg] = ndgrid(globalMin(1):voxelSize:globalMax(1), ...
    globalMin(2):voxelSize:globalMax(2), ...
    globalMin(3):voxelSize:globalMax(3));
gridPoints = [xg(:) yg(:) zg(:)];
approx.(typeName).gridSpheres  = gridPoints;
approx.(typeName).boundSpheres = [globalMin; globalMax];

% Loop over prevalence levels
for k = 1:numel(pctNames)
    pctName = pctNames{k};
    if isfield(approx.(typeName), pctName) && isfield(approx.(typeName).(pctName), 'models')
        nModels = numel(approx.(typeName).(pctName).models);
        for modelNr = 1:nModels
            model = approx.(typeName).(pctName).models(modelNr);
            spheres   = model.spheres;

            t1 = tic;
            % Points in sphere model
            inSpheres = false(size(gridPoints,1),1);
            for s = 1:size(spheres,2)
                center = spheres(1:3,s)';
                radius = spheres(4,s);
                distances = sqrt(sum((gridPoints - center).^2, 2));
                inSpheres = inSpheres | (distances <= radius);
            end
            pointsInSpheres = gridPoints(inSpheres,:);
            approx.(typeName).(pctName).modelVol(modelNr).pointsInSpheres = pointsInSpheres;

            % Points in reference
            t2 = tic;
            inRefLocal = inpolyhedron(refF,refV,pointsInSpheres);
            toc(t2)
            intersectMask = false(size(gridPoints,1),1);
            intersectMask(inSpheres) = inRefLocal;
            intersectIdx = find(intersectMask);
            %approx.(typeName).(pctName).modelVol(modelNr).intersectMaskGlobal = intersectMask; % save memory
            approx.(typeName).(pctName).modelVol(modelNr).intersectIdxGlobal = intersectIdx;

            % Volume
            approx.(typeName).(pctName).modelVol(modelNr).intersectVolume = sum(intersectMask) * voxelSize^3;

            % Intersection % save memory
            %intersectionPoints = gridPoints(intersectIdx, :);
            %approx.(typeName).(pctName).modelVol(modelNr).intersectPointsPelvis = intersectionPoints;

            % Time and memory
            approx.(typeName).(pctName).modelVol(modelNr).intersectTime = toc(t1);
            tmp = approx.(typeName).(pctName).modelVol(modelNr);
            mem = whos('tmp');
            approx.(typeName).(pctName).modelVol(modelNr).memoryBytes = mem.bytes;
        end
    end
end
disp(['Volume intersection calculated: i=' num2str(i) ', typeName=' typeName ', pctName=' pctName]);

end