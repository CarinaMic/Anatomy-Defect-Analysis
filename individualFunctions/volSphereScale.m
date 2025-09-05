%% Volume approximation with scaling: intersection of sphere model and reference pelvis

% Input:    typeName: 'alpha' | 'refinedAlpha' | 'inside' | 'grid'
%           approx: struct to be updated (output container)
%           gridPoints: array of global voxel-grid point coordinates
%           pctName: subfield label (e.g., 'p50', 'p75')
%           voxelVol: Scalar voxel volume (e.g., voxelSize^3), units consistent with targetVolume
%           targetVolume: Desired reference-intersection volume (same units as voxelVol, e.g., mm^3)
%           tolerancePct: tolerance (e.g., 0.01 for ±1%) for target volume
%           maxIter: Maximum number of scaling iterations
%           initialScale: Initial multiplicative scale factor applied to centers and radii
%           scalingExp: Exponent/gain for multiplicative update: scaleFactor <- scaleFactor * (targetVolume/volNow)^scalingExp

% Output:   approx: updated struct

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [approx] = volSphereScale(approx, i, gridPoints, typeName, pctName, voxelVol, targetVolume, tolerancePct, maxIter, initialScale, scalingExp)

if isfield(approx.(typeName), pctName) && isfield(approx.(typeName).(pctName), 'modelVol')

    nModels = numel(approx.(typeName).(pctName).modelVol);
    for modelNr = 1:nModels
        t1  = tic;
        model     = approx.(typeName).(pctName).models(modelNr);
        spheresOrig    = model.spheres;
        if isfield(approx.(typeName).pairs, "modelVol")
            idxInit = approx.(typeName).pairs.modelVol(1).intersectIdxGlobal; % pairs with the largest volume of grid points
        else
            warning(['No "pairs"-model (no upscaling possible): i=' num2str(i) ', typeName=' typeName]);
            idxInit = approx.(typeName).(pctName).modelVol(modelNr).intersectIdxGlobal;  % no upscaling possible
        end

        % CoM of model
        CoM = mean(spheresOrig(1:3,:), 2);
        scaleFactor = initialScale;
        tolerance = tolerancePct * targetVolume;
        % Points in sphere model and reference pelvis
        intersectPoints = gridPoints(idxInit,:);

        % Iterative scaling
        for iter = 1:maxIter
            % Scale sphere model around CoM
            scaledSpheres = spheresOrig;
            scaledSpheres(1:3,:) = CoM + scaleFactor*(spheresOrig(1:3,:) - CoM);
            scaledSpheres(4,:)   = scaleFactor*spheresOrig(4,:);

            % Check which points of the first volume intersection lie within the scaled model
            inScaled = false(size(intersectPoints, 1), 1);
            for s = 1:size(scaledSpheres,2)
                center = scaledSpheres(1:3,s)';
                radius = scaledSpheres(4,s);
                distances = sqrt(sum((intersectPoints - center).^2, 2));
                inScaled = inScaled | (distances <= radius);
            end
            % New volume of intersection
            volNow = sum(inScaled) * voxelVol;
            errorVol = volNow - targetVolume;

            % Cancellation condition
            if abs(errorVol) <= tolerance
                disp(['Target volume achieved: ' num2str(volNow) ' mm^3 in iteration ' num2str(iter)]);
                break;
            end

            % Update scaleFactor
            scaleFactor = scaleFactor * (targetVolume / volNow)^scalingExp;
        end

        % GlobalGrid (Index & Maske)
        intersectIdxGlobal = idxInit(inScaled);
        intersectMaskGlobal = false(size(gridPoints,1), 1);
        intersectMaskGlobal(intersectIdxGlobal) = true;

        % Sphere model volume (without reference intersection)
        inModel = false(size(gridPoints,1),1);
        for s = 1:size(scaledSpheres,2)
            ctr = scaledSpheres(1:3,s)';
            rad = scaledSpheres(4,s);
            d   = sqrt(sum((gridPoints - ctr).^2, 2));
            inModel = inModel | (d <= rad);
        end
        modelIdxGlobal = find(inModel);
        modelVolume    = numel(modelIdxGlobal) * voxelVol;
        volRatio = volNow / modelVolume;

        % Save
        approx.(typeName).(pctName).modelScale(modelNr).intersectIdxGlobal = intersectIdxGlobal;
        approx.(typeName).(pctName).modelScale(modelNr).intersectMaskGlobal = intersectMaskGlobal;
        approx.(typeName).(pctName).modelScale(modelNr).scaledSpheres = scaledSpheres;
        approx.(typeName).(pctName).modelScale(modelNr).intersectVolumeScaled = volNow;
        approx.(typeName).(pctName).modelScale(modelNr).scaleFactorToRef = scaleFactor;
        approx.(typeName).(pctName).modelScale(modelNr).volScaledSpheres = modelVolume;
        approx.(typeName).(pctName).modelScale(modelNr).volRatio = volRatio;
        % Time and memory
        approx.(typeName).(pctName).modelScale(modelNr).intersectTime = toc(t1);
        tmp = approx.(typeName).(pctName).modelScale(modelNr);
        memInfo = whos('tmp');
        approx.(typeName).(pctName).modelScale(modelNr).memoryBytes = memInfo.bytes;
    end
end
disp(['Volume intersection with scaling: i=' num2str(i) ', typeName=' typeName ', pctName=' pctName]);

end