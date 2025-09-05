 classdef Approxi
    % Approxi class for approximation of the defect model (intersection)
    
    % Developed by C.Micheler,
    % Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
    % Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich
    
    properties 
        alpha = struct();           % Approximation type alpha
        refinedAlpha = struct();    % Approximation type refined alpha
        inside = struct();          % Approximation type inside
        grid = struct();            % Approximation type grid 
        comp = struct()             % Comparison of types
    end

    methods
        %% Constructer: generate object
        function obj = Approxi()
            disp('class approxi initialized')
        end

        %% Approximation with the clumped sphere algorithmn of Sihaeri
        function obj = clumpedSphere(obj,typeName,pctName,fileName,nNodes,scaleFactor,smoothFactor)

            % sihaeri/DEM-ClumpedSphere: https://de.mathworks.com/matlabcentral/fileexchange/67754-sihaeri-dem-clumpedsphere
            % Paramter for pupulateSpheres:
            % fileName, nNodes, scaleFactor, smoothFact, writeTec, writeLammpsTemp, lammpsType, writeStlScaled, rho, findSTLmI
            [assembly, ~, volSpheres, ~] = populateSpheres(fileName, nNodes, scaleFactor, smoothFactor, 0, 0, 1, 0, 1, 0);
            obj.(typeName).(pctName).spheres = assembly;
            obj.(typeName).(pctName).volSpheres = volSpheres;

        end
        %% Volume approximation: intersection of sphere model and reference
        function obj = volSphereIntersect(obj, i, reference, typeName, pctNames, voxelSize)            

            % Reference
            refF = reference.faces;
            refV = reference.vertices;

            % Global bounding box 
            globalMin = [ Inf  Inf  Inf];
            globalMax = [-Inf -Inf -Inf];
            for k = 1:numel(pctNames)
                pctName = pctNames{k};
                if isfield(obj.(typeName).(pctName), 'models')
                    model = obj.(typeName).(pctName).models;
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
            obj.(typeName).gridSpheres  = gridPoints;
            obj.(typeName).boundSpheres = [globalMin; globalMax];

            % Loop over prevalence levels
            for k = 1:numel(pctNames)
                pctName = pctNames{k};
                if isfield(obj.(typeName), pctName) && isfield(obj.(typeName).(pctName), 'models')
                    nModels = numel(obj.(typeName).(pctName).models);
                    for modelNr = 1:nModels
                        model = obj.(typeName).(pctName).models(modelNr);
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
                        obj.(typeName).(pctName).modelVol(modelNr).pointsInSpheres = pointsInSpheres;

                        % Points in reference
                        t2 = tic;
                        inRefLocal = inpolyhedron(refF,refV,pointsInSpheres);
                        toc(t2)
                        intersectMask = false(size(gridPoints,1),1);
                        intersectMask(inSpheres) = inRefLocal;
                        intersectIdx = find(intersectMask);
                        %obj.(typeName).(pctName).modelVol(modelNr).intersectMaskGlobal = intersectMask; % save memory 
                        obj.(typeName).(pctName).modelVol(modelNr).intersectIdxGlobal = intersectIdx;

                        % Volume
                        obj.(typeName).(pctName).modelVol(modelNr).intersectVolume = sum(intersectMask) * voxelSize^3;

                        % Intersection % save memory
                        %intersectionPoints = gridPoints(intersectIdx, :);
                        %obj.(typeName).(pctName).modelVol(modelNr).intersectPointsPelvis = intersectionPoints;

                        % Time and memory
                        obj.(typeName).(pctName).modelVol(modelNr).intersectTime = toc(t1);
                        tmp = obj.(typeName).(pctName).modelVol(modelNr); 
                        mem = whos('tmp');
                        obj.(typeName).(pctName).modelVol(modelNr).memoryBytes = mem.bytes;
                    end
                end
            end
            disp(['Volume intersection calculated: i=' num2str(i) ', typeName=' typeName ', pctName=' pctName]);
        
        end 
        %% Volume approximation with scaling: intersection of sphere model and reference pelvis
        function obj = volSphereScale(obj, i, gridPoints, typeName, pctName, voxelVol, targetVolume, tolerancePct, maxIter, initialScale, scalingExp)

            if isfield(obj.(typeName), pctName) && isfield(obj.(typeName).(pctName), 'modelVol')

                nModels = numel(obj.(typeName).(pctName).modelVol);
                for modelNr = 1:nModels
                    t1  = tic;
                    model     = obj.(typeName).(pctName).models(modelNr);
                    spheresOrig    = model.spheres;                         
                    if isfield(obj.(typeName).pairs, "modelVol")
                        idxInit = obj.(typeName).pairs.modelVol(1).intersectIdxGlobal; % pairs with the largest volume of grid points
                    else
                        warning(['No "pairs"-model (no upscaling possible): i=' num2str(i) ', typeName=' typeName]);
                        idxInit = obj.(typeName).(pctName).modelVol(modelNr).intersectIdxGlobal;  % no upscaling possible
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
                    obj.(typeName).(pctName).modelScale(modelNr).intersectIdxGlobal = intersectIdxGlobal;
                    %obj.(typeName).(pctName).modelScale(modelNr).intersectMaskGlobal = intersectMaskGlobal;
                    obj.(typeName).(pctName).modelScale(modelNr).scaledSpheres = scaledSpheres;
                    obj.(typeName).(pctName).modelScale(modelNr).intersectVolumeScaled = volNow;
                    obj.(typeName).(pctName).modelScale(modelNr).scaleFactorToRef = scaleFactor;
                    obj.(typeName).(pctName).modelScale(modelNr).volScaledSpheres = modelVolume;
                    obj.(typeName).(pctName).modelScale(modelNr).volRatio = volRatio;
                    % Time and memory
                    obj.(typeName).(pctName).modelScale(modelNr).intersectTime = toc(t1);
                    tmp = obj.(typeName).(pctName).modelScale(modelNr);
                    memInfo = whos('tmp');
                    obj.(typeName).(pctName).modelScale(modelNr).memoryBytes = memInfo.bytes;
                end
            end
            disp(['Volume intersection with scaling: i=' num2str(i) ', typeName=' typeName ', pctName=' pctName]);

        end
         %% Comparison approxi models (point frequency, voxel-based)
        function obj = compSpheres(obj, i, typeName, pctNames)

            % Global grid points
            gridPoints = obj.(typeName).gridSpheres;
            nPoints = size(gridPoints,1);

            % Point frequency
            pointFreq = zeros(nPoints,1);
            modelList  = struct('pctName',{},'modelNr',{},'pointFreq',{},'scoreSum',{},'scoreMean',{},'scoreGlobal',{});

            dupMap  = containers.Map('KeyType','char','ValueType','logical'); % for duplicates
            for k = 1:numel(pctNames)
                pctName = pctNames{k};
                if isfield(obj.(typeName), pctName) && isfield(obj.(typeName).(pctName), 'modelScale')

                    modelScale = obj.(typeName).(pctName).modelScale;
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
                idx = obj.(typeName).(pctName).modelScale(modelNr).intersectIdxGlobal(:);
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
            obj.(typeName).modelFit = modelList;

            disp(['Compared sphere models: i=' num2str(i) ', typeName=' typeName]);

        end
        %% Create mesh für sphere model (Fibonacci)
        function obj = sphereMesh(obj, i, typeName,  bestModelType, spheres, spacing)

            nSpheres = size(spheres,2);
            outerPts = [];

            % Fibonacci points on unit sphere
            goldAng = pi * (3 - sqrt(5));   % golden angle

            for s = 1:nSpheres
                centre = spheres(1:3,s)';
                radius = spheres(4,s);
                % Required number of points for desired spacing
                approxN = max(round(4*pi*radius^2 / spacing^2), 50);   % min. 50 points
                k = (0:approxN-1)';       % Fibonacci Index
                z = 1 - 2*(k+0.5)/approxN;
                rad = sqrt(1 - z.^2);
                phi = k * goldAng;
                ptsLocal = [rad .* cos(phi), rad .* sin(phi), z];  % unit sphere
                ptsWorld = radius * ptsLocal + centre;  % scale and transform

                % Filter points inside
                inOther = false(size(ptsWorld,1),1);
                for j = 1:nSpheres
                    if j == s 
                        continue; 
                    end
                    cj = spheres(1:3,j)';
                    rj = spheres(4,j);
                    d  = vecnorm(ptsWorld - cj, 2, 2);            % Distance to sphere j
                    inOther = inOther | (d < rj - 1e-6);          % Tolerance
                end
                outerPts = [outerPts; ptsWorld(~inOther,:)];                                  
            end

            obj.(typeName).bestModelDist.(['sphere' bestModelType 'Mesh']) = outerPts; 

            disp(['Created sphere mesh: i=' num2str(i) ', typeName=' typeName]);

        end
        %% Vertex-to-Nearest-Neighbour Distance
        function obj = sphereDistance(obj, i, defIdx, verticesDefect, facesDefect, typeName, bestModelType)

            % Spheres
            spherePts = obj.(typeName).bestModelDist.(['sphere' bestModelType 'Mesh']);
            bestModelNr = obj.(typeName).bestModel.(['score' bestModelType 'ModelNr']);
            bestModelPct = obj.(typeName).bestModel.(['score' bestModelType 'Pct']);
            scaledSpheres = obj.(typeName).(bestModelPct).modelScale(bestModelNr).scaledSpheres; % scaled spheres

            % Nearest Neighbour
            Mdl = KDTreeSearcher(spherePts);
            [~, dSurface] = knnsearch(Mdl, verticesDefect);

            % Points in sphere model
            centers = scaledSpheres(1:3,:)'; 
            radii   = scaledSpheres(4,:)';    
            D = pdist2(verticesDefect, centers);  
            insideMask = any(D <= radii', 2);      % Points inside sphere model?
            %Inside = 0, Outside = vertex-to-nearest neighbour distance
            minDistances = dSurface;
            minDistances(insideMask) = 0;

            % Colourised nearest neighbour
            nearVertexMatrix = minDistances(facesDefect);
            faceMean = mean(nearVertexMatrix, 2);
            minVal = min(faceMean); maxVal = max(faceMean);
            normVals = (faceMean - minVal) / (maxVal - minVal);
            normVals = max(0, min(1, normVals)); 
            colourNum = 32768;
            colourMap = viridis(colourNum);
            colourIdx = round(normVals * (colourNum-1)) + 1;
            rgbColour = colourMap(colourIdx,:);

            % Save
            obj.(typeName).bestModelDist.(['nearDist' bestModelType]){defIdx,1} = minDistances;
            obj.(typeName).bestModelDist.(['nearDist' bestModelType 'Face']){defIdx,1} = faceMean;
            obj.(typeName).bestModelDist.(['nearDist' bestModelType 'Norm']){defIdx,1} = normVals;
            obj.(typeName).bestModelDist.(['nearDist' bestModelType 'FaceColour']){defIdx,1} = rgbColour;
            obj.(typeName).bestModelDist.(['nearDist' bestModelType 'Mean'])(defIdx,1) = mean(minDistances);
            obj.(typeName).bestModelDist.(['nearDist' bestModelType 'Max'])(defIdx,1) = max(minDistances);
            obj.(typeName).bestModelDist.(['nearDist' bestModelType 'Std'])(defIdx,1) = std(minDistances);

            disp(['Vertex-to-Nearest-Neighbour Distance: i=' num2str(i) ', typeName=' typeName ', defectNr=' num2str(defIdx)]);

        end   
        %% Vertex-to-Nearest-Neighbour Distance
        % function obj = sphereDistance(obj, i, defIdx, verticesDefect, facesDefect, typeName, bestModelType)      
        % 
        %     % Sphere Mesh
        %     spherePts = obj.(typeName).bestModelDist.(['sphere' bestModelType 'Mesh']);
        %     % KD-Tree on sphere mesh
        %     Mdl = KDTreeSearcher(spherePts);
        % 
        %     % Nearest neighbour
        %     idx = knnsearch(Mdl, verticesDefect);
        %     minDistances = sqrt(sum((verticesDefect - spherePts(idx,:)).^2, 2));
        % 
        %     % Colourised nearest neighbour
        %     nearVertexMatrix = minDistances(facesDefect);
        %     faceMean = mean(nearVertexMatrix, 2);
        %     minVal = min(faceMean); maxVal = max(faceMean);
        %     normVals = (faceMean - minVal) / (maxVal - minVal);
        %     normVals = max(0, min(1, normVals)); % Clamping
        %     colourNum = 32768;
        %     colourMap = viridis(colourNum);
        %     colourIdx = round(normVals * (colourNum-1)) + 1;
        %     rgbColour = colourMap(colourIdx,:);
        % 
        %     % Save
        %     obj.(typeName).bestModelDist.(['nearDist' bestModelType]){defIdx,1} = minDistances;
        %     obj.(typeName).bestModelDist.(['nearDist' bestModelType 'Face']){defIdx,1} = faceMean;
        %     obj.(typeName).bestModelDist.(['nearDist' bestModelType 'Norm']){defIdx,1} = normVals;
        %     obj.(typeName).bestModelDist.(['nearDist' bestModelType 'FaceColour']){defIdx,1} = rgbColour;
        %     obj.(typeName).bestModelDist.(['nearDist' bestModelType 'Mean'])(defIdx,1) = mean(minDistances);
        %     obj.(typeName).bestModelDist.(['nearDist' bestModelType 'Max'])(defIdx,1) = max(minDistances);
        %     obj.(typeName).bestModelDist.(['nearDist' bestModelType 'Std'])(defIdx,1) = std(minDistances);
        % 
        %     disp(['Vertex-to-Nearest-Neighbour Distance: i=' num2str(i) ', typeName=' typeName ', defectNr=' num2str(defIdx)]);
        % 
        % end   
    end
 end