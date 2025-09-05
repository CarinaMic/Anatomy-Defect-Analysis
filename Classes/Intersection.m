 classdef Intersection
    % Intersect class for distribution of the defect models
    
    % Developed by C.Micheler,
    % Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
    % Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich
    

    properties 
        alpha = struct();           % Intersection type alpha
        refinedAlpha = struct();    % Intersection type refined alpha
        inside = struct();          % Intersection type inside
        grid = struct();            % Intersection type grid 
        comp = struct()             % Comparison of types
        
    end
    
    methods
        %% Constructer: generate object
        function obj = Intersection()
            disp('class intersect initialized')
        end
        
        %% Intersection sequence (function inpolyhedron): with defect volume and defect surface points; new defect volume with shrink wrap 
        function obj = intersectSeq(obj, typeName, intersectPct, baseIntersect, ratio, n, overlapCount, importVertices, importFaces)
            % Calculate intersection combinations: binomial coefficient
            defectCombis = nchoosek(n, overlapCount);
            obj.(typeName).(intersectPct).(baseIntersect).combis = defectCombis;
            defectCombiCount = size(defectCombis, 1);

            % Preallocate for overlap results
            overlapVertices = cell(defectCombiCount, 1);
            overlapFaces = cell(defectCombiCount, 1);

            for i = 1:defectCombiCount
                % Initial vertices and faces
                overlapVerticesSet = importVertices{defectCombis(i, 1)};
                overlapFacesSet = importFaces{defectCombis(i, 1)};

                % Perform overlap computation for the current combination
                for j = 1:(overlapCount - 1)
                    % Overlap 1: Check points of current set in next defect
                    overlap1 = inpolyhedron(overlapFacesSet, overlapVerticesSet, importVertices{defectCombis(i, j + 1)});
                    overlap1Vertices = importVertices{defectCombis(i, j + 1)}(overlap1, :);

                    % Overlap 2: Check points of next defect in current set
                    overlap2 = inpolyhedron(importFaces{defectCombis(i, j + 1)}, importVertices{defectCombis(i, j + 1)}, overlapVerticesSet);
                    overlap2Vertices = overlapVerticesSet(overlap2, :);

                    % Combine overlaps
                    if ~isempty(overlap1Vertices) || ~isempty(overlap2Vertices)
                        % Update overlap vertices and shrink wrap
                        overlapVerticesSet = [overlap1Vertices; overlap2Vertices];
                        overlapFacesSet = boundary(overlapVerticesSet, 0.8); % Shrink wrap with alpha shape
                        if isempty(overlapFacesSet)
                            break; % Exit loop if no faces generated
                        end
                    else
                        break; % Exit loop if no intersection
                    end
                end

                % Save results for current combination
                overlapVertices{i, 1} = overlapVerticesSet;
                overlapFaces{i, 1} = overlapFacesSet;
                disp(['Processed defect combination ', num2str(i)]);
            end

            % Store results in the object
            obj.(typeName).(intersectPct).(baseIntersect).vertices = overlapVertices;
            obj.(typeName).(intersectPct).(baseIntersect).faces = overlapFaces;

            % Log the intersection completion
            disp(['Defect intersection for overlap: ', num2str(ratio), ' % completed.']);
        end

        %% Intersection of boundary box and defect (Pre-filter)
        function obj = intersectBox(obj,typeName,pointsNum,volNum,box,inputPoints,inputPointNums,inputMasksIdx,inputMaskNums)

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
                    %obj.(typeName).pairs.boxDefectInsideMask = isInsideAll; % logical mask     %%% to save memory
                    obj.(typeName).pairs.boxDefectInsideMaskIdx = find(isInsideAll); % idx
                case {'inside', 'grid'}
                    % Logical mask base:
                    % inside: defect Mesh + refined patches - reference pelvis vertices - grid points
                    isDefect = isInsideAll(1:inputPointNums(1),1);
                    isVertices = isInsideAll((inputPointNums(1)+1):(inputPointNums(1)+inputPointNums(2)), 1);
                    isPoints = isInsideAll((inputPointNums(1)+inputPointNums(2)+1):end, 1);
                    inputMaskDefect = inputMasksIdx(1:inputPointNums(1),1);
                    inputMaskVertices = inputMasksIdx((inputPointNums(1)+1):(inputPointNums(1)+inputPointNums(2)), 1);
                    inputMaskPoints = inputMasksIdx((inputPointNums(1)+inputPointNums(2)+1):end, 1);
                    obj.(typeName).pairs.boxDefectInsideMaskIdx = inputMaskDefect(isDefect,1); % idx
                    obj.(typeName).pairs.boxVerticesInsideMaskIdx = inputMaskVertices(isVertices,1);
                    obj.(typeName).pairs.boxPointsInsideMaskIdx = inputMaskPoints(isPoints,1);
                    % Logical mask                                                              %%% to save memory
                    %obj.(typeName).pairs.boxDefectInsideMask = false(inputMaskNums(1),1); 
                    %obj.(typeName).pairs.boxVerticesInsideMask = false(inputMaskNums(2),1);
                    %obj.(typeName).pairs.boxPointsInsideMask = false(inputMaskNums(3),1);
                    %obj.(typeName).pairs.boxDefectInsideMask(obj.(typeName).pairs.boxDefectInsideMaskIdx) = true; 
                    %obj.(typeName).pairs.boxVerticesInsideMask(obj.(typeName).pairs.boxVerticesInsideMaskIdx) = true;
                    %obj.(typeName).pairs.boxPointsInsideMask(obj.(typeName).pairs.boxPointsInsideMaskIdx) = true;
            end
              
            % Points in global cosy                                                             %%% to save memory
            %pointsInside = inputPoints(isInsideAll,:);
            %obj.(typeName).pairs.boxInside = pointsInside;                                     
            
            disp(['intersection "', (typeName), '" defect - box: pelvis pair ', num2str(pointsNum), '-', num2str(volNum)]);

        end

        %% Intersection of two defects
        function obj = intersectPairs(obj,typeName,pointsNum,volNum,inputIdx,inputNums,inputPointBase,inputBaseNums,inputVolume)
            
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
            insideVol = inpolyhedron(inputVolume.faces,inputVolume.vertices,... 
                inputAll);

            % % Plot for control
            % figure;
            % hold on;
            % grid on;
            % axis equal;
            % title('Inside and Outside Points Test');
            % xlabel('X');
            % ylabel('Y');
            % zlabel('Z');
            % % Plot the volume (triangulated surface)
            % trisurf(inputVolume.faces, inputVolume.vertices(:,1), inputVolume.vertices(:,2), inputVolume.vertices(:,3), ...
            %     'FaceColor', 'cyan', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            % % Plot inside points (in blue)
            % scatter3(inputAll(insideVol,1), inputAll(insideVol,2), inputAll(insideVol,3), ...
            %     50, 'b', 'filled', 'DisplayName', 'Inside Points');
            % % Plot outside points (in red)
            % scatter3(inputAll(~insideVol,1), inputAll(~insideVol,2), inputAll(~insideVol,3), ...
            %     50, 'r', 'filled', 'DisplayName', 'Outside Points');
            % view(3)

            % Extract the points that are inside 
            switch typeName
                case {'alpha', 'refinedAlpha'}
                    % Logical mask base:
                    % alpha: defect Mesh
                    % refinedAlpha: defect Mesh + refined patches
                    obj.(typeName).pairs.defectInsideMaskIdx = inputIdx(insideVol,1); % idx
                    % Logical mask                                                              %%% to save memory
                    %obj.(typeName).pairs.defectInsideMask = false(inputMaskNums,1);
                    %obj.(typeName).pairs.defectInsideMask(obj.(typeName).pairs.defectInsideMaskIdx) = true; 
                case {'inside', 'grid'}
                    % Logical mask base:
                    % inside: defect Mesh + refined patches - reference pelvis vertices - grid points
                    isDefect = insideVol(1:inputNums(1),1);
                    isVertices = insideVol((inputNums(1)+1):(inputNums(1)+inputNums(2)), 1);
                    isPoints = insideVol((inputNums(1)+inputNums(2)+1):end, 1);
                    inputMaskDefect = inputIdx(1:inputNums(1),1);
                    inputMaskVertices = inputIdx((inputNums(1)+1):(inputNums(1)+inputNums(2)), 1);
                    inputMaskPoints = inputIdx((inputNums(1)+inputNums(2)+1):end, 1);
                    obj.(typeName).pairs.defectInsideMaskIdx = inputMaskDefect(isDefect,1); % idx
                    obj.(typeName).pairs.verticesInsideMaskIdx = inputMaskVertices(isVertices,1);
                    obj.(typeName).pairs.pointsInsideMaskIdx = inputMaskPoints(isPoints,1);
                    % Logical mask                                                              %%% to save memory
                    %obj.(typeName).pairs.defectInsideMask = false(inputMaskNums(1),1); 
                    %obj.(typeName).pairs.verticesInsideMask = false(inputMaskNums(2),1);
                    %obj.(typeName).pairs.pointsInsideMask = false(inputMaskNums(3),1);
                    %obj.(typeName).pairs.defectInsideMask(obj.(typeName).pairs.defectInsideMaskIdx) = true; 
                    %obj.(typeName).pairs.verticesInsideMask(obj.(typeName).pairs.verticesInsideMaskIdx) = true;
                    %obj.(typeName).pairs.pointsInsideMask(obj.(typeName).pairs.pointsInsideMaskIdx) = true;    
            end

            
            % Points
            %pointsInside = inputAll(insideVol,:);                                              %%% to save memory
            %obj.(typeName).pairs.inside = pointsInside;            

            disp(['intersection "', (typeName), '" defect: pelvis pair ', num2str(pointsNum), '-', num2str(volNum)]);

        end

       %% Intersection combinations (counts)
        function obj = intersectCount(obj, baseIntersect, typeName, intersectPctNames,...
                referenceDefectIdx, referenceVerticesIdx, referencePointsIdx)

            inputIDs = obj.(typeName).IDs;
            matchingIdx = obj.(typeName).pairs.idxOfAll;

            % Preallocation
            numPointBases = size(inputIDs, 1); % num of point bases = num IDs
            numMatching = numPointBases - 1;
            obj.(typeName).pairs.frequencyDefect = cell(numPointBases, 1);
            if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
                obj.(typeName).pairs.frequencyVertices = cell(numPointBases, 1);
                obj.(typeName).pairs.frequencyPoints = cell(numPointBases, 1);
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
                obj.(typeName).pairs.frequencyDefect{i,1} = table(uniqueIndicesDefect, pointBaseFreqDefect, ...
                    'VariableNames', {'PointBaseIndex', 'Frequency'});
                if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
                    % Vertices
                    [uniqueIndicesVertices, ~, idxCountsVertices] = unique(currentPointBaseVerticesIdx);
                    pointBaseFreqVertices = accumarray(idxCountsVertices, 1);
                    obj.(typeName).pairs.frequencyVertices{i,1} = table(uniqueIndicesVertices, pointBaseFreqVertices, ...
                        'VariableNames', {'PointBaseIndex', 'Frequency'});
                    % Points
                    [uniqueIndicesPoints, ~, idxCountsPoints] = unique(currentPointBasePointsIdx);
                    pointBaseFreqPoints = accumarray(idxCountsPoints, 1);
                    obj.(typeName).pairs.frequencyPoints{i,1} = table(uniqueIndicesPoints, pointBaseFreqPoints, ...
                        'VariableNames', {'PointBaseIndex', 'Frequency'});
                end
            end


            % Filter frequencies based on intersection percentages
            timeIntersectFilter = NaN(length(intersectPctNames), 1);
            for k = 1:length(intersectPctNames)
                pctCount = obj.(typeName).(intersectPctNames{k}).count;
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
                        filteredFreqDefect{i} = obj.(typeName).pairs.frequencyDefect{i}( ...
                            obj.(typeName).pairs.frequencyDefect{i}.Frequency >= filterThreshold, :);
                        % Filter frequencyVertices and frequencyPoints
                        if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
                            filteredFreqVertices{i} = obj.(typeName).pairs.frequencyVertices{i}( ...
                                obj.(typeName).pairs.frequencyVertices{i}.Frequency >= filterThreshold, :);
                            filteredFreqPoints{i} = obj.(typeName).pairs.frequencyPoints{i}( ...
                                obj.(typeName).pairs.frequencyPoints{i}.Frequency >= filterThreshold, :);
                        end
                    end
                    % Store filtered tables
                    obj.(typeName).(intersectPctNames{k}).filteredFreqDefect = filteredFreqDefect;
                    if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
                        obj.(typeName).(intersectPctNames{k}).filteredFreqVertices = filteredFreqVertices;
                        obj.(typeName).(intersectPctNames{k}).filteredFreqPoints = filteredFreqPoints;
                    end
                    timeIntersectFilter(k,1) = toc(t2);
                end
            end         
            obj.(typeName).timeIntersectFilter = timeIntersectFilter;

            disp(['defect intersection "', (typeName), '": frequency tables for pelvis defect categorie "', (baseIntersect), '"']);

        end

        %% Function to allocate coordinates to defects, vertices, and points
        function obj = allocateCoords(obj, baseIntersect, typeName, inputIDs, pelvisDefect, refData, ...
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
            obj.(typeName).(currentPctName) = currentPct;
            % Data storage of intersection data
            sizePct = whos('currentPct');
            obj.(typeName).(currentPctName).sizeIntersectAllo = sizePct.bytes;

            disp(['defect intersection "', (typeName), '": coordinates for "', (baseIntersect), '" and "', (currentPctName), '"']);
        end

        %% Shrink Wrap with alphaShape (without user-input)
        function obj = shrinkAlpha(obj, baseIntersect, typeName, intersectName, aRadius, numDefect, numVertices, ...
                inputPoints, addPoints)
            
            % Combine input points and additional points
            allVerticesPoints = [inputPoints; addPoints];

            % Initialize
            isClosed = false;
            maxIterations = 150;
            iterationCount = 0;
            while ~isClosed && iterationCount < maxIterations
                try
                    % Create alphaShape and get the boundary facets
                    shrinkWrap = alphaShape(allVerticesPoints, aRadius);
                    shrinkFaces = boundaryFacets(shrinkWrap);
                    % Check if the mesh is closed
                    edgesList = sort([shrinkFaces(:,1) shrinkFaces(:,2);
                        shrinkFaces(:,2) shrinkFaces(:,3);
                        shrinkFaces(:,3) shrinkFaces(:,1)], 2);
                    [~, ~, edgeOccurrences] = unique(edgesList, 'rows', 'stable');
                    edgeCounts = accumarray(edgeOccurrences, 1);
                    % Check if any edge is not shared by exactly two triangles
                    isClosed = all(edgeCounts == 2);
                    % If the mesh is closed, break out of the loop (no increase of alpha radius)
                    if isClosed
                        break;
                    end
                catch
                    % If an error occurs, increase aRadius and try again
                    isClosed = false;
                end

                % Increment alpha radius and iteration count
                aRadius = aRadius + 1;
                iterationCount = iterationCount + 1;
            end

            % Store results
            obj.(typeName).(intersectName).radius = aRadius;
            %obj.(typeName).(intersectName).alpha = shrinkWrap;
            obj.(typeName).(intersectName).faces = shrinkFaces;
            obj.(typeName).(intersectName).allVerticesPoints = allVerticesPoints; % all vertices 
            % Calculate used vertices and defect count
            idxUsed = unique(shrinkFaces(:));
            idxDefect    = idxUsed(idxUsed <= numDefect);
            idxVertices  = idxUsed(idxUsed > numDefect & idxUsed <= numDefect+numVertices);
            idxPoints    = idxUsed(idxUsed > numDefect+numVertices);
            obj.(typeName).(intersectName).numOutputDefect   = numel(idxDefect);
            obj.(typeName).(intersectName).numOutputVertices = numel(idxVertices);
            obj.(typeName).(intersectName).numOutputPoints   = numel(idxPoints);
            obj.(typeName).(intersectName).usedDefectCount   = numel(idxDefect) / numDefect;
            obj.(typeName).(intersectName).usedCount = numel(idxUsed) / length(allVerticesPoints);

             % Volume (divergence theorem)
            volumeAlpha = 0; % Initialize volume
            % Loop through each face
            for i = 1:size(shrinkFaces, 1) 
                % Get the vertices of the face
                v1 = allVerticesPoints(shrinkFaces(i, 1), :); 
                v2 = allVerticesPoints(shrinkFaces(i, 2), :); 
                v3 = allVerticesPoints(shrinkFaces(i, 3), :); 
                % Calculate the signed volume of the tetrahedron
                tetraVolume = dot(v1, cross(v2, v3)) / 6;
                % Add to total volume
                volumeAlpha = volumeAlpha + tetraVolume;
            end
            % The volume should be positive
            obj.(typeName).(intersectName).volume = abs(volumeAlpha);

            disp(['alphaShape calculated: intersect type "', (typeName),'", pelvis defect categorie ', ...
                (baseIntersect), ' and ' ,(intersectName)])
           
        end
        %% Comparison of intersection types
         function obj = compareIntersect(obj, type1, type2, baseName, pctName, numIDs, type1Data, type2Data)
            
            combi = [type1 '_' type2];
            
            % Initialize storage
            obj.comp.(combi).(baseName).(pctName).uniqueDefectType1 = cell(numIDs, 1);
            obj.comp.(combi).(baseName).(pctName).uniqueDefectType2 = cell(numIDs, 1);
            obj.comp.(combi).(baseName).(pctName).sharedDefect = cell(numIDs, 1);
            obj.comp.(combi).(baseName).(pctName).lengthUniqueDefectType1 = zeros(numIDs, 1);
            obj.comp.(combi).(baseName).(pctName).lengthUniqueDefectType2 = zeros(numIDs, 1);
            obj.comp.(combi).(baseName).(pctName).lengthSharedDefect = zeros(numIDs, 1);
            if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
                % Vertices
                obj.comp.(combi).(baseName).(pctName).uniqueVerticesType1 = cell(numIDs, 1);
                obj.comp.(combi).(baseName).(pctName).uniqueVerticesType2 = cell(numIDs, 1);
                obj.comp.(combi).(baseName).(pctName).sharedVertices = cell(numIDs, 1);
                obj.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType1 = zeros(numIDs, 1);
                obj.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType2 = zeros(numIDs, 1);
                obj.comp.(combi).(baseName).(pctName).lengthSharedVertices = zeros(numIDs, 1);
                % Points
                obj.comp.(combi).(baseName).(pctName).uniquePointsType1 = cell(numIDs, 1);
                obj.comp.(combi).(baseName).(pctName).uniquePointsType2 = cell(numIDs, 1);
                obj.comp.(combi).(baseName).(pctName).sharedPoints = cell(numIDs, 1);
                obj.comp.(combi).(baseName).(pctName).lengthUniquePointsType1 = zeros(numIDs, 1);
                obj.comp.(combi).(baseName).(pctName).lengthUniquePointsType2 = zeros(numIDs, 1);
                obj.comp.(combi).(baseName).(pctName).lengthSharedPoints = zeros(numIDs, 1);
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
                obj.comp.(combi).(baseName).(pctName).uniqueDefectType1{b} = setdiff(defect1, defect2);
                obj.comp.(combi).(baseName).(pctName).uniqueDefectType2{b} = setdiff(defect2, defect1);
                obj.comp.(combi).(baseName).(pctName).sharedDefect{b} = intersect(defect1, defect2);
                obj.comp.(combi).(baseName).(pctName).lengthUniqueDefectType1(b) = numel(obj.comp.(combi).(baseName).(pctName).uniqueDefectType1{b});
                obj.comp.(combi).(baseName).(pctName).lengthUniqueDefectType2(b) = numel(obj.comp.(combi).(baseName).(pctName).uniqueDefectType2{b});
                obj.comp.(combi).(baseName).(pctName).lengthSharedDefect(b) = numel(obj.comp.(combi).(baseName).(pctName).sharedDefect{b});

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
                        obj.comp.(combi).(baseName).(pctName).sharedVertices{b} = [];
                        obj.comp.(combi).(baseName).(pctName).sharedPoints{b} = [];
                        obj.comp.(combi).(baseName).(pctName).lengthSharedVertices(b) = 0;
                        obj.comp.(combi).(baseName).(pctName).lengthSharedPoints(b) = 0;

                        obj.comp.(combi).(baseName).(pctName).uniqueVerticesType1{b} = vertices1;
                        obj.comp.(combi).(baseName).(pctName).uniqueVerticesType2{b} = vertices2;
                        obj.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType1(b) = numel(vertices1);
                        obj.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType2(b) = numel(vertices2);

                        obj.comp.(combi).(baseName).(pctName).uniquePointsType1{b} = points1;
                        obj.comp.(combi).(baseName).(pctName).uniquePointsType2{b} = points2;
                        obj.comp.(combi).(baseName).(pctName).lengthUniquePointsType1(b) = numel(points1);
                        obj.comp.(combi).(baseName).(pctName).lengthUniquePointsType2(b) = numel(points2);
                    else
                        % Compare Vertices
                        obj.comp.(combi).(baseName).(pctName).uniqueVerticesType1{b} = setdiff(vertices1, vertices2);
                        obj.comp.(combi).(baseName).(pctName).uniqueVerticesType2{b} = setdiff(vertices2, vertices1);
                        obj.comp.(combi).(baseName).(pctName).sharedVertices{b} = intersect(vertices1, vertices2);
                        obj.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType1(b) = numel(obj.comp.(combi).(baseName).(pctName).uniqueVerticesType1{b});
                        obj.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType2(b) = numel(obj.comp.(combi).(baseName).(pctName).uniqueVerticesType2{b});
                        obj.comp.(combi).(baseName).(pctName).lengthSharedVertices(b) = numel(obj.comp.(combi).(baseName).(pctName).sharedVertices{b});

                        % Compare Points
                        obj.comp.(combi).(baseName).(pctName).uniquePointsType1{b} = setdiff(points1, points2);
                        obj.comp.(combi).(baseName).(pctName).uniquePointsType2{b} = setdiff(points2, points1);
                        obj.comp.(combi).(baseName).(pctName).sharedPoints{b} = intersect(points1, points2);
                        obj.comp.(combi).(baseName).(pctName).lengthUniquePointsType1(b) = numel(obj.comp.(combi).(baseName).(pctName).uniquePointsType1{b});
                        obj.comp.(combi).(baseName).(pctName).lengthUniquePointsType2(b) = numel(obj.comp.(combi).(baseName).(pctName).uniquePointsType2{b});
                        obj.comp.(combi).(baseName).(pctName).lengthSharedPoints(b) = numel(obj.comp.(combi).(baseName).(pctName).sharedPoints{b});
                    end
                end
            end

            % Summarize results
            obj.comp.(combi).(baseName).(pctName).sumUniqueDefectType1 = sum(obj.comp.(combi).(baseName).(pctName).lengthUniqueDefectType1);
            obj.comp.(combi).(baseName).(pctName).sumUniqueDefectType2 = sum(obj.comp.(combi).(baseName).(pctName).lengthUniqueDefectType2);
            obj.comp.(combi).(baseName).(pctName).sumSharedDefect = sum(obj.comp.(combi).(baseName).(pctName).lengthSharedDefect);
            if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
                obj.comp.(combi).(baseName).(pctName).sumUniqueVerticesType1 = sum(obj.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType1);
                obj.comp.(combi).(baseName).(pctName).sumUniqueVerticesType2 = sum(obj.comp.(combi).(baseName).(pctName).lengthUniqueVerticesType2);
                obj.comp.(combi).(baseName).(pctName).sumSharedVertices = sum(obj.comp.(combi).(baseName).(pctName).lengthSharedVertices);
                obj.comp.(combi).(baseName).(pctName).sumUniquePointsType1 = sum(obj.comp.(combi).(baseName).(pctName).lengthUniquePointsType1);
                obj.comp.(combi).(baseName).(pctName).sumUniquePointsType2 = sum(obj.comp.(combi).(baseName).(pctName).lengthUniquePointsType2);
                obj.comp.(combi).(baseName).(pctName).sumSharedPoints = sum(obj.comp.(combi).(baseName).(pctName).lengthSharedPoints);
            end

            disp(['comparison intersect types "', (type1),'" - "',(type2), '": pelvis defect categorie ', ...
                (baseName), ' and ' ,(pctName)])
         end

         %% Comparison of intersection types: allocate coordinates
         function obj = allocateCoordsComp(obj, combi, baseName, pctName, inputIDs, pelvisDefect, pelvis, allPelvis, type1, type2, weightKabsch)
             numIDs = length(inputIDs);

             % Initialize storage in allDefect
             obj.comp.(combi).(baseName).(pctName).insideDefectType1 = cell(numIDs, 1);
             obj.comp.(combi).(baseName).(pctName).insideDefectType2 = cell(numIDs, 1);
             obj.comp.(combi).(baseName).(pctName).insideDefectShared = cell(numIDs, 1);
             if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
                 obj.comp.(combi).(baseName).(pctName).insideVerticesType1 = cell(numIDs, 1);
                 obj.comp.(combi).(baseName).(pctName).insideVerticesType2 = cell(numIDs, 1);
                 obj.comp.(combi).(baseName).(pctName).insideVerticesShared = cell(numIDs, 1);
                 obj.comp.(combi).(baseName).(pctName).insidePointsType1 = cell(numIDs, 1);
                 obj.comp.(combi).(baseName).(pctName).insidePointsType2 = cell(numIDs, 1);
                 obj.comp.(combi).(baseName).(pctName).insidePointsShared = cell(numIDs, 1);
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
                 defectIdx1 = obj.comp.(combi).(baseName).(pctName).uniqueDefectType1{l};
                 defectIdx2 = obj.comp.(combi).(baseName).(pctName).uniqueDefectType2{l};
                 defectShared = obj.comp.(combi).(baseName).(pctName).sharedDefect{l};
                 if ~isempty(defectIdx1)
                     obj.comp.(combi).(baseName).(pctName).insideDefectType1{l} = pointBaseDefectType1(defectIdx1, :);
                 else
                     obj.comp.(combi).(baseName).(pctName).insideDefectType1{l} = [];
                 end
                 if ~isempty(defectIdx2)
                     obj.comp.(combi).(baseName).(pctName).insideDefectType2{l} = pointBaseDefectType2(defectIdx2, :);
                 else
                     obj.comp.(combi).(baseName).(pctName).insideDefectType2{l} = [];
                 end
                 if ~isempty(defectShared) % Assuming shared defects are same in both (points of refined patches at the end)
                     obj.comp.(combi).(baseName).(pctName).insideDefectShared{l} = pointBaseDefectType1(defectShared, :); % Assuming shared indices are valid for type1
                 else
                     obj.comp.(combi).(baseName).(pctName).insideDefectShared{l} = [];
                 end
                 
                 % Allocate vertices and points for inside/grid
                 if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
                     verticesIdx1 = obj.comp.(combi).(baseName).(pctName).uniqueVerticesType1{l};
                     verticesIdx2 = obj.comp.(combi).(baseName).(pctName).uniqueVerticesType2{l};
                     verticesShared = obj.comp.(combi).(baseName).(pctName).sharedVertices{l};
                     pointsIdx1 = obj.comp.(combi).(baseName).(pctName).uniquePointsType1{l};
                     pointsIdx2 = obj.comp.(combi).(baseName).(pctName).uniquePointsType2{l};
                     pointsShared = obj.comp.(combi).(baseName).(pctName).sharedPoints{l};

                     if ~isempty(verticesIdx1)
                         obj.comp.(combi).(baseName).(pctName).insideVerticesType1{l} = pointBaseVerticesType1(verticesIdx1, :);
                     else
                         obj.comp.(combi).(baseName).(pctName).insideVerticesType1{l} = [];
                     end
                     if ~isempty(verticesIdx2)
                         obj.comp.(combi).(baseName).(pctName).insideVerticesType2{l} = pointBaseVerticesType2(verticesIdx2, :);
                     else
                         obj.comp.(combi).(baseName).(pctName).insideVerticesType2{l} = [];
                     end
                     if ~isempty(verticesShared)
                         obj.comp.(combi).(baseName).(pctName).insideVerticesShared{l} = pointBaseVerticesType1(verticesShared, :);
                     else
                         obj.comp.(combi).(baseName).(pctName).insideVerticesShared{l} = [];
                     end
                     if ~isempty(pointsIdx1)
                         obj.comp.(combi).(baseName).(pctName).insidePointsType1{l} = pointBasePointsType1(pointsIdx1, :);
                     else
                         obj.comp.(combi).(baseName).(pctName).insidePointsType1{l} = [];
                     end
                     if ~isempty(pointsIdx2)
                         obj.comp.(combi).(baseName).(pctName).insidePointsType2{l} = pointBasePointsType2(pointsIdx2, :);
                     else
                         obj.comp.(combi).(baseName).(pctName).insidePointsType2{l} = [];
                     end
                     if ~isempty(pointsShared)
                         obj.comp.(combi).(baseName).(pctName).insidePointsShared{l} = pointBasePointsType1(pointsShared, :);
                     else
                         obj.comp.(combi).(baseName).(pctName).insidePointsShared{l} = [];
                     end

                 end
             end
         end
    end
 end