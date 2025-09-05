 classdef Transform
    % Transform class for the transformation of a model (stl)
    
    % Developed by C.Micheler,
    % Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
    % Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich
    
    
    properties 
        scaled = struct();  % Scaled model
        trans = struct();   % transformed model (translation)
        trafo = struct();   % transformed model (translation + rotation)   
    end
    
    methods
        %% Constructer: generate object
        function obj = Transform()
            disp('class transform initialized')
        end
        
        %% Scale
        function obj = scale(obj,pelvisNum,import,scaleFactor)
            
            % Scale Mesh: https://de.mathworks.com/help/driving/ref/extendedobjectmesh.scale.html
            % Attention: Scaling leads to new position of the geometry
            mesh = extendedObjectMesh(import.vertices,import.faces);
            scaledMesh = scale(mesh,scaleFactor);
            obj.scaled.vertices = scaledMesh.Vertices;
            obj.scaled.faces = scaledMesh.Faces;
           
            % List of all landmarks
            allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', ...
                'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};

            % Scale available landmarks (for transformation)
            for idx = 1:length(allLandmarks)
                landmark = allLandmarks{idx};
                if isfield(import, landmark)
                    obj.scaled.(landmark) = import.(landmark) * scaleFactor;
                end
            end

            % Normals
            obj.scaled.normals = import.normals;

            % Calculate the centroid of the triangular faces in the mesh
            V1 = scaledMesh.Vertices(scaledMesh.Faces(:,1),:);
            V2 = scaledMesh.Vertices(scaledMesh.Faces(:,2),:);
            V3 = scaledMesh.Vertices(scaledMesh.Faces(:,3),:);
            faceVertices = cat(3, V1, V2, V3);
            obj.scaled.centreFaces = mean(faceVertices, 3);

            disp(['model scaled: pelvis ',num2str(pelvisNum)])

        end

        %% Translation to reference point
        function obj = translation(obj,pelvisNum,refPoint,importData)

            % Translation
            obj.trans.transVector = refPoint - importData.acentre;
            obj.trans.vertices = importData.vertices + obj.trans.transVector;

            % List of all landmarks
            allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre',... 
                'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};

            % Translate available landmarks
            for idx = 1:length(allLandmarks)
                landmark = allLandmarks{idx};
                if isfield(importData, landmark)
                    if strcmp(landmark, 'acentre')
                        obj.trans.(landmark) = refPoint;
                    else
                        obj.trans.(landmark) = importData.(landmark) + obj.trans.transVector;
                    end
                end
            end

            % Faces / Normals
            obj.trans.faces = importData.faces;
            obj.trans.normals = importData.normals;
            % Calculate the centroid of the triangular faces in the mesh
            V1 = obj.trans.vertices(importData.faces(:,1),:);
            V2 = obj.trans.vertices(importData.faces(:,2),:);
            V3 = obj.trans.vertices(importData.faces(:,3),:);
            faceVertices = cat(3, V1, V2, V3);
            obj.trans.centreFaces = mean(faceVertices, 3);

            disp(['vertices shifted: pelvis ',num2str(pelvisNum)])
        end   

        %% Transformation: Translation and Rotation / Orientation
        function obj = trafoKabsch(obj, pelvisNum, refLandmarks, weight, weightName, importData)

            allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', ...
                'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};

            % Collect landmarks
            Landmarks = zeros(length(weight), 3); % weight defines number of landmarks
            for idx = 1:length(weight)
                if isfield(importData, allLandmarks{idx})
                    Landmarks(idx, :) = importData.(allLandmarks{idx});
                else
                    Landmarks(idx, :) = NaN;
                end
            end
            % Filter out rows with NaN values
            validRows = ~any(isnan(Landmarks), 2);
            validIndices = find(validRows); % Get the valid indices (those not NaN)
            Landmarks = Landmarks(validRows, :);
            refLandmarks = refLandmarks(validRows, :);
            weight = weight(validRows);

            % Kabsch
            [obj.trafo.(weightName).R, obj.trafo.(weightName).t, obj.trafo.(weightName).lrms] = ...
                Kabsch(Landmarks', refLandmarks', weight);
            rotated_Kabsch = obj.trafo.(weightName).R * Landmarks' + obj.trafo.(weightName).t;

            % Euler angles (z,y,x)
            obj.trafo.(weightName).euler = rad2deg(rotm2eul(obj.trafo.(weightName).R));

            % Assign to obj.rotation using the valid landmarks list
            for k = 1:length(validIndices)
                idx = validIndices(k);  % Retrieve the actual landmark index
                obj.trafo.(weightName).(allLandmarks{idx}) = rotated_Kabsch(:, k)';  % Assign the rotated value
            end
            % For the landmarks that weren't valid, assign NaN
            invalidIndices = setdiff(1:length(allLandmarks), validIndices);
            for idx = invalidIndices
                obj.trafo.(weightName).(allLandmarks{idx}) = NaN(1, 3);
            end

            % Transform vertices
            rotatedVertices = obj.trafo.(weightName).R * importData.vertices' + obj.trafo.(weightName).t;
            obj.trafo.(weightName).vertices = rotatedVertices';
            % Faces / Normals
            obj.trafo.(weightName).faces = importData.faces;
            obj.trafo.(weightName).normals = importData.normals;
            % Calculate the centroid of the triangular faces in the mesh
            V1 = obj.trafo.(weightName).vertices(importData.faces(:,1),:);
            V2 = obj.trafo.(weightName).vertices(importData.faces(:,2),:);
            V3 = obj.trafo.(weightName).vertices(importData.faces(:,3),:);
            faceVertices = cat(3, V1, V2, V3);
            obj.trafo.(weightName).centreFaces = mean(faceVertices, 3);

            % Landmark errors
            for idx = 1:length(allLandmarks)
                landmark = allLandmarks{idx};
                if any(validIndices == idx)
                    % Calculate Euclidean distance between transformed and reference landmarks
                    transformedLandmark = obj.trafo.(weightName).(landmark);
                    referenceLandmark = refLandmarks(validIndices == idx, :);
                    obj.trafo.(weightName).(['error_' landmark]) = norm(transformedLandmark - referenceLandmark);
                else
                    % Landmark was filtered out, assign NaN and note it
                    obj.trafo.(weightName).(['error_' landmark]) = NaN;
                end
            end

            disp(['vertices transformed: pelvis ', num2str(pelvisNum)]);
        end

        %% Vertex-to-Nearest-Neighbour Distance
        function obj = vertexNeighbour(obj, pelvisNum, verticesRef, weightName, vertices, faces)

            % Vertex-to-Nearest-Neighbour calculates the distance from each vertex to its
            % nearest neighbor in the vertices list.
            % Using a KD-tree searcher for improved performance
            % Minimum distances from each vertex to the nearest vertex in reference

            % Create a KD-tree searcher for verticesB
            MdlB = KDTreeSearcher(verticesRef);

            % Find the nearest neighbor in verticesB for each vertex in verticesA
            idx = knnsearch(MdlB,vertices);

            % Calculate the distances based on the nearest neighbor indices
            minDistances = sqrt(sum((vertices - verticesRef(idx,:)).^2, 2));

            % Pelvis with vertex-to-nearest-neighbour color
            % Retrieve nearVertex for each vertex of each face
            nearVertexMatrix = minDistances(faces);
            % Calculate mean for each face
            obj.trafo.(weightName).nearVertexFace = mean(nearVertexMatrix,2);
            % Normalisation
            lowerBound = min(obj.trafo.(weightName).nearVertexFace);
            upperBound = max(obj.trafo.(weightName).nearVertexFace);
            % Normalize data within this range
            normalData = (obj.trafo.(weightName).nearVertexFace - lowerBound) / (upperBound - lowerBound);
            normalData(normalData < 0) = 0;   % Clamp values below the range to 0
            normalData(normalData > 1) = 1;   % Clamp values above the range to 1
            % Apply colormap (recommended: viridis)
            % 16-bit colour: alpha: 1bit R: 5bit G: 5bit B: 5bit -> 32768 colors
            colourNum = 32768;
            colourMap = viridis(colourNum);
            colourIdx = round(normalData * (colourNum-1)) + 1;
            rgbColour = colourMap(colourIdx,:);

            obj.trafo.(weightName).nearVertexFaceNorm = rgbColour;
            obj.trafo.(weightName).nearVertex = minDistances;
            obj.trafo.(weightName).nearVertexMax = max(minDistances);
            obj.trafo.(weightName).nearVertexMean = mean(minDistances);
            obj.trafo.(weightName).nearVertexStd = std(minDistances);

            disp(['Vertex-to-Nearest-Neighbour: pelvis ', num2str(pelvisNum)]);

        end

        %% Hausdorff Distance
        function obj = hausdorffDistance(obj, pelvisNum, verticesRef, weightName, vertices)

            % Hausdorff-Distance calculates the Hausdorff distance between two point clouds.
            % This is done by calculating the maximum of the minimum distances from each vertex in
            % vertices to the nearest vertex in verticesRef.

            % Create KD-tree searchers for both vertex sets
            MdlRef = KDTreeSearcher(verticesRef);
            MdlTest = KDTreeSearcher(vertices);

            % First, calculate distance from each vertex in 'vertices' to its nearest neighbor in 'verticesRef'
            idxRef = knnsearch(MdlRef, vertices);
            distancesRef = sqrt(sum((vertices - verticesRef(idxRef,:)).^2, 2));

            % Second, calculate distance from each vertex in 'verticesRef' to its nearest neighbor in 'vertices'
            idxTest = knnsearch(MdlTest, verticesRef);
            distancesTest = sqrt(sum((verticesRef - vertices(idxTest,:)).^2, 2));

            % Calculate the Hausdorff Distance
            % Maximum of the minimal distances in both directions
            hausdorffDist = max(max(distancesRef), max(distancesTest));

            % Store results in obj structure
            obj.trafo.(weightName).hausdorffDist = hausdorffDist;

            disp(['Hausdorff Distance: pelvis ', num2str(pelvisNum)]);

        end

    end
    
 end