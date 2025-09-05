classdef Curvature 
    % Curvature class to calculate the curvature of a stl area

    % Developed by C.Micheler,
    % Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
    % Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

    properties      

        normal = struct();      % Curvature normal
        cMean = struct();       % Curvature mean
        gauss = struct();       % Curvature gauss 

        boundFaces = [];
        boundFacesNr = [];
        boundVertices = []; 
        boundVerticesNr = []; 

        allVertexFaceMap = {[]};  
        vertexFaceMap = {[]}; 
        vertexMap = {[]}; 

        vertexAdjMap = {[]}; 
        faceAdjMap = [];
        edgeAdjMap = [];

        vertexComPoint = {[]};  
        comEdge = {[]};
        comEdgeNorm = {[]};

        baryAreaFaces = [];
        areaFaceMap = [];

        angleFaceNormal = {[]};
     
    end
    
    methods
        %% Constructer: generate object
        function obj = Curvature()
            disp('class curvature initialized')
        end

        %% Identify boundary vertices / faces
        function obj = bound(obj,pelvisNum,importData) % derived from function edging in Volume

            % Identify boundary vertices / faces
            % Boundary edge: triangular side (combination of two vertices) only once in row of faces list

            faces = importData.ConnectivityList;
            vertices = importData.Points;

            % Extract all edges from the faces
            edges = [faces(:, [1, 2]); faces(:, [2, 3]); faces(:, [3, 1])];
            sortedEdges = sort(edges, 2);

            % Identify unique edges and their counts
            [uniqueEdges, ~, edgeIndices] = unique(sortedEdges, 'rows');
            edgeCounts = accumarray(edgeIndices, 1);

            % Boundary edges appear only once
            isBoundaryEdge = edgeCounts == 1;
            boundaryEdges = uniqueEdges(isBoundaryEdge, :);

            % Precompute face edges for all faces
            faceEdges = [sort(faces(:, [1, 2]), 2), (1:size(faces, 1))'; ...
                sort(faces(:, [2, 3]), 2), (1:size(faces, 1))'; ...
                sort(faces(:, [3, 1]), 2), (1:size(faces, 1))'];

            % Find faces that include boundary edges
            boundaryFaces = ismember(faceEdges(:, 1:2), boundaryEdges, 'rows');
            boundaryFaceIndices = faceEdges(boundaryFaces, 3);

            % Extract boundary faces
            obj.boundFacesNr = unique(boundaryFaceIndices);
            obj.boundFaces = faces(obj.boundFacesNr,:);

            % Store boundary vertices
            obj.boundVerticesNr = unique(boundaryEdges(:)); % Flatten the array and then find unique vertices
            obj.boundVertices = vertices(obj.boundVerticesNr,:);

            disp(['Boundary identified: pelvis ',num2str(pelvisNum)])

        end

        %% Identify adjacent vertices of a centre vertex (vertex face map)
        function obj = adjVertex(obj,pelvisNum,verticesNrROI,facesROI)

            % Precompute mapping from vertices to faces
            tempAllVertexFaceMap = cell(max(facesROI(:)),1);
            for f = 1:size(facesROI,1)
                for v = facesROI(f,:)
                    tempAllVertexFaceMap{v} = [tempAllVertexFaceMap{v},f];
                end
            end
            obj.allVertexFaceMap = tempAllVertexFaceMap; 

            % Preallocation cell array (cells with different size)
            verticesCount = size(verticesNrROI,1);
            tempVertexFaceMap = cell(verticesCount,1);
            tempVertexMap = cell(verticesCount,1);

            % Map faces to vertices
            for i = 1:verticesCount
                % Face numbers (row) with corresponding vertice in it
                % Cell row = row number in verticesNrROI -> corresponding vertice
                vertex = verticesNrROI(i,1);
                if ismember(vertex, obj.boundVerticesNr) % for open mesh
                    % Skip if vertex is a boundary vertex
                    continue;
                end
                tempVertexFaceMap{i} = tempAllVertexFaceMap{vertex};
                % Faces with the vertices number
                vertexFaces = facesROI(tempVertexFaceMap{i},:)';
                % Adjacent Vertices (Remove the central vertice)
                vertexFaces(vertexFaces == vertex) = [];
                tempVertexMap{i} = reshape(vertexFaces,2,[]).';
            end

            % Store in obj (creates overhead in loop)
            obj.vertexFaceMap = tempVertexFaceMap; 
            obj.vertexMap = tempVertexMap;

            disp(['Vertex face map: pelvis ', num2str(pelvisNum)]);

       end
      
       %% Identify sorted adjacent faces around centre vertice (circle) 
       function obj = adjVertexSeq(obj,pelvisNum,verticesNrROI)     

            verticesCount = size(verticesNrROI,1);

            % Preallocation cell array (cells with different size)
            obj.vertexAdjMap = cell(verticesCount,1); % Sorted faces structure of the centre vertex (circle of adjacent faces)
            obj.vertexComPoint = cell(verticesCount,1); % Second common point (sorted)
            % Temporary variables for storing results
            tempVertexAdjMap = cell(verticesCount,1);
            tempVertexComPoint = cell(verticesCount,1);
            
            for i = 1:verticesCount

                % Skip if vertex is a boundary vertex
                if ismember(verticesNrROI(i,1), obj.boundVerticesNr) % for open mesh
                    continue;
                end

                % Add row number of the faces / row number of the normals
                % Column 1-2: row numbers of the vertices of the faces; without centre vertice
                % Column 3: row number of the faces = row number of the normals
                currVertexFaces = obj.vertexMap{i,1};
                currVertexFaces(:,3) = obj.vertexFaceMap{i,1};
                % Initialization
                sortedFaces = zeros(size(currVertexFaces));
                comPoints = zeros(size(currVertexFaces,1),1);
                
                % Start with the first face
                sortedFaces(1,:) = currVertexFaces(1,:);
                comPoints(1) = currVertexFaces(1,1);
                nextVertex = currVertexFaces(1,1);

                % Use a logical array to keep track of remaining faces
                remaining = true(size(currVertexFaces,1),1);
                remaining(1) = false;  % First face is already used

                % Reorder remaining faces based on adjacency
                for j = 2:size(currVertexFaces,1)
                    nextRow = find((currVertexFaces(:,1) == nextVertex | currVertexFaces(:,2) == nextVertex) & remaining, 1);
                    if isempty(nextRow)
                        break; % Break the loop if no next face is found
                    end
                    sortedFaces(j,:) = currVertexFaces(nextRow,:);
                    nextVertex = setdiff(currVertexFaces(nextRow, 1:2), nextVertex);
                    comPoints(j,1) = nextVertex;
                    remaining(nextRow) = false;  % Mark as used
                end

                % Store
                tempVertexAdjMap{i,1} = sortedFaces;
                tempVertexComPoint{i,1} = comPoints;
            end
            
            % Store in obj (creates overhead in loop)
            obj.vertexAdjMap = tempVertexAdjMap;
            obj.vertexComPoint = tempVertexComPoint;

            disp(['Vertex face map sorted: pelvis ', num2str(pelvisNum)]);

       end
       
       %% Calculate vector of every edge in ROI and the area of the faces
       % For mean (2) and gauss (3) curvature
       function obj = edgesArea(obj,pelvisNum,verticesNrROI,verticesROI)

           verticesCount = size(verticesNrROI, 1);

           % Coordinates of the centre vertice
           vertices = verticesROI(verticesNrROI(:, 1), :);

           % Preallocate arrays
           % Coordinates of the adjacent vertices
           adjVertices = cell(verticesCount,1);
           % Vector of the edges
           edges = cell(verticesCount,1);
           % Common edge
           tempComEdge = cell(verticesCount,1);
           tempComEdgeNorm = cell(verticesCount,1);
           % Calculate the area of the faces
           crossProduct = cell(verticesCount,1);
           areaFaces = cell(verticesCount,1);
           tempSumAreaFaces = zeros(verticesCount,1);

           for i = 1:verticesCount

               % Skip if vertex is a boundary vertex
               if ismember(verticesNrROI(i,1), obj.boundVerticesNr) % for open mesh
                   continue;
               end

               % Coordinates of the adjacent vertices
               % Vertice for edge1 of the face
               adjVertices{i,1}(:,:,1) = verticesROI(obj.vertexAdjMap{i,1}(:,1),:);
               % Vertice for edge2 of the face
               adjVertices{i,1}(:,:,2) = verticesROI(obj.vertexAdjMap{i,1}(:,2),:);
               % Calculate vector of the adjacent edges
               % edge{:,:}(:,:,1) -> edge1 of the face
               % edge{:,:}(:,:,2) -> edge2 of the face
               edges{i,1} = adjVertices{i,1} - vertices(i,:);

               % Common edge and its norm 
               % Circle back to the first element
               tempComEdge{i,1} = verticesROI(obj.vertexComPoint{i,1}(:,1),:) - vertices(i,:);
               tempComEdgeNorm{i,1} = vecnorm(tempComEdge{i,1},2,2);

               % Calculate the area of the faces
               % Norm of the cross product of two vectors correspond to the area of the parallelogram
               % -> The area of the triangle is the half
               % |a x b| = |a||b|sin(phi)
               % normal vector = cross product
               % normals of the stl-file: length of the normals normalized to length 1
               % -> new calculation with the original edges
               crossProduct{i,1} = cross(edges{i,1}(:,:,1),edges{i,1}(:,:,2));
               areaFaces{i,1} = 0.5*vecnorm(crossProduct{i,1},2,2);
               tempSumAreaFaces(i,1) = sum(areaFaces{i,1});
           end

           % Store in obj
           obj.comEdge = tempComEdge; 
           obj.comEdgeNorm = tempComEdgeNorm; 
           obj.baryAreaFaces = tempSumAreaFaces / 3;
           
           disp(['Edges and face areas: pelvis ', num2str(pelvisNum)]);

      end

      %% Calculate angle between the faces (normal vector) 
      % For (1) normal and (2) mean curvature
      function obj = angleFaces(obj,pelvisNum,verticesNrROI,normalsROI)

          verticesCount = size(verticesNrROI,1);

          % Preallocation cell array (cells with different size)
          tempAngleFaceNormal = cell(verticesCount,1);
          vertexNormalMap = cell(verticesCount,1);
          for i = 1:verticesCount

              % Skip if vertex is a boundary vertex
              if ismember(verticesNrROI(i,1), obj.boundVerticesNr) % for open mesh
                  continue;
              end

              % Adjacent normals
              vertexNormal = normalsROI(obj.vertexAdjMap{i,1}(:,3),:); % sorted
              % Store adjusted normals
              vertexNormalMap{i,1} = vertexNormal;
              % Circle: first normal at the end for circular calculation (angles between the normal vectors)
              vertexNormal = [vertexNormal; vertexNormal(1,:)];

              % Calculate angles between adjacent normals/faces (in degree); cutting angle
              % angle = acosd(dot(u,v)/(norm(u)*norm(v)))
              % cos(phi)=dot(u,v)/(norm(u)*norm(v))
              dotProducts = sum(vertexNormal(1:end-1,:) .* vertexNormal(2:end,:),2);
              norms = vecnorm(vertexNormal(1:end-1,:),2,2) .* vecnorm(vertexNormal(2:end,:),2,2);
              fraction = dotProducts ./ norms;
              % Due to rounding, the value may lie outside the value range of cos [-1, 1]
              % hence the limitation of the value
              limitFraction = min(1, max(-1, fraction)); 
              angles = acosd(limitFraction);

              % Store angles
              tempAngleFaceNormal{i,1} = angles;
          end

          obj.angleFaceNormal = tempAngleFaceNormal; 

          disp(['Face angles: pelvis ', num2str(pelvisNum)]);

      end

      %% Adjacent faces: Find faces that share two vertices with the current face 
      % For (1) normal and (2) mean curvature
      function obj = adjFaces(obj,pelvisNum,facesROI,verticesROI)
          
          % Preallocation
          tempFaceAdjMap = zeros(size(facesROI,1),3);
          tempEdgeAdjMap = zeros(size(facesROI,1),3);
          tempAreaFace = zeros(size(facesROI,1),1);
          crossProduct = zeros(size(facesROI,1),3);
          facesNrROI = (1:size(facesROI,1))';

          for i = 1:size(facesROI,1)

              % Skip if face is a boundary face
              if ismember(facesNrROI(i,1), obj.boundFacesNr) % for open mesh
                  continue;
              end

              currentRow = facesROI(i,:);
              % 3 cells for the 3 vertices of currentRow (one cell per vertex -> cell: faceNr which contains the vertex)
              vertexCells = obj.allVertexFaceMap(currentRow); 
              % Concatenate cell contents
              allNumbers = [vertexCells{:}];
              % Count occurrences
              uniqueNumbers = unique(allNumbers);
              counts = histc(allNumbers, uniqueNumbers);
              % Adjacent faces: two common vertices -> two occurrences of the FaceNr
              faceNrTwice = uniqueNumbers(counts == 2);
              % Store adjacent faces
              % row N -> 3 facesNrs of adjacent faces to faceNr N
              tempFaceAdjMap(i,:) = faceNrTwice;

              % Map of common edges (norm) of adjacent faces 
              tempEdges = zeros(3,3);
              for j = 1:3
                  allNumbersE = [currentRow facesROI(tempFaceAdjMap(i,j),:)];
                  uniqueNumbersE = unique(allNumbersE);
                  countsE = histc(allNumbersE, uniqueNumbersE);
                  % Common edge: two common vertices -> two occurrences of the verticeNr
                  verticeNrTwice = uniqueNumbersE(countsE == 2);
                  tempEdges(j,:) = verticesROI(verticeNrTwice(2),:) - verticesROI(verticeNrTwice(1),:);
                  tempEdgeAdjMap(i,j) = vecnorm(tempEdges(j,:));                  
              end
              crossProduct(i,:) = cross(tempEdges(1,:),tempEdges(2,:));
              tempAreaFace(i,1) = 0.5*vecnorm(crossProduct(i,:),2,2);
          end

          % Sum of the adjacent faces (without centre face)
          for i = 1:size(facesROI, 1)
              if ismember(facesNrROI(i, 1), obj.boundFacesNr)
                  tempSumAreaFaces(i, 1) = 0; % boundary faces (for open mesh)
              else
                  tempSumAreaFaces(i, 1) = sum(tempAreaFace(tempFaceAdjMap(i, :)));
              end
          end

          obj.faceAdjMap = tempFaceAdjMap;
          obj.edgeAdjMap = tempEdgeAdjMap;
          obj.areaFaceMap = tempSumAreaFaces / 3;

          disp(['Adjacent faces: pelvis ', num2str(pelvisNum)]);

      end
      
      %% (1) Curvature: normal curvature (input: angleFaceNormal)
      function obj = curveNormal(obj,pelvisNum,pelvisID,normalsROI,facesROI,verticesROI,topCurveRange)

          % stl
          TR = triangulation(facesROI,verticesROI);

          %%%% Curvature per vertex (allocated faces in previous calculations)); also for open mesh %%%%
          % Mean angle of every vertice (angle between the faces / normal vector)
          curveVertex = nan(size(obj.angleFaceNormal,1),1); % Preallocation
          for i=1:size(obj.angleFaceNormal,1)
              curveVertex(i,1) = mean(obj.angleFaceNormal{i,1}); 
          end
          obj.normal.vertex = curveVertex;
          % Convert to faces (for stlwrite)
          % Retrieve curvature for each vertex of each face
          vertexCurvatures = curveVertex(facesROI);
          % Calculate mean curvature for each face
          obj.normal.vertexFace = mean(vertexCurvatures,2);
          % Curvature range (vertices)
          [ascendCurve,ascendCurveIdx] = sort(curveVertex);
          curveRange(2,1) = length(ascendCurve);
          curveRange(1,1) = curveRange(2,1) - round((topCurveRange/100)*curveRange(2,1));
          obj.normal.topCurveVertexIdx = sort(ascendCurveIdx(curveRange(1):curveRange(2),1));
          
          % Coloured stl (curvature): Write curvature in stl (stl binary) as attribute byte count
          stlName = 'curveNormalVertexFace';
          curveName = 'normVertexFace';
          obj = obj.colourSTL(TR,obj.normal.vertexFace,pelvisID,stlName,'normal',curveName);

          %%%% Curvature per face (3 edges); also for open mesh %%%%
          % Mean angle of every face (angle between the faces / normal vector)
          angles = nan(size(obj.faceAdjMap,1),3); % Preallocate matrix for angles
          facesNrROI = (1:size(facesROI,1))';
          for i = 1:size(obj.faceAdjMap,1)
              % Skip if face is a boundary face
              if ismember(facesNrROI(i,1), obj.boundFacesNr) % for open mesh
                  continue;
              end
              centerNormal = normalsROI(i,:);
              adjNormals = normalsROI(obj.faceAdjMap(i,:),:);
              % Vectorized angle calculation
              cosTheta = dot(repmat(centerNormal,3,1), adjNormals, 2) ./ ...
                  (vecnorm(centerNormal,2,2) .* vecnorm(adjNormals,2,2));
              % Due to rounding, the value may lie outside the value range of cos [-1, 1]
              % hence the limitation of the value
              limitCosTheta = min(1, max(-1, cosTheta));
              angles(i,:) = acosd(limitCosTheta); % Convert to degrees
          end
          obj.normal.face = mean(angles,2); 
          % Curvature range (faces)
          [ascendCurve,ascendCurveIdx] = sort(obj.normal.face);
          curveRange(2,1) = length(ascendCurve);
          curveRange(1,1) = curveRange(2,1) - round((topCurveRange/100)*curveRange(2,1));
          faceRangeIdx = sort(ascendCurveIdx(curveRange(1):curveRange(2),1));
          vertexIndices = facesROI(faceRangeIdx,:);
          flatVertexIndices = vertexIndices(:);
          obj.normal.topCurveFaceIdx= unique(flatVertexIndices);

          % Coloured stl (curvature): Write curvature in stl (stl binary) as attribute byte count
          stlName = 'curveNormalFace';
          curveName = 'normFace';
          obj = obj.colourSTL(TR,obj.normal.face,pelvisID,stlName,'normal',curveName);

          disp(['Normal curvature: pelvis ', num2str(pelvisNum)]);

      end

      %% (2) Curvature: mean curvature 
      % Mean curvature (per vertex) by Subburaj
      function obj = curveMean(obj,pelvisNum,pelvisID,normalsROI,facesROI,verticesROI,topCurveRange)
          
          % stl
          TR = triangulation(facesROI,verticesROI);

          %%%% Curvature per vertex (allocated faces in previous calculations); also for open mesh %%%%
          curveVertex = nan(size(obj.angleFaceNormal,1),1); % Preallocation
          for i=1:size(obj.angleFaceNormal,1)
              curveVertex(i,1) = (sum(obj.comEdgeNorm{i,1}.*abs(obj.angleFaceNormal{i,1})))/...
                  obj.baryAreaFaces(i,1); % Area around vertex: baryAreaFaces
          end
          obj.cMean.vertex = curveVertex;
          % Convert to faces (for stlwrite) -> Curvature method to identify defect area
          % Retrieve curvature for each vertex of each face
          vertexCurvatures = curveVertex(facesROI);
          % Calculate mean curvature for each face
          obj.cMean.vertexFace = mean(vertexCurvatures, 2);
          % Curvature range (vertices)
          [ascendCurve,ascendCurveIdx] = sort(curveVertex);
          curveRange(2,1) = length(ascendCurve);
          curveRange(1,1) = curveRange(2,1) - round((topCurveRange/100)*curveRange(2,1));
          obj.cMean.topCurveVertexIdx = sort(ascendCurveIdx(curveRange(1):curveRange(2),1));
          
          % Coloured stl (curvature): Write curvature in stl (stl binary) as attribute byte count
          stlName = 'curveMeanVertexFace';
          curveName = 'normVertexFace';
          obj = obj.colourSTL(TR,obj.cMean.vertexFace,pelvisID,stlName,'cMean',curveName);


          %%%% Curvature per face (3 edges); modified mean curvature; also for open mesh %%%%
          angles = nan(size(obj.faceAdjMap,1),3); % Preallocate matrix for angles
          curveFace = nan(size(obj.faceAdjMap,1),1);
          facesNrROI = (1:size(facesROI,1))';
          for i = 1:size(obj.faceAdjMap,1)
              % Skip if face is a boundary face
              if ismember(facesNrROI(i,1), obj.boundFacesNr) % for open mesh
                  continue;
              end
              centerNormal = normalsROI(i,:);
              adjNormals = normalsROI(obj.faceAdjMap(i,:),:);
              % Vectorized angle calculation
              cosTheta = dot(repmat(centerNormal,3,1), adjNormals, 2) ./ ...
                  (vecnorm(centerNormal,2,2) .* vecnorm(adjNormals,2,2));
              % Due to rounding, the value may lie outside the value range of cos [-1, 1]
              % hence the limitation of the value
              limitCosTheta = min(1, max(-1, cosTheta));
              angles(i,:) = acosd(limitCosTheta); % Convert to degrees
              % Curvature calculation
              curveFace(i,1) = sum(obj.edgeAdjMap(i, :) .* abs(angles(i, :))) / obj.areaFaceMap(i);
          end
          obj.cMean.face = curveFace;       
          % Curvature range (faces)
          [ascendCurve,ascendCurveIdx] = sort(obj.cMean.face);
          curveRange(2,1) = length(ascendCurve);
          curveRange(1,1) = curveRange(2,1) - round((topCurveRange/100)*curveRange(2,1));
          faceRangeIdx = sort(ascendCurveIdx(curveRange(1):curveRange(2),1));
          vertexIndices = facesROI(faceRangeIdx,:);
          flatVertexIndices = vertexIndices(:);
          obj.cMean.topCurveFaceIdx = unique(flatVertexIndices);

          % Coloured stl (curvature): Write curvature in stl (stl binary) as attribute byte count
          stlName = 'curveMeanFace';
          curveName = 'normFace';
          obj = obj.colourSTL(TR,curveFace,pelvisID,stlName,'cMean',curveName);

          disp(['Mean curvature: pelvis ', num2str(pelvisNum)]);

      end

      %% (3) Curvature: gauss curvature 
      % Gauss curvature (per vertex) by Subburaj
      function obj = curveGauss(obj,pelvisNum,pelvisID,verticesNrROI,facesROI,verticesROI,topCurveRange)

          % stl
          TR = triangulation(facesROI,verticesROI);

          %%%% Curvature per vertex (allocated faces in previous calculations) %%%%
          % Angle between edges
          % angle = acosd(dot(u,v)/(norm(u)*norm(v))); % cos(phi)=dot(u,v)/(norm(u)/*norm(v))
          angleEdge_ROI = cell(size(obj.comEdge, 1), 1);
          for i = 1:size(obj.comEdge, 1)
              % Skip if vertex is a boundary vertex
              if ismember(verticesNrROI(i,1), obj.boundVerticesNr) % for open mesh
                  continue;
              end
              u = obj.comEdge{i}(1:end,:);
              v = [obj.comEdge{i}(2:end,:); obj.comEdge{i}(1,:)];
              dotProduct = sum(u .* v, 2);
              normsProduct = vecnorm(u, 2, 2) .* vecnorm(v, 2, 2);
              fraction = dotProduct ./ normsProduct;
              % Due to rounding, the value may lie outside the value range of cos [-1, 1]
              % hence the limitation of the value
              limitFraction = min(1, max(-1, fraction));
              angleEdge_ROI{i} = acosd(limitFraction);
          end

          % Gaussian curvature by Subburaj
          % Curvature per vertex (allocated faces in previous calculations)
          totalAngleSum = cellfun(@sum, angleEdge_ROI);
          % Adapted for open mesh
          curveVertex = nan(size(totalAngleSum)); % Initialize with NaN values
          nonZeroIndices = obj.baryAreaFaces ~= 0; % Find indices where baryAreaFaces is not zero
          curveVertex(nonZeroIndices,1) = (360 - totalAngleSum(nonZeroIndices)) ./ obj.baryAreaFaces(nonZeroIndices);
          obj.gauss.vertex = curveVertex;
          % Curvature range (vertices)
          [ascendCurve,ascendCurveIdx] = sort(curveVertex);
          curveRange(2,1) = length(ascendCurve);
          curveRange(1,1) = curveRange(2,1) - round((topCurveRange/100)*curveRange(2,1));
          obj.gauss.topCurveVertexIdx = sort(ascendCurveIdx(curveRange(1):curveRange(2),1));

          % Convert to faces (for stlwrite)
          % Retrieve curvature for each vertex of each face
          vertexCurvatures = curveVertex(facesROI);
          % Calculate mean curvature for each face
          obj.gauss.vertexFace = mean(vertexCurvatures, 2);

          % Coloured stl (curvature): Write curvature in stl (stl binary) as attribute byte count
          stlName = 'curveGaussVertexFace';
          curveName = 'normVertexFace';
          obj = obj.colourSTL(TR,obj.gauss.vertexFace,pelvisID,stlName,'gauss',curveName);

          %%%% Only curvature per vertex %%%%

          disp(['Gauss curvature: pelvis ', num2str(pelvisNum)]);

      end

      %% Coloured stl
      function obj = colourSTL(obj, TR, attribute, pelvisID, stlName, curveType, curveName)

          % Coloured stl (curvature): Write curvature in stl (stl binary) as attribute byte count (attribute per face)
          % Data distribution: calculate the 2.5th and 97.5th percentiles (95% of the data)
          % Filter out NaN values for percentile calculation (for open mesh)
          validAttribute = attribute(~isnan(attribute));
          lowerBound = prctile(validAttribute, 2.5);
          upperBound = prctile(validAttribute, 97.5);
          % Normalize data within this range
          normalData = (attribute - lowerBound) / (upperBound - lowerBound);
          normalData(normalData < 0) = 0;   % Clamp values below the range to 0
          normalData(normalData > 1) = 1;   % Clamp values above the range to 1
          obj.(curveType).(curveName) = normalData; 
          % Set NaN values to white color (for open mesh)
          nanIndices = isnan(normalData);
          normalData(nanIndices) = 1; % Set to max value for mapping to white color

          % Apply colormap (recommended: viridis)
          % 16-bit colour: alpha: 1bit R: 5bit G: 5bit B: 5bit -> 32768 colors
          colourNum = 32768;
          colourMap = viridis(colourNum);
          colourIdx = round(normalData * (colourNum-1)) + 1;
          rgbColour = colourMap(colourIdx,:);
          % Set NaN values to white color in RGB (for open mesh)
          rgbColour(nanIndices, :) = repmat([1, 1, 1], sum(nanIndices), 1);
          % Curvature colour
          rgbName = ['RGB',curveName];
          obj.(curveType).(rgbName) = rgbColour;

          % Convert to 16-bit colour
          alpha = bitshift(uint16(1),15); % don't forget!
          red = bitshift(uint16(round(rgbColour(:,1) * 31)), 10);
          green = bitshift(uint16(round(rgbColour(:,2) * 31)), 5);
          blue = uint16(round(rgbColour(:,3) * 31));
          colour = bitor(bitor(bitor(alpha, red), green), blue);
          % Set NaN values to white color in 16-bit format (for open mesh)
          colour(nanIndices) = 65535; % White color in 16-bit format

          % Stl write
          stlwrite(TR, ['.\GeometriesCurves\pelvis', pelvisID, stlName, '.stl'], ...
              'binary','Attribute', colour)
          
          disp(['write stl: pelvis ', pelvisID]);

      end

      %% Calculation of the distance from the boundary of the defect area to the top curvature values (curve range)
      function obj = curveRangeBound(obj, pelvisNum, vertices, topCurveVerticesIdx, curveType, curveName)

          % Coordinates for top curvature range
          topCurveVertices = vertices(topCurveVerticesIdx, :); % whole pelvis
          % Boundary vertices
          boundaryVertices = obj.boundVertices; % defect area (section of the entire pelvis)

          % Find the nearest top curvature vertex for each boundary vertex
          [idx, distances] = knnsearch(topCurveVertices, boundaryVertices);
          obj.(curveType).([curveName,'BoundDistance']) = distances;

          % Calculate the average of the minimum distances
          obj.(curveType).([curveName,'BoundDistanceMean']) = mean(distances);
          obj.(curveType).([curveName,'BoundDistanceStd']) = std(distances);
          obj.(curveType).([curveName,'BoundDistanceMax']) = max(distances);
          obj.(curveType).([curveName,'BoundDistanceMin']) = min(distances);

          % Display the result
          disp(['Validation defect area: pelvis ', num2str(pelvisNum)]);

      end

    end
end