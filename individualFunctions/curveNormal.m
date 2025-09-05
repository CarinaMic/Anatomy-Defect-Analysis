% (1) Curvature: normal curvature (input: angleFaceNormal)
%               angles between the normals of adjacent faces 
% (1a) curveVertex: curvature related to the centre vertices (any number of triangles around centre vertices)
%                   -> curvature averaged to the faces (curveVertexFace)
% (1b) curveFace: curvature related to the centre faces (3 adjacent faces -> 3 angles)
% Curvature related to the faces in order to be able to write the curvature values to the stl (attribute byte count)
% For open and closed meshes

% Input:    boundFacesNr (for open meshes): faces structure of boundary faces
%                                           (parameter data from function boundSTL)
%           faceAdjMap: row number of faces which are adjacent to the certain face
%                       row number of faceAdjMap is the centre face
%                       3 adjacent faces to every face
%           angleFaceNormal (cell format): angles between the adjacent faces around centre vertice 
%                                           sorted as in vertexAdjMap (circle)
%                                           without boundary vertices (for open meshes)
%                                           (parameter data from function anglesFaces)
%           pelvisID: ID of the pelvis (for naming)
%           normalsROI: coordinates of the normals in ROI
%           facesROI: faces structure of the vertices in ROI
%           verticesROI: coordinates of the vertices in ROI 
%           topCurveRange: range with the maximum curvature values (in percent)

% Output:   curveVertex: curvature related to the centre vertices (any number of triangles around centre vertices)
%           curveVertexFace: curveVertex averaged to the faces 
%           topCurveVertexIdx: indices of the vertices with the maximum curvature values (curveVertex)
%           normVertexFace: normalized curveVertexFace
%                           non-linear normalization for non-linear colour distribution, 
%                           otherwise colour differentiation of the clustered value range is hardly possible
%                           colour distribution adjusted to 95% of the data (percentiles 2.5% and 97.5%)
%                           (parameter data from function colourSTL)
%           RGBnormVertexFace: normVertexFace as the corresponding RGB colour value
%                               (parameter data from function colourSTL)
% 
%           curveFace: curvature related to the centre faces (3 adjacent faces -> 3 angles)
%           topCurveFaceIdx: indices of the vertices with the maximum curvature values (curveFace)
%           normFace: normalized curveFace
%                       non-linear normalization for non-linear colour distribution, 
%                       otherwise colour differentiation of the clustered value range is hardly possible
%                       colour distribution adjusted to 95% of the data (percentiles 2.5% and 97.5%)
%                       (parameter data from function colourSTL)
%           RGBnormFace: normFace as the corresponding RGB colour value
%                       (parameter data from function colourSTL)

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function [curveVertex,curveVertexFace,topCurveVertexIdx,normVertexFace,RGBnormVertexFace, ...
    curveFace,topCurveFaceIdx,normFace,RGBnormFace] = curveNormal(boundFacesNr,faceAdjMap,angleFaceNormal, ...
    pelvisID,normalsROI,facesROI,verticesROI,topCurveRange)

    % stl
    TR = triangulation(facesROI,verticesROI);
    
    %%%% Curvature per vertex (allocated faces in previous calculations)); also for open mesh %%%%
    % Mean angle of every vertice (angle between the faces / normal vector)
    curveVertex = nan(size(angleFaceNormal,1),1); % Preallocation
    for i=1:size(angleFaceNormal,1)
        curveVertex(i,1) = mean(angleFaceNormal{i,1});
    end
    % Convert to faces (for stlwrite)
    % Retrieve curvature for each vertex of each face
    vertexCurvatures = curveVertex(facesROI);
    % Calculate mean curvature for each face
    curveVertexFace = mean(vertexCurvatures,2);
    % Curvature range (vertices)
    [ascendCurve,ascendCurveIdx] = sort(curveVertex);
    [~,curveRange(2,1)] = max(ascendCurve);
    curveRange(1,1) = curveRange(2,1) - round((topCurveRange/100)*size(curveVertex,1));
    topCurveVertexIdx = sort(ascendCurveIdx(curveRange(1):curveRange(2),1));
    
    % Coloured stl (curvature): Write curvature in stl (stl binary) as attribute byte count
    stlName = 'curveNormalVertexFace';
    [normVertexFace,RGBnormVertexFace] = colourSTL(TR,curveVertexFace,pelvisID,stlName); % function colourSTL


    %%%% Curvature per face (3 edges); also for open mesh %%%%
    % Mean angle of every face (angle between the faces / normal vector)
    angles = nan(size(faceAdjMap,1),3); % Preallocate matrix for angles
    facesNrROI = (1:size(facesROI,1))';
    for i = 1:size(faceAdjMap,1)
        % Skip if face is a boundary face
        if ismember(facesNrROI(i,1), boundFacesNr) % for open mesh
            continue;
        end
        centerNormal = normalsROI(i,:);
        adjNormals = normalsROI(faceAdjMap(i,:),:);
        % Vectorized angle calculation
        cosTheta = dot(repmat(centerNormal,3,1), adjNormals, 2) ./ ...
            (vecnorm(centerNormal,2,2) .* vecnorm(adjNormals,2,2));
        % Due to rounding, the value may lie outside the value range of cos [-1, 1]
        % hence the limitation of the value
        limitCosTheta = min(1, max(-1, cosTheta));
        angles(i,:) = acosd(limitCosTheta); % Convert to degrees
    end
    curveFace = mean(angles,2);
    % Curvature range (faces)
    [ascendCurve,ascendCurveIdx] = sort(curveFace);
    [~,curveRange(2,1)] = max(ascendCurve);
    curveRange(1,1) = curveRange(2,1) - round((topCurveRange/100)*size(curveFace,1));
    faceRangeIdx = sort(ascendCurveIdx(curveRange(1):curveRange(2),1));
    vertexIndices = facesROI(faceRangeIdx,:);
    flatVertexIndices = vertexIndices(:);
    topCurveFaceIdx = unique(flatVertexIndices);

    % Coloured stl (curvature): Write curvature in stl (stl binary) as attribute byte count
    stlName = 'curveNormalFace';
    [normFace,RGBnormFace] = colourSTL(TR,curveFace,pelvisID,stlNamee); % function colourSTL
    
    disp('Normal curvature calculated');

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
          % Convert to faces (for stlwrite)
          % Retrieve curvature for each vertex of each face
          vertexCurvatures = curveVertex(facesROI);
          % Calculate mean curvature for each face
          vertexFace = mean(vertexCurvatures,2);
          % Curvature range (vertices)
          [ascendCurve,ascendCurveIdx] = sort(curveVertex);
          curveRange(2,1) = length(ascendCurve);
          curveRange(1,1) = curveRange(2,1) - round((topCurveRange/100)*curveRange(2,1));
          topCurveVertexIdx = sort(ascendCurveIdx(curveRange(1):curveRange(2),1));
          
          % Coloured stl (curvature): Write curvature in stl (stl binary) as attribute byte count
          stlName = 'curveNormalVertexFace';
          curveName = 'normVertexFace';
          obj = obj.colourSTL(TR,vertexFace,pelvisID,stlName,'normal',curveName);

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