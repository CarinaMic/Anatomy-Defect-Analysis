% (3) Curvature: gauss curvature 
%                   Gauss curvature (per vertex) from Subburaj el al. (2009) and adjusted
%                   Subburaj, K., Ravi, B., & Agarwal, M. (2009). 
%                   Automated identification of anatomical landmarks on 3D bone models reconstructed from CT scan images. 
%                   Computerized Medical Imaging and Graphics, 33(5), 359-368.
% (3a) curveVertex: curvature related to the centre vertices (any number of triangles around centre vertices)
%                   -> curvature averaged to the faces (curveVertexFace)
% Curvature related to the faces in order to be able to write the curvature values to the stl (attribute byte count)
% For open and closed meshes

% Input:    comEdge (cell format): coordinates of the common edges of the adjacent faces (coordinates)
%                                   sorted as in vertexAdjMap (circle)
%                                   without boundary vertices (for open meshes)
%                                   (parameter data from function edgesArea)
%           baryAreaFaces: barycentric area which is one third of the area of the triangles around the centre vertice
%                           without boundary vertices (for open meshes)
%                           (parameter data from function edgesArea)
%           pelvisID: ID of the pelvis (for naming)
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

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function [curveVertex,curveVertexFace,topCurveVertexIdx,normVertexFace,RGBnormVertexFace] = curveGauss(...
    comEdge,baryAreaFaces, pelvisID,facesROI,verticesROI,topCurveRange)

    % stl
    TR = triangulation(facesROI,verticesROI);
    
    %%%% Curvature per vertex (allocated faces in previous calculations) %%%%
    % Angle between edges
    % angle = acosd(dot(u,v)/(norm(u)*norm(v))); % cos(phi)=dot(u,v)/(norm(u)/*norm(v))
    angleEdge_ROI = cell(size(comEdge, 1), 1);
    for i = 1:size(comEdge, 1)
        u = comEdge{i}(1:end-1, :);
        v = comEdge{i}(2:end, :);
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
    nonZeroIndices = baryAreaFaces ~= 0; % Find indices where baryAreaFaces is not zero
    curveVertex(nonZeroIndices,1) = (360 - totalAngleSum(nonZeroIndices)) ./ baryAreaFaces(nonZeroIndices);
    % Curvature range (vertices)
    [ascendCurve,ascendCurveIdx] = sort(curveVertex);
    [~,curveRange(2,1)] = max(ascendCurve);
    curveRange(1,1) = curveRange(2,1) - round((topCurveRange/100)*size(curveVertex,1));
    topCurveVertexIdx = sort(ascendCurveIdx(curveRange(1):curveRange(2),1));
    
    % Convert to faces (for stlwrite)
    % Retrieve curvature for each vertex of each face
    vertexCurvatures = curveVertex(facesROI);
    % Calculate mean curvature for each face
    curveVertexFace = mean(vertexCurvatures, 2);
    
    % Coloured stl (curvature): Write curvature in stl (stl binary) as attribute byte count
    stlName = 'curveGaussVertexFace';
    [normVertexFace,RGBnormVertexFace] = colourSTL(TR,curveVertexFace,pelvisID,stlName);
    
    %%%% Only curvature per vertex %%%%
    
    disp('Gauss curvature calculated');

end