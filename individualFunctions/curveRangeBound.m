% Calculation of the distance from the boundary of the defect area to the top curvature values (curve range)

% Input:    vertices: vertices of the whole pelvis
%           topCurveVerticesIdx: indices of the vertices with the maximum curvature values 
%           boundVertices: Coordinates of boundary vertices 
%                           (parameter data from function boundSTL)

% Output:   distances: distances between boundary vertices and nearest point of top curvature values         
%           meanDistance: mean of distances (boundary of defect area)
%           maxDistance: max distance between a boundary vertice and nearest point of top curvature values 
%           minDistance: min distance between a boundary vertice and nearest point of top curvature values 

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function [distances, meanDistance, maxDistance, minDistance] = curveRangeBound(vertices, topCurveVerticesIdx, boundVertices)
    
    % Coordinates for top curvature range
    topCurveVertices = vertices(topCurveVerticesIdx, :); % whole pelvis
    % Boundary vertices
    boundaryVertices = boundVertices; % defect area (section of the entire pelvis)
    
    % Find the nearest top curvature vertex for each boundary vertex
    [idx, distances] = knnsearch(topCurveVertices, boundaryVertices);
    
    % Calculate the average of the minimum distances
    meanDistance = mean(distances);
    maxDistance = max(distances);
    minDistance = min(distances);
    
    % Display the result
    disp('Validation defect area');

end