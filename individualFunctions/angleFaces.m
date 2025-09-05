% Calculate angle between the faces (normal vector)
% Function for (1) normal and (2) mean curvature
% For open and closed meshes

% Input:    boundVerticesNr (for open meshes): row number of boundary vertices 
%                                               (parameter data from function boundSTL)
%           vertexAdjMap (cell format): Sorted adjacent faces around centre vertice (circle); without centre vertice
%                                        Column 1-2: row numbers of the vertices of the faces; without centre vertice
%                                        Column 3: row number of the faces = row number of the normals
%                                       without boundary vertices (for open meshes)
%                                       (parameter data from function adjVertexSeq)
%           verticesNrROI: numbering of the vertices in ROI (region of interest)
%           normalsROI: coordinates of the normals in ROI

% Output:   angleFaceNormal (cell format): angles between the adjacent faces around centre vertice 
%                                           sorted as in vertexAdjMap (circle)
%                                           without boundary vertices (for open meshes)

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function [angleFaceNormal] = angleFaces(boundVerticesNr,vertexAdjMap,verticesNrROI,normalsROI)

    verticesCount = size(verticesNrROI,1);
    
    % Preallocation cell array (cells with different size)
    angleFaceNormal = cell(verticesCount,1);
    vertexNormalMap = cell(verticesCount,1);
    for i = 1:verticesCount
    
        % Skip if vertex is a boundary vertex
        if ismember(verticesNrROI(i,1), boundVerticesNr) % for open mesh
            continue;
        end
    
        % Adjacent normals
        vertexNormal = normalsROI(vertexAdjMap{i,1}(:,3),:); % sorted
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
        angleFaceNormal{i,1} = angles;
    end
        
    disp('Face angles calculated');

end
