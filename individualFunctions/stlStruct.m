% Restructuring of vertices, faces and normals (for File Exchange stlread function)

% Input:    vertices 
%           faces strcuture of stl

% Output:   vertices 
%           faces structure of the stl
% Mirroring, if necessary

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [verticesStruct,facesStruct] = stlStruct(vertices,faces,mirror)
    % For File Exchange stlread function: vertices can have duplicates
    % Structure of vertices and faces after stlread-function:
    % Vertices: Array with three columns(x,y,z-coordinates);
    %           Vertices are duplicated for the face structure;
    % Faces:    Array with three columns;
    %           Sequentially numbered in ascending order [1 2 3; 4 5 6; ...]
    %           -> Duplicated vertices (f.e. vertices 3 & 4 are identic);
    %           Each row presents one face;
    %           Each line containts the vertices numbers (row number) for one face;
    
    % New structure:
    % Vertices: Array with three columns(x,y,z-coordinates);
    %           Vertices are not duplicated for the face structure;
    %           -> Reduces computing time later;
    % Faces:    Array with three columns;
    %           Each row presents one face;
    %           Each line containts the vertices numbers (row number) for one face
    %           -> New order to link the correct vertices;
    
    % Find duplicates in vertices, remove them and ...
    % modify the array for the faces so that the correct vertices are linked
    
    % Find duplicates in vertices (unique)
    [vertices_woDup,~,vertices_numDup] = unique(vertices,'rows','stable');
    
    faces_in = faces'; % Transpose for following calculation
    faces_in(:) = vertices_numDup(:); % Replace the duplicates with the first occurence
    facesStruct = faces_in'; % Output: new faces structure
    
    verticesStruct = vertices_woDup; % Output: new vertices structure
    
    % Mirroring 
    if mirror==1
        % Mirroring at yz-plane -> x-Coordinate * (-1)
        verticesStruct(:,1) = (-1) * verticesStruct(:,1);
        % Mesh: normal orientation to the center -> defines the outside of the mesh
        % The orientation of the normal is defined by the ordering of the vertices in a triangle (faces strcuture).
        % The ordering follows the right-hand rule.
        facesStruct = faces_in([1;3;2],:)';
    end
    
    % Normals are the same as before
    
    disp('STL data restructured')

end