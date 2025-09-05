%% Intersection sequence (function inpolyhedron): with defect volume and defect surface points; new defect volume with shrink wrap

% Input:    n: number of defects; nchoosek(n, overlapCount) combinations are formed
%           overlapCount: number of defects to overlap per combination (k aus n)
%           importVertices: list of the vertices of the import mesh
%           importFaces: list of the faces of the import mesh

% Output:   overlapCombis: all chosen combinations (nchoosek)
%           overlapVertices: intersection vertices per combination
%           overlapFaces: intersection faces per combination

% Developed by C.Micheler, 
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [overlapCombis, overlapVertices, overlapFaces] = intersectSeq(n, overlapCount, importVertices, importFaces)

% Calculate intersection combinations: binomial coefficient
defectCombis = nchoosek(n, overlapCount);
overlapCombis = defectCombis;
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

% Log the intersection completion
disp('Defect intersection completed.');
end