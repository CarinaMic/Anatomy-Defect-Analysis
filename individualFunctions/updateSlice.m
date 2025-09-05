%% Function to update the slice with the reference pelvis

% Input:    slider: UI control with a numeric value (z-position of the slice)
%           vertices: point coordinates to display (same order as rgbColour)
%           rgbColour: per-point RGB colors in [0,1], aligned with `vertices`
%           pelvis: struct array with reference pelvis data
%           highlightIndices: vector of global point indices to highlight (relative to
%                   `verticesMaskIdx`); pass [] for no highlighting
%           verticesMaskIdx  mapping from rows of `vertices` to global grid indices

% Output:   (none)  The current axes are cleared and redrawn. The function
%                   updates the visualization only and returns no variables.

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [] = updateSlice(slider, vertices, rgbColour, pelvis, highlightIndices, verticesMaskIdx)

z = slider.Value;
tolerance = 0.5;
% Clear figure, but keep pelvis and light visible
cla;
hold on;  % Keep elements visible during updates

% Replot the reference pelvis
patch('Faces', pelvis(1).import.processed.faces, ...
    'Vertices', pelvis(1).import.processed.vertices, ...
    'FaceColor', [0.8, 0.8, 0.8], ... % TUMcolors.grey20
    'FaceAlpha', 0.5, ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.5);

% Replot the reference acetabulum center
plot3(pelvis(1).import.processed.acentre(1), ...
    pelvis(1).import.processed.acentre(2), ...
    pelvis(1).import.processed.acentre(3), ...
    '.', 'Color', [0 1 0], 'MarkerSize', 50);  % pelvisDefect(1).import.processed.acentre

% Add the light for 3D effect
light('Position', [1 1 5], 'Style', 'infinite');

% Plot the current slice based on Z-value
idx_slice = abs(vertices(:,3) - z) < tolerance;  % Select points near the Z-value (logical)
% Highlighting
if isempty(highlightIndices)
    % no highlighting
    % idx Base: vertices (allPelvis.refPoints.gridPoints.inside)
    scatter3(vertices(idx_slice, 1), vertices(idx_slice, 2), vertices(idx_slice, 3), ...
        30, rgbColour(idx_slice, :), 'filled');
else
    % Highlighting
    % idx Base: grid (allPelvis.refPoints.gridPoints.pointsBox)
    pointIndices = verticesMaskIdx(idx_slice,:);
    isHighlighted = ismember(pointIndices, highlightIndices); % result: logical to pointIndices;
    idx_highlighted = pointIndices(isHighlighted);
    idx_normal = pointIndices(~isHighlighted);
    % transform to idx Base: vertices
    [~, idx_highlighted_in_vertices] = ismember(idx_highlighted, verticesMaskIdx);
    [~, idx_normal_in_vertices] = ismember(idx_normal, verticesMaskIdx);
    % not in pct
    scatter3(vertices(idx_normal_in_vertices,1), vertices(idx_normal_in_vertices,2), vertices(idx_normal_in_vertices,3), ...
        30, rgbColour(idx_normal_in_vertices,:), 'filled');
    % Highlighted points (in pct)
    scatter3(vertices(idx_highlighted_in_vertices,1), vertices(idx_highlighted_in_vertices,2), vertices(idx_highlighted_in_vertices,3), ...
        30, rgbColour(idx_highlighted_in_vertices,:), 'filled', ...
        'MarkerEdgeColor', 'magenta', 'MarkerEdgeAlpha', 0.25, 'LineWidth', 0.1);
end

% Plot all points below the current slice
idx_below = vertices(:,3) < z;  % Select points below the Z-value
scatter3(vertices(idx_below, 1), vertices(idx_below, 2), vertices(idx_below, 3), ...
    2, rgbColour(idx_below, :), 'filled');

% Set display properties
daspect([1 1 1]);
view(3);
hold off;
end