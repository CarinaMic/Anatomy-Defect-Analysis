 classdef Distribution
    % Intersect class for distribution of the defect models
    
    % Developed by C.Micheler,
    % Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
    % Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich
    

    properties 
        inside = struct();          % Intersection type inside
    end
    
    methods
        %% Constructer: generate object
        function obj = Distribution()
            disp('class Distribution initialized')
        end

        %% Function to update the slice with the reference pelvis
        function obj = updateSlice(obj, slider, vertices, rgbColour, pelvis, highlightIndices, verticesMaskIdx)
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

        %% Display pelvis defect distribution in slices: function to update the view
        function obj = updateViews(obj, src, ax1, ax2, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, viewType, ...
                vertices, rgbColour, useHighlight, highlightIndices, verticesMaskIdx)
            handles = guidata(src);

            % Active views
            sagittalActive = isfield(handles, 'sagittalSlider') && isvalid(handles.sagittalSlider);
            coronalActive = isfield(handles, 'coronalSlider') && isvalid(handles.coronalSlider);
            axialActive = isfield(handles, 'axialSlider') && isvalid(handles.axialSlider);
            if strcmp(viewType, 'sagittal')
                slider = handles.sagittalSlider;
                axisSelection = [2, 3]; % Y-Z-Plane
                xlabelText = 'Y'; ylabelText = 'Z';
                planeAxis = 1; % X-Position
                titleText = 'Sagittal View';
            elseif strcmp(viewType, 'coronal')
                slider = handles.coronalSlider;
                axisSelection = [1, 3]; % X-Z-Plane
                xlabelText = 'X'; ylabelText = 'Z';
                planeAxis = 2; % Y-Position
                titleText = 'Coronal View';
            elseif strcmp(viewType, 'axial')
                slider = handles.axialSlider;
                axisSelection = [1, 2]; % X-Y-Plane
                xlabelText = 'X'; ylabelText = 'Y';
                planeAxis = 3; % Z-Position
                titleText = 'Axial View';
            end
            %D = handles.D;
            %rgbColour = handles.rgbColour;

            % Actual value of slider
            currentValue = get(slider, 'Value');
            tolerance = 0.5;  % Tolerance for shift selection
            % Selection of points near the current X value
            idx = abs(vertices(:, planeAxis) - currentValue) < tolerance;

            % Update 2D view
            axes(ax1); cla(ax1);
            hold(ax1, 'on');
            % Highlighting
            if useHighlight && ~isempty(highlightIndices)
                % Highlighting
                % idx Base: grid (allPelvis.refPoints.gridPoints.pointsBox)
                pointIndices = verticesMaskIdx(idx);
                isHighlighted = ismember(pointIndices, highlightIndices);
                idx_highlighted = pointIndices(isHighlighted);
                idx_normal = pointIndices(~isHighlighted);
                % transform to idx Base: vertices
                [~, idx_highlighted_in_vertices] = ismember(idx_highlighted, verticesMaskIdx);
                [~, idx_normal_in_vertices] = ismember(idx_normal, verticesMaskIdx);

                % not in pct
                scatter(ax1,vertices(idx_normal_in_vertices, axisSelection(1)), vertices(idx_normal_in_vertices, axisSelection(2)),  ...
                    30, rgbColour(idx_normal_in_vertices,:), 'filled');
                % Highlighted points (in pct)
                scatter(ax1, vertices(idx_highlighted_in_vertices, axisSelection(1)), vertices(idx_highlighted_in_vertices, axisSelection(2)), ...
                    30, rgbColour(idx_highlighted_in_vertices,:), 'filled', ...
                    'MarkerEdgeColor', 'magenta', 'MarkerEdgeAlpha', 0.25, 'LineWidth', 0.1);
            else
                % no highlighting
                % idx Base: vertices (allPelvis.refPoints.gridPoints.inside)
                scatter(ax1, vertices(idx, axisSelection(1)), vertices(idx, axisSelection(2)), ...
                    30, rgbColour(idx, :), 'filled');
            end

            % Acetabulum Centre in 2D view
            hold(ax1, 'on');
            % Remove previous acetabulum centre
            oldCentre = findobj(ax1, 'Tag', 'AcetabulumCentre');
            if ~isempty(oldCentre)
                delete(oldCentre);
            end
            % Acetabulum-Center
            centre = pelvis(1).import.processed.acentre;
            if abs(centre(planeAxis) - currentValue) < (tolerance + 0.5)
                plot(ax1, centre(axisSelection(1)), centre(axisSelection(2)), '.', ...
                    'Color', [0 0.3961 0.7412] , 'MarkerSize', 75); 
                %plot(ax1, centre(axisSelection(1)), centre(axisSelection(2)), 'x', ...
                %    'Color', 'm', 'MarkerSize', 32, 'LineWidth',11);
            end

            title([titleText ' - Value = ' num2str(currentValue, '%.2f')]);
            xlabel(ax1, xlabelText);
            ylabel(ax1, ylabelText);
            daspect(ax1, [1 1 1]);

            switch viewType
                case 'sagittal'
                    axis(ax1, [y_min y_max z_min z_max]);
                case 'coronal'
                    axis(ax1, [x_min x_max z_min z_max]);
                case 'axial'
                    axis(ax1, [x_min x_max y_min y_max]);
            end
            hold(ax1, 'off');


            % Update 3D view
            axes(ax2);
            hSurfaces = findobj(ax2, 'Type', 'Surface');
            delete(hSurfaces);
            hPatch = findobj(ax2, 'Type', 'patch');
            if isempty(hPatch)
                patch('Faces', pelvis(1).import.processed.faces, ...
                    'Vertices', pelvis(1).import.processed.vertices, ...
                    'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.5, ...
                    'EdgeColor', 'none', 'Parent', ax2);
            else
                set(hPatch, 'Vertices', pelvis(1).import.processed.vertices);
            end
            view(ax2, 37.5, 30);
            daspect(ax2, [1 1 1]);
            title('3D View');
            xlabel(ax2, 'X');
            ylabel(ax2, 'Y');
            zlabel(ax2, 'Z');

            % Light
            if ~isfield(handles, 'light') || ~isvalid(handles.light)
                handles.light = camlight(ax2);
                guidata(src, handles);
            end
            lighting(ax2, 'phong');
            hold(ax2, 'on');

            % 3D-Plane (Slice)
            [X, Y] = meshgrid(linspace(x_min, x_max, 10), linspace(y_min, y_max, 10)); % Axial
            [X2, Z] = meshgrid(linspace(x_min, x_max, 10), linspace(z_min, z_max, 10)); % Sagittal
            [Y2, Z2] = meshgrid(linspace(y_min, y_max, 10), linspace(z_min, z_max, 10)); % Coronal

            % Plane in 3D view
            if sagittalActive && coronalActive && axialActive
                % Show all planes
                surf(X, Y, ones(size(X)) * get(handles.axialSlider, 'Value'), ...
                    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'red', 'Parent', ax2); % Axial (Z)
                surf(ones(size(Y2)) * get(handles.sagittalSlider, 'Value'), Y2, Z2, ...
                    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'blue', 'Parent', ax2); % Sagittal (X)
                surf(X2, ones(size(X2)) * get(handles.coronalSlider, 'Value'), Z, ...
                    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'green', 'Parent', ax2); % Coronal (Y)
            else
                % Show only one plane
                switch viewType
                    case 'axial'
                        surf(X, Y, ones(size(X)) * currentValue, ...
                            'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'red', 'Parent', ax2);
                    case 'sagittal'
                        surf(ones(size(Y2)) * currentValue, Y2, Z2, ...
                            'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'blue', 'Parent', ax2);
                    case 'coronal'
                        surf(X2, ones(size(X2)) * currentValue, Z, ...
                            'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'green', 'Parent', ax2);
                end
            end

            % Acetabulum Centre in 3D view
            plot3(ax2, centre(1), centre(2), centre(3), '.', ...
                'Color', [0 0.3961 0.7412], 'MarkerSize', 40); 
            %plot3(ax2, centre(1), centre(2), centre(3), 'x', ...
            %    'Color', 'm', 'MarkerSize', 15, 'LineWidth',5);
            hold(ax2, 'off');

            % Reduce unnecessary calls to drawnow
            drawnow expose;
        end
    end
 end