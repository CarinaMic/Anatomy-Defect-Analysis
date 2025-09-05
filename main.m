%% Matlab initialisation script

% Project for the analysis of pelvises and acetabular defects

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


% Clear environment and close figures
clc;
clear variables;
close all;

% Add necessary paths
addpath('./Functions/'); % Add Function folder

% Parallel computing setup
if isempty(gcp('nocreate'))
    parpool('local'); % Start parallel pool if not already running
end

% Define TUM color palette
TUMcolors = struct( ...
    'blue300', [0, 101/255, 189/255], ...
    'blue540', [0, 51/255, 89/255], ...
    'blue301', [0, 82/255, 147/255], ...
    'blue542', [100/255, 160/255, 200/255], ...
    'blue283', [152/255, 198/255, 234/255], ...
    'grey80',  [51/255, 51/255, 51/255], ...
    'grey50',  [128/255, 128/255, 128/255], ...
    'grey20',  [204/255, 204/255, 204/255], ...
    'cream',   [218/255, 215/255, 203/255], ...
    'orange',  [227/255, 114/255, 34/255], ...
    'green',   [162/255, 173/255, 0/255] ...
);

%% %%%%%%%%%%%% Class Pelvis %%%%%%%%%%%%%%%

% Folder of the geometries
filedir = './GeometriesPelves/';
% Lists of stl and xlsx files
stlfiles = dir(fullfile(filedir, '*.stl'));
xlsfiles = dir(fullfile(filedir, '*.xlsx'));
% Use stl/xls files count to determine data count
dataCount = length(stlfiles);
%dataCount = length(xlsfiles);

% Check for mismatch or absence of files
if isempty(stlfiles)
    error('No stl files found in the specified directory.');
elseif isempty(xlsfiles)
    error('No xlsx files found in the specified directory.');
elseif length(stlfiles) ~= length(xlsfiles)
    error('Number of stl files does not match number of xlsx files.');
end

% Initalization: Generate object array
pelvis(dataCount,1) = Pelvis;

%% %%%%%%%%%%%% Class Import %%%%%%%%%%%%%%%
% Import pelvis data

% File paths
filePath_pelvisData = fullfile(filedir, {xlsfiles.name}); % Pelvis data/infos
filePath_pelvis = fullfile(filedir, {stlfiles.name}); % Pelvis geometry data

% Parallel loop for data import and processing
parfor i = 1:dataCount

   % Initalization: Generate object
   pelvisImport = Import(i); % cache (to ensure that the parallel calculations are independent) 
      
   % Import the pelvis data/infos
   pelvisImport = pelvisImport.importData(i,filePath_pelvisData{i});
  
   % Import the pelvis geometry data
   pelvisImport = pelvisImport.importSTL(i, filePath_pelvis{i}) ...
                               .centroid(i) ... % Centroid of the faces (for the normal vector)
                               .edgeLength(i); % Edge length of  mesh

   % stl properties
   verticesPelvis(i,1) = size(pelvisImport.loaded.TR.Points,1);
   facesPelvis(i,1) = size(pelvisImport.loaded.TR.ConnectivityList,1);

   % Mean edge length
   edgeLength(i,1) = pelvisImport.processed.meanEdgeLength;

   % Paprosky
   paproskyType(i,1) = pelvisImport.loaded.paprosky;

   % Store the result back in pelvis
   pelvis(i).import = pelvisImport;

end

% Paprosky classes
allPelvis.paprosky.type = paproskyType;
paproskyTypes = {'2a', '2b', '2c', '3a', '3b', 'NaN'};
paproskyNames = {'IIa', 'IIb', 'IIc', 'IIIa', 'IIIb', 'NaN'};
allPelvis.paprosky.counts = struct();
for i = 1:length(paproskyNames)
    allPelvis.paprosky.counts.(paproskyNames{i}) = 0; 
end
for i = 1:length(paproskyType)
    idx = find(strcmp(paproskyType{i}, paproskyTypes));
    if ~isempty(idx)
        allPelvis.paprosky.counts.(paproskyNames{idx}) = allPelvis.paprosky.counts.(paproskyNames{idx}) + 1;
    end
end
disp(allPelvis.paprosky.counts);
% with reference pelvis
% STL properties
allPelvis.stl.vertices = verticesPelvis;
allPelvis.stl.verticesMean = mean(verticesPelvis);
allPelvis.stl.verticesStd = std(verticesPelvis);
allPelvis.stl.verticesMax = max(verticesPelvis);
allPelvis.stl.verticesMin = min(verticesPelvis);
allPelvis.stl.faces = facesPelvis;
allPelvis.stl.facesMean = mean(facesPelvis);
allPelvis.stl.facesStd = std(facesPelvis);
allPelvis.stl.facesMax = max(facesPelvis);
allPelvis.stl.facesMin = min(facesPelvis);
% Mean edge Length
allPelvis.stl.edgeLength = edgeLength;
allPelvis.stl.edgeMean = mean(edgeLength);
allPelvis.stl.edgeStd = std(edgeLength);
allPelvis.stl.edgeMax = max(edgeLength); % max of mean
allPelvis.stl.edgeMin = min(edgeLength); % min of mean
% without reference pelvis
% % STL properties
% allPelvis.stl.verticesMean = mean(verticesPelvis(2:end));
% allPelvis.stl.verticesStd = std(verticesPelvis(2:end));
% allPelvis.stl.verticesMax = max(verticesPelvis(2:end));
% allPelvis.stl.verticesMin = min(verticesPelvis(2:end));
% allPelvis.stl.facesMean = mean(facesPelvis(2:end));
% allPelvis.stl.facesStd = std(facesPelvis(2:end));
% allPelvis.stl.facesMax = max(facesPelvis(2:end));
% allPelvis.stl.facesMin = min(facesPelvis(2:end));
% % Mean edge Length
% allPelvis.stl.edgeMean = mean(edgeLength(2:end));
% allPelvis.stl.edgeStd = std(edgeLength(2:end));
% allPelvis.stl.edgeMax = max(edgeLength(2:end)); % max of mean
% allPelvis.stl.edgeMin = min(edgeLength(2:end)); % min of mean

clear pelvisImport verticesPelvis facesPelvis edgeLength paproskyType

%% Save and load properties of Class Import

% Save class import data (properties)
savePelvisDataImport = struct();
% Meta information of curve
metaImport = metaclass(pelvis(1).import);
propertiesImport = {metaImport.PropertyList.Name};
for i = 1:dataCount
    for j = 1:length(propertiesImport)
        propertyName = propertiesImport{j};
        % Cache
        savePelvisDataImport(i).(propertyName) = pelvis(i).import.(propertyName);
    end
end
% Save data
save('.\pelvisDataImport.mat', 'savePelvisDataImport', '-v7.3'); % Adapt storage location
clear savePelvisDataImport metaImport propertiesImport propertyName

% Load class import data (properties)
loadPelvisDataImport = load('.\pelvisDataImport.mat', 'savePelvisDataImport'); % Adapt storage location
%pelvis(dataCount) = Pelvis; % Initialisation
metaImport = metaclass(pelvis(1).import);
%propertiesImport = {metaImport.PropertyList.Name}; % load all properties
propertiesImport = {'processed'}; % load selected properties: {'loaded', 'processed'}
for i = 1:dataCount
    for j = 1:length(propertiesImport)
        propertyName = propertiesImport{j};
        if isfield(loadPelvisDataImport.savePelvisDataImport(i), propertyName)
            pelvis(i).import.(propertyName) = loadPelvisDataImport.savePelvisDataImport(i).(propertyName);
        end
    end
end
clear loadPelvisDataImport metaImport propertiesImport propertyName

%% Display the pelvis with the face normals (for control)

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount 
    i = 1;  % Pelvis number
    %figure('Visible','off')
    figure
    hold on
    % Pelvis
    patch('Faces',pelvis(i).import.processed.faces,...
        'Vertices',pelvis(i).import.processed.vertices,...
        ... % 'FaceVertexCData' Number of colors must be equal to number of verticex; 
        ... % 'FaceVertexCData',(1:size(pelvis(i).import.processed.vertices,1))',...
        ... % Color interpolated: 'FaceColor','interp'
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',TUMcolors.grey50,...    % Edge color
        'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis 3Dview (lightning)
    % patch('Faces',pelvis(i).import.processed.faces,...
    %     'Vertices',pelvis(i).import.processed.vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor','none',...              % Edge color
    %     'EdgeAlpha',0.25,...                % Transparency of the edges
    %     ... % Ligthing for 3d effect
    %     'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
    %     'AmbientStrength', 0.5);
    % light('Position', [1 1 5], 'Style', 'infinite');
    
    % Normals of the faces (optional)
    % quiver3(pelvis(i).import.processed.centreFaces(:,1),...
    %     pelvis(i).import.processed.centreFaces(:,2),...
    %     pelvis(i).import.processed.centreFaces(:,3),...
    %     pelvis(i).import.processed.normals(:,1),...
    %     pelvis(i).import.processed.normals(:,2),...
    %     pelvis(i).import.processed.normals(:,3),2,...
    %     'Color',TUMcolors.orange)

    % Display acetabulum centre
    plot3(pelvis(i).import.processed.acentre(1), ...
        pelvis(i).import.processed.acentre(2), ...
        pelvis(i).import.processed.acentre(3), ...
        '*', 'MarkerEdgeColor', TUMcolors.blue300, ...
        'MarkerSize', 10);

    % Set view and axis properties
    title(['Pelvis ' num2str(i)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z');
    % grid minor
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')3D.fig'])
%end

%% Display of the distribution of the edge length in a histogram (for control)

% Distribution
% Define histogram edges
histoEdges = 0.025:0.05:1.025;

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount 
    i = 1;  % Pelvis number
    %figure('Visible','off')
    figure
    histogram(pelvis(i).import.processed.edgeLength, histoEdges, 'Normalization', 'probability'); % Histogram
    % Setting up labels
    title(['Edge Length: Pelvis ' num2str(i)]); 
    xAxisLabels = arrayfun(@(a, b) sprintf('%.3f-%.3f', a, b), ...
        histoEdges(1:end-1), histoEdges(2:end), 'UniformOutput', false);
    xticks(histoEdges(1:end-1) + 0.025); % Set x-ticks at the center of each bin
    xticklabels(xAxisLabels);
    xlabel('Edge Length');
    xlim([0.025, 1.025]);
    ylim([0 1]);
    ylabel('Frequency');

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')edgeLengthHisto.fig'])
%end

clear histoEdges xAxisLabels

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Class Curvature %%%%%%%%%%%%

% Curvature of the whole pelvis 
% (1) Normal Curvature
% (2) Mean Curvature
% (3) Gauss Curvature

% Determine the curvature
% with orginial data (not mirrored or transformed)
% Preallocate array for storing times
executionTimesCurve = zeros(dataCount, 1); 
parfor i = 1:dataCount 
    tic
    % Determine adjacent vertices / faces (vertex face map)
    verticesNr = (1:size(pelvis(i).import.loaded.TR.Points,1))';
    pelvis(i).curve = pelvis(i).curve.adjVertex(i,verticesNr,pelvis(i).import.loaded.TR.ConnectivityList);
    % Sort adjacent vertices / faces (vertex face map)
    pelvis(i).curve = pelvis(i).curve.adjVertexSeq(i,verticesNr);
    % Calculate the edges and face areas
    % For mean (2) and gauss (3) curvature
    pelvis(i).curve = pelvis(i).curve.edgesArea(i,verticesNr,pelvis(i).import.loaded.TR.Points);
    % Calculate the angles between the faces (normal vector)
    % For normal (1) and mean (2) curvature
    pelvis(i).curve = pelvis(i).curve.angleFaces(i,verticesNr,pelvis(i).import.loaded.normals);
    % Determine adjacent faces
    pelvis(i).curve = pelvis(i).curve.adjFaces(i,pelvis(i).import.loaded.TR.ConnectivityList,...
        pelvis(i).import.loaded.TR.Points);  

    % Extract pelvis name
    underscoreIndex = strfind(stlfiles(i).name, '_');
    pelvisID = stlfiles(i).name(1:underscoreIndex(1)-1);
 
    % Determine the curvature (with stl write)
    topCurveRange = 2.5; % in percent
    % (1) Normal Curvature
    pelvis(i).curve = pelvis(i).curve.curveNormal(i,pelvisID,pelvis(i).import.loaded.normals,...
        pelvis(i).import.loaded.TR.ConnectivityList,pelvis(i).import.loaded.TR.Points,topCurveRange);
    % (2) Mean Curvature
    pelvis(i).curve = pelvis(i).curve.curveMean(i,pelvisID,pelvis(i).import.loaded.normals,...
        pelvis(i).import.loaded.TR.ConnectivityList,pelvis(i).import.loaded.TR.Points,topCurveRange);
    % (3) Gauss Curvature
    pelvis(i).curve = pelvis(i).curve.curveGauss(i,pelvisID,verticesNr,...
        pelvis(i).import.loaded.TR.ConnectivityList,pelvis(i).import.loaded.TR.Points,topCurveRange);
    executionTimesCurve(i) = toc; % End timer and store time in array
end

% Execution time of curvature
allPelvis.time.curve = executionTimesCurve;
allPelvis.time.curveMean = mean(executionTimesCurve);
allPelvis.time.curveStd = std(executionTimesCurve);
allPelvis.time.curveMax = max(executionTimesCurve);
allPelvis.time.curveMin = min(executionTimesCurve);

% Data storage of curvature data
% Preallocation
allPelvis.curve.size.curveInfo = cell(length(pelvis), 1);
allPelvis.curve.size.curve = zeros(length(pelvis),1); 
% Size of pelvis(i).curve variable (complete pelvis, without boundary) 
for i = 1:dataCount
    sizeCurve = pelvis(i).curve; % cache
    allPelvis.curve.size.curveInfo{i} = whos('sizeCurve');
    allPelvis.curve.size.curve(i) = allPelvis.curve.size.curveInfo{i}.bytes;
end
allPelvis.curve.size.curveMean = mean(allPelvis.curve.size.curve);
allPelvis.curve.size.curveStd = std(allPelvis.curve.size.curve);
allPelvis.curve.size.curveMax = max(allPelvis.curve.size.curve);
allPelvis.curve.size.curveMin = min(allPelvis.curve.size.curve);

clear executionTimesCurve verticesNr underscoreIndex pelvisID topCurveRange

%% Comparison of the curvature types (normalized data)

pelvisCurve = cell(dataCount, 1);
for i = 1:dataCount 
    % Normalized curvature data
    % Retrieve all curvature data for the current dataset
    curveNormalVertexFace = pelvis(i).curve.normal.normVertexFace;
    curveNormalFace = pelvis(i).curve.normal.normFace;
    curveMeanVertexFace = pelvis(i).curve.cMean.normVertexFace;
    curveMeanFace = pelvis(i).curve.cMean.normFace;
    curveGaussVertexFace = pelvis(i).curve.gauss.normVertexFace;

     % Pair the datasets for comparison
    curvePairs = {
        {'cNormal', curveNormalVertexFace, curveNormalFace},
        {'cMean', curveMeanVertexFace, curveMeanFace},
        {'verNormalGauss', curveNormalVertexFace, curveGaussVertexFace},
        {'verMeanGauss', curveMeanVertexFace, curveGaussVertexFace},
        {'normalMean', curveNormalFace, curveMeanFace},
        {'verNormalMean', curveNormalVertexFace, curveMeanVertexFace}
 
        {'normalVerMean', curveNormalVertexFace, curveMeanFace},
        {'normalMeanVer', curveNormalFace, curveMeanVertexFace}, 
        {'normalGaussVer', curveNormalFace, curveGaussVertexFace},
        {'meanGaussVer', curveMeanFace, curveGaussVertexFace}
    };
   
   % Loop through each pair and perform the comparison
    for pairIdx = 1:length(curvePairs)
        pairName = curvePairs{pairIdx}{1};
        curveData1 = curvePairs{pairIdx}{2};
        curveData2 = curvePairs{pairIdx}{3};

        % Calculate difference (mean absolute error)
        diff = curveData1 - curveData2;
        pelvisCurve{i}.(pairName).meanAbsError = mean(abs(diff));
        pelvisCurve{i}.(pairName).meanSquaredError = mean(diff.^2);
        
        % Calculate the correlation between the two curves
        pelvisCurve{i}.(pairName).correlation = corr(curveData1, curveData2);

        % Store values for similarity in allPelvis structure
        allPelvis.curve.meanAbsError.(pairName)(i,1) = pelvisCurve{i}.(pairName).meanAbsError;
        allPelvis.curve.meanSquaredError.(pairName)(i,1) = pelvisCurve{i}.(pairName).meanSquaredError;
        allPelvis.curve.corr.(pairName)(i,1) = pelvisCurve{i}.(pairName).correlation;

        % Define common bin edges for histograms
        pelvisCurve{i}.(pairName).binEdges = linspace(0, 1, 20);

        % Calculate histograms 
        pelvisCurve{i}.(pairName).counts1 = histcounts(curveData1, pelvisCurve{i}.(pairName).binEdges);
        pelvisCurve{i}.(pairName).counts2 = histcounts(curveData2, pelvisCurve{i}.(pairName).binEdges);

        % Normalize counts to create a probability distribution
        pelvisCurve{i}.(pairName).pdf1 = pelvisCurve{i}.(pairName).counts1 / sum(pelvisCurve{i}.(pairName).counts1);
        pelvisCurve{i}.(pairName).pdf2 = pelvisCurve{i}.(pairName).counts2 / sum(pelvisCurve{i}.(pairName).counts2);

        % Calculate Intersection of histograms (sum of minimums)
        pelvisCurve{i}.(pairName).intersection = sum(min(pelvisCurve{i}.(pairName).pdf1, ...
            pelvisCurve{i}.(pairName).pdf2));
        % Calculate Bhattacharyya distance
        pelvisCurve{i}.(pairName).bhattacharyya = -log(sum(sqrt(pelvisCurve{i}.(pairName).pdf1 .* ...
            pelvisCurve{i}.(pairName).pdf2)));
 
        % Store intersections and distances in allPelvis structure
        allPelvis.curve.inter.(pairName)(i,1) = pelvisCurve{i}.(pairName).intersection;
        allPelvis.curve.bhat.(pairName)(i,1) = pelvisCurve{i}.(pairName).bhattacharyya;
    end
end

% Mean and standard deviation of intersection and Bhattacharyya distance (all pelvises)
for i = 1:length(curvePairs)
    allPelvis.curve.inter.([curvePairs{i}{1} '_mean']) = mean(allPelvis.curve.inter.([curvePairs{i}{1}]));
    allPelvis.curve.inter.([curvePairs{i}{1} '_std']) = std(allPelvis.curve.inter.([curvePairs{i}{1}]));
    allPelvis.curve.bhat.([curvePairs{i}{1} '_mean']) = mean(allPelvis.curve.bhat.([curvePairs{i}{1}]));
    allPelvis.curve.bhat.([curvePairs{i}{1} '_std']) = std(allPelvis.curve.bhat.([curvePairs{i}{1}]));
    allPelvis.curve.meanAbsError.([curvePairs{i}{1} '_mean']) = mean(allPelvis.curve.meanAbsError.([curvePairs{i}{1}]));
    allPelvis.curve.meanAbsError.([curvePairs{i}{1} '_std']) = std(allPelvis.curve.meanAbsError.([curvePairs{i}{1}]));
    allPelvis.curve.meanSquaredError.([curvePairs{i}{1} '_mean']) = mean(allPelvis.curve.meanSquaredError.([curvePairs{i}{1}]));
    allPelvis.curve.meanSquaredError.([curvePairs{i}{1} '_std']) = std(allPelvis.curve.meanSquaredError.([curvePairs{i}{1}]));
    allPelvis.curve.corr.([curvePairs{i}{1} '_mean']) = mean(allPelvis.curve.corr.([curvePairs{i}{1}]));
    allPelvis.curve.corr.([curvePairs{i}{1} '_std']) = std(allPelvis.curve.corr.([curvePairs{i}{1}]));
    % without reference pelvis
    % allPelvis.curve.inter.([curvePairs{i}{1} '_mean']) = mean(allPelvis.curve.inter.([curvePairs{i}{1}])(2:end));
    % allPelvis.curve.inter.([curvePairs{i}{1} '_std']) = std(allPelvis.curve.inter.([curvePairs{i}{1}])(2:end));
    % allPelvis.curve.bhat.([curvePairs{i}{1} '_mean']) = mean(allPelvis.curve.bhat.([curvePairs{i}{1}])(2:end));
    % allPelvis.curve.bhat.([curvePairs{i}{1} '_std']) = std(allPelvis.curve.bhat.([curvePairs{i}{1}])(2:end));
    % allPelvis.curve.meanAbsError.([curvePairs{i}{1} '_mean']) = mean(allPelvis.curve.meanAbsError.([curvePairs{i}{1}])(2:end));
    % allPelvis.curve.meanAbsError.([curvePairs{i}{1} '_std']) = std(allPelvis.curve.meanAbsError.([curvePairs{i}{1}])(2:end));
    % allPelvis.curve.meanSquaredError.([curvePairs{i}{1} '_mean']) = mean(allPelvis.curve.meanSquaredError.([curvePairs{i}{1}])(2:end));
    % allPelvis.curve.meanSquaredError.([curvePairs{i}{1} '_std']) = std(allPelvis.curve.meanSquaredError.([curvePairs{i}{1}])(2:end));
    % allPelvis.curve.corr.([curvePairs{i}{1} '_mean']) = mean(allPelvis.curve.corr.([curvePairs{i}{1}])(2:end));
    % allPelvis.curve.corr.([curvePairs{i}{1} '_std']) = std(allPelvis.curve.corr.([curvePairs{i}{1}])(2:end));
end

clear curveData1 curveData2 curveNormalVertexFace curveGaussVertexFace curveNormalFace curveMeanVertexFace curveMeanFace ...
curveGaussVertexFace curvePairs diff pairIdx pairName %pelvisCurve

%% Save and load properties of Class Curvature
% Class curvature generates a huge amount of data for detailed pelvis meshes -> data are outsourced

% Save class curve data (properties)
savePelvisDataCurve = struct();
% Meta information of curve
metaCurve = metaclass(pelvis(1).curve);
propertiesCurve = {metaCurve.PropertyList.Name};
for i = 1:dataCount
    for j = 1:length(propertiesCurve)
        propertyName = propertiesCurve{j};
        % Cache
        savePelvisDataCurve(i).(propertyName) = pelvis(i).curve.(propertyName);
    end
end
% Save data
save('.\pelvisDataCurveReduced2.mat', 'savePelvisDataCurve', '-v7.3'); % Adapt storage location
clear savePelvisDataCurve metaCurve propertiesCurve propertyName

% Clear class curve (selected)
propertiesCurveToClear = {'boundFaces', 'boundFacesNr', 'boundVertices', 'boundVerticesNr',...
    'allVertexFaceMap', 'vertexFaceMap', 'vertexMap', 'vertexAdjMap', 'faceAdjMap', 'edgeAdjMap',...
    'vertexComPoint', 'comEdge', 'comEdgeNorm', 'baryAreaFaces', 'areaFaceMap', 'angleFaceNormal'}; % Adapt
for i = 1:dataCount
    % Empty fields in pelvis(i).curve
    for j = 1:length(propertiesCurveToClear)
        propertyName = propertiesCurveToClear{j};
        fieldType = class(pelvis(i).curve.(propertyName));
        switch fieldType
            case 'double'
                pelvis(i).curve.(propertyName) = []; % Empty double
            case 'struct'
                pelvis(i).curve.(propertyName) = struct(); % Empty struct
            case 'cell'
                pelvis(i).curve.(propertyName) = {}; % Empty cell
            otherwise
                pelvis(i).curve.(propertyName) = []; % Default to empty
        end
    end
end
clear propertiesCurveToClear propertyName fieldType

% Load class curve data (properties)
loadPelvisDataCurve = load('.\pelvisDataCurveReduced.mat', 'savePelvisDataCurve'); % Adapt storage location
%pelvis(dataCount) = Pelvis; % Initialisation
metaCurve = metaclass(pelvis(1).curve); 
propertiesCurve = {metaCurve.PropertyList.Name}; % load all properties
%propertiesCurve = {'normal', 'cMean', 'gauss'}; % load selected properties
for i = 1:dataCount
    for j = 1:length(propertiesCurve)
        propertyName = propertiesCurve{j};
        if isfield(loadPelvisDataCurve.savePelvisDataCurve(i), propertyName)
            pelvis(i).curve.(propertyName) = loadPelvisDataCurve.savePelvisDataCurve(i).(propertyName);
        end
    end
end
clear loadPelvisDataCurve metaCurve propertiesCurve propertyName

%% Display curvature for the whole pelvis (for control)

% Loop with save figure to save the data/figures
%parfor i = 1 : dataCount
    i = 1;  % Pelvis number
    %figure('Visible','off')  
    figure
    hold on
    curve = 'cMeanVertexFace'; % Choose 'normalVertexFace', 'normalFace', 'cMeanVertexFace', 'cMeanFace', or 'gauss'
    
    switch curve
        case 'normalVertexFace'
            curveType = 'normal';
            rgbType = 'RGBnormVertexFace';
            rangeType = 'topCurveVertexIdx';
        case 'normalFace'
            curveType = 'normal';
            rgbType = 'RGBnormFace';
            rangeType = 'topCurveFaceIdx'; 
            
        case 'cMeanVertexFace'
            curveType = 'cMean';
            rgbType = 'RGBnormVertexFace'; 
            rangeType = 'topCurveVertexIdx'; 
        case 'cMeanFace'
            curveType = 'cMean';
            rgbType = 'RGBnormFace';
            rangeType = 'topCurveFaceIdx'; 
            
        case 'gauss'
            curveType = 'gauss';
            rgbType = 'RGBnormVertexFace';
            rangeType = 'topCurveVertexIdx';
    end
    
    % Fetch the appropriate RGB data and range indices
    rgb = pelvis(i).curve.(curveType).(rgbType);
    rangeIdx = pelvis(i).curve.(curveType).(rangeType);
    
    % Plot the mesh with curvature color mapping
    patch('Faces', pelvis(i).import.processed.faces, ...
          'Vertices', pelvis(i).import.processed.vertices, ...
          'FaceVertexCData', rgb, ...
          'FaceColor', 'flat', ...
          'FaceAlpha', 1, ...
          'EdgeColor', 'none', ...
          'EdgeAlpha', 0.25);
    % Pelvis 3Dview (lightning)
    % patch('Faces',pelvis(i).import.processed.faces,...
    %     'Vertices',pelvis(i).import.processed.vertices,...
    %     'FaceVertexCData', rgb, ...
    %     'FaceColor', 'flat', ...
    %     'FaceAlpha', 1, ...
    %     'EdgeColor','none',...              % Edge color
    %     'EdgeAlpha',0.25,...                % Transparency of the edges
    %     ... % Ligthing for 3d effect
    %     'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
    %     'AmbientStrength', 0.5);
    % light('Position', [1 1 5], 'Style', 'infinite');
    
    % Display acetabulum centre
    plot3(pelvis(i).import.processed.acentre(1), ...
          pelvis(i).import.processed.acentre(2), ...
          pelvis(i).import.processed.acentre(3), ...
          '*', 'MarkerEdgeColor', TUMcolors.blue300, 'MarkerSize', 10);
    % Display curvature range (optional): topCurveRange
    l1 = plot3(pelvis(i).import.processed.vertices(rangeIdx, 1), ...
          pelvis(i).import.processed.vertices(rangeIdx, 2), ...
          pelvis(i).import.processed.vertices(rangeIdx, 3), ...
          '.', 'MarkerEdgeColor', TUMcolors.orange, 'MarkerSize', 5);

    % Legend: colormap
    values = linspace(0, 1, 256);
    colormap('viridis');
    c = colorbar('Ticks', [0, 1], 'TickLabels', {'0', '1'});
    ylabel(c, 'Normalised Curvature', 'Rotation', 90, 'HorizontalAlignment', 'center');
    set(c, 'Position', [0.8, 0.1, 0.05, 0.8]); % [left, bottom, width, height]
    % Set view and axis properties
    title(['Curvature ',curve,': Pelvis ' num2str(i)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z');
    %legend(l1,'Highest Curve Values'); % for topCurveRange
    daspect([1, 1, 1]);
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')',curve,'.fig'])
%end

% Legend: colormap
% values = linspace(0, 1, 256);
% figure%('Position', [100, 100, 400, 1600]);
% colormap('viridis');
% c = colorbar('Ticks', [0, 1], 'TickLabels', {'0', '1'});
% ylabel(c, 'Normalised Curvature', 'Rotation', 90, 'HorizontalAlignment', 'center');
% set(c, 'Position', [0.8, 0.1, 0.05, 0.8]); % [left, bottom, width, height]
% axis off;

clear curve curveType rgbType rangeType rgb rangeIdx l1 values c

%% Display histograms for curvature comparisons

pairName = {
        {'cNormal'},
        {'cMean'},
        {'verNormalGauss'},
        {'verMeanGauss'},
        {'normalMean'},
        {'verNormalMean'},

        {'normalVerMean'},
        {'normalMeanVer'},
        {'normalGaussVer'},
        {'meanGaussVer'}
    };

% Loop with save figure to save the data/figures
%parfor i = 1 : dataCount
    i = 1;  % Pelvis number
    %figure('Visible','off')
    figure
    for pairIdx = 1:length(pairName)
        maxY = max([max(pelvisCurve{i}.(pairName{pairIdx}{1}).pdf1), max(pelvisCurve{i}.(pairName{pairIdx}{1}).pdf2)]);
        figure;
        subplot(1,2,1);
        bar(pelvisCurve{i}.(pairName{pairIdx}{1}).binEdges(1:end-1), pelvisCurve{i}.(pairName{pairIdx}{1}).pdf1, 'BarWidth', 1);
        ylim([0, maxY]);
        title('Curvature Data Set 1');
        subplot(1,2,2);
        bar(pelvisCurve{i}.(pairName{pairIdx}{1}).binEdges(1:end-1), pelvisCurve{i}.(pairName{pairIdx}{1}).pdf2, 'BarWidth', 1);
        ylim([0, maxY]);
        title('Curvature Data Set 2');
        sgtitle(['Comparison ', pairName{pairIdx}{1},': Pelvis ' num2str(i)]);

        % Save figure (figure unvisible, but saved visible)
        %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
        %savefig(['./Figures/Pelvis(',num2str(i),')',char(pairName{pairIdx}),'Histo.fig'])
    end
%end

clear pairName pairIdx maxY

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%% Class Boundary %%%%%%%%%%%%%%
% Scaling Methods (different boundary methods)

% User-defined method %%%
% Options: 'boundBoxHull', 'boundBoxSVD', 'boundSphere', 'distanceLM', 'distanceCentreLM'
methodBound = 'distanceLM'; 

% Boundary: scaling methods
for i = 1:dataCount   
    switch methodBound
        case 'boundBoxHull'
            % Boundary box / cuboid: hull
            level = 1; % User-defined %%%
            pelvis(i).boundaries = pelvis(i).boundaries.boundBoxHull(i, pelvis(i).import.processed.vertices, level); 
        case 'boundBoxSVD'
            % Boundary box / cuboid 2: SVD of the points on the convex hull 
            pelvis(i).boundaries = pelvis(i).boundaries.boundBoxSVD(i, pelvis(i).import.processed.vertices);
        case 'boundSphere'
            % Boundary sphere
            pelvis(i).boundaries = pelvis(i).boundaries.boundSphere(i, pelvis(i).import.processed.vertices);
        case 'distanceLM'
            % Distances landmarks
            pelvis(i).boundaries = pelvis(i).boundaries.distanceLM(i, pelvis(i).import.processed);
        case 'distanceCentreLM'
            % Distances landmarks to acentre
            pelvis(i).boundaries = pelvis(i).boundaries.distanceCentreLM(i, pelvis(i).import.processed);
        otherwise
            % Default: Distances landmarks
            pelvis(i).boundaries = pelvis(i).boundaries.distanceLM(i, pelvis(i).import.processed);
    end
end

% Counts for closest landmarks to acentre
allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', 'aiis', 'iim', ...
    'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};
for j = 1:length(allLandmarks)
    allPelvis.scale.landmarks.centre.LMcounts.(allLandmarks{j}) = 0;
end
% methodBound: 'distanceCentreLM'
if strcmp(methodBound, 'distanceCentreLM')
    for i = 1:dataCount
        % Counts for closest landmarks to acentre
        allPelvis.scale.landmarks.centre.closestLM = pelvis(i).boundaries.landmarks.centre.closestLM;
        for j = 1:length(allPelvis.scale.landmarks.centre.closestLM)
            landmark = allPelvis.scale.landmarks.centre.closestLM{j};
            if isfield(allPelvis.scale.landmarks.centre.LMcounts, landmark)
                allPelvis.scale.landmarks.centre.LMcounts.(landmark) = allPelvis.scale.landmarks.centre.LMcounts.(landmark) + 1;
            end
        end
    end
end 

clear level landmark

%% Scaling pelvis (different boundary methods)

% Reference: first pelvis data
switch methodBound
    % Boundary box / cuboid: hull
    case 'boundBoxHull'
        boundMethod = 'cuboid';
        boundType = 'hull';
        refVariable = 'diagonal'; % Reference variable: diagonal (bounding box)
    % Boundary box / cuboid 2: SVD of the points on the convex hull 
    case 'boundBoxSVD'
        boundMethod = 'cuboid';
        boundType = 'SVD';
        refVariable = 'diagonal'; % Reference variable: diagonal (bounding box)
    % Boundary sphere
    case 'boundSphere'
        boundMethod = 'sphere';
        boundType = 'hull';
        refVariable = 'R';        % Reference variable: radius (sphere)
    % Distances landmarks
    case 'distanceLM'
        boundMethod = 'landmarks';
        boundType = 'all';
        refVariable = 'distance'; % Reference variable: distances (landmarks)
    % Distances landmarks to acentre
    case 'distanceCentreLM'
        boundMethod = 'landmarks';
        boundType = 'centre';
        refVariable = 'distance'; % Reference variable: distances (landmarks)
    otherwise
        boundMethod = 'landmarks';
        boundType = 'all';
        refVariable = 'distance'; % Reference variable: distances (landmarks)
end
refSize = pelvis(1).boundaries.(boundMethod).(boundType).(refVariable);

% Size / scale
for i = 1:dataCount
    pelvis(i).boundaries.(boundMethod).(boundType).sizeAll = ...
        pelvis(i).boundaries.(boundMethod).(boundType).(refVariable) ./ refSize; % Size to reference (percent)
    pelvis(i).boundaries.(boundMethod).(boundType).size = nanmean(pelvis(i).boundaries.(boundMethod).(boundType).sizeAll); % For boundary methods with lanbdmark distances
    pelvis(i).boundaries.(boundMethod).(boundType).scale = 1 / pelvis(i).boundaries.(boundMethod).(boundType).size; % Scaling factor
    sizeAll(i) = pelvis(i).boundaries.(boundMethod).(boundType).size;
    scaleAll(i) = pelvis(i).boundaries.(boundMethod).(boundType).scale;
    if strcmp(methodBound, 'distanceLM') || strcmp(methodBound, 'distanceCentreLM')
        sizePelvis(i,1) = max(pelvis(i).boundaries.(boundMethod).(boundType).(refVariable));
    else
        sizePelvis(i,1) = pelvis(i).boundaries.(boundMethod).(boundType).(refVariable);
    end
end

% Size / Scale (all pelvises)
% without reference pelvis
allPelvis.scale.(boundMethod).(boundType).size(:,1) = sizeAll;
allPelvis.scale.(boundMethod).(boundType).sizeMean = mean(sizeAll(2:end));
allPelvis.scale.(boundMethod).(boundType).sizeStd = std(sizeAll(2:end));
allPelvis.scale.(boundMethod).(boundType).sizeMax = max(sizeAll(2:end));
allPelvis.scale.(boundMethod).(boundType).sizeMin = min(sizeAll(2:end));
allPelvis.scale.(boundMethod).(boundType).scale(:,1) = scaleAll;
allPelvis.scale.(boundMethod).(boundType).scaleMean = mean(scaleAll(2:end));
allPelvis.scale.(boundMethod).(boundType).scaleStd = std(scaleAll(2:end));
allPelvis.scale.(boundMethod).(boundType).scaleMax = max(scaleAll(2:end));
allPelvis.scale.(boundMethod).(boundType).scaleMin = min(scaleAll(2:end));
% Pelvis size
allPelvis.scale.(boundMethod).(boundType).(refVariable) = sizePelvis;
allPelvis.scale.(boundMethod).(boundType).([refVariable 'Mean']) = mean(sizePelvis);
allPelvis.scale.(boundMethod).(boundType).([refVariable 'Std']) = std(sizePelvis);
% without reference pelvis
% allPelvis.scale.(boundMethod).(boundType).([refVariable 'Mean']) = mean(sizePelvis(2:end));
% allPelvis.scale.(boundMethod).(boundType).([refVariable 'Std']) = std(sizePelvis(2:end));

clear refVariable refSize sizeAll scaleAll sizePelvis

%% Save and load properties of Class Boundary
% Save class boundary data (properties)
savePelvisDataBound = struct();
% Meta information of boundary
metaBound = metaclass(pelvis(1).boundaries);
propertiesBound = {metaBound.PropertyList.Name};
for i = 1:dataCount
    for j = 1:length(propertiesBound)
        propertyName = propertiesBound{j};
        % Cache
        savePelvisDataBound(i).(propertyName) = pelvis(i).boundaries.(propertyName);
    end
end
% Save data
save('.\pelvisDataBound.mat', 'savePelvisDataBound', '-v7.3'); % Adapt storage location
clear savePelvisDataBound metaBound propertiesBound propertyName

% Load class boundary data (properties)
loadPelvisDataBound = load('.\pelvisDataBound.mat', 'savePelvisDataBound'); % Adapt storage location
%pelvis(dataCount) = Pelvis; % Initialisation
metaBound = metaclass(pelvis(1).boundaries);
propertiesBound = {metaBound.PropertyList.Name}; % load all properties
%propertiesBound = {'cuboid', 'sphere'}; % load selected properties
for i = 1:dataCount
    for j = 1:length(propertiesBound)
        propertyName = propertiesBound{j};
        if isfield(loadPelvisDataBound.savePelvisDataBound(i), propertyName)
            pelvis(i).boundaries.(propertyName) = loadPelvisDataBound.savePelvisDataBound(i).(propertyName);
        end
    end
end
clear loadPelvisDataBound metaBound propertiesBound propertyName

%% Display pelvis with the boundary box: based on convex hull and its faces or edges (for control)
% Box / cuboid

% Loop with save figure to save the data/figures
%parfor i = 1 : dataCount
    i = 1;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on

    % Pelvis
    % patch('Faces', pelvis(i).import.processed.faces, ...
    %       'Vertices', pelvis(i).import.processed.vertices, ...
    %       'FaceColor', [0.9 0.75 0.68], ...
    %       'FaceAlpha', 1, ...
    %       'EdgeColor', TUMcolors.grey50, ...
    %       'EdgeAlpha', 0.25);
    % Pelvis 3Dview (lightning)
    patch('Faces',pelvis(i).import.processed.faces,...
        'Vertices',pelvis(i).import.processed.vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Display acetabulum centre
    plot3(pelvis(i).import.processed.acentre(1), ...
        pelvis(i).import.processed.acentre(2), ...
    pelvis(i).import.processed.acentre(3), ... 
    '*', 'MarkerEdgeColor', TUMcolors.blue300, ...
        'MarkerSize', 10);

    % Display bounding box
    trisurf(pelvis(i).boundaries.cuboid.hull.tri,...
        pelvis(i).boundaries.cuboid.hull.cornerpoints(:,1),...
        pelvis(i).boundaries.cuboid.hull.cornerpoints(:,2),...
        pelvis(i).boundaries.cuboid.hull.cornerpoints(:,3),...
        'FaceColor',TUMcolors.blue300,'EdgeColor',TUMcolors.blue300,'FaceAlpha',0.25);

    % Display cuboid edges
    for j=1:3
        quiver3(pelvis(i).boundaries.cuboid.hull.cornerpoints(1,1),...
            pelvis(i).boundaries.cuboid.hull.cornerpoints(1,2),...
            pelvis(i).boundaries.cuboid.hull.cornerpoints(1,3),... % start point
            pelvis(i).boundaries.cuboid.hull.edgeVector(j,1),...
            pelvis(i).boundaries.cuboid.hull.edgeVector(j,2),...
            pelvis(i).boundaries.cuboid.hull.edgeVector(j,3),0,'LineWidth',2)
    end

    % Format and display properties
    title(['Pelvis ' num2str(i)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid off
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')CuboidMinHull.fig'])
%end

%% Display pelvis with the boundary box: SVD (for control) 
% Box / cuboid

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount
    i = 1;  % Pelvis number    
    %figure('Visible','off')
    figure
    hold on

    % Pelvis
    % patch('Faces',pelvis(i).import.processed.faces,...
    %     'Vertices',pelvis(i).import.processed.vertices,...
    %     'FaceColor', [0.9 0.75 0.68], ...
    %     'FaceAlpha', 1, ...
    %     'EdgeColor', TUMcolors.grey50, ...
    %     'EdgeAlpha', 0.25);
    % Pelvis 3Dview (lightning)
    patch('Faces',pelvis(i).import.processed.faces,...
        'Vertices',pelvis(i).import.processed.vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Display acetabulum centre
    plot3(pelvis(i).import.processed.acentre(1), ...
        pelvis(i).import.processed.acentre(2), ...
    pelvis(i).import.processed.acentre(3), ... 
    '*', 'MarkerEdgeColor', TUMcolors.blue300, ...
        'MarkerSize', 10);

    % Display bounding box
    trisurf(pelvis(i).boundaries.cuboid.SVD.tri,...
        pelvis(i).boundaries.cuboid.SVD.cornerpoints(:,1),...
        pelvis(i).boundaries.cuboid.SVD.cornerpoints(:,2),...
        pelvis(i).boundaries.cuboid.SVD.cornerpoints(:,3),...
        'FaceColor',TUMcolors.blue300,'EdgeColor',TUMcolors.blue300,'FaceAlpha',0.25);

    % Display cuboid edges
    for j=1:3
        quiver3(pelvis(i).boundaries.cuboid.SVD.cornerpoints(1,1),...
            pelvis(i).boundaries.cuboid.SVD.cornerpoints(1,2),...
            pelvis(i).boundaries.cuboid.SVD.cornerpoints(1,3),... % start point
            pelvis(i).boundaries.cuboid.SVD.edgeVector(j,1),...
            pelvis(i).boundaries.cuboid.SVD.edgeVector(j,2),...
            pelvis(i).boundaries.cuboid.SVD.edgeVector(j,3),0,'LineWidth',2)
    end
    
    % Format and display properties
    title(['Pelvis ' num2str(i)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid off
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')CuboidSVD.fig'])
%end

%% Display pelvis with the boundary sphere (for control)
% Sphere

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount
    i = 1;  % Pelvis number
    %figure('Visible','off')
    figure
    hold on
    
    % Pelvis
    % patch('Faces',pelvis(i).import.processed.faces,...
    %     'Vertices',pelvis(i).import.processed.vertices,...
    %     'FaceColor',[0.8 0.8 0.8], ...      % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',[0.6 0.6 0.6],...       % Edge color
    %     'EdgeAlpha',0.5);                   % Transparency of the edges
    % Pelvis 3Dview (lightning)
    patch('Faces',pelvis(i).import.processed.faces,...
        'Vertices',pelvis(i).import.processed.vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    [x,y,z] = sphere;
    % Scale to desire radius
    x = x * pelvis(i).boundaries.sphere.hull.R;
    y = y * pelvis(i).boundaries.sphere.hull.R;
    z = z * pelvis(i).boundaries.sphere.hull.R;
    % C: Translate sphere to new location
    % Figure: Plot as surface
    surf(x+pelvis(i).boundaries.sphere.hull.centre(1),y+pelvis(i).boundaries.sphere.hull.centre(2),...
        z+pelvis(i).boundaries.sphere.hull.centre(3), ...
        'FaceColor',TUMcolors.blue300,'EdgeColor',TUMcolors.blue300,'FaceAlpha',0.25);
    % scatter3(pelvis(i).import.processed.vertices(:,1),...
    %     pelvis(i).import.processed.vertices(:,2),...
    %     pelvis(i).import.processed.vertices(:,3));

    % Display acetabulum centre
    plot3(pelvis(i).import.processed.acentre(1), ...
        pelvis(i).import.processed.acentre(2), ...
    pelvis(i).import.processed.acentre(3), ... 
    '*', 'MarkerEdgeColor', TUMcolors.blue300, ...
    'MarkerSize', 10);
    
    % Format and display properties
    title(['Pelvis ' num2str(i)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid off
    axis equal;
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')Sphere.fig'])
%end

clear x y z

%% Display landmark distances in bar graph

% Depending on the number of available landmarks
% If there are less than 14 landmarks available, then there are not less landmark distance combinations.
num = nchoosek(length(allLandmarks),2);

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount
    i = 1;  % Pelvis number
    %figure('Visible','off')
    figure
    
    bar(1:num,pelvis(i).boundaries.landmarks.all.distance(:,1))

    % Format and display properties
    title(['Pelvis ' num2str(i)]); 
    xlabel('Distance Nr.')
    ylabel('Distance in mm')
    grid on

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')LM.fig'])
%end

% All Landmark distances
for i = 1:dataCount
    allPelvis.scale.landmarks.all.allDistances(:, i) = pelvis(i).boundaries.landmarks.all.distance(:,1);
end

%figure('Visible','off')
figure
bar(1:num,allPelvis.scale.landmarks.all.allDistances)
% Format and display properties
title('Distances for each Pelvis') % Titel des Diagramms
xlabel('Distance Nr.')
ylabel('Distance in mm')
grid on

% Save figure (figure unvisible, but saved visible)
%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%savefig('./Figures/PelvisLM.fig')

clear num

%% Display landmark distance to acetabulum centre in bar graph

% Depending on the number of available landmarks
% If there are less than 14 landmarks available, then there are not less landmark distance combinations.
num = length(allLandmarks) - 1;

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount
    i = 1;  % Pelvis number
    %figure('Visible','off')
    figure
    
    bar(1:num,pelvis(i).boundaries.landmarks.centre.distance(:,1))

    % Format and display properties
    title(['Pelvis ' num2str(i)]); 
    xlabel('Distance Nr.')
    ylabel('Distance in mm')
    grid on

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')LMcentre.fig'])
%end

% All landmark distances to acetabulum centre
for i = 1:dataCount
    allPelvis.scale.landmarks.centre.centreDistances(:,i) = pelvis(i).boundaries.landmarks.centre.distance(:,1);
end
%figure('Visible','off')
figure
bar(1:num,allPelvis.scale.landmarks.centre.centreDistances)
% Format and display properties
xlabel('Distance Nr.')
ylabel('Distance in mm')
grid on

% Save figure (figure unvisible, but saved visible)
%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%savefig('./Figures/PelvisLMcentre.fig')

clear num

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Class Transform %%%%%%%%%%%%

% Reference: acetabulum centre of the reference pelvis (first pelvis data)
refPoint = pelvis(1).import.processed.acentre; %refPoint acentre 
allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', ...
    'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};
% Initialize refLandmarks with NaNs
refLandmarks = NaN(length(allLandmarks), 3);
% Extract available landmarks
for idx = 1:length(allLandmarks)
    if isfield(pelvis(1).import.processed, allLandmarks{idx})

        refLandmarks(idx, :) = pelvis(1).import.processed.(allLandmarks{idx});
    end
end

% Transformation: Translation/Rotation (SVD/Kabsch)
% User-defined weighting (depends on the number of landmarks) %%%
wKabsch{1} = [1,1,1,1,1,1];
wKabsch{2} = [1,1,1,1,1,1000000000];
wKabsch{3} = [1,1,1,1,1,1,1,1,1,1,1,1,1,1];
wKabsch{4} = [1,1,1,1,1,1000000000,1,1,1,1,1,1,1,1];
wKabsch{5} = [1,1,1,1,1000,1000000000,1000,1000,1000,1,1,1,1000,1]; %%%
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};

parfor i = 1:dataCount

    % Scale
    pelvis(i).transform = pelvis(i).transform.scale(i, pelvis(i).import.processed, ...
        pelvis(i).boundaries.(boundMethod).(boundType).scale); % Scaling method defined in class boundary (see section)

    % Translation to refPoint (acentre) %%% With SCALED data
    % Attention: Scaling leads to new position of the geometry
    pelvis(i).transform = pelvis(i).transform.translation(i,refPoint, ...
        pelvis(i).transform.scaled);

    % Transformation: Translation/Rotation (SVD/Kabsch) %%% With SCALED and TRANSLATED data
    % for w = 1:length(weightKabsch)
    %     pelvis(i).transform = pelvis(i).transform.trafoKabsch(i, refLandmarks, ...
    %         wKabsch{w}, weightKabsch{w}, pelvis(i).transform.trans);
    % end
    
    % Transformation: Translation/Rotation (SVD/Kabsch) %%% With SCALED data
    for w = 1:length(weightKabsch)
        pelvis(i).transform = pelvis(i).transform.trafoKabsch(i, refLandmarks, ...
            wKabsch{w}, weightKabsch{w}, pelvis(i).transform.scaled);
    end
    
end

% Least Root Mean Square of transformation
for i = 1:length(weightKabsch) % Weighting
    allPelvis.trafo.(weightKabsch{i}).lrms = NaN;
    for j = 2:dataCount % without reference pelvis    
        allPelvis.trafo.(weightKabsch{i}).lrms(j,1) = pelvis(j).transform.trafo.(weightKabsch{i}).lrms;
    end
    % LRMS across all pelvises without reference pelvis
    allPelvis.trafo.(weightKabsch{i}).lrmsMean = mean(allPelvis.trafo.(weightKabsch{i}).lrms(2:dataCount));
    allPelvis.trafo.(weightKabsch{i}).lrmsStd = std(allPelvis.trafo.(weightKabsch{i}).lrms(2:dataCount));
end

%% Evaluate transformation: landmark error

% Define all landmarks
allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', ...
    'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};

% Collect errors for each landmark from each pelvis
for i = 1:length(weightKabsch) % Weighting
    for j = 1:length(allLandmarks)
        landmark = allLandmarks{j};
        allPelvis.trafo.(weightKabsch{i}).errorLM.(landmark) = NaN; % Initialize
        for k = 2:numel(pelvis) % First pelvis: reference pelvis
            if isfield(pelvis(k).transform.trafo.(weightKabsch{i}), ['error_' landmark])
                allPelvis.trafo.(weightKabsch{i}).errorLM.(landmark)(k,1) = ...
                    pelvis(k).transform.trafo.(weightKabsch{i}).(['error_' landmark]);
            end
        end
    end
end

% Calculate mean and standard deviation
selectedLandmarks = [allPelvis.scale.landmarks.centre.closestLM, 'acentre'];
for i = 1:length(weightKabsch)
    for j = 1:length(allLandmarks)
        landmark = allLandmarks{j};
        % Calculate mean and standard deviation per landmark
        allPelvis.trafo.(weightKabsch{i}).errorLMmean.(landmark) = ...
            mean(allPelvis.trafo.(weightKabsch{i}).errorLM.(landmark), 'omitnan');
        allPelvis.trafo.(weightKabsch{i}).errorLMstd.(landmark) = ...
            std(allPelvis.trafo.(weightKabsch{i}).errorLM.(landmark), 'omitnan');
    end
    % Calculate mean and standard deviation across all landmarks
    allPelvis.trafo.(weightKabsch{i}).errorLMmeanAll = ... 
        mean(cell2mat(struct2cell(allPelvis.trafo.(weightKabsch{i}).errorLMmean)), 'omitnan');
    allPelvis.trafo.(weightKabsch{i}).errorLMstdAll = ...
        std(cell2mat(struct2cell(allPelvis.trafo.(weightKabsch{i}).errorLMmean)), 'omitnan');
    % Calculate mean and standard deviation for the closest landmarks to acentre (selected landmarks)
    allPelvis.trafo.(weightKabsch{i}).errorLMmeanSelected = [];
    allPelvis.trafo.(weightKabsch{i}).LMselected = selectedLandmarks;
    for j = 1:length(selectedLandmarks)
        lm = selectedLandmarks{j};
        allPelvis.trafo.(weightKabsch{i}).errorLMmeanSelected(end+1,1) = allPelvis.trafo.(weightKabsch{i}).errorLMmean.(lm); 
    end
    % Calculate mean and standard deviation across selected landmarks
    allPelvis.trafo.(weightKabsch{i}).errorLMmeanAllselected = mean(allPelvis.trafo.(weightKabsch{i}).errorLMmeanSelected);
    allPelvis.trafo.(weightKabsch{i}).errorLMstdAllselected = std(allPelvis.trafo.(weightKabsch{i}).errorLMmeanSelected);
end

%% Evaluate transformation: Vertex-to-Nearest-Neighbour Distance

refVertices = pelvis(1).import.processed.vertices;
for i = 1:length(weightKabsch) % Weighting
    allPelvis.trafo.(weightKabsch{i}).vertexError = NaN;
    allPelvis.trafo.(weightKabsch{i}).vertexErrorMax = NaN;
    for j = 2:dataCount % without reference pelvis
        % First: Kabsch transformation of all pelvis vertices (in function trafoKabsch)
        % Vertex-to-Nearest-Neighbour distance
        pelvis(j).transform = pelvis(j).transform.vertexNeighbour(j, refVertices, weightKabsch{i}, ...
            pelvis(j).transform.trafo.(weightKabsch{i}).vertices,pelvis(j).transform.trafo.(weightKabsch{i}).faces);
        % Store Vertex-to-Nearest-Neighbour distance
        allPelvis.trafo.(weightKabsch{i}).vertexError(j,1) = pelvis(j).transform.trafo.(weightKabsch{i}).nearVertexMean;
        allPelvis.trafo.(weightKabsch{i}).vertexErrorMax(j,1) = pelvis(j).transform.trafo.(weightKabsch{i}).nearVertexMax;
    end
    % Error across all pelvises without reference pelvis
    allPelvis.trafo.(weightKabsch{i}).vertexErrorMean = mean(allPelvis.trafo.(weightKabsch{i}).vertexError(2:dataCount));
    allPelvis.trafo.(weightKabsch{i}).vertexErrorStd = std(allPelvis.trafo.(weightKabsch{i}).vertexError(2:dataCount));
    allPelvis.trafo.(weightKabsch{i}).vertexErrorMaxMean = mean(allPelvis.trafo.(weightKabsch{i}).vertexErrorMax(2:dataCount));
    allPelvis.trafo.(weightKabsch{i}).vertexErrorMaxStd = std(allPelvis.trafo.(weightKabsch{i}).vertexErrorMax(2:dataCount));
end

% Global colour coding of verntex-to-nearest-neighbor
for i = 1:length(weightKabsch) % Weighting
    allFaceValues = [];
    for j = 2:dataCount
        allFaceValues = [allFaceValues; pelvis(j).transform.trafo.(weightKabsch{i}).nearVertexFace];
    end
    globalMin = min(allFaceValues);
    globalMax = max(allFaceValues);
    colourNum = 32768;
    colourMap = viridis(colourNum);
    for j = 2:dataCount
        normalData = (pelvis(j).transform.trafo.(weightKabsch{i}).nearVertexFace - globalMin) / (globalMax - globalMin);
        normalData = max(0, min(1, normalData));
        colourIdx = round(normalData * (colourNum - 1)) + 1;
        pelvis(j).transform.trafo.(weightKabsch{i}).nearVertexFaceNormAll = colourMap(colourIdx, :); % rgbColour
    end
end
clear allFacesValues globalMin globalMax colourNum colourMap normalData colourIdx 

%% Evaluate transformation: hounsdorff distance

refVertices = pelvis(1).import.processed.vertices;
for i = 1:length(weightKabsch) % Weighting
    allPelvis.trafo.(weightKabsch{i}).hausdorffDist = NaN;
    for j = 2:dataCount % without reference pelvis        
        % Hausdorff-Distance
        pelvis(j).transform = pelvis(j).transform.hausdorffDistance(j, refVertices, weightKabsch{i}, ...
            pelvis(j).transform.trafo.(weightKabsch{i}).vertices);
        % Store Hausdorff-Distance
        allPelvis.trafo.(weightKabsch{i}).hausdorffDist(j,1) = pelvis(j).transform.trafo.(weightKabsch{i}).hausdorffDist;
    end
    % Distance across all pelvises without reference pelvis
    allPelvis.trafo.(weightKabsch{i}).hausdorffMean = mean(allPelvis.trafo.(weightKabsch{i}).hausdorffDist(2:dataCount));
    allPelvis.trafo.(weightKabsch{i}).hausdorffStd = std(allPelvis.trafo.(weightKabsch{i}).hausdorffDist(2:dataCount));
end

%% Save and load properties of Class Transform

% Save class transform data (properties)
savePelvisDataTransform = struct();
% Meta information of curve
metaTransform = metaclass(pelvis(1).transform);
propertiesTransform = {metaTransform.PropertyList.Name};
for i = 1:dataCount
    for j = 1:length(propertiesTransform)
        propertyName = propertiesTransform{j};
        % Cache
        savePelvisDataTransform(i).(propertyName) = pelvis(i).transform.(propertyName);
    end
end
% Save data
save('.\pelvisDataTransform.mat', 'savePelvisDataTransform', '-v7.3'); % Adapt storage location
clear savePelvisDataTransform metaTransform propertiesTransform propertyName

% Clear class transform (selected)
propertiesTransformToClear = {'scaled', 'trans'}; % Adapt
for i = 1:dataCount
    % Empty fields in pelvis(i).transform
    for j = 1:length(propertiesTransformToClear)
        propertyName = propertiesTransformToClear{j};
        fieldType = class(pelvis(i).transform.(propertyName));
        switch fieldType
            case 'double'
                pelvis(i).transform.(propertyName) = []; % Empty double
            case 'struct'
                pelvis(i).transform.(propertyName) = struct(); % Empty struct
            case 'cell'
                pelvis(i).transform.(propertyName) = {}; % Empty cell
            otherwise
                pelvis(i).transform.(propertyName) = []; % Default to empty
        end
    end
end

% Clear class transform.trafo (selected)
propertiesTransformToClear = {'w1', 'w2', 'w3', 'w4'}; % Adapt
for i = 1:dataCount
    % Empty fields in pelvis(i).transform.trafo
    for j = 1:length(propertiesTransformToClear)
        propertyName = propertiesTransformToClear{j};
        fieldType = class(pelvis(i).transform.trafo.(propertyName));
        switch fieldType
            case 'double'
                pelvis(i).transform.trafo.(propertyName) = []; % Empty double
            case 'struct'
                pelvis(i).transform.trafo.(propertyName) = struct(); % Empty struct
            case 'cell'
                pelvis(i).transform.trafo.(propertyName) = {}; % Empty cell
            otherwise
                pelvis(i).transform.trafo.(propertyName) = []; % Default to empty
        end
    end
end
clear propertiesTransformToClear propertyName fieldType

% Load class transform data (properties)
loadPelvisDataTransform = load('.\pelvisDataTransformReduced.mat', 'savePelvisDataTransform'); % Adapt storage location
%pelvis(dataCount) = Pelvis; % Initialisation
metaTransform = metaclass(pelvis(1).transform);
%propertiesTransform = {metaTransform.PropertyList.Name}; % load all properties
propertiesTransform = {'trafo'}; % load selected properties
for i = 1:dataCount
    for j = 1:length(propertiesTransform)
        propertyName = propertiesTransform{j};
        if isfield(loadPelvisDataTransform.savePelvisDataTransform(i), propertyName)
            pelvis(i).transform.(propertyName) = loadPelvisDataTransform.savePelvisDataTransform(i).(propertyName);
        end
    end
end
clear loadPelvisDataTransform metaTransform propertiesTransform propertyName

%% Display scaled pelvis (for control)

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount
    i = 1;  % Pelvis number
    %figure('Visible','off')
    figure
    hold on
    
    % Pelvis original (optional)
    % patch('Faces',pelvis(i).import.processed.faces,...
    %     'Vertices',pelvis(i).import.processed.vertices,...
    %     'FaceColor',TUMcolors.grey20, ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis original 3Dview (lightning)
    patch('Faces',pelvis(i).import.processed.faces,...
        'Vertices',pelvis(i).import.processed.vertices,...
        'FaceColor',TUMcolors.grey20, ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Pelvis scaled 
    % patch('Faces',pelvis(i).transform.scaled.faces,...
    %     'Vertices',pelvis(i).transform.scaled.vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis scaled 3Dview (lightning)
    patch('Faces',pelvis(i).transform.scaled.faces,...
        'Vertices',pelvis(i).transform.scaled.vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Format and display properties
    title(['Pelvis ' num2str(i)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z');
    legend('Pelvis Original', 'Pelvis Scaled');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')scaled.fig'])
%end

%% Display scaled and shifted pelvises (for control)

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount
    i = 12;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    handles = [];
    labels = {};

    % Pelvis scaled and shifted 
    % patch('Faces',pelvis(i).transform.trans.faces,...
    %     'Vertices',pelvis(i).transform.trans.vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis scaled and shifted 3Dview (lightning)
    patch('Faces',pelvis(i).transform.trans.faces,...
        'Vertices',pelvis(i).transform.trans.vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Landmarks
    for idx = 1:length(allLandmarks)
        landmark = allLandmarks{idx};
        % Check if the landmark is available for the transformed pelvis
        if isfield(pelvis(i).transform.trans, landmark)
            coords = pelvis(i).transform.trans.(landmark);
            % Plot Landmarks
            h1 = plot3(coords(1), coords(2), coords(3), '*', 'Color', TUMcolors.orange, 'MarkerSize', 10);
            handles(1) = h1;
            labels{1} = 'Transformed Landmark';
        end
        % Check if the landmark is available for the reference
        if isfield(pelvis(1).import.processed, landmark)
            ref_coords = pelvis(1).import.processed.(landmark);
            % Plot Reference Landmarks
            h2 = plot3(ref_coords(1), ref_coords(2), ref_coords(3), '*', 'Color', TUMcolors.blue300, 'MarkerSize', 10);
            handles(2) = h2;
            labels{2} = 'Reference Landmark';
        end
    end
    
    % Display acetabulum centre (reference)
    h3 = plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
        pelvis(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10); 
    handles(end+1) = h3;
    labels{end+1} = 'Acetabulum Centre';
    
    % Format and display properties
    title(['Pelvis ' num2str(i)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    legend(handles, labels);
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')trans.fig'])
%end

clear h1 h2 h3 handles labels landmark

%% Display scaled and transformed pelvises (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount
    i = 2;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    handles = [];
    labels = {};
    
    % Pelvis scaled and transformed
    % patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis scaled and transformed 3Dview (lightning)
    patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');
    
    % Landmarks
    for idx = 1:length(allLandmarks)
        landmark = allLandmarks{idx};
        % Check if the landmark is available for the transformed pelvis
        if isfield(pelvis(i).transform.trafo.(weightKabsch{w}), landmark)
            coords = pelvis(i).transform.trafo.(weightKabsch{w}).(landmark);
            % Plot Landmarks
            h1 = plot3(coords(1), coords(2), coords(3), '*', 'Color', TUMcolors.orange, 'MarkerSize', 10);
            handles(1) = h1;
            labels{1} = 'Transformed Landmark';
        end
        % Check if the landmark is available for the reference
        if isfield(pelvis(1).import.processed, landmark)
            ref_coords = pelvis(1).import.processed.(landmark);
            % Plot Reference Landmarks
            h2 = plot3(ref_coords(1), ref_coords(2), ref_coords(3), '*', 'Color', TUMcolors.blue300, 'MarkerSize', 10);
            handles(2) = h2;
            labels{2} = 'Reference Landmark';
        end
    end 
    
    % Reference acetabulum centre
    h3 = plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
        pelvis(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10); 
    handles(end+1) = h3;
    labels{end+1} = 'Acetabulum Centre';
    
    % Format and display properties
    title(['Pelvis ' num2str(i)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    legend(handles, labels);
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')trafo.fig'])
%end

clear h1 h2 h3 handles labels landmark

%% Display all scaled and transformed pelvises (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

%figure('Visible','off')
figure
hold on
% Predefined color lists
patchColors = lines(dataCount); % Using the MATLAB 'lines' colormap
markerColors = jet(dataCount);  % Using the MATLAB 'jet' colormap

for i = 1:dataCount
    % Plot with changing FaceColor
    % Pelvis scaled and transformed
    % patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
    %     'FaceColor',patchColors(i,:), ...   % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis scaled and transformed 3Dview (lightning)
    patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',patchColors(i,:), ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    currentMarkerColor = markerColors(i,:); % Single marker color for all landmarks in this iteration

    % Plot each available landmark with the currentMarkerColor
    for idx = 1:length(allLandmarks)
        landmark = allLandmarks{idx};
        % Check if the landmark is available for the transformed pelvis
        if isfield(pelvis(i).transform.trafo.(weightKabsch{w}), landmark)
            coords = pelvis(i).transform.trafo.(weightKabsch{w}).(landmark);
            % Plot Landmarks
            plot3(coords(1), coords(2), coords(3), '*', 'Color', currentMarkerColor, 'MarkerSize', 10);
        end
    end
end

% % without color change
% for i = 1:dataCount
%     % Pelvis scaled and transformed
%     % patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
%     %     'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
%     %     'FaceColor',[0.9 0.75 0.68], ...   % Face color
%     %     'FaceAlpha',1,...                   % Transparency of the faces
%     %     'EdgeColor',TUMcolors.grey50,...    % Edge color
%     %     'EdgeAlpha',0.25);                  % Transparency of the edges
%     % Pelvis scaled and transformed 3Dview (lightning)
%     patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
%         'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
%         'FaceColor',[0.9 0.75 0.68], ...   % Face color
%         'FaceAlpha',1,...                   % Transparency of the faces
%         'EdgeColor','none',...              % Edge color
%         'EdgeAlpha',0.25,...                % Transparency of the edges
%         ... % Ligthing for 3d effect
%         'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
%         'AmbientStrength', 0.5);
%     light('Position', [1 1 5], 'Style', 'infinite');
%
%     for idx = 1:length(allLandmarks)
%         landmark = allLandmarks{idx};
%         % Check if the landmark is available for the transformed pelvis
%         if isfield(pelvis(i).transform.trafo.(weightKabsch{w}), landmark)
%             coords = pelvis(i).transform.trafo.(weightKabsch{w}).(landmark);
%             % Plot Landmarks
%             plot3(coords(1), coords(2), coords(3), '*', 'Color', TUMcolors.orange, 'MarkerSize', 10);
%         end
%     end
% end

for idx = 1:length(allLandmarks)
    landmark = allLandmarks{idx};
    % Check if the landmark is available for the reference
    if isfield(pelvis(1).import.processed, landmark)
        ref_coords = pelvis(1).import.processed.(landmark);
        % Plot Reference Landmarks
        plot3(ref_coords(1), ref_coords(2), ref_coords(3), '*', 'Color', TUMcolors.blue300, 'MarkerSize', 10);
    end
end

% Reference acetabulum centre 
plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
    pelvis(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10);

% Format and display properties
title(['Pelvis ' num2str(i)]); 
xlabel('X'); ylabel('Y'); zlabel('Z');
daspect([1, 1, 1]); % Equal aspect ratio for the axes
view(3);
hold off;

% Save figure (figure unvisible, but saved visible)
%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%savefig(['./Figures/PelvisTrafo.fig'])

clear landmark currentMarkerColor coords markerColors patchColors

%% Display scaled and transformed pelvis with reference pelvis (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount
    i = 1;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    handles = [];
    labels = {};
    
    % Pelvis scaled and transformed
    % h1 = patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...   % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis scaled and transformed 3Dview (lightning)
    h1 = patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');
    handles(1) = h1;
    labels{1} = 'Transformed Pelvis';
    % Reference pelvis
    h2 = patch('Faces',pelvis(1).import.processed.faces,...
        'Vertices',pelvis(1).import.processed.vertices,...
        'FaceColor',TUMcolors.grey20, ...    % Face color
        'FaceAlpha',0.5,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');
    handles(2) = h2;
    labels{2} = 'Reference Pelvis';

    % Landmarks
    for idx = 1:length(allLandmarks)
        landmark = allLandmarks{idx};
        % Check if the landmark is available for the transformed pelvis
        if isfield(pelvis(i).transform.trafo.(weightKabsch{w}), landmark)
            coords = pelvis(i).transform.trafo.(weightKabsch{w}).(landmark);
            % Plot Landmarks
            h3 = plot3(coords(1), coords(2), coords(3), '*', 'Color', TUMcolors.orange, 'MarkerSize', 10);
            handles(3) = h3;
            labels{3} = 'Transformed Landmark';
        end
        % Check if the landmark is available for the reference
        if isfield(pelvis(1).import.processed, landmark)
            ref_coords = pelvis(1).import.processed.(landmark);
            % Plot Reference Landmarks
            h4 = plot3(ref_coords(1), ref_coords(2), ref_coords(3), '*', 'Color', TUMcolors.blue300, 'MarkerSize', 10);
            handles(4) = h4;
            labels{4} = 'Reference Landmark';
        end
    end 
    
    % Reference acetabulum centre
    h5 = plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
        pelvis(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10); 
    handles(5) = h5;
    labels{5} = 'Acetabulum Centre';
    
    % Format and display properties
    title(['Pelvis ' num2str(i)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    xlim([-50 150]); % for same scaling in different plots
    ylim([-20 150]);
    zlim([-250 20]);
    %view(3);
    view(80,-12.5) % View for images
    legend(handles, labels);
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')RefTrafo.fig'])
%end

clear h1 h2 h3 h4 h5 handles labels landmark ref_coords

%% Display scaled and transformed pelvis (incl. vertex-to-neares-neighbour color) with reference pelvis (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 2:dataCount
    i = 1;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    handles = [];
    labels = {};
    
    % Pelvis scaled and transformed
    % patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
    %     'FaceVertexCData',pelvis(i).transform.trafo.(weightKabsch{w}).nearVertexFaceNorm, ...   % Face color
    %     'FaceColor', 'flat', ...
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis scaled and transformed 3Dview (lightning)
    patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
        ... % Face color: % local: .nearVertexFaceNorm % global: .nearVertexFaceNormAll
        'FaceVertexCData',pelvis(i).transform.trafo.(weightKabsch{w}).nearVertexFaceNorm, ...  
        'FaceColor', 'flat', ...
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                   % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 1);
   light('Position', [1 1 -5], 'Style', 'infinite');
    material dull; % Alternativen: shiny, dull, metal
    % Reference pelvis
    patch('Faces',pelvis(1).import.processed.faces,...
        'Vertices',pelvis(1).import.processed.vertices,...
        'FaceColor',TUMcolors.grey20, ...    % Face color
        'FaceAlpha',0.5,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Landmarks
    for idx = 1:length(allLandmarks)
        landmark = allLandmarks{idx};
        % Check if the landmark is available for the transformed pelvis
        if isfield(pelvis(i).transform.trafo.(weightKabsch{w}), landmark)
            coords = pelvis(i).transform.trafo.(weightKabsch{w}).(landmark);
            % Plot Landmarks
            h1 = plot3(coords(1), coords(2), coords(3), '.', 'Color', 'r', 'MarkerSize', 40); % TUMcolors.orange
            handles(1) = h1;
            labels{1} = 'Transformed Landmark';
        end
        % Check if the landmark is available for the reference
        if isfield(pelvis(1).import.processed, landmark)
            ref_coords = pelvis(1).import.processed.(landmark);
            % Plot Reference Landmarks
            h2 = plot3(ref_coords(1), ref_coords(2), ref_coords(3), '.', 'Color', 'g', 'MarkerSize', 40); % TUMcolors.blue300
            handles(2) = h2;
            labels{2} = 'Reference Landmark';
        end
    end 
    
    % Reference acetabulum centre 
    h3 = plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
        pelvis(1).import.processed.acentre(3), ...
        '.','Color','m', 'MarkerSize', 45); %'r'
    plot3(pelvis(1).import.processed.acentre(1), ...
          pelvis(1).import.processed.acentre(2), ...
          pelvis(1).import.processed.acentre(3), ...
          'x', 'MarkerEdgeColor', 'm', 'MarkerSize', 22, 'LineWidth', 5);
    handles(end+1) = h3;
    labels{end+1} = 'Acetabulum Centre';
    
    % Legend: colormap
    colormap('viridis');
    c = colorbar('Ticks', [0, 1], 'TickLabels', {'0', '1'});
    ylabel(c, 'Normalised Curvature', 'Rotation', 90, 'HorizontalAlignment', 'center');
    set(c, 'Position', [0.8, 0.1, 0.02, 0.8]); % [left, bottom, width, height]
    % Format and display properties
    title(['Pelvis ' num2str(i)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    xlim([-50 160]); % for same scaling in different plots
    ylim([-20 150]);
    zlim([-250 20]);
    %view(3);
    view(80,-12.5) % View for images
    legend(handles, labels);
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')trafoNearVertex.fig'])
%end

clear h1 h2 h3 handles labels landmark

%% Display scaled and transformed pelvis with curvature (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount
    i = 1;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    handles = [];
    labels = {};

    curve = 'cMeanVertexFace'; % Choose 'normalVertexFace', 'normalFace', 'cMeanVertexFace', 'cMeanFace', or 'gauss'
    
    switch curve
        case 'normalVertexFace'
            curveType = 'normal';
            rgbType = 'RGBnormVertexFace';
            rangeType = 'topCurveVertexIdx';
        case 'normalFace'
            curveType = 'normal';
            rgbType = 'RGBnormFace';
            rangeType = 'topCurveFaceIdx'; 
            
        case 'cMeanVertexFace'
            curveType = 'cMean';
            rgbType = 'RGBnormVertexFace'; 
            rangeType = 'topCurveVertexIdx'; 
        case 'cMeanFace'
            curveType = 'cMean';
            rgbType = 'RGBnormFace';
            rangeType = 'topCurveFaceIdx'; 
            
        case 'gauss'
            curveType = 'gauss';
            rgbType = 'RGBnormVertexFace';
            rangeType = 'topCurveVertexIdx';
    end
    % Fetch the appropriate RGB data and range indices
    rgb = pelvis(i).curve.(curveType).(rgbType);
    
    % Pelvis scaled and transformed
    % patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
    %     'FaceVertexCData',rgb, ...   % Face color
    %     'FaceColor', 'flat', ...
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor','none',...              % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis scaled and transformed 3Dview (lightning)
    patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceVertexCData',rgb, ...   % Face color
        'FaceColor', 'flat', ...
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 1);
    light('Position', [1 1 -5], 'Style', 'infinite');
    material dull; % Alternativen: shiny, dull, metal
    % Reference pelvis
    patch('Faces',pelvis(1).import.processed.faces,...
        'Vertices',pelvis(1).import.processed.vertices,...
        'FaceColor',TUMcolors.grey20, ...    % Face color
        'FaceAlpha',0.5,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Landmarks (optional)
    % for idx = 1:length(allLandmarks)
    %     landmark = allLandmarks{idx};
    %     % Check if the landmark is available for the transformed pelvis
    %     if isfield(pelvis(i).transform.trafo.(weightKabsch{w}), landmark)
    %         coords = pelvis(i).transform.trafo.(weightKabsch{w}).(landmark);
    %         % Plot Landmarks
    %         h1 = plot3(coords(1), coords(2), coords(3), '.', 'Color', TUMcolors.orange, 'MarkerSize', 30);
    %         handles(1) = h1;
    %         labels{1} = 'Transformed Landmark';
    %     end
    %     % Check if the landmark is available for the reference
    %     if isfield(pelvis(1).import.processed, landmark)
    %         ref_coords = pelvis(1).import.processed.(landmark);
    %         % Plot Reference Landmarks
    %         h2 = plot3(ref_coords(1), ref_coords(2), ref_coords(3), '.', 'Color', TUMcolors.blue300, 'MarkerSize', 30);
    %         handles(2) = h2;
    %         labels{2} = 'Reference Landmark';
    %     end
    % end 

    % Reference acetabulum centre 
    h3 = plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
        pelvis(1).import.processed.acentre(3), '.','Color','m', 'MarkerSize', 50); %'r'
    handles(end+1) = h3;
    labels{end+1} = 'Acetabulum Centre';
    
    % Legend: colormap
    colormap('viridis');
    c = colorbar('Ticks', [0, 1], 'TickLabels', {'0', '1'});
    ylabel(c, 'Normalised Curvature', 'Rotation', 90, 'HorizontalAlignment', 'center');
    set(c, 'Position', [0.8, 0.1, 0.02, 0.8]); % [left, bottom, width, height]
    % Format and display properties
    title(['Pelvis ' num2str(i) ' - ' curve]); 
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    xlim([-50 160]); % for same scaling in different plots
    ylim([-20 150]);
    zlim([-250 20]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(3);
    view(80,-12.5) % View for images
    legend(handles, labels);
    grid off;
    hold off;
        
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')trafoCurve.fig'])
%end

clear h1 h2 h3 handles labels landmark rangeType rgb rgbType curveType curve

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Class Defect %%%%%%%%%%%%%%%

% Folder geometries defects
filedirDefect = './GeometriesDefects/';
% Lists of stl and xlsx files
stlfilesDefect = dir(fullfile(filedirDefect, '*.stl'));
xlsfilesDefect = dir(fullfile(filedirDefect, '*.xlsx'));
% Use stl/xls files count to determine data count
dataCountDefect = length(stlfilesDefect);
%dataCountDefect = length(xlsfilesDefect);

% Check for mismatch or absence of files
if isempty(stlfilesDefect)
    error('No stl files for defect found in the specified directory.');
elseif isempty(xlsfilesDefect)
    error('No xlsx files for defect found in the specified directory.');
elseif length(stlfilesDefect) ~= length(xlsfilesDefect)
    error('Number of stl files for defect does not match number of xlsx files.');
end

% Initalization: Generate object array
pelvisDefect(dataCountDefect,1) = Defect;

%% %%%%%%%%%%%% Class Import %%%%%%%%%%%%%%%
% Import pelvis defect data

% File paths
filePath_defectData = fullfile(filedirDefect, {xlsfilesDefect.name});   % Pelvis defect data/infos
filePath_defect = fullfile(filedirDefect, {stlfilesDefect.name});       % Pelvis defect geometry data

% Parallel loop for data import and processing
parfor i = 1:dataCountDefect

   % Initalization: Generate object
   pelvisDefectImport = Import(i); % cache
      
   % Import the defect data/infos
   pelvisDefectImport = pelvisDefectImport.importData(i,filePath_defectData{i});
  
   % Import the defect geometry data
   pelvisDefectImport = pelvisDefectImport.importSTL(i, filePath_defect{i}) ...
                               .centroid(i) ... % Centroid of the faces (for the normal vector)
                               .edgeLength(i); % Edge length of  mesh

   % stl properties
   verticesDefect(i,1) = size(pelvisDefectImport.loaded.TR.Points,1);
   facesDefect(i,1) = size(pelvisDefectImport.loaded.TR.ConnectivityList,1);

   % Mean edge length
   edgeLengthDefect(i,1) = pelvisDefectImport.processed.meanEdgeLength;

   % Paprosky
   paproskyTypeDefect(i,1) = pelvisDefectImport.loaded.paprosky;

   % Store the result back in pelvis
   pelvisDefect(i).import = pelvisDefectImport;

end

% Paprosky classes
allDefect.paprosky.type = paproskyTypeDefect;
paproskyTypes = {'2a', '2b', '2c', '3a', '3b', 'NaN'};
paproskyNames = {'IIa', 'IIb', 'IIc', 'IIIa', 'IIIb', 'NaN'};
allDefect.paprosky.counts = struct();
for i = 1:length(paproskyNames)
    allDefect.paprosky.counts.(paproskyNames{i}) = 0; 
end
for i = 1:length(paproskyTypeDefect)
    idx = find(strcmp(paproskyTypeDefect{i}, paproskyTypes));
    if ~isempty(idx)
        allDefect.paprosky.counts.(paproskyNames{idx}) = allDefect.paprosky.counts.(paproskyNames{idx}) + 1;
    end
end
disp(allDefect.paprosky.counts);

% STL properties
allDefect.stl.vertices = verticesDefect;
allDefect.stl.verticesMean = mean(verticesDefect);
allDefect.stl.verticesStd = std(verticesDefect);
allDefect.stl.verticesMax = max(verticesDefect);
allDefect.stl.verticesMin = min(verticesDefect);
allDefect.stl.faces = facesDefect;
allDefect.stl.facesMean = mean(facesDefect);
allDefect.stl.facesStd = std(facesDefect);
allDefect.stl.facesMax = max(facesDefect);
allDefect.stl.facesMin = min(facesDefect);
% Mean edge Length
allDefect.stl.edgeLength = edgeLengthDefect;
allDefect.stl.edgeMean = mean(edgeLengthDefect);
allDefect.stl.edgeStd = std(edgeLengthDefect);
allDefect.stl.edgeMax = max(edgeLengthDefect); % max of mean
allDefect.stl.edgeMin = min(edgeLengthDefect); % min of mean
% without reference pelvis
% allDefect.stl.verticesMean = mean(verticesDefect(2:end));
% allDefect.stl.verticesStd = std(verticesDefect(2:end));
% allDefect.stl.verticesMax = max(verticesDefect(2:end));
% allDefect.stl.verticesMin = min(verticesDefect(2:end));
% allDefect.stl.facesMean = mean(facesDefect(2:end));
% allDefect.stl.facesStd = std(facesDefect(2:end));
% allDefect.stl.facesMax = max(facesDefect(2:end));
% allDefect.stl.facesMin = min(facesDefect(2:end));
% % Mean edge Length
% allDefect.stl.edgeMean = mean(edgeLengthDefect(2:end));
% allDefect.stl.edgeStd = std(edgeLengthDefect(2:end));
% allDefect.stl.edgeMax = max(edgeLengthDefect(2:end)); % max of mean
% allDefect.stl.edgeMin = min(edgeLengthDefect(2:end)); % min of mean

clear pelvisDefectImport

%% Save and load properties of Class Import (Defect)

% Save class import data (properties)
savePelvisDefectDataImport = struct();
% Meta information of import
metaImportDefect = metaclass(pelvisDefect(1).import);
propertiesImportDefect = {metaImportDefect.PropertyList.Name};
parfor i = 1:dataCountDefect
    for j = 1:length(propertiesImportDefect)
        propertyName = propertiesImportDefect{j};
        % Cache
        savePelvisDefectDataImport(i).(propertyName) = pelvisDefect(i).import.(propertyName);
    end
end
% Save data
save('.\pelvisDefectDataImport.mat', 'savePelvisDefectDataImport', '-v7.3'); % Adapt storage location
clear savePelvisDefectDataImport metaImportDefect propertiesImportDefect propertyName

% Load class import data (properties)
loadPelvisDefectDataImport = load('.\pelvisDefectDataImport.mat', 'savePelvisDefectDataImport'); % Adapt storage location
pelvisDefect(dataCountDefect) = Defect; % Initialisation
metaImportDefect = metaclass(pelvisDefect(1).import);
%propertiesImportDefect = {metaImportDefect.PropertyList.Name}; % load all properties
propertiesImportDefect = {'processed'}; % load selected properties
for i = 1:dataCountDefect
    for j = 1:length(propertiesImportDefect)
        propertyName = propertiesImportDefect{j};
        if isfield(loadPelvisDefectDataImport.savePelvisDefectDataImport(i), propertyName)
            pelvisDefect(i).import.(propertyName) = loadPelvisDefectDataImport.savePelvisDefectDataImport(i).(propertyName);
        end
    end
end
clear loadPelvisDefectDataImport metaImportDefect propertiesImportDefect propertyName

%% Display defect with the face normals (for control)

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 1;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on

    % Pelvis defect
    patch('Faces',pelvisDefect(i).import.processed.faces,...
        'Vertices',pelvisDefect(i).import.processed.vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',TUMcolors.grey50,...    % Edge color
        'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect 3Dview (lightning)
    % patch('Faces',pelvisDefect(i).import.processed.faces,...
    %     'Vertices',pelvisDefect(i).import.processed.vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor','none',...              % Edge color
    %     'EdgeAlpha',0.25,...                % Transparency of the edges
    %     ... % Ligthing for 3d effect
    %     'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
    %     'AmbientStrength', 0.5);
    % light('Position', [1 1 5], 'Style', 'infinite');

    % Normals of the faces (optional)
    % quiver3(pelvisDefect(i).import.processed.centreFaces(:,1),...
    %     pelvisDefect(i).import.processed.centreFaces(:,2),...
    %     pelvisDefect(i).import.processed.centreFaces(:,3),...
    %     pelvisDefect(i).import.processed.normals(:,1),...
    %     pelvisDefect(i).import.processed.normals(:,2),...
    %     pelvisDefect(i).import.processed.normals(:,3),2,...
    %     'Color',TUMcolors.orange)

    % Display acetabulum centre
    plot3(pelvisDefect(i).import.processed.acentre(1), ...
        pelvisDefect(i).import.processed.acentre(2), ...
        pelvisDefect(i).import.processed.acentre(3), '*', ... 
        'MarkerEdgeColor', TUMcolors.blue300, ...
        'MarkerSize', 10);

    % Set view and axis properties
    % grid minor
    title(['Pelvis ' num2str(i)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),').fig'])
%end

%% Display defect with the whole pelvis (for control)

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 1;  % Pelvis number
    %figure('Visible','off')  
    figure
    hold on

    % Pelvis 3Dview (lightning)
    patch('Faces', pelvis(i).import.processed.faces, ...
        'Vertices', pelvis(i).import.processed.vertices, ...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Defect area
    patch('Faces', pelvisDefect(i).import.processed.faces, ...
          'Vertices', pelvisDefect(i).import.processed.vertices, ...
          'FaceColor', TUMcolors.green, ...
          'FaceAlpha', 1, ...
          'EdgeColor', 'none', ...
          'EdgeAlpha', 0.25);
    view(3)
    % Display acetabulum centre
    plot3(pelvisDefect(i).import.processed.acentre(1), ...
          pelvisDefect(i).import.processed.acentre(2), ...
          pelvisDefect(i).import.processed.acentre(3), ...
          '*', 'MarkerEdgeColor', TUMcolors.blue300, 'MarkerSize', 10);

    % Set view and axis properties
    title(['Pelvis defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    legend('Pelvis', 'Pelvis Defect', 'Acetabulum Centre');
    daspect([1, 1, 1]); % equal axis
    view(3);
    %view(80,-12.5) % View for images
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/Pelvis(',num2str(i),')PelvisDefect.fig'])

%end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Class Curvature %%%%%%%%%%%%

% Curvature of defect area (open mesh)
% (1) Normal Curvature
% (2) Mean Curvature
% (3) Gauss Curvature

% Determine the curvature
% with orginial data (not mirrored or transformed)
% Preallocate array for storing times
executionTimesCurve = zeros(dataCountDefect, 1);
for i = 1:dataCountDefect %parfor
    tic
    % Open mesh: boundary without curvature (faces and vertices)
    pelvisDefect(i).curve = pelvisDefect(i).curve.bound(i,pelvisDefect(i).import.loaded.TR);

    % Determine adjacent vertices / faces (vertex face map)
    verticesNr = (1:size(pelvisDefect(i).import.loaded.TR.Points,1))';
    pelvisDefect(i).curve = pelvisDefect(i).curve.adjVertex(i,verticesNr,pelvisDefect(i).import.loaded.TR.ConnectivityList);
    % Sort adjacent vertices / faces (vertex face map)
    pelvisDefect(i).curve = pelvisDefect(i).curve.adjVertexSeq(i,verticesNr);
    % Calculate the edges and face areas
    % For mean (2) and gauss (3) curvature
    pelvisDefect(i).curve = pelvisDefect(i).curve.edgesArea(i,verticesNr,pelvisDefect(i).import.loaded.TR.Points);
    % Calculate the angles between the faces (normal vector)
    % For normal (1) and mean (2) curvature
    pelvisDefect(i).curve = pelvisDefect(i).curve.angleFaces(i,verticesNr,pelvisDefect(i).import.loaded.normals);
    % Determine adjacent faces
    pelvisDefect(i).curve = pelvisDefect(i).curve.adjFaces(i,pelvisDefect(i).import.loaded.TR.ConnectivityList,...
        pelvisDefect(i).import.loaded.TR.Points);

    % Extract pelvis name
    underscoreIndex = strfind(stlfilesDefect(i).name, '_');
    pelvisID = stlfilesDefect(i).name(1:underscoreIndex(1)-1);

    % Determine the curvature (with stl write)
    topCurveRange = 2.5; % in percent
    % (1) Normal Curvature
    pelvisDefect(i).curve = pelvisDefect(i).curve.curveNormal(i,pelvisID,pelvisDefect(i).import.loaded.normals,...
        pelvisDefect(i).import.loaded.TR.ConnectivityList,pelvisDefect(i).import.loaded.TR.Points,topCurveRange);
    % (2) Mean Curvature
    pelvisDefect(i).curve = pelvisDefect(i).curve.curveMean(i,pelvisID,pelvisDefect(i).import.loaded.normals,...
        pelvisDefect(i).import.loaded.TR.ConnectivityList,pelvisDefect(i).import.loaded.TR.Points,topCurveRange);
    % (3) Gauss Curvature
    pelvisDefect(i).curve = pelvisDefect(i).curve.curveGauss(i,pelvisID,verticesNr,...
        pelvisDefect(i).import.loaded.TR.ConnectivityList,pelvisDefect(i).import.loaded.TR.Points,topCurveRange);
    executionTimesCurve(i) = toc; % End timer and store time in array

end

% Execution time of curvature
allDefect.time.curve = executionTimesCurve;
allDefect.time.curveMean = mean(executionTimesCurve);
allDefect.time.curveStd = std(executionTimesCurve);
allDefect.time.curveMax = max(executionTimesCurve);
allDefect.time.curveMin = min(executionTimesCurve);

% Data storage of curvature data
% Preallocation
allDefect.size.curveInfo = cell(length(pelvisDefect), 1);
allDefect.size.curve = zeros(length(pelvisDefect),1); 
% Size of pelvisDefect(i).curve variable 
for i = 1:dataCount
    sizeCurve = pelvisDefect(i).curve; % cache
    allDefect.size.curveInfo{i} = whos('sizeCurve');
    allDefect.size.curve(i) = allDefect.size.curveInfo{i}.bytes;
end
allDefect.size.curveMean = mean(allDefect.size.curve);
allDefect.size.curveStd = std(allDefect.size.curve);
allDefect.size.curveMax = max(allDefect.size.curve);
allDefect.size.curveMin = min(allDefect.size.curve);
clearvars sizeCurve

clear executionTimesCurve verticesNr underscoreIndex pelvisID topCurveRange

%% Validation defect area (with curve range) 5

curveMethods = {'normalVertexFace', 'normalFace', 'cMeanVertexFace', 'cMeanFace', 'gauss'};
% Choose 'normalVertexFace', 'normalFace', 'cMeanVertexFace', 'cMeanFace', or 'gauss'

for methodIdx = 1:length(curveMethods)
    curve = curveMethods{methodIdx};

    switch curve
        case 'normalVertexFace'
            curveType = 'normal';
            rangeType = 'topCurveVertexIdx';
            curveName = 'curveVertex';
        case 'normalFace'
            curveType = 'normal';
            rangeType = 'topCurveFaceIdx';
            curveName = 'curveFace';

        case 'cMeanVertexFace'
            curveType = 'cMean';
            rangeType = 'topCurveVertexIdx';
            curveName = 'curveVertex';
        case 'cMeanFace'
            curveType = 'cMean';
            rangeType = 'topCurveFaceIdx';
            curveName = 'curveFace';

        case 'gauss'
            curveType = 'gauss';
            rangeType = 'topCurveVertexIdx';
            curveName = 'curveVertex';
    end

    for i = 1:dataCountDefect
        % Top curvature values
        topCurveVerticesIdx = pelvis(i).curve.(curveType).(rangeType);

        % Validation
        % (previously: calculated boundary vertices in class curvature)
        pelvisDefect(i).curve = pelvisDefect(i).curve.curveRangeBound(i,pelvis(i).import.loaded.TR.Points,...
            topCurveVerticesIdx,curveType,curveName);
        % Cache
        curveBoundDistance(i,1) = pelvisDefect(i).curve.(curveType).([curveName, 'BoundDistanceMean']);
    end

    clear topCurveVerticesIdx

    % Validation of defect area
    % Without reference pelvis; reference pelvis without defect
    allDefect.curve.(curveType).([curveName, 'BoundDistanceMean']) = mean(curveBoundDistance(2:end));
    allDefect.curve.(curveType).([curveName, 'BoundDistanceStd']) = std(curveBoundDistance(2:end));
    allDefect.curve.(curveType).([curveName, 'BoundDistanceMax']) = max(curveBoundDistance(2:end));
    allDefect.curve.(curveType).([curveName, 'BoundDistanceMin']) = min(curveBoundDistance(2:end));

end

clear curveMethods methodIdx curveType rangeType curveName topCurveVerticesIdx curveBoundDistance

%% Save and load properties of Class Curvature (Defect)
% Class curvature generates a huge amount of data for detailed pelvis meshes -> data are outsourced

% Save class curve data (properties)
savePelvisDefectDataCurve = struct();
% Meta information of curve
metaCurveDefect = metaclass(pelvisDefect(1).curve);
propertiesCurveDefect = {metaCurveDefect.PropertyList.Name};
for i = 1:dataCountDefect
    for j = 1:length(propertiesCurveDefect)
        propertyName = propertiesCurveDefect{j};
        % Cache
        savePelvisDefectDataCurve(i).(propertyName) = pelvisDefect(i).curve.(propertyName);
    end
end
% Save data
save('.\pelvisDefectDataCurve.mat', 'savePelvisDefectDataCurve', '-v7.3'); % Adapt storage location
clear savePelvisDefectDataCurve metaCurveDefect propertiesCurveDefect propertyName

% Clear class curve (selected)
propertiesCurveDefectToClear = {'boundFaces', 'boundFacesNr', ...
    'allVertexFaceMap', 'vertexFaceMap', 'vertexMap', 'vertexAdjMap', 'faceAdjMap', 'edgeAdjMap',...
    'vertexComPoint', 'comEdge', 'comEdgeNorm', 'baryAreaFaces', 'areaFaceMap', 'angleFaceNormal'}; % Adapt
for i = 1:dataCountDefect
    % Empty fields in pelvis(i).curve
    for j = 1:length(propertiesCurveDefectToClear)
        propertyName = propertiesCurveDefectToClear{j};
        fieldType = class(pelvisDefect(i).curve.(propertyName));
        switch fieldType
            case 'double'
                pelvisDefect(i).curve.(propertyName) = []; % Empty double
            case 'struct'
                pelvisDefect(i).curve.(propertyName) = struct(); % Empty struct
            case 'cell'
                pelvisDefect(i).curve.(propertyName) = {}; % Empty cell
            otherwise
                pelvisDefect(i).curve.(propertyName) = []; % Default to empty
        end
    end
end
clear propertiesCurveDefectToClear propertyName

% Load class curve data (properties)
loadPelvisDefectDataCurve = load('.\pelvisDefectDataCurveReduced.mat', 'savePelvisDefectDataCurve'); % Adapt storage location
%pelvisDefect(dataCountDefect) = Defect; % Initialisation
metaCurveDefect = metaclass(pelvisDefect(1).curve);
propertiesCurveDefect = {metaCurveDefect.PropertyList.Name};
%propertiesCurveDefect = {'normal', 'cMean'}; % load selected properties
for i = 1:dataCountDefect
    for j = 1:length(propertiesCurveDefect)
        propertyName = propertiesCurveDefect{j};
        if isfield(loadPelvisDefectDataCurve.savePelvisDefectDataCurve(i), propertyName)
            pelvisDefect(i).curve.(propertyName) = loadPelvisDefectDataCurve.savePelvisDefectDataCurve(i).(propertyName);
        end
    end
end
clear loadPelvisDefectDataCurve metaCurveDefect propertiesCurveDefect propertyName

%% Display curvature for defect area (for control)

% Loop with save figure to save the data/figures
%parfor i = 1 : dataCount
    i = 2;  % Pelvis number
    %figure('Visible','off')
    figure
    hold on
    curve = 'cMeanVertexFace'; % Choose 'normalVertexFace', 'normalFace', 'cMeanVertexFace', 'cMeanFace', or 'gauss'
    
    switch curve
        case 'normalVertexFace'
            curveType = 'normal';
            rgbType = 'RGBnormVertexFace';
            rangeType = 'topCurveVertexIdx';
        case 'normalFace'
            curveType = 'normal';
            rgbType = 'RGBnormFace';
            rangeType = 'topCurveFaceIdx'; 
            
        case 'cMeanVertexFace'
            curveType = 'cMean';
            rgbType = 'RGBnormVertexFace'; 
            rangeType = 'topCurveVertexIdx'; 
        case 'cMeanFace'
            curveType = 'cMean';
            rgbType = 'RGBnormFace';
            rangeType = 'topCurveFaceIdx'; 
            
        case 'gauss'
            curveType = 'gauss';
            rgbType = 'RGBnormVertexFace';
            rangeType = 'topCurveVertexIdx';
    end
    
    % Fetch the appropriate RGB data and range indices
    rgb = pelvisDefect(i).curve.(curveType).(rgbType);
    rangeIdx = pelvisDefect(i).curve.(curveType).(rangeType);
    
    % Plot the mesh with curvature color mapping
    patch('Faces', pelvisDefect(i).import.processed.faces, ...
          'Vertices', pelvisDefect(i).import.processed.vertices, ...
          'FaceVertexCData', rgb, ...
          'FaceColor', 'flat', ...
          'FaceAlpha', 1, ...
          'EdgeColor', TUMcolors.grey50, ...
          'EdgeAlpha', 0.25);
    % Pelvis 3Dview (lightning)
    % patch('Faces',pelvisDefect(i).import.processed.faces,...
    %     'Vertices',pelvisDefect(i).import.processed.vertices,...
    %     'FaceVertexCData', rgb, ...
    %     'FaceColor', 'flat', ...
    %     'FaceAlpha', 1, ...
    %     'EdgeColor','none',...              % Edge color
    %     'EdgeAlpha',0.25,...                % Transparency of the edges
    %     ... % Ligthing for 3d effect
    %     'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
    %     'AmbientStrength', 0.5);
    % light('Position', [1 1 5], 'Style', 'infinite');
    
    % Display acetabulum centre
    plot3(pelvisDefect(i).import.processed.acentre(1), ...
          pelvisDefect(i).import.processed.acentre(2), ...
          pelvisDefect(i).import.processed.acentre(3), ...
          '.', 'MarkerEdgeColor', 'm', 'MarkerSize', 50); %TUMcolors.blue300
    % Display curvature range (optional)
    % plot3(pelvisDefect(i).import.processed.vertices(rangeIdx, 1), ...
    %       pelvisDefect(i).import.processed.vertices(rangeIdx, 2), ...
    %       pelvisDefect(i).import.processed.vertices(rangeIdx, 3), ...
    %       '*', 'MarkerEdgeColor', TUMcolors.orange, 'MarkerSize', 10);
    
    % Legend: colormap
    values = linspace(0, 1, 256);
    colormap('viridis');
    c = colorbar('Ticks', [0, 1], 'TickLabels', {'0', '1'});
    ylabel(c, 'Normalised Curvature', 'Rotation', 90, 'HorizontalAlignment', 'center');
    set(c, 'Position', [0.8, 0.1, 0.05, 0.8]); % [left, bottom, width, height]
    % Set view and axis properties
    title(['Curvature ',curve,': Pelvis ' num2str(i)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z');
    %legend(l1,'Highest Curve Values'); % for topCurveRange
    daspect([1, 1, 1]);
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')',curve,'.fig'])

%end

clear curve curveType rgbType rangeType rgb rangeIdx values c

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Class Transform %%%%%%%%%%%%

% User-defined method (see class transform of pelvis)
% Options: 'boundBoxHull', 'boundBoxSVD', 'boundSphere', 'distanceLM', 'distanceCentreLM'
% methodBound = 'distanceLM'; 
boundMethod = 'landmarks';
boundType = 'all';

% Reference: first pelvis data
% Reference: acetabulum centre (first defect data)
refPoint = pelvisDefect(1).import.processed.acentre;
% Reference landmarks 
allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', ...
    'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'}; 
% allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre'};
% Initialize refLandmarks with NaNs
refLandmarks = NaN(length(allLandmarks), 3);
% Extract available landmarks
for i = 1:length(allLandmarks)
    if isfield(pelvisDefect(1).import.processed, allLandmarks{i})
        refLandmarks(i,:) = pelvisDefect(1).import.processed.(allLandmarks{i});
    end
end

% Transformation: Translation/Rotation (SVD/Kabsch)
% User-defined weighting (depends on the number of landmarks) %%%
wKabsch{1} = [1,1,1,1,1,1];
wKabsch{2} = [1,1,1,1,1,1000000000];
wKabsch{3} = [1,1,1,1,1,1,1,1,1,1,1,1,1,1];
wKabsch{4} = [1,1,1,1,1,1000000000,1,1,1,1,1,1,1,1];
wKabsch{5} = [1,1,1,1,1000,1000000000,1000,1000,1000,1,1,1,1000,1]; %%%
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};

parfor i = 1:dataCountDefect 
    
    % Scale
    pelvisDefect(i).transform = pelvisDefect(i).transform.scale(i, pelvisDefect(i).import.processed, ...
        pelvis(i).boundaries.(boundMethod).(boundType).scale); % Scaling method defined in class boundary (see section)

    % Translation to refPoint (acentre) %%% With SCALED data
    % Attention: Scaling leads to new position of the geometry
    pelvisDefect(i).transform = pelvisDefect(i).transform.translation(i,refPoint, ...
        pelvisDefect(i).transform.scaled);

    % Transformation: Translation/Rotation (SVD/Kabsch) %%% With SCALED and TRANSLATED data
    % for w = 1:length(weightKabsch)
    %     pelvisDefect(i).transform = pelvisDefect(i).transform.trafoKabsch(i, refLandmarks, ...
    %         wKabsch{w}, weightKabsch{w}, pelvisDefect(i).transform.trans);
    % end
    
    % Transformation: Translation/Rotation (SVD/Kabsch) %%% With SCALED data   
    for w = 1:length(weightKabsch)
        pelvisDefect(i).transform = pelvisDefect(i).transform.trafoKabsch(i, refLandmarks, ...
            wKabsch{w}, weightKabsch{w},pelvisDefect(i).transform.scaled);
    end
    
end

% Least Root Mean Square of transformation
for i = 1:length(weightKabsch) % Weighting
    allDefect.trafo.(weightKabsch{i}).lrms = NaN;
    for j = 2:dataCount % without reference pelvis    
        allDefect.trafo.(weightKabsch{i}).lrms(j,1) = pelvisDefect(j).transform.trafo.(weightKabsch{i}).lrms;
    end
    % LRMS across all pelvises without reference pelvis
    allDefect.trafo.(weightKabsch{i}).lrmsMean = mean(allDefect.trafo.(weightKabsch{i}).lrms(2:dataCount));
    allDefect.trafo.(weightKabsch{i}).lrmsStd = std(allDefect.trafo.(weightKabsch{i}).lrms(2:dataCount));
end

clear methodBound boundMethod boundType refPoint refLandmarks wKabsch

%% Evaluate transformation: landmark error

% Define all landmarks
allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', ...
    'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};

% Collect errors for each landmark from each pelvis
for i = 1:length(weightKabsch) % Weighting
    for j = 1:length(allLandmarks)
        landmark = allLandmarks{j};
        allDefect.trafo.(weightKabsch{i}).errorLM.(landmark) = NaN; % Initialize
        for k = 2:numel(pelvis) % First pelvis: reference pelvis
            if isfield(pelvisDefect(k).transform.trafo.(weightKabsch{i}), ['error_' landmark])
                allDefect.trafo.(weightKabsch{i}).errorLM.(landmark)(k,1) = ...
                    pelvisDefect(k).transform.trafo.(weightKabsch{i}).(['error_' landmark]);
            end
        end
    end
end

% Calculate mean and standard deviation
selectedLandmarks = [allPelvis.scale.landmarks.centre.closestLM, 'acentre'];
for i = 1:length(weightKabsch)
    for j = 1:length(allLandmarks)
        landmark = allLandmarks{j};
        % Calculate mean and standard deviation per landmark
        allDefect.trafo.(weightKabsch{i}).errorLMmean.(landmark) = ...
            mean(allDefect.trafo.(weightKabsch{i}).errorLM.(landmark), 'omitnan');
        allDefect.trafo.(weightKabsch{i}).errorLMstd.(landmark) = ...
            std(allDefect.trafo.(weightKabsch{i}).errorLM.(landmark), 'omitnan');
    end
    % Calculate mean and standard deviation across all landmarks
    allDefect.trafo.(weightKabsch{i}).errorLMmeanAll = ... 
        mean(cell2mat(struct2cell(allDefect.trafo.(weightKabsch{i}).errorLMmean)), 'omitnan');
    allDefect.trafo.(weightKabsch{i}).errorLMstdAll = ...
        std(cell2mat(struct2cell(allDefect.trafo.(weightKabsch{i}).errorLMmean)), 'omitnan');
    % Calculate mean and standard deviation for the closest landmarks to acentre (selected landmarks)
    allDefect.trafo.(weightKabsch{i}).errorLMmeanSelected = [];
    allDefect.trafo.(weightKabsch{i}).LMselected = selectedLandmarks;
    for j = 1:length(selectedLandmarks)
        lm = selectedLandmarks{j};
        allDefect.trafo.(weightKabsch{i}).errorLMmeanSelected(end+1,1) = allDefect.trafo.(weightKabsch{i}).errorLMmean.(lm); 
    end
    % Calculate mean and standard deviation across selected landmarks
    allDefect.trafo.(weightKabsch{i}).errorLMmeanAllselected = mean(allDefect.trafo.(weightKabsch{i}).errorLMmeanSelected);
    allDefect.trafo.(weightKabsch{i}).errorLMstdAllselected = std(allDefect.trafo.(weightKabsch{i}).errorLMmeanSelected);
end

clear selectedLandmarks 

%% Evaluate transformation: Vertex-to-Nearest-Neighbour Distance

refVertices = pelvisDefect(1).import.processed.vertices;
for i = 1:length(weightKabsch) % Weighting
    allDefect.trafo.(weightKabsch{i}).vertexError = NaN;
    allDefect.trafo.(weightKabsch{i}).vertexErrorMax = NaN;
    for j = 2:dataCountDefect % without reference pelvis
        % First: Kabsch transformation of all pelvis vertices (in function trafoKabsch)
        % Vertex-to-Nearest-Neighbour distance
        pelvisDefect(j).transform = pelvisDefect(j).transform.vertexNeighbour(j, refVertices, weightKabsch{i}, ...
            pelvisDefect(j).transform.trafo.(weightKabsch{i}).vertices,pelvisDefect(j).transform.trafo.(weightKabsch{i}).faces);
        % Store Vertex-to-Nearest-Neighbour distance
        allDefect.trafo.(weightKabsch{i}).vertexError(j,1) = pelvisDefect(j).transform.trafo.(weightKabsch{i}).nearVertexMean;
        allDefect.trafo.(weightKabsch{i}).vertexErrorMax(j,1) = pelvisDefect(j).transform.trafo.(weightKabsch{i}).nearVertexMax;
    end
    % Error across all pelvises without reference pelvis
    allDefect.trafo.(weightKabsch{i}).vertexErrorMean = mean(allDefect.trafo.(weightKabsch{i}).vertexError(2:end));
    allDefect.trafo.(weightKabsch{i}).vertexErrorStd = std(allDefect.trafo.(weightKabsch{i}).vertexError(2:end));
    allDefect.trafo.(weightKabsch{i}).vertexErrorMaxMean = mean(allDefect.trafo.(weightKabsch{i}).vertexErrorMax(2:end));
    allDefect.trafo.(weightKabsch{i}).vertexErrorMaxStd = std(allDefect.trafo.(weightKabsch{i}).vertexErrorMax(2:end));
end

%% Evaluate transformation: hounsdorff distance

refVertices = pelvisDefect(1).import.processed.vertices;
for i = 1:length(weightKabsch) % Weighting
    allDefect.trafo.(weightKabsch{i}).hausdorffDist = NaN;
    for j = 2:dataCount % without reference pelvis        
        % Hausdorff-Distance
        pelvisDefect(j).transform = pelvisDefect(j).transform.hausdorffDistance(j, refVertices, weightKabsch{i}, ...
            pelvisDefect(j).transform.trafo.(weightKabsch{i}).vertices);
        % Store Hausdorff-Distance
        allDefect.trafo.(weightKabsch{i}).hausdorffDist(j,1) = pelvisDefect(j).transform.trafo.(weightKabsch{i}).hausdorffDist;
    end
    % Distance across all pelvises without reference pelvis
    allDefect.trafo.(weightKabsch{i}).hausdorffMean = mean(allDefect.trafo.(weightKabsch{i}).hausdorffDist(2:dataCount));
    allDefect.trafo.(weightKabsch{i}).hausdorffStd = std(allDefect.trafo.(weightKabsch{i}).hausdorffDist(2:dataCount));
end

clear refVertices

%% Save and load properties of Class Transform

% Save class transform data (properties)
savePelvisDefectDataTransform = struct();
% Meta information of curve
metaTransformDefect = metaclass(pelvisDefect(1).transform);
propertiesTransformDefect = {metaTransformDefect.PropertyList.Name};
for i = 1:dataCountDefect
    for j = 1:length(propertiesTransformDefect)
        propertyName = propertiesTransformDefect{j};
        % Cache
        savePelvisDefectDataTransform(i).(propertyName) = pelvisDefect(i).transform.(propertyName);
    end
end
% Save data
save('.\pelvisDefectDataTransform.mat', 'savePelvisDefectDataTransform', '-v7.3'); % Adapt storage location
clear savePelvisDefectDataTransform metaTransformDefect propertiesTransformDefect propertyName

% Clear class transform (selected)
propertiesTransformDefectToClear = {'scaled', 'trans'}; % Adapt
for i = 1:dataCountDefect
    % Empty fields in pelvisDefect(i).transform
    for j = 1:length(propertiesTransformDefectToClear)
        propertyName = propertiesTransformDefectToClear{j};
        fieldType = class(pelvisDefect(i).transform.(propertyName));
        switch fieldType
            case 'double'
                pelvisDefect(i).transform.(propertyName) = []; % Empty double
            case 'struct'
                pelvisDefect(i).transform.(propertyName) = struct(); % Empty struct
            case 'cell'
                pelvisDefect(i).transform.(propertyName) = {}; % Empty cell
            otherwise
                pelvisDefect(i).transform.(propertyName) = []; % Default to empty
        end
    end
end

% Clear class transform.trafo (selected)
propertiesTransformDefectToClear = {'w1', 'w2', 'w3', 'w4'}; % Adapt
for i = 1:dataCount
    % Empty fields in pelvis(i).transform.trafo
    for j = 1:length(propertiesTransformDefectToClear)
        propertyName = propertiesTransformDefectToClear{j};
        fieldType = class(pelvisDefect(i).transform.trafo.(propertyName));
        switch fieldType
            case 'double'
                pelvisDefect(i).transform.trafo.(propertyName) = []; % Empty double
            case 'struct'
                pelvisDefect(i).transform.trafo.(propertyName) = struct(); % Empty struct
            case 'cell'
                pelvisDefect(i).transform.trafo.(propertyName) = {}; % Empty cell
            otherwise
                pelvisDefect(i).transform.trafo.(propertyName) = []; % Default to empty
        end
    end
end
clear propertiesTransformDefectToClear propertyName fieldType

% Load class transform data (properties)
loadPelvisDefectDataTransform = load('.\pelvisDefectDataTransformReduced.mat', 'savePelvisDefectDataTransform'); % Adapt storage location
%pelvisDefect(dataCountDefect) = Defect; % Initialisation
metaTransformDefect = metaclass(pelvisDefect(1).transform);
%propertiesTransformDefect = {metaTransformDefect.PropertyList.Name};
propertiesTransformDefect = {'trafo'}; % load selected properties
for i = 1:dataCountDefect
    for j = 1:length(propertiesTransformDefect)
        propertyName = propertiesTransformDefect{j};
        if isfield(loadPelvisDefectDataTransform.savePelvisDefectDataTransform(i), propertyName)
            pelvisDefect(i).transform.(propertyName) = ...
                loadPelvisDefectDataTransform.savePelvisDefectDataTransform(i).(propertyName);
        end
    end
end
clear loadPelvisDefectDataTransform metaTransformDefect propertiesTransformDefect propertyName

%% Display scaled defect (for control)

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 2;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on

    % Pelvis defect original (optional)
    % patch('Faces',pelvisDefect(i).import.processed.faces,...
    %     'Vertices',pelvisDefect(i).import.processed.vertices,...
    %     'FaceColor',TUMcolors.grey20, ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect original 3Dview (lightning)
    patch('Faces',pelvisDefect(i).import.processed.faces,...
        'Vertices',pelvisDefect(i).import.processed.vertices,...
        'FaceColor',TUMcolors.grey20, ...   % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Pelvis defect scaled 
    % patch('Faces',pelvisDefect(i).transform.scaled.faces,...
    %     'Vertices',pelvisDefect(i).transform.scaled.vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect scaled 3Dview (lightning)
    patch('Faces',pelvisDefect(i).transform.scaled.faces,...
        'Vertices',pelvisDefect(i).transform.scaled.vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Format and display properties
    title(['Pelvis defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    legend('Pelvis Original', 'Pelvis Scaled');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')scaled.fig'])
%end

%% Display scaled and shifted defects (for control)

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 2;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    handles = [];
    labels = {};
    
    % Pelvis defect scaled and shifted 
    % patch('Faces',pelvisDefect(i).transform.trans.faces,...
    %     'Vertices',pelvisDefect(i).transform.trans.vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect scaled and shifted 3Dview (lightning)
    patch('Faces',pelvisDefect(i).transform.trans.faces,...
        'Vertices',pelvisDefect(i).transform.trans.vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Landmarks
    for idx = 1:length(allLandmarks)
        landmark = allLandmarks{idx};
        % Check if the landmark is available for the transformed pelvis
        if isfield(pelvisDefect(i).transform.trans, landmark)
            coords = pelvisDefect(i).transform.trans.(landmark);
            % Plot Landmarks
            h1 = plot3(coords(1), coords(2), coords(3), '*', 'Color', TUMcolors.orange, 'MarkerSize', 10);
            handles(1) = h1;
            labels{1} = 'Transformed Landmark';
        end
        % Check if the landmark is available for the reference
        if isfield(pelvisDefect(1).import.processed, landmark)
            ref_coords = pelvisDefect(1).import.processed.(landmark);
            % Plot Reference Landmarks
            h2 = plot3(ref_coords(1), ref_coords(2), ref_coords(3), '*', 'Color', TUMcolors.blue300, 'MarkerSize', 10);
            handles(2) = h2;
            labels{2} = 'Reference Landmark';
        end
    end
    
    % Display acetabulum centre (reference)
    h3 =  plot3(pelvisDefect(1).import.processed.acentre(1), pelvisDefect(1).import.processed.acentre(2), ...
        pelvisDefect(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10); 
    handles(end+1) = h3;
    labels{end+1} = 'Acetabulum Centre';
    
    % Format and display properties
    title(['Pelvis defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;
    legend(handles, labels);

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')trans.fig'])
%end
clear h1 h2 h3 handles labels landmark idx

%% Display scaled and transformed defects (for control)

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 2;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    handles = [];
    labels = {};
    
    % Pelvis defect scaled and transformed
    % patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect scaled and transformed 3Dview (lightning)
    patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');
    
    % Landmarks
    for idx = 1:length(allLandmarks)
        landmark = allLandmarks{idx};
        % Check if the landmark is available for the transformed pelvis
        if isfield(pelvisDefect(i).transform.trafo.(weightKabsch{w}), landmark)
            coords = pelvisDefect(i).transform.trafo.(weightKabsch{w}).(landmark);
            % Plot Landmarks
            h1 = plot3(coords(1), coords(2), coords(3), '*', 'Color', TUMcolors.orange, 'MarkerSize', 10);
            handles(1) = h1;
            labels{1} = 'Transformed Landmark';
        end
        % Check if the landmark is available for the reference
        if isfield(pelvisDefect(1).import.processed, landmark)
            ref_coords = pelvisDefect(1).import.processed.(landmark);
            % Plot Reference Landmarks
            h2 = plot3(ref_coords(1), ref_coords(2), ref_coords(3), '*', 'Color', TUMcolors.blue300, 'MarkerSize', 10);
            handles(2) = h2;
            labels{2} = 'Reference Landmark';
        end
    end 
    
    % Reference acetabulum centre
    h3 = plot3(pelvisDefect(1).import.processed.acentre(1), pelvisDefect(1).import.processed.acentre(2), ...
        pelvisDefect(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10); 
    handles(end+1) = h3;
    labels{end+1} = 'Acetabulum Centre';

    % Format and display properties
    title(['Pelvis defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    legend(handles, labels);
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')trafo.fig'])
%end

clear h1 h2 h3 handles labels landmark idx

%% Display all scaled and transformed defects (for control)

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

%figure('Visible','off')
figure
hold on

% Reference pelvis
patch('Faces',pelvis(1).import.processed.faces,...
    'Vertices',pelvis(1).import.processed.vertices,...
    'FaceColor',TUMcolors.grey20, ...    % Face color
    'FaceAlpha',0.5,...                   % Transparency of the faces
    'EdgeColor','none',...              % Edge color
    'EdgeAlpha',0.25,...                % Transparency of the edges
    ... % Ligthing for 3d effect
    'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
    'AmbientStrength', 0.5);
light('Position', [1 1 5], 'Style', 'infinite');

% Predefined color lists
patchColors = lines(dataCountDefect); % Using the MATLAB 'lines' colormap
markerColors = jet(dataCountDefect);  % Using the MATLAB 'jet' colormap

for i = 1:dataCountDefect

    % Plot with changing FaceColor
    % Pelvis defect scaled and transformed
    % patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
    %     'FaceColor',patchColors(i,:), ...   % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect scaled and transformed 3Dview (lightning)
    patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',patchColors(i,:), ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    currentMarkerColor = markerColors(i,:); % Single marker color for all landmarks in this iteration

    % Plot each available landmark with the currentMarkerColor
    for idx = 1:length(allLandmarks)
        landmark = allLandmarks{idx};
        % Check if the landmark is available for the transformed pelvis
        if isfield(pelvisDefect(i).transform.trafo.(weightKabsch{w}), landmark)
            coords = pelvisDefect(i).transform.trafo.(weightKabsch{w}).(landmark);
            % Plot Landmarks
            plot3(coords(1), coords(2), coords(3), '*', 'Color', currentMarkerColor, 'MarkerSize', 10);
        end
    end
end

% % without color change
% for i = 1:dataCount
%     % Pelvis defect scaled and transformed
%     % patch('Faces',pelvisdefect(i).transform.trafo.(weightKabsch{w}).faces,...
%     %     'Vertices',pelvisdefect(i).transform.trafo.(weightKabsch{w}).vertices,...
%     %     'FaceColor',[0.9 0.75 0.68], ...   % Face color
%     %     'FaceAlpha',1,...                   % Transparency of the faces
%     %     'EdgeColor',TUMcolors.grey50,...    % Edge color
%     %     'EdgeAlpha',0.25);                  % Transparency of the edges
%     % Pelvis defect scaled and transformed 3Dview (lightning)
%     patch('Faces',pelvisdefect(i).transform.trafo.(weightKabsch{w}).faces,...
%         'Vertices',pelvisdefect(i).transform.trafo.(weightKabsch{w}).vertices,...
%         'FaceColor',[0.9 0.75 0.68], ...   % Face color
%         'FaceAlpha',1,...                   % Transparency of the faces
%         'EdgeColor','none',...              % Edge color
%         'EdgeAlpha',0.25,...                % Transparency of the edges
%         ... % Ligthing for 3d effect
%         'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
%         'AmbientStrength', 0.5);
%     light('Position', [1 1 5], 'Style', 'infinite');
%
%     for idx = 1:length(allLandmarks)
%         landmark = allLandmarks{idx};
%         % Check if the landmark is available for the transformed pelvis
%         if isfield(pelvisdefect(i).transform.trafo.(weightKabsch{w}), landmark)
%             coords = pelvisdefect(i).transform.trafo.(weightKabsch{w}).(landmark);
%             % Plot Landmarks
%             plot3(coords(1), coords(2), coords(3), '*', 'Color', TUMcolors.orange, 'MarkerSize', 10);
%         end
%     end
% end

for idx = 1:length(allLandmarks)
    landmark = allLandmarks{idx};
    % Check if the landmark is available for the reference
    if isfield(pelvisDefect(1).import.processed, landmark)
        ref_coords = pelvisDefect(1).import.processed.(landmark);
        % Plot Reference Landmarks
        plot3(ref_coords(1), ref_coords(2), ref_coords(3), '*', 'Color', TUMcolors.blue300, 'MarkerSize', 10);
    end
end

% Reference acetabulum centre
plot3(pelvisDefect(1).import.processed.acentre(1), pelvisDefect(1).import.processed.acentre(2), ...
    pelvisDefect(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10);

% Format and display properties
title(['Pelvis defect ' num2str(i)]);
xlabel('X'); ylabel('Y'); zlabel('Z');
daspect([1, 1, 1]); % Equal aspect ratio for the axes
view(3);
hold off;

% Save figure (figure unvisible, but saved visible)
%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%savefig('./Figures/PelvisDefectTrafo.fig')

clear patchColors markerColors currentMarkerColor landmark rec_coords idx

%% Display scaled and transformed defect with reference pelvis (for control)

% Reference landmarks 
allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', ...
    'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'}; 
% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weightingw

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount
    i = 16;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    handles = [];
    labels = {};
        
    % Pelvis defect scaled and transformed
    % patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...   % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect scaled and transformed 3Dview (lightning)
    patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');
    % Reference pelvis
    patch('Faces',pelvis(1).import.processed.faces,...
        'Vertices',pelvis(1).import.processed.vertices,...
        'FaceColor',TUMcolors.grey20, ...    % Face color
        'FaceAlpha',0.5,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Landmarks
    for idx = 1:length(allLandmarks)
        landmark = allLandmarks{idx};
        % Check if the landmark is available for the transformed pelvis
        if isfield(pelvisDefect(i).transform.trafo.(weightKabsch{w}), landmark)
            coords = pelvisDefect(i).transform.trafo.(weightKabsch{w}).(landmark);
            % Plot Landmarks
            h1 = plot3(coords(1), coords(2), coords(3), '*', 'Color', TUMcolors.orange, 'MarkerSize', 10);
            handles(1) = h1;
            labels{1} = 'Transformed Landmark';
        end
        % Check if the landmark is available for the reference
        if isfield(pelvisDefect(1).import.processed, landmark)
            ref_coords = pelvisDefect(1).import.processed.(landmark);
            % Plot Reference Landmarks
            h2 = plot3(ref_coords(1), ref_coords(2), ref_coords(3), '*', 'Color', TUMcolors.blue300, 'MarkerSize', 10);
            handles(2) = h2;
            labels{2} = 'Reference Landmark';
        end
    end 
    
    % Reference acetabulum centre 
    h3 = plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
        pelvis(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10); 
    handles(end+1) = h3;
    labels{end+1} = 'Acetabulum Centre';
    
    % Format and display properties
    title(['Pelvis defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    xlim([-50 150]); % for same scaling in different plots
    ylim([-20 150]);
    zlim([-250 20]);
    %view(3);
    view(80,-12.5) % View for images
    legend(handles, labels);
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')RefTrafo.fig'])
%end

clear h1 h2 h3 handles labels ref_coords coords landmark idx

%% Display scaled and transformed pelvis+defect with reference pelvis (for control)

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting
% Landmarks
allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', ...
    'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'}; 

% Loop with save figure to save the data/figures
%parfor i = 1:dataCount
    i = 1;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    handles = [];
    labels = {};

    % Pelvis scaled and transformed
    % patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...   % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis scaled and transformed 3Dview (lightning)
    patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Pelvis defect scaled and transformed
    % patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
    %     'FaceColor',TUMcolors.green, ...   % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect scaled and transformed 3Dview (lightning)
    patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',TUMcolors.green, ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Reference pelvis
    patch('Faces',pelvis(1).import.processed.faces,...
        'Vertices',pelvis(1).import.processed.vertices,...
        'FaceColor',TUMcolors.grey20, ...    % Face color
        'FaceAlpha',0.5,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Landmarks
    for idx = 1:length(allLandmarks)
        landmark = allLandmarks{idx};
        % Check if the landmark is available for the transformed pelvis
        if isfield(pelvisDefect(i).transform.trafo.(weightKabsch{w}), landmark)
            coords = pelvisDefect(i).transform.trafo.(weightKabsch{w}).(landmark);
            % Plot Landmarks
            h1 = plot3(coords(1), coords(2), coords(3), '*', 'Color', TUMcolors.orange, 'MarkerSize', 10);
            handles(1) = h1;
            labels{1} = 'Transformed Landmark';
        end
        % Check if the landmark is available for the reference
        if isfield(pelvisDefect(1).import.processed, landmark)
            ref_coords = pelvisDefect(1).import.processed.(landmark);
            % Plot Reference Landmarks
            h2 = plot3(ref_coords(1), ref_coords(2), ref_coords(3), '*', 'Color', TUMcolors.blue300, 'MarkerSize', 10);
            handles(2) = h2;
            labels{2} = 'Reference Landmark';
        end
    end 
    
    % Reference acetabulum centre
    h3 = plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
        pelvis(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10); 
    handles(end+1) = h3;
    labels{end+1} = 'Acetabulum Centre';
    
    % Format and display properties
    title(['Pelvis ' num2str(i)]);
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    xlim([-50 150]); % for same scaling in different plots
    ylim([-20 150]);
    zlim([-250 20]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    %view(3);
    view(80,-12.5) % View for images
    legend(handles, labels);
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')RefTrafo(2).fig'])
%end

clear h1 h2 h3 handles labels ref_coords coords landmark idx

%% Display curvature for the whole pelvis (scaled&transformed) with defect area and topCurveRange (for control)

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%for i = 1:dataCountDefect
    i = 1;  % Pelvis number
    %figure('Visible','off')  
    figure
    hold on
    curve = 'cMeanVertexFace'; % Choose 'normalVertexFace', 'normalFace', 'cMeanVertexFace', 'cMeanFace', or 'gauss'
    
    switch curve
        case 'normalVertexFace'
            curveType = 'normal';
            rgbType = 'RGBnormVertexFace';
            rangeType = 'topCurveVertexIdx';
        case 'normalFace'
            curveType = 'normal';
            rgbType = 'RGBnormFace';
            rangeType = 'topCurveFaceIdx'; 
            
        case 'cMeanVertexFace'
            curveType = 'cMean';
            rgbType = 'RGBnormVertexFace'; 
            rangeType = 'topCurveVertexIdx'; 
        case 'cMeanFace'
            curveType = 'cMean';
            rgbType = 'RGBnormFace';
            rangeType = 'topCurveFaceIdx'; 
            
        case 'gauss'
            curveType = 'gauss';
            rgbType = 'RGBnormVertexFace';
            rangeType = 'topCurveVertexIdx';
    end
    
    % Fetch the appropriate RGB data and range indices
    rgb = pelvis(i).curve.(curveType).(rgbType);
    rangeIdx = pelvis(i).curve.(curveType).(rangeType);
    
    % Pelvis
    % patch('Faces', pelvis(i).transform.trafo.(weightKabsch{w}).faces, ...
    %       'Vertices', pelvis(i).transform.trafo.(weightKabsch{w}).vertices, ...
    %       'FaceVertexCData', rgb, ...
    %       'FaceColor', 'flat', ...
    %       'FaceAlpha', 1, ...
    %       'EdgeColor', 'none', ...
    %       'EdgeAlpha', 0.25);
    % Pelvis 3Dview (lightning)
    patch('Faces',pelvis(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvis(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceVertexCData', rgb, ...
        'FaceColor', 'flat', ...
        'FaceAlpha', 1, ...
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 1);
    light('Position', [1 1 -5], 'Style', 'infinite');
    material dull; % Alternativen: shiny, dull, metal

    % Defect area
    patch('Faces', pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces, ...
          'Vertices', pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices, ...
          'FaceVertexCData', rgb, ...
          'FaceColor', 'flat', ...
          'FaceAlpha', 1, ...
          'EdgeColor', TUMcolors.grey50, ...
          'EdgeAlpha', 0.25);
    
    % Display acetabulum centre
    plot3(pelvis(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
          pelvis(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
          pelvis(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
          '.', 'MarkerEdgeColor', 'm', 'MarkerSize', 45); %TUMcolors.blue300
    plot3(pelvis(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
          pelvis(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
          pelvis(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
          'x', 'MarkerEdgeColor', 'm', 'MarkerSize', 22, 'LineWidth', 5);
    % % Display curvature range (optional): topCurveRange
    % l1 = plot3(pelvis(i).transform.trafo.(weightKabsch{w}).vertices(rangeIdx, 1), ...
    %       pelvis(i).transform.trafo.(weightKabsch{w}).vertices(rangeIdx, 2), ...
    %       pelvis(i).transform.trafo.(weightKabsch{w}).vertices(rangeIdx, 3), ...
    %       '.', 'MarkerEdgeColor', TUMcolors.orange, 'MarkerSize', 5);
    % Display defect area boundary
    l2 = plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices(pelvisDefect(i).curve.boundVerticesNr, 1), ...
          pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices(pelvisDefect(i).curve.boundVerticesNr, 2), ...
          pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices(pelvisDefect(i).curve.boundVerticesNr, 3), ...
          '.', 'MarkerEdgeColor', 'r', 'MarkerSize', 20);  %TUMcolors.grey20

    % Legend: colormap
    colormap('viridis');
    c = colorbar('Ticks', [0, 1], 'TickLabels', {'0', '1'});
    ylabel(c, 'Normalised Curvature', 'Rotation', 90, 'HorizontalAlignment', 'center');
    set(c, 'Position', [0.8, 0.1, 0.02, 0.8]); % [left, bottom, width, height]
    % Set view and axis properties
    title(['Curvature ',curve,': Pelvis ' num2str(i)]); 
    %legend([l1, l2], 'Highest Curve Values', 'Pelvis Defect'); % for topCurveRange
    %legend(l2, 'Pelvis Defect');
    xlim([-50 150]); % for same scaling in different plots
    ylim([-20 150]);
    zlim([-250 20]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]);
    %view(3);
    view(80,-12.5) % View for images
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')',curve,'CurveRange.fig'])

%end

clear rgb rangeIdx pelvisDefectBound rgbType curve curveType rangeType l1 l2 c values

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%% Class Volume %%%%%%%%%%%%%%%%
% Class Defect: Fill reference pelvis with points (for volume/intersection)

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% For reference model filled with points
% Boundary box / cuboid around reference pelvis and all defects 
tic
allDefect.stl.verticesAll = pelvis(1).import.processed.vertices;
% For reference model filled with points and only fill the relevant defect area with points (time-saving)
% Boundary box / cuboid around all defects
% allDefect.stl.verticesAll = [];
for i = 1:dataCountDefect  
    allDefect.stl.verticesAll = [allDefect.stl.verticesAll; pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices];
end
% Initalization: Generate object
allDefect.boundaries = Boundary(); 
% Boundary box / cuboid: hull
level = 1;
allDefect.boundaries = allDefect.boundaries.boundBoxHull(0,allDefect.stl.verticesAll,level);
allDefect.time.boundAll = toc;

% Fill reference pelvis with points 
% Attention: check COSY of boundary box
tic
allPelvis.refPoints = Volume(); % Initalization
minDist = 0.5; % Distance for points in boundary box (grid)
allPelvis.refPoints = allPelvis.refPoints.fillPoints(...
    0, pelvis(1).import.processed, minDist, allDefect.boundaries.cuboid.hull);
allPelvis.time.refPointsInside = toc;

% Save data
savePelvisDataRefPointsInside = allPelvis.refPoints.gridPoints;
save('.\PelvisDataRefPointsInside.mat', 'savePelvisDataRefPointsInside', '-v7.3'); % Adapt storage location
clear savePelvisDataRefPointsInside
% Load data
allPelvis.refPoints = Volume(); 
loadPelvisDataRefPointsInside = load('.\pelvisDataRefPointsInside.mat'); % Adapt storage location
allPelvis.refPoints.gridPoints = loadPelvisDataRefPointsInside.savePelvisDataRefPointsInside;
clear loadPelvisDataRefPointsInside

clear minDist level

%% Display reference pelvis with the boundary box (convex hull) and its points inside (for control)

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

%figure('Visible','off')
figure
hold on

% Predefined color lists
patchColors = lines(dataCountDefect); % Using the MATLAB 'lines' colormap
markerColors = jet(dataCountDefect);  % Using the MATLAB 'jet' colormap

% Reference pelvis
patch('Faces',pelvis(1).import.processed.faces,...
    'Vertices',pelvis(1).import.processed.vertices,...
    'FaceColor',TUMcolors.grey20, ...    % Face color
    'FaceAlpha',0.5,...                   % Transparency of the faces
    'EdgeColor','none',...              % Edge color
    'EdgeAlpha',0.25,...                % Transparency of the edges
    ... % Ligthing for 3d effect
    'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
    'AmbientStrength', 0.5);
light('Position', [1 1 5], 'Style', 'infinite');

% Reference acetabulum centre
% plot3(pelvisDefect(1).import.processed.acentre(1), pelvisDefect(1).import.processed.acentre(2), ...
%     pelvisDefect(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10);

% Plot the pelvis defects
for i = 1:dataCountDefect
    patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',patchColors(i,:), ...
        'FaceAlpha',1,...                   % Transparency of the faces
            'EdgeColor','none',...              % Edge color
            'EdgeAlpha',0.25,...                % Transparency of the edges
            ... % Ligthing for 3d effect
            'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
            'AmbientStrength', 0.5);

    currentMarkerColor = markerColors(i,:); % Single marker color
end

% Display bounding box
% Define the 8 corners of the bounding box
trisurf(allDefect.boundaries.cuboid.hull.tri,...
    allDefect.boundaries.cuboid.hull.cornerpoints(:,1),...
    allDefect.boundaries.cuboid.hull.cornerpoints(:,2),...
    allDefect.boundaries.cuboid.hull.cornerpoints(:,3),...
    'FaceColor',TUMcolors.blue300,'EdgeColor',TUMcolors.blue300,'FaceAlpha',0.25);

% Display cuboid edges
for j=1:3
    quiver3(allDefect.boundaries.cuboid.hull.cornerpoints(1,1),...
        allDefect.boundaries.cuboid.hull.cornerpoints(1,2),...
        allDefect.boundaries.cuboid.hull.cornerpoints(1,3),... % start point
        allDefect.boundaries.cuboid.hull.edgeVector(j,1),...
        allDefect.boundaries.cuboid.hull.edgeVector(j,2),...
        allDefect.boundaries.cuboid.hull.edgeVector(j,3),0,'LineWidth',2)
end

% Plot points inside relevant area (max bounding box)
plot3(allPelvis.refPoints.gridPoints.inside(:,1), ...
    allPelvis.refPoints.gridPoints.inside(:,2), ...
    allPelvis.refPoints.gridPoints.inside(:,3), '*','Color',TUMcolors.orange, 'MarkerSize', 2);

% Format and display properties
title('Reference pelvis with defects');
xlabel('X'); ylabel('Y'); zlabel('Z');
daspect([1, 1, 1]); % Equal aspect ratio for the axes
view(3);
hold off;

% Save figure (figure unvisible, but saved visible)
%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%savefig('./Figures/PelvisDefectFillPoints.fig')

clear patchColors markerColors currentMarkerColor

%% Identifiy the edge of the defects

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Edge loops of defect
parfor i = 1:dataCountDefect 
    tic
    % Identify the boundary points at the edge
    pelvisDefect(i).volume = pelvisDefect(i).volume.edging(i,pelvisDefect(i).transform.trafo.(weightKabsch{w})); % Transformed data

    % Identify the edge loops
    % Acetabulum edgeLoop assumption: acetabulum hole with the most vertices (max length) -> first edge loop
    pelvisDefect(i).volume = pelvisDefect(i).volume.edgeLoop(i,pelvisDefect(i).volume.edge);

    % Vertices of the edge loops
    for j = 1:size(pelvisDefect(i).volume.edge.edgeLoops,1)
        pelvisDefect(i).volume.edge.verticesLoops{j,1} = pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices(...
            pelvisDefect(i).volume.edge.edgeLoops{j},:);
    end
    timeEdge(i,1) = toc;
end  

allDefect.time.edge = timeEdge;
allDefect.time.edgeMean = mean(timeEdge);
allDefect.time.edgeStd = std(timeEdge);

clear timeEdge

%% Display scaled, transformed defects with the edge (for control)

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 1;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on

    % Pelvis defect scaled and transformed
    % patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect scaled and transformed 3Dview (lightning)
    patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Edge
    plot3(pelvisDefect(i).volume.edge.vertices(:,1),pelvisDefect(i).volume.edge.vertices(:,2),...
        pelvisDefect(i).volume.edge.vertices(:,3),...
        '*', 'Color', TUMcolors.blue300, 'MarkerSize', 10)
    
    % Format and display properties
    title(['Edge of Pelvis Defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')Edge.fig'])
%end

%% Display scaled, transformed defects with the edge loops (for control)

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 1;  % Pelvis number
    %figure('Visible','off')  
    figure
    hold on
    
    % Use the MATLAB jet colormap
    numLoops = size(pelvisDefect(i).volume.edge.verticesLoops, 1);
    colors = jet(numLoops);  % Get distinct colors for each loop

    % Pelvis defect scaled and transformed
    % patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
    %     'FaceColor',[0.9 0.75 0.68], ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect scaled and transformed 3Dview (lightning)
    patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');
    
    % Edge loops
    for k = 1:numLoops 
        currentColor = colors(k, :);  % Get color from the colormap
        plot3(pelvisDefect(i).volume.edge.verticesLoops{k}(:,1),pelvisDefect(i).volume.edge.verticesLoops{k}(:,2),...
            pelvisDefect(i).volume.edge.verticesLoops{k}(:,3),...
            '*', 'Color', currentColor, 'MarkerSize', 10); 
    end

    % Format and display properties
    title(['Edge of Pelvis Defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
    %savefig(['./Figures/PelvisDefect(',num2str(i),')EdgeLoops.fig'])
%end

clear numLoops colors currentColor

%% Patches - fill holes in pelvis defect mesh

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Mesh edge loops -> patches
% Uses comPoint in function (sorted edge loop points)
% types = 'all', 'acentre, 'acetabulum'
parfor i = 1:dataCountDefect 
    tic
    % Patch all edge loops and refine patches (refined mesh not with acetabulum patch)
    pelvisDefect(i).volume = pelvisDefect(i).volume.meshPatch(i,pelvisDefect(i).transform.trafo.(weightKabsch{w}),'all');
    timePatch(i,1) = toc;

    % Patch of acetabulum edge loop + acentre -> doesn't work very well!!
    % pelvisDefect(i).volume = pelvisDefect(i).volume.meshPatch(i,pelvisDefect(i).transform.trafo.(weightKabsch{w}),...
    %     pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre,'acentre');

    % Patch of acetabulum edge loop + acetabulum edge loop reference -> doesn't work very well!!
    % Acetabulum edge loop reference
    % refEdgeLoop = pelvisDefect(1).transform.trafo.(weightKabsch{w}).vertices(pelvisDefect(1).volume.edge.comPoint{1,1},:);
    % pelvisDefect(i).volume = pelvisDefect(i).volume.meshPatch(i,pelvisDefect(i).transform.trafo.(weightKabsch{w}),...
    %     refEdgeLoop,'acetabulum');
end
allDefect.time.patch = timePatch;
allDefect.time.patchMean = mean(timePatch);
allDefect.time.patchStd = std(timePatch);

% Defects with more than one edge loop (one edge loop = acetabulum loop)
for i = 1:dataCountDefect
    % without acetabulum loop and without patches of more components
    allDefect.volume.refinedAlpha.patches(i,1) = length(pelvisDefect(i).volume.patches.all.vertices) - ...
        pelvisDefect(i).volume.patches.all.numComponents;
end
allDefect.volume.refinedAlpha.patchesIdx = find(allDefect.volume.refinedAlpha.patches > 0); % idx of pelvises

clear timePatch

%% Patches - remesh pelvis defect mesh with patches

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Check orientation of normals (target state: outwards); only for type 'all'
allDefect.time.patchDir(1:dataCountDefect,1) = NaN;
for i = allDefect.volume.refinedAlpha.patchesIdx' % for-Loop because of plot and user input     
    tic
    % number of edge loops; not acetabulum patch
    for k = (pelvisDefect(i).volume.patches.all.numComponents+1):(size(pelvisDefect(i).volume.edge.edgeLoops,1)) 
        % Check inside / outside of meshed and refined patches
        pelvisDefect(i).volume = pelvisDefect(i).volume.flipNormals(i,k,pelvisDefect(i).transform.trafo.(weightKabsch{w}),...
            pelvisDefect(i).volume.patches.all.vertices{k}, pelvisDefect(i).volume.patches.all.faces{k},...
            pelvisDefect(i).volume.patches.all.normals{k}, pelvisDefect(i).volume.patches.all.centreFaces{k},'all');
        % Patch and their normals without refined clearer in plot
    end
    allDefect.time.patchDir(i,1) = toc;
end
allDefect.time.patchDirMean = mean(allDefect.time.patchDir,'omitnan');
allDefect.time.patchDirStd = std(allDefect.time.patchDir,'omitnan');

% Save user-input
savePelvisDefectDataFlipNormals = struct();
for i = allDefect.volume.refinedAlpha.patchesIdx'
    % Cache 
    savePelvisDefectDataFlipNormals(i).flip = pelvisDefect(i).volume.patches.all.flip;
end
% Save data
save('.\pelvisDefectDataFlipNormals.mat', 'savePelvisDefectDataFlipNormals', '-v7.3'); % Adapt storage location
clear savePelvisDefectDataFlipNormals

% Load user-input
loadPelvisDefectDataFlipNormals = load('.\pelvisDefectDataFlipNormals.mat'); % Adapt storage location
savePelvisDefectDataFlipNormals = loadPelvisDefectDataFlipNormals.savePelvisDefectDataFlipNormals;
for i = allDefect.volume.refinedAlpha.patchesIdx'
    pelvisDefect(i).volume.patches.all.flip = savePelvisDefectDataFlipNormals(i).flip;
end
clear savePelvisDefectDataFlipNormals

% Combine meshes (with refined patches)
parfor i = 1:dataCountDefect
    tic
    if ismember(i, allDefect.volume.refinedAlpha.patchesIdx)
        % Combine meshes (defect mesh + refined patches); only for type 'all' and without acetabulum patch
        pelvisDefect(i).volume = pelvisDefect(i).volume.remeshPatch(i,pelvisDefect(i).transform.trafo.(weightKabsch{w}),...
            pelvisDefect(i).volume.patches.all.refinedFaces, pelvisDefect(i).volume.patches.all.refinedVertices, 'all');
    else
        % Defects without patches
        pelvisDefect(i).volume.patches.all.comVertices = pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices;
        pelvisDefect(i).volume.patches.all.comFaces = pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces(:, [3, 2, 1]);
        pelvisDefect(i).volume.patches.all.comNormals = -1 .* pelvisDefect(i).transform.trafo.(weightKabsch{w}).normals; 
        pelvisDefect(i).volume.patches.all.comCentreFaces = pelvisDefect(i).transform.trafo.(weightKabsch{w}).centreFaces;
        disp(['Combined mesh with patches: pelvis ',num2str(i)])
    end
    comVerticesDefect(i,1) = size(pelvisDefect(i).volume.patches.all.comVertices,1);
    comFacesDefect(i,1) = size(pelvisDefect(i).volume.patches.all.comFaces,1);
    timePatchMesh(i,1) = toc;
end
% No change of direction of the normals for defect of reference pelvis with backside 
pelvisDefect(1).volume.patches.all.comNormals = pelvisDefect(1).transform.trafo.(weightKabsch{w}).normals; 
pelvisDefect(1).volume.patches.all.comFaces = pelvisDefect(1).transform.trafo.(weightKabsch{w}).faces;
allDefect.time.patchMesh = timePatchMesh;
allDefect.time.patchMeshMean = mean(timePatchMesh);
allDefect.time.patchMeshStd = std(timePatchMesh);
% stl properties
allDefect.volume.refinedAlpha.comVertices = comVerticesDefect;
allDefect.volume.refinedAlpha.comVerticesMean = mean(comVerticesDefect);
allDefect.volume.refinedAlpha.comVerticesStd = std(comVerticesDefect);
allDefect.volume.refinedAlpha.comVerticesMax = max(comVerticesDefect);
allDefect.volume.refinedAlpha.comVerticesMin = min(comVerticesDefect);
allDefect.volume.refinedAlpha.comFaces = comFacesDefect;
allDefect.volume.refinedAlpha.comFacesMean = mean(comFacesDefect);
allDefect.volume.refinedAlpha.comFacesStd = std(comFacesDefect);
allDefect.volume.refinedAlpha.comFacesMax = max(comFacesDefect);
allDefect.volume.refinedAlpha.comFacesMin = min(comFacesDefect);

% Check stl
for i = 1:dataCountDefect
    TR = triangulation(pelvisDefect(i).volume.patches.all.comFaces, pelvisDefect(i).volume.patches.all.comVertices);
    stlwrite(TR,['.\GeometriesDefectsRemesh\PelvisDefect(',num2str(i),')Remesh.stl']); 
    disp(['write stl: pelvis defect ', num2str(i)]);
end

clear timePatchMesh comVerticesDefect comFacesDefect TR

%% Display scaled, transformed defects with meshed patches - all (for control)

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect 
%for i = allDefect.volume.refinedAlpha.patchesIdx'
    i = 1;  % Pelvis number
    %figure('Visible','off')  
    figure
    hold on
   
    % Use the MATLAB jet colormap
    numLoops = size(pelvisDefect(i).volume.edge.edgeLoops,1);
    numComponents = pelvisDefect(i).volume.patches.all.numComponents;

    colors = jet(numLoops);  % Get distinct colors for each loop

    % Pelvis defect scaled and transformed
    patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',TUMcolors.grey50,...    % Edge color
        'EdgeAlpha',0.25);                  % Transparency of the edges

    % Patches
    if size(pelvisDefect(i).volume.edge.edgeLoops,1) > 1
        for k = (numComponents+1):numLoops % Other patches
            patch('Faces',pelvisDefect(i).volume.patches.all.faces{k,1},... % pelvisDefect(i).volume.patches.all.refinedFaces{k,1}
                'Vertices',pelvisDefect(i).volume.patches.all.vertices{k,1},... % pelvisDefect(i).volume.patches.all.refinedVertices{k,1}
                'FaceColor',[0.9 0.75 0.68], ...    % Face color
                'FaceAlpha',1,...                   % Transparency of the faces
                'EdgeColor',TUMcolors.grey50,...    % Edge color
                'EdgeAlpha',0.25);                  % Transparency of the edges
            currentColor = colors(k, :);  % Get color from the colormap
            plot3(pelvisDefect(i).volume.edge.verticesLoops{k}(:,1),pelvisDefect(i).volume.edge.verticesLoops{k}(:,2),...
                pelvisDefect(i).volume.edge.verticesLoops{k}(:,3),...
                '*', 'Color', currentColor, 'MarkerSize', 10);
        end
    end

    % Display acetabulum centre 
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10);

    % Format and display properties
    title(['Patches of Pelvis Defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
    %savefig(['./Figures/PelvisDefect(',num2str(i),')Patches.fig'])
%end

clear numLoops numComponents colors currentColor

%% Display scaled, transformed defects with meshed and refined patches - combined mesh (for control) 

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
%for i = allDefect.volume.refinedAlpha.patchesIdx'
    i = 1; % Pelvis number
    %figure('Visible','off')  
    figure
    hold on
    
    % Pelvis defect remeshed (scaled and transformed)
    patch('Faces',pelvisDefect(i).volume.patches.all.comFaces,...
        'Vertices',pelvisDefect(i).volume.patches.all.comVertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',TUMcolors.grey50,...    % Edge color
        'EdgeAlpha',0.25);                  % Transparency of the edges
    % Plot the mesh with curvature color mapping (optional)
    % curveType = 'cMean';
    % rgbType = 'RGBnormVertexFace';
    % rgb = pelvisDefect(i).curve.(curveType).(rgbType);
    % patch('Faces', pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces, ...
    %     'Vertices', pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices, ...
    %     'FaceVertexCData', rgb, ...
    %     'FaceColor', 'flat', ...
    %     'FaceAlpha', 1, ...
    %     'EdgeColor', 'none', ...
    %     'EdgeAlpha', 0.1);    

    % Use the MATLAB jet colormap
    numLoops = size(pelvisDefect(i).volume.edge.edgeLoops,1);
    colors = jet(numLoops);  % Get distinct colors for each loop
    if size(pelvisDefect(i).volume.edge.edgeLoops,1) > 1
        for k = (pelvisDefect(i).volume.patches.all.numComponents+1):numLoops % Other patches
            patch('Faces',pelvisDefect(i).volume.patches.all.refinedFaces{k,1},...
                'Vertices',pelvisDefect(i).volume.patches.all.refinedVertices{k,1},...
                'FaceColor',[0.9 0.75 0.68], ...    % Face color [0.9 0.75 0.68][1 0.6 0.2]
                'FaceAlpha',1,...                   % Transparency of the faces
                'EdgeColor','none',...    % Edge color % Figure Paper 'none' TUMcolors.grey50
                'EdgeAlpha',0.25);                  % Transparency of the edges
            hold on
            currentColor = colors(k,:);  % Get color from the colormap
            plot3(pelvisDefect(i).volume.edge.verticesLoops{k}(:,1),pelvisDefect(i).volume.edge.verticesLoops{k}(:,2),...
                pelvisDefect(i).volume.edge.verticesLoops{k}(:,3),...
                '.', 'Color', currentColor, 'MarkerSize', 30); % Figure Paper 'k'
        end
        plot3(pelvisDefect(i).volume.edge.verticesLoops{1}(:,1),pelvisDefect(i).volume.edge.verticesLoops{1}(:,2),...
                pelvisDefect(i).volume.edge.verticesLoops{1}(:,3),...
                '.', 'Color', 'r', 'MarkerSize', 30); 
    end

    % % Normals of the faces (optional)
    % quiver3(pelvisDefect(i).volume.patches.all.comCentreFaces(:,1),...
    %     pelvisDefect(i).volume.patches.all.comCentreFaces(:,2),...
    %     pelvisDefect(i).volume.patches.all.comCentreFaces(:,3),...
    %     pelvisDefect(i).volume.patches.all.comNormals(:,1),...
    %     pelvisDefect(i).volume.patches.all.comNormals(:,2),...
    %     pelvisDefect(i).volume.patches.all.comNormals(:,3),2,...
    %     'Color',TUMcolors.orange)

    % Display acetabulum centre 
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
        '.','Color','m', 'MarkerSize', 80); % pelvis(i).transform.trafo.(weightKabsch{w}).acentre
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
        'x', 'MarkerEdgeColor', 'm', 'MarkerSize', 35, 'LineWidth', 10); % pelvis(i).transform.trafo.(weightKabsch{w}).acentre

    % Format and display properties
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['Patches of Pelvis Defect ' num2str(i)]);
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    view(80, -10) % Figure Paper
    grid off;
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
    %savefig(['./Figures/PelvisDefect(',num2str(i),')PatchesRemesh.fig'])
%end

clear numLoops colors currentColor curveType curve rgb rgbType

%% Simple meshing / triangulation with acentre - generate volume 

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% ShrinkWrap (alphaShape) 
% Alpha radius:     
% Inf, where alphaShape produces the convex hull
% 0, where alphaShape produces an empty alpha shape

% ShrinkWrap (alphaShape) without user interface
% Alpha radius: first time mesh is closed and manifold
shrinkType = 'alpha'; % alpha shape with centre of acetabulum
aRadius = 1; % Start value for alpha radius
parfor  i = 1:dataCountDefect 
    tic
    numData = size(pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,1);
    pelvisDefect(i).volume = pelvisDefect(i).volume.shrinkAlpha(i,shrinkType,aRadius,numData,...
        pelvisDefect(i).transform.trafo.(weightKabsch{w})); 
    usedCount(i,1) = pelvisDefect(i).volume.shrink.alpha.usedDefectVerticesCount;
    volumeAlpha(i,1) = pelvisDefect(i).volume.shrink.alpha.volume;
    alphaRadius(i,1) = pelvisDefect(i).volume.shrink.alpha.radius;
    timeAlpha(i,1) = toc;
end 
allDefect.time.alpha = timeAlpha;
allDefect.time.alphaMean = mean(timeAlpha);
allDefect.time.alphaStd = std(timeAlpha);

% ShrinkWrap (alphaShape) with user interface
shrinkType = 'alpha';
%aRadius = 1; % Start value for alpha radius
timeAlphaUI(1:dataCountDefect,1) = NaN;
for  i = 1:dataCountDefect % for-Loop because of plot and user input  
    tic
    aRadius = pelvisDefect(i).volume.shrink.alpha.radius; % Use alpha radius from previous calculation without user input
    numData = size(pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,1);
    pelvisDefect(i).volume = pelvisDefect(i).volume.shrinkAlphaUI(i,shrinkType,aRadius,numData,...
        pelvisDefect(i).transform.trafo.(weightKabsch{w})); 
    usedCount(i,1) = pelvisDefect(i).volume.shrink.alpha.usedDefectVerticesCount;
    volumeAlpha(i,1) = pelvisDefect(i).volume.shrink.alpha.volume;
    alphaRadius(i,1) = pelvisDefect(i).volume.shrink.alpha.radius;
    timeAlphaUI(i,1) = toc;
end 
allDefect.time.alphaUI = timeAlphaUI;
allDefect.time.alphaUIMean = mean(timeAlphaUI,'omitnan');
allDefect.time.alphaUIStd = std(timeAlphaUI,'omitnan');

% AlphaRadius, used vertices for alphaShape and generated volume
allDefect.volume.alpha.radius = alphaRadius;
allDefect.volume.alpha.radiusMean = mean(alphaRadius);
allDefect.volume.alpha.radiusStd = std(alphaRadius);
allDefect.volume.alpha.radiusMax = max(alphaRadius);
allDefect.volume.alpha.radiusMin = min(alphaRadius);
allDefect.volume.alpha.usedDefectVerticesCount = usedCount;
allDefect.volume.alpha.usedDefectVerticesCountMean = mean(usedCount);
allDefect.volume.alpha.usedDefectVerticesCountStd = std(usedCount);
allDefect.volume.alpha.usedDefectVerticesCountMax = max(usedCount);
allDefect.volume.alpha.usedDefectVerticesCountMin = min(usedCount);
allDefect.volume.alpha.volume = volumeAlpha;
allDefect.volume.alpha.volumeMean = mean(volumeAlpha);
allDefect.volume.alpha.volumeStd = std(volumeAlpha);
allDefect.volume.alpha.volumeMax = max(volumeAlpha);
allDefect.volume.alpha.volumeMin = min(volumeAlpha);

% Define the different Paprosky types and their corresponding variable names
paproskyTypes = {'2a', '2b', '2c', '3a', '3b', 'NaN'};
paproskyNames = {'IIa', 'IIb', 'IIc', 'IIIa', 'IIIb', 'NaN'};
% Initialize a structure to store the volumes for each Paprosky type
allDefect.volume.alpha.paproskyVol = struct();
for i = 1:length(paproskyNames)
    allDefect.volume.alpha.paproskyVol.(paproskyNames{i}) = [];
end

% Collect the volumes for each Paprosky type
for i = 1:dataCountDefect
    paproskyType = pelvisDefect(i).import.processed.paprosky;
    idx = find(strcmp(paproskyTypes, paproskyType));
    if ~isempty(idx)
        varName = paproskyNames{idx};
        allDefect.volume.alpha.paproskyVol.(varName) = [allDefect.volume.alpha.paproskyVol.(varName); ...
            pelvisDefect(i).volume.shrink.alpha.volume];
    end
end
% Calculate the mean volume for each Paprosky type
for i = 1:length(paproskyNames)
    varName = paproskyNames{i};
    allDefect.volume.alpha.paproskyVolMean.(varName) = mean(allDefect.volume.alpha.paproskyVol.(varName), 'omitnan');
    allDefect.volume.alpha.paproskyVolStd.(varName) = std(allDefect.volume.alpha.paproskyVol.(varName), 'omitnan');
    allDefect.volume.alpha.paproskyVolMax.(varName) = max(allDefect.volume.alpha.paproskyVol.(varName));
    allDefect.volume.alpha.paproskyVolMin.(varName) = min(allDefect.volume.alpha.paproskyVol.(varName));
end

clear aRadius numData usedCount volumeAlpha alphaRadius timeAlpha timeAlphaUI idx varName shrinkType

%% Reference pelvis: acetabulum sphere - sphere with Fibonacci spiral (points)

% Sphere parameter
allPelvis.acentre.sphere.radius = 24.107;  % Radius of the sphere
allPelvis.acentre.sphere.centre = pelvisDefect(1).import.processed.acentre;  % centre
allPelvis.acentre.sphere.pointDistance = 0.5;  % Distance between the points of the sphere
% Number of points
allPelvis.acentre.sphere.surfaceArea = 4 * pi * allPelvis.acentre.sphere.radius^2;
allPelvis.acentre.sphere.numPoints = round(allPelvis.acentre.sphere.surfaceArea / (allPelvis.acentre.sphere.pointDistance^2));  
% Fibonacci-Spiral on sphere
allPelvis.acentre.sphere.indices = (0:allPelvis.acentre.sphere.numPoints-1)';
allPelvis.acentre.sphere.phi = pi * (3 - sqrt(5));  % golden angle in radiant
% z-Coordinate and radius of the circles in different heights
allPelvis.acentre.sphere.z = 2 * allPelvis.acentre.sphere.indices / allPelvis.acentre.sphere.numPoints - 1;
allPelvis.acentre.sphere.radCircles = sqrt(1 - allPelvis.acentre.sphere.z.^2);
% Azimuth angle
allPelvis.acentre.sphere.theta = allPelvis.acentre.sphere.phi * allPelvis.acentre.sphere.indices;
% Convert to cartesian coordinates
allPelvis.acentre.sphere.pointsX = allPelvis.acentre.sphere.radius * allPelvis.acentre.sphere.radCircles .* cos(allPelvis.acentre.sphere.theta) + ...
    allPelvis.acentre.sphere.centre(1);
allPelvis.acentre.sphere.pointsY = allPelvis.acentre.sphere.radius * allPelvis.acentre.sphere.radCircles .* sin(allPelvis.acentre.sphere.theta) + ...
    allPelvis.acentre.sphere.centre(2);
allPelvis.acentre.sphere.pointsZ = allPelvis.acentre.sphere.radius * allPelvis.acentre.sphere.z + allPelvis.acentre.sphere.centre(3);
allPelvis.acentre.sphere.points = [allPelvis.acentre.sphere.pointsX, allPelvis.acentre.sphere.pointsY, allPelvis.acentre.sphere.pointsZ];  

% Distance to next neighbour with knnsearch
% 'K', 2: next neighbour incl. point itself
[allPelvis.acentre.sphere.nextIdx, allPelvis.acentre.sphere.nextDist] = knnsearch(allPelvis.acentre.sphere.points, allPelvis.acentre.sphere.points, 'K', 2); 
allPelvis.acentre.sphere.nextDistance = allPelvis.acentre.sphere.nextDist(:, 2); 
% Mean distance (mean edge length)
allPelvis.acentre.sphere.meanEdgeLength = mean(allPelvis.acentre.sphere.nextDistance);
% Percentage of points in range
allPelvis.acentre.sphere.lowerBound = 0.45;  
allPelvis.acentre.sphere.upperBound = 0.55;  
allPelvis.acentre.sphere.inRange = allPelvis.acentre.sphere.nextDistance >= ...
    allPelvis.acentre.sphere.lowerBound & allPelvis.acentre.sphere.nextDistance <= allPelvis.acentre.sphere.upperBound;
allPelvis.acentre.sphere.PercentinRange = (sum(allPelvis.acentre.sphere.inRange) / length(allPelvis.acentre.sphere.nextDistance)) * 100;

% Visualize
figure;
scatter3(allPelvis.acentre.sphere.pointsX, allPelvis.acentre.sphere.pointsY, allPelvis.acentre.sphere.pointsZ, 5, 'filled');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Points on the sphere');

%% Simple meshing / triangulation with acetabulum geometry - generate volume 

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Find vertices / indices of reference geometrie 
% (acetabulum incl. acetabular rim and backside of acetabulum)
% Which vertices of the reference pelvis are used for the generate volume method
% Distance to next neighbour with knnsearch
[allPelvis.acentre.refGeometry.nextIdx, allPelvis.acentre.refGeometry.nextDist] = knnsearch(...
    pelvis(1).import.processed.vertices, pelvisDefect(1).import.processed.vertices, 'K',1); 
allPelvis.acentre.refGeometry.threshold = 0.01; % for identical points
allPelvis.acentre.refGeometry.inThreshold = find(allPelvis.acentre.refGeometry.nextDist < allPelvis.acentre.refGeometry.threshold);
% Points/Vertices in pelvis point quantity
allPelvis.acentre.refGeometry.pelvisIdx = allPelvis.acentre.refGeometry.nextIdx(allPelvis.acentre.refGeometry.inThreshold);

% ShrinkWrap (alphaShape) for reference pelvis same as alpha
pelvisDefect(1).volume.shrink.alphaCupBack = pelvisDefect(1).volume.shrink.alpha; 
usedCount(1,1) = pelvisDefect(1).volume.shrink.alpha.usedDefectVerticesCount;
volumeAlpha(1,1) = pelvisDefect(1).volume.shrink.alpha.volume;
alphaRadius(1,1) = pelvisDefect(1).volume.shrink.alpha.radius;
timeAlpha(1,1) = NaN;

% ShrinkWrap (alphaShape) without user interface
% Alpha radius: first time mesh is closed and manifold
shrinkType = 'alphaCupBack'; % addPoints: generously cut acetabulum incl. acetabular rim + backside of acetabulum
addPoints = pelvisDefect(1).import.processed.vertices;
aRadius = 1; % Start value for alpha radius
parfor i = 2:dataCountDefect 
    tic
    numData = size(pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,1);
    pelvisDefect(i).volume = pelvisDefect(i).volume.shrinkAlpha(i,shrinkType,aRadius,numData,...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}), ...
        addPoints);
    usedCount(i,1) = pelvisDefect(i).volume.shrink.alphaCupBack.usedDefectVerticesCount;
    volumeAlpha(i,1) = pelvisDefect(i).volume.shrink.alphaCupBack.volume;
    alphaRadius(i,1) = pelvisDefect(i).volume.shrink.alphaCupBack.radius;
    timeAlpha(i,1) = toc;
end 
allDefect.time.alphaCupBack = timeAlpha;
allDefect.time.alphaCupBackMean = mean(timeAlpha,'omitnan');
allDefect.time.alphaCupBackStd = std(timeAlpha,'omitnan');

% ShrinkWrap (alphaShape) with user interface
shrinkType = 'alphaCupBack';
addPoints = pelvisDefect(1).import.processed.vertices;
%aRadius = 1; % Start value for alpha radius
timeAlphaUI(1:dataCountDefect,1) = NaN;
for i = 2:dataCountDefect % for-Loop because of plot and user input  % Start value for alpha radius
    tic
    aRadius = pelvisDefect(i).volume.shrink.alphaSphere.radius; % Use alpha radius from previous calculation without user input
    numData = size(pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,1);
    pelvisDefect(i).volume = pelvisDefect(i).volume.shrinkAlphaUI(i,shrinkType,aRadius,numData,...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}), addPoints); 
    usedCount(i,1) = pelvisDefect(i).volume.shrink.alphaCupBack.usedDefectVerticesCount;
    volumeAlpha(i,1) = pelvisDefect(i).volume.shrink.alphaCupBack.volume;
    alphaRadius(i,1) = pelvisDefect(i).volume.shrink.alphaCupBack.radius;
    timeAlphaUI(i,1) = toc;
end 
allDefect.time.alphaCupBackUI = timeAlphaUI;
allDefect.time.alphaCupBackUIMean = mean(timeAlphaUI,'omitnan');
allDefect.time.alphaCupBackUIStd = std(timeAlphaUI,'omitnan');

% AlphaRadius, used vertices for alphaShape and generated volume
allDefect.volume.alphaCupBack.radius = alphaRadius;
allDefect.volume.alphaCupBack.radiusMean = mean(alphaRadius);
allDefect.volume.alphaCupBack.radiusStd = std(alphaRadius);
allDefect.volume.alphaCupBack.radiusMax = max(alphaRadius);
allDefect.volume.alphaCupBack.radiusMin = min(alphaRadius);
allDefect.volume.alphaCupBack.usedDefectVerticesCount = usedCount;
allDefect.volume.alphaCupBack.usedDefectVerticesCountMean = mean(usedCount);
allDefect.volume.alphaCupBack.usedDefectVerticesCountStd = std(usedCount);
allDefect.volume.alphaCupBack.usedDefectVerticesCountMax = max(usedCount);
allDefect.volume.alphaCupBack.usedDefectVerticesCountMin = min(usedCount);
allDefect.volume.alphaCupBack.volume = volumeAlpha;
allDefect.volume.alphaCupBack.volumeMean = mean(volumeAlpha);
allDefect.volume.alphaCupBack.volumeStd = std(volumeAlpha);
allDefect.volume.alphaCupBack.volumeMax = max(volumeAlpha);
allDefect.volume.alphaCupBack.volumeMin = min(volumeAlpha);

% Define the different Paprosky types and their corresponding variable names
paproskyTypes = {'2a', '2b', '2c', '3a', '3b', 'NaN'};
paproskyNames = {'IIa', 'IIb', 'IIc', 'IIIa', 'IIIb', 'NaN'};
% Initialize a structure to store the volumes for each Paprosky type
allDefect.volume.alphaCupBack.paproskyVol = struct();
for i = 1:length(paproskyNames)
    allDefect.volume.alphaCupBack.paproskyVol.(paproskyNames{i}) = [];
end

% Collect the volumes for each Paprosky type
for i = 1:dataCountDefect
    paproskyType = pelvisDefect(i).import.processed.paprosky;
    idx = find(strcmp(paproskyTypes, paproskyType));
    if ~isempty(idx)
        varName = paproskyNames{idx};
        allDefect.volume.alphaCupBack.paproskyVol.(varName) = [allDefect.volume.alphaCupBack.paproskyVol.(varName); ...
            pelvisDefect(i).volume.shrink.alpha.volume];
    end
end
% Calculate the mean volume for each Paprosky type
for i = 1:length(paproskyNames)
    varName = paproskyNames{i};
    allDefect.volume.alphaCupBack.paproskyVolMean.(varName) = mean(allDefect.volume.alphaCupBack.paproskyVol.(varName), 'omitnan');
    allDefect.volume.alphaCupBack.paproskyVolStd.(varName) = std(allDefect.volume.alphaCupBack.paproskyVol.(varName), 'omitnan');
    allDefect.volume.alphaCupBack.paproskyVolMax.(varName) = max(allDefect.volume.alphaCupBack.paproskyVol.(varName));
    allDefect.volume.alphaCupBack.paproskyVolMin.(varName) = min(allDefect.volume.alphaCupBack.paproskyVol.(varName));
end

clear shrinkType addPoints aRadius numData usedCount volumeAlpha alphaRadius timeAlpha timeAlphaUI idx varName

%% Display scaled, transformed and alphaShaped defect (for control) 

% Transformation
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 1;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    
    % Pelvis defect alphaShape
    shrinkType = 'alphaCupBack'; % type: 'alpha', 'alphaCupBack' ... 
    % Pelvis defect alphaShape
    % patch('Faces',pelvisDefect(i).volume.shrink.(shrinkType).faces,...
    %     'Vertices',pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,...
    %     'FaceColor',[0.8 0.8 0.8], ...      % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect alphaShape (lightning)
    patch('Faces',pelvisDefect(i).volume.shrink.(shrinkType).faces,...
        'Vertices',pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,...
        'FaceColor',[0.8 0.8 0.8], ...      % Face color
        'FaceAlpha',1,...                   % Transparency of the faces;    % Figure Paper 0.5 
        'EdgeColor',TUMcolors.grey80  ,...  % Edge color;                    
        'EdgeAlpha',0.5,...                % Transparency of the edges         
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Pelvis defect
    patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...    % Edge color;                   % Figure Paper 'none'
        'EdgeAlpha',0.25,...                % Transparency of the edges      
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Edges of patches
    if shrinkType == "alphaCupBack"
        numLoops = size(pelvisDefect(i).volume.edge.edgeLoops,1);
        colors = jet(numLoops);  % Get distinct colors for each loop
        if size(pelvisDefect(i).volume.edge.edgeLoops,1) > 1
            for k = (pelvisDefect(i).volume.patches.all.numComponents+1):numLoops % Other patches
                currentColor = colors(k,:);  % Get color from the colormap colors(k,:)
                plot3(pelvisDefect(i).volume.edge.verticesLoops{k}(:,1),pelvisDefect(i).volume.edge.verticesLoops{k}(:,2),...
                    pelvisDefect(i).volume.edge.verticesLoops{k}(:,3),...
                    '.', 'Color', currentColor, 'MarkerSize', 30); % Figure Paper 'k' 
            end
        end
    end

    % % Pelvis defect vertices (+acentre)
    % plot3(pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints(:,1), ...
    %     pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints(:,2), ...
    %     pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints(:,3),...
    %     '.','Color',[0.2 0.2 0.2],'MarkerSize', 1)

    % Display edge (acetabulum edge)
    for k = 1:pelvisDefect(i).volume.patches.all.numComponents
        plot3(pelvisDefect(i).volume.edge.verticesLoops{k}(:,1),...
            pelvisDefect(i).volume.edge.verticesLoops{k}(:,2),...
            pelvisDefect(i).volume.edge.verticesLoops{k}(:,3),...
            '.','Color','r','MarkerSize',30) % Figure Paper 'r' 
    end

    % Used vertices of original defect for alphaShape (optional)
    % idxUsed = unique(pelvisDefect(i).volume.shrink.(shrinkType).faces(:));
    % Pused   = pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints(idxUsed, :);
    % plot3(Pused(:,1), Pused(:,2), Pused(:,3), '.','MarkerSize',30, 'Color', 'g');

    % Reference acetabulum centre
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
        '.','Color','m', 'MarkerSize', 80); % pelvis(i).transform.trafo.(weightKabsch{w}).acentre
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
        'x', 'MarkerEdgeColor', 'm', 'MarkerSize', 35, 'LineWidth', 10); % pelvis(i).transform.trafo.(weightKabsch{w}).acentre

    % Format and display properties
    title(['Alpha Shape: Pelvis Defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    view(80, -10) % Figure Paper
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')OriginalAlphaShapeCupBack.fig'])
%end

clear shrinkType idxUsed Pused
 
%% Display scaled, transformed and alphaShaped defect with reference pelvis (for control)

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 1;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    
    % Pelvis defect alphaShape
    shrinkType = 'alphaCupBack'; % type: 'alpha, 'alphaCupBack' ... %%%
    % Pelvis defect scaled and transformed
    % patch('Faces',pelvisDefect(i).volume.shrink.(shrinkType).faces,...
    %     'Vertices',pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,...
    %     'FaceColor',[0.9 0.75 0.68], ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect alphaShape (lightning)
    patch('Faces',pelvisDefect(i).volume.shrink.(shrinkType).faces,...
        'Vertices',pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Reference pelvis
    patch('Faces',pelvis(1).import.processed.faces,...
        'Vertices',pelvis(1).import.processed.vertices,...
        'FaceColor',TUMcolors.grey20, ...    % Face color
        'FaceAlpha',0.5,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');
    
    % Reference acetabulum centre
    plot3(pelvisDefect(1).import.processed.acentre(1), pelvisDefect(1).import.processed.acentre(2), ...
        pelvisDefect(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10); 
    
    % Format and display properties
    title(['Alpha Shape: Pelvis Defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')RefAlphaShapeCupBack.fig'])
%end
clear shrinkType

%% Meshing / triangulation - preprocess - reference inside 
% Combine methods (hybrid) - called inside, because of points inside

% Points inside alphaShape
shrinkType = 'alphaCupBack';
parfor  i = 1:dataCountDefect
    tic
    pelvisDefect(i).volume = pelvisDefect(i).volume.shrinkInside(i, shrinkType, ... % shrinkType
        pelvis(1).import.processed,... % pelvis(1).import.processed = pelvis(1).transform.trafo.(weightKabsch{w})
        allPelvis.refPoints.gridPoints.inside, allPelvis.acentre.refGeometry.pelvisIdx,...
        allPelvis.refPoints.gridPoints.pointsBoxNum, allPelvis.refPoints.gridPoints.insideMaskIdx); 
    verticesPointsNum(i,1) = size(pelvisDefect(i).volume.shrink.inside.refVerticesInside,1) + ...
        size(pelvisDefect(i).volume.shrink.inside.refPointsInside,1);
    timeRefInsideAlpha(i,1) = toc;
end 
allDefect.time.refInsideAlpha = timeRefInsideAlpha;
allDefect.time.refInsideAlphaMean = mean(timeRefInsideAlpha);
allDefect.time.refInsideAlphaStd = std(timeRefInsideAlpha);
allDefect.volume.inside.verticesPointsNum = verticesPointsNum; 

% Thicken mesh
thickness = 0.5;
step_size = 0.5;
parfor i = 1:dataCountDefect
    tic
    pelvisDefect(i).volume = pelvisDefect(i).volume.thickenMesh(i, ...
        pelvisDefect(i).volume.patches.all, thickness, step_size); % with combined mesh (defect mesh incl. refined patches) 
    thickenMeshPoints(i,1) = size(pelvisDefect(i).volume.shrink.inside.verticesThicken,1);
    timeThickenMesh(i,1) = toc;
end
allDefect.time.thickenMesh = timeThickenMesh;
allDefect.time.thickenMeshMean = mean(timeThickenMesh);
allDefect.time.thickenMeshStd = std(timeThickenMesh);
% Additional points
allDefect.volume.inside.thickenMeshPoints = thickenMeshPoints;
allDefect.volume.inside.thickenMeshMean = mean(thickenMeshPoints);
allDefect.volume.inside.thickenMeshStd = std(thickenMeshPoints);
allDefect.volume.inside.thickenMeshMax = max(thickenMeshPoints);
allDefect.volume.inside.thickenMeshMin = min(thickenMeshPoints);

clear verticesPointsNum timeRefInsideAlpha shrinkType thickenMeshPoints timeThickenMesh thickness step_size

%% Meshing / triangulation - preprocess - reference inside (filter)
% Combine methods (hybrid) - called inside, because of points inside

% Filter out points: 
% Identify points (reference vertices and points inside) which are in-/outside of the pelvis defect 
parfor  i = 1:dataCountDefect 
    tic
    pelvisDefect(i).volume = pelvisDefect(i).volume.inoutPoints(i,pelvisDefect(i).volume.patches.all,... % with combined mesh (defect mesh incl. refined patches)
        pelvis(1).import.processed); % pelvis(1).import.processed = pelvis(1).transform.trafo.(weightKabsch{w})
    verticesPointsNum(i,1) = size(pelvisDefect(i).volume.shrink.inside.inVertices,1) + ...
        size(pelvisDefect(i).volume.shrink.inside.inPoints,1);
    timeInoutPoints(i,1) = toc;
end
allDefect.time.inoutPoints = timeInoutPoints;
allDefect.time.inoutPointsMean = mean(timeInoutPoints);
allDefect.time.inoutPointsStd = std(timeInoutPoints);
allDefect.volume.inside.filteredInOutNum = verticesPointsNum; 

% Filter out points: 
% Identify areas with low point density (outliers) - also for clustering
thresholdLowDensity = 2.5; % threshold 2.5% percentile
parfor  i = 1:dataCountDefect
    tic
    pelvisDefect(i).volume = pelvisDefect(i).volume.densityPoints(i,thresholdLowDensity,...
        pelvis(1).import.processed); % pelvis(1).import.processed = pelvis(1).transform.trafo.(weightKabsch{w})
    executionTimesPointDensity(i,1) = pelvisDefect(i).volume.shrink.inside.executionTimeDensity; % only function mvksdensity
    executionTimesNumPointDensity(i,1) = size(pelvisDefect(i).volume.shrink.inside.pointCloudDensity,1);
    verticesPointsNum(i,1) = size(pelvisDefect(i).volume.shrink.inside.filteredInAll,1);
    timePointDensityAll(i,1) = toc;
end
allDefect.time.pointDensityAll = timePointDensityAll;
allDefect.time.pointDensityAllMean = mean(timePointDensityAll);
allDefect.time.pointDensityAllStd = std(timePointDensityAll);
allDefect.volume.inside.filteredDensityNum = verticesPointsNum; 

% Execution time of densityPoints (density function mvksdensity)
allDefect.time.pointDensity = executionTimesPointDensity;
allDefect.time.pointDensityMean = mean(executionTimesPointDensity);
allDefect.time.pointDensityStd = std(executionTimesPointDensity);
allDefect.time.pointDensityMax = max(executionTimesPointDensity);
allDefect.time.pointDensityMin = min(executionTimesPointDensity);
allDefect.volume.inside.pointDensityNum = executionTimesNumPointDensity;
allDefect.volume.inside.pointDensityNumMean = mean(executionTimesNumPointDensity);
allDefect.volume.inside.pointDensityNumStd = std(executionTimesNumPointDensity);
allDefect.volume.inside.pointDensityNumMax = max(executionTimesNumPointDensity);
allDefect.volume.inside.pointDensityNumMin = min(executionTimesNumPointDensity);

% Filter out points: cluster
% Identify clusters and outliers which are seperated (outside of pelvis defect) - keep cluster with most points
% DBSCAN Parameter (cluster method)
epsilon = 0.6;  % Radius for neighbours (depend on data set; grid of points 0.5 mm)
minPts = 5;     % Minimal Points to create cluster
for  i = 1:dataCountDefect 
    tic
    pelvisDefect(i).volume = pelvisDefect(i).volume.clusterPoints(i,epsilon,minPts,...
        pelvisDefect(i).volume.patches.all, ... % with combined mesh (defect mesh incl. refined patches)
        pelvis(1).import.processed); % pelvis(1).import.processed / pelvis(1).transform.trafo.(weightKabsch{w})
    executionTimesPointCluster(i,1) = pelvisDefect(i).volume.shrink.inside.executionTimeCluster; % only function dbscan
    executionTimesNumPointCluster(i,1) = size(pelvisDefect(i).volume.shrink.inside.clusterLabels,1);  
    filteredVerticesPointsNum(i,1) = size(pelvisDefect(i).volume.shrink.inside.mainClusterVertices,1) + ... 
        size(pelvisDefect(i).volume.shrink.inside.mainClusterPoints,1);
    timePointClusterAll(i,1) = toc;
end 
allDefect.time.pointClusterAll = timePointClusterAll;
allDefect.time.pointClusterAllMean = mean(timePointClusterAll);
allDefect.time.pointClusterAllStd = std(timePointClusterAll);
allDefect.volume.inside.filteredVerticesPointsNum = filteredVerticesPointsNum;

% Execution time of clusterPoints (cluster function dbscan)
allDefect.time.pointCluster = executionTimesPointCluster;
allDefect.time.pointClusterMean = mean(executionTimesPointCluster);
allDefect.time.pointClusterStd = std(executionTimesPointCluster);
allDefect.time.pointClusterMax = max(executionTimesPointCluster);
allDefect.time.pointClusterMin = min(executionTimesPointCluster);
allDefect.volume.inside.clusterNum = executionTimesNumPointCluster;
allDefect.volume.inside.clusterNumMean = mean(executionTimesNumPointCluster);
allDefect.volume.inside.clusterNumStd = std(executionTimesNumPointCluster);
allDefect.volume.inside.clusterNumMax = max(executionTimesNumPointCluster);
allDefect.volume.inside.clusterNumMin = min(executionTimesNumPointCluster); 

% Cluster parameter 
allDefect.volume.inside.detectedOtherClusters = cell(dataCountDefect, 1); 
allDefect.volume.inside.keptOtherClusters = cell(dataCountDefect, 1);     
for i = 1:dataCountDefect
    allDefect.volume.inside.detectedOtherClusters{i} = []; 
    allDefect.volume.inside.keptOtherClusters{i} = [];     
    % Detected "otherClusters"
    if ~isempty(pelvisDefect(i).volume.shrink.inside.otherClustersVertices) && ...
       numel(pelvisDefect(i).volume.shrink.inside.otherClustersVertices) >= 1
       for j = 1:length(pelvisDefect(i).volume.shrink.inside.otherClustersVertices)
            % Detected "otherClusters", but not kept
            if size(pelvisDefect(i).volume.shrink.inside.otherClustersVertices, 1) == 1 && ...
               isempty(pelvisDefect(i).volume.shrink.inside.otherClustersVertices{j}) && ...
               isempty(pelvisDefect(i).volume.shrink.inside.otherClustersPoints{j})
                allDefect.volume.inside.detectedOtherClusters{i} = [allDefect.volume.inside.detectedOtherClusters{i}, j];
            end
            % Kept "otherClusters"
            if ~isempty(pelvisDefect(i).volume.shrink.inside.otherClustersVertices{j}) || ...
               ~isempty(pelvisDefect(i).volume.shrink.inside.otherClustersPoints{j})
                allDefect.volume.inside.detectedOtherClusters{i} = [allDefect.volume.inside.detectedOtherClusters{i}, j];
                allDefect.volume.inside.keptOtherClusters{i} = [allDefect.volume.inside.keptOtherClusters{i}, j];
            end
       end
    end
    allDefect.volume.inside.has_outliers(i,1) = ~isempty(pelvisDefect(i).volume.shrink.inside.outliersVertices) || ...
                   ~isempty(pelvisDefect(i).volume.shrink.inside.outliersPoints);
end
% Cluster Idx
allDefect.volume.inside.detectedOtherClustersIdx = [];
allDefect.volume.inside.keptOtherClustersIdx = [];
for i = 1:dataCountDefect
    if ~isempty(allDefect.volume.inside.detectedOtherClusters{i})
        allDefect.volume.inside.detectedOtherClustersIdx = [allDefect.volume.inside.detectedOtherClustersIdx; i];
    end
    if ~isempty(allDefect.volume.inside.keptOtherClusters{i})
        allDefect.volume.inside.keptOtherClustersIdx = [allDefect.volume.inside.keptOtherClustersIdx; i];
    end
end
allDefect.volume.inside.outliersClusters = find(allDefect.volume.inside.has_outliers);

% Logical mask for reference vertices/points inside
for i = 1:dataCountDefect
    % Vertices
    pelvisDefect(i).volume.shrink.inside.isVertexInside = pelvisDefect(i).volume.shrink.inside.mainClusterVerticesMask;
    % Points
    pelvisDefect(i).volume.shrink.inside.isPointInside = pelvisDefect(i).volume.shrink.inside.mainClusterPointsMask;
end

clear timeInoutPoints verticesPointsNum thresholdLowDensity executionTimesPointDensity executionTimesNumPointDensity timePointDensityAll ...
    epsilon minPts executionTimesPointCluster executionTimesNumPointCluster filteredVerticesPointsNum timePointClusterAl

%% Display scaled, transformed defect with reference pelvis points/vertices inside - shrinkInside (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 1; % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    
    % Pelvis defect remeshed (scaled and transformed)
    % For Figure Paper
    % patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
    %      'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
    patch('Faces',pelvisDefect(i).volume.patches.all.comFaces,...
         'Vertices',pelvisDefect(i).volume.patches.all.comVertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',TUMcolors.grey50,...    % Edge color % Figure Paper 'none'
        'EdgeAlpha',0.25);                  % Transparency of the edges
    light('Position', [1 1 5], 'Style', 'infinite');

    % % Reference pelvis mesh inside
    patch('Faces',pelvisDefect(i).volume.shrink.inside.refFacesInside,...
        'Vertices',pelvis(1).transform.trafo.(weightKabsch{w}).vertices,... % pelvis(1).import.processed 
        'FaceColor',[0.8 0.8 0.8], ...      % Face color
        'FaceAlpha',0.5,...                   % Transparency of the faces % Figure Paper 0.25
        'EdgeColor',TUMcolors.grey50,...    % Edge color
        'EdgeAlpha',0.25);                  % Transparency of the edges
    % Reference pelvis mesh inside: vertices
    plot3(pelvisDefect(i).volume.shrink.inside.refVerticesInside(:,1), ...
        pelvisDefect(i).volume.shrink.inside.refVerticesInside(:,2), ...
        pelvisDefect(i).volume.shrink.inside.refVerticesInside(:,3), ...
        '.','Color',TUMcolors.orange, 'MarkerSize', 1);  % Figure Paper 'c' 

    % Reference pelvis points inside
    plot3(pelvisDefect(i).volume.shrink.inside.refPointsInside(:,1), ...
        pelvisDefect(i).volume.shrink.inside.refPointsInside(:,2), ...
        pelvisDefect(i).volume.shrink.inside.refPointsInside(:,3), ...
        '.','Color',TUMcolors.blue300, 'MarkerSize', 1); % Figure Paper 'c' 

    % Reference acetabulum centre
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
        '.','Color','m', 'MarkerSize', 80); % pelvis(i).transform.trafo.(weightKabsch{w}).acentre
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
        'x', 'MarkerEdgeColor', 'm', 'MarkerSize', 35, 'LineWidth', 10); % pelvis(i).transform.trafo.(weightKabsch{w}).acentre

    % Display edge (acetabulum edge) (optional)
    % for k = 1:pelvisDefect(i).volume.patches.all.numComponents
    %     plot3(pelvisDefect(i).volume.edge.verticesLoops{k}(:,1),...
    %         pelvisDefect(i).volume.edge.verticesLoops{k}(:,2),...
    %         pelvisDefect(i).volume.edge.verticesLoops{k}(:,3),...
    %         '.','Color', TUMcolors.orange,'MarkerSize',30) % Figure Paper 'r' 
    % end

    % Edges of patches (optional)
    % numLoops = size(pelvisDefect(i).volume.edge.edgeLoops,1);
    % colors = jet(numLoops);  % Get distinct colors for each loop
    % if size(pelvisDefect(i).volume.edge.edgeLoops,1) > 1
    %     for k = (pelvisDefect(i).volume.patches.all.numComponents+1):numLoops % Other patches
    %         currentColor = 'k';  % Get color from the colormap colors(k,:)
    %         plot3(pelvisDefect(i).volume.edge.verticesLoops{k}(:,1),pelvisDefect(i).volume.edge.verticesLoops{k}(:,2),...
    %             pelvisDefect(i).volume.edge.verticesLoops{k}(:,3),...
    %             '.', 'Color', currentColor, 'MarkerSize', 30); % Figure Paper 'm'
    %     end
    % end

    % Format and display properties
    title(['Vetices/Points Inside (shrink wrap): Pelvis Defect ' num2str(i)]);
    legend('Pelvis defect', 'Reference pelvis inside', 'Reference vertices inside', 'Reference points inside', 'Acetabulum centre');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    view(80, -10) % Figure Paper
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')shrinkInside.fig'])
%end

%% Display scaled, transformed and thicken defect with filtered reference pelvis points/vertices inside (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect 
    i = 2; % Pelvis number
    %figure('Visible','off')   
    figure
    hold on

    % Pelvis defect remeshed (scaled and transformed)
    patch('Faces',pelvisDefect(i).volume.patches.all.comFaces,...
        'Vertices',pelvisDefect(i).volume.patches.all.comVertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...    % Edge color TUMcolors.grey50
        'EdgeAlpha',0.25);                  % Transparency of the edges
    light('Position', [1 1 5], 'Style', 'infinite');

    % Reference pelvis mesh inside
    patch('Faces',pelvisDefect(i).volume.shrink.inside.mainClusterFaces,...
        'Vertices',pelvis(1).transform.trafo.(weightKabsch{w}).vertices,... % pelvis(1).import.processed
        'FaceColor',[0.8 0.8 0.8], ...      % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',TUMcolors.grey50,...    % Edge color
        'EdgeAlpha',0.25);                  % Transparency of the edges

    % Main cluster: vertices
    plot3(pelvisDefect(i).volume.shrink.inside.mainClusterVertices(:,1), ...
        pelvisDefect(i).volume.shrink.inside.mainClusterVertices(:,2), ...
        pelvisDefect(i).volume.shrink.inside.mainClusterVertices(:,3), '.','Color', TUMcolors.orange, 'MarkerSize', 2);
    % Main cluster: points
    plot3(pelvisDefect(i).volume.shrink.inside.mainClusterPoints(:,1), ...
        pelvisDefect(i).volume.shrink.inside.mainClusterPoints(:,2), ...
        pelvisDefect(i).volume.shrink.inside.mainClusterPoints(:,3), '.','Color', TUMcolors.blue300, 'MarkerSize', 2);

    % Thicken mesh
    plot3(pelvisDefect(i).volume.shrink.inside.verticesThicken(:,1), ...
        pelvisDefect(i).volume.shrink.inside.verticesThicken(:,2), ...
        pelvisDefect(i).volume.shrink.inside.verticesThicken(:,3), '.','Color', TUMcolors.grey20, 'MarkerSize', 2);

    % Format and display properties
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['Thicken Pelvis Defect ' num2str(i)]);
    legend('Pelvis Defect','Reference pelvis','Reference vertices inside', 'Reference points inside','Thicken defect mesh');
    %legend('Pelvis Defect','Reference pelvis','Thicken defect mesh');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')thicken.fig'])
%end

%% Display scaled, transformed defect with filtered out (inoutPoints) reference pelvis points/vertices inside (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 1; % Pelvis number
    %figure('Visible','off')   
    figure
    hold on

    % Pelvis defect remeshed (scaled and transformed)
    patch('Faces',pelvisDefect(i).volume.patches.all.comFaces,...
        'Vertices',pelvisDefect(i).volume.patches.all.comVertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',TUMcolors.grey50,...    % Edge color
        'EdgeAlpha',0.25);                  % Transparency of the edges
    light('Position', [1 1 5], 'Style', 'infinite');

    % Reference pelvis mesh inside
    patch('Faces',pelvisDefect(i).volume.shrink.inside.inFaces,...
        'Vertices',pelvis(1).transform.trafo.(weightKabsch{w}).vertices,... % pelvis(1).import.processed
        'FaceColor',[0.8 0.8 0.8], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',TUMcolors.grey50,...    % Edge color
        'EdgeAlpha',0.25);                  % Transparency of the edges

    % Reference pelvis points/vertices
    % Reference pelvis mesh outside: vertices (optional)
    % plot3(pelvisDefect(i).volume.shrink.inside.outVertices(:,1), ...
    %     pelvisDefect(i).volume.shrink.inside.outVertices(:,2), ...
    %     pelvisDefect(i).volume.shrink.inside.outVertices(:,3), '.','Color', TUMcolors.orange, 'MarkerSize', 2); % green
    % Reference pelvis mesh inside: vertices
    plot3(pelvisDefect(i).volume.shrink.inside.inVertices(:,1), ...
        pelvisDefect(i).volume.shrink.inside.inVertices(:,2), ...
        pelvisDefect(i).volume.shrink.inside.inVertices(:,3), '.','Color', TUMcolors.green, 'MarkerSize', 5); % orange
    % Reference pelvis points outside (optional)
    % plot3(pelvisDefect(i).volume.shrink.inside.outPoints(:,1), ...
    %     pelvisDefect(i).volume.shrink.inside.outPoints(:,2), ...
    %     pelvisDefect(i).volume.shrink.inside.outPoints(:,3), '.', 'Color', 'k', 'MarkerSize', 2); % black
    % Reference pelvis points inside
    plot3(pelvisDefect(i).volume.shrink.inside.inPoints(:,1), ...
        pelvisDefect(i).volume.shrink.inside.inPoints(:,2), ...
        pelvisDefect(i).volume.shrink.inside.inPoints(:,3), '.','Color', TUMcolors.blue300, 'MarkerSize', 5); % blue

    % Reference acetabulum centre
    plot3(pelvisDefect(1).import.processed.acentre(1), pelvisDefect(1).import.processed.acentre(2), ...
        pelvisDefect(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10); 

    % Format and display properties
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['Filtered Points/vetices (in/out): Pelvis Defect ' num2str(i)]);
    %legend('Pelvis defect','Reference pelvis','Reference vertices outside','Reference vertices inside', 'Reference points outside', 'Reference points inside', 'Acetabulum centre');
    legend('Pelvis defect','Reference pelvis','Reference vertices inside',  'Reference points inside', 'Acetabulum centre');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')inoutPointsInside.fig'])
%end

%% Display scaled, transformed defect with filtered out (densityPoints) reference pelvis points/vertices inside (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 1; % Pelvis number
    %figure('Visible','off')   
    figure
    hold on

    % Pelvis defect remeshed (scaled and transformed)
    patch('Faces',pelvisDefect(i).volume.patches.all.comFaces,...
        'Vertices',pelvisDefect(i).volume.patches.all.comVertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...    % Edge color TUMcolors.grey50
        'EdgeAlpha',0.25);                  % Transparency of the edges
    light('Position', [1 1 5], 'Style', 'infinite');

    % Reference pelvis mesh inside
    patch('Faces',pelvisDefect(i).volume.shrink.inside.filteredInFaces,...
        'Vertices',pelvis(1).transform.trafo.(weightKabsch{w}).vertices,... % pelvis(1).import.processed
        'FaceColor',[0.8 0.8 0.8], ...      % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',TUMcolors.grey50,...    % Edge color
        'EdgeAlpha',0.25);                  % Transparency of the edges

    % Point cloud with density (colorbar); all points/vertices (optional)
    insideVerticesPoints = [pelvisDefect(i).volume.shrink.inside.inVertices; pelvisDefect(i).volume.shrink.inside.inPoints];
    scatter3(insideVerticesPoints(:,1), insideVerticesPoints(:,2), insideVerticesPoints(:,3), 10, ...
        pelvisDefect(i).volume.shrink.inside.densityRGB, 'filled');

    % Reference pelvis mesh inside: vertices
    plot3(pelvisDefect(i).volume.shrink.inside.filteredInVertices(:,1), ...
        pelvisDefect(i).volume.shrink.inside.filteredInVertices(:,2), ...
        pelvisDefect(i).volume.shrink.inside.filteredInVertices(:,3), '.','Color', TUMcolors.orange, 'MarkerSize', 2); % orange
    % Reference pelvis points inside
    plot3(pelvisDefect(i).volume.shrink.inside.filteredInPoints(:,1), pelvisDefect(i).volume.shrink.inside.filteredInPoints(:,2), ...
        pelvisDefect(i).volume.shrink.inside.filteredInPoints(:,3), '.','Color', TUMcolors.blue300, 'MarkerSize', 2); % blue  
    % Filtered out reference vertices/points (optional)
    % plot3(pelvisDefect(i).volume.shrink.inside.filteredOutAll(:,1), pelvisDefect(i).volume.shrink.inside.filteredOutAll(:,2), ...
    %     pelvisDefect(i).volume.shrink.inside.filteredOutAll(:,3), '.','Color', 'r', 'MarkerSize', 5);

    % Format and display properties
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['Density point cloud: pelvis defect ' num2str(i)]);
    legend('Pelvis Defect','Reference pelvis','Reference vertices inside', 'Reference points inside', 'Filtered out vertices/points');
    %legend('Pelvis Defect','Reference pelvis','Reference vertices inside', 'Reference points inside');
    colormap('viridis');
    colorbar
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')densityPointsRGB.fig'])
%end

clear insideVerticesPoints

%% Display scaled, transformed defect with filtered out (clusterPoints) reference pelvis points/vertices inside (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect 
    i = 1; % Pelvis number
    %figure('Visible','off')   
    figure
    hold on

    % Pelvis defect remeshed (scaled and transformed)
    % For Figure Paper
    % patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
    %     'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
    p1 = patch('Faces',pelvisDefect(i).volume.patches.all.comFaces,...
        'Vertices',pelvisDefect(i).volume.patches.all.comVertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',TUMcolors.grey50,...    % Edge color % Figure Paper 'none' 
        'EdgeAlpha',0.25);                  % Transparency of the edges
    light('Position', [1 1 5], 'Style', 'infinite');

    % Reference pelvis mesh inside
    p2 = patch('Faces',pelvisDefect(i).volume.shrink.inside.mainClusterFaces,...
        'Vertices',pelvis(1).transform.trafo.(weightKabsch{w}).vertices,... % pelvis(1).import.processed
        'FaceColor',[0.8 0.8 0.8], ...      % Face color
        'FaceAlpha',1,...                   % Transparency of the faces % Figure Paper 0.5 
        'EdgeColor',TUMcolors.grey50,...    % Edge color
        'EdgeAlpha',0.25);                  % Transparency of the edges

    % Main cluster: vertices
    hMainVertices = plot3(pelvisDefect(i).volume.shrink.inside.mainClusterVertices(:,1), ...
        pelvisDefect(i).volume.shrink.inside.mainClusterVertices(:,2), ...
        pelvisDefect(i).volume.shrink.inside.mainClusterVertices(:,3), '.','Color', TUMcolors.orange, 'MarkerSize', 2); % Figure Paper 'c' 
    % Main cluster: points
    hMainPoints = plot3(pelvisDefect(i).volume.shrink.inside.mainClusterPoints(:,1), ...
        pelvisDefect(i).volume.shrink.inside.mainClusterPoints(:,2), ...
        pelvisDefect(i).volume.shrink.inside.mainClusterPoints(:,3), '.','Color', TUMcolors.blue300, 'MarkerSize', 2); % Figure Paper 'c' 
    
    % Other clusters with different colors
    % Create a color map for the clusters
    hOtherClusters = [];
    num_clusters = length(pelvisDefect(i).volume.shrink.inside.clusterLabels);
    colors = lines(num_clusters);
    for k = 1:length(pelvisDefect(i).volume.shrink.inside.otherClustersVertices)
        % Other clusters vertices 
        if ~isempty(pelvisDefect(i).volume.shrink.inside.otherClustersVertices{k})
            hOtherClustersVertices = plot3(pelvisDefect(i).volume.shrink.inside.otherClustersVertices{k}(:,1), ...
                pelvisDefect(i).volume.shrink.inside.otherClustersVertices{k}(:,2), ...
                pelvisDefect(i).volume.shrink.inside.otherClustersVertices{k}(:,3), '.', ...
                'Color', colors(k + 1,:), 'MarkerSize', 2); 
            hOtherClusters = [hOtherClusters, hOtherClustersVertices];
        end
        % Other clusters points 
        if ~isempty(pelvisDefect(i).volume.shrink.inside.otherClustersPoints{k})
            hOtherClustersPoints = plot3(pelvisDefect(i).volume.shrink.inside.otherClustersPoints{k}(:,1), ...
                pelvisDefect(i).volume.shrink.inside.otherClustersPoints{k}(:,2), ...
                pelvisDefect(i).volume.shrink.inside.otherClustersPoints{k}(:,3), '.', ...
                'Color', colors(k + 1,:), 'MarkerSize', 2); 
            hOtherClusters = [hOtherClusters, hOtherClustersPoints];
        end
    end
    % Outliers
    hOutliers = [];
    if ~isempty(pelvisDefect(i).volume.shrink.inside.outliersVertices)
        hOutliersVertices = plot3(pelvisDefect(i).volume.shrink.inside.outliersVertices(:,1), ...
            pelvisDefect(i).volume.shrink.inside.outliersVertices(:,2), ...
            pelvisDefect(i).volume.shrink.inside.outliersVertices(:,3), ...
            '.', 'Color', 'r', 'MarkerSize', 2); 
        hOutliers = [hOutliers, hOutliersVertices];
    end
    if ~isempty(pelvisDefect(i).volume.shrink.inside.outliersPoints)
        hOutliersPoints = plot3(pelvisDefect(i).volume.shrink.inside.outliersPoints(:,1), ...
            pelvisDefect(i).volume.shrink.inside.outliersPoints(:,2), ...
            pelvisDefect(i).volume.shrink.inside.outliersPoints(:,3), ...
            'x', 'Color', 'r', 'MarkerSize', 5);
        hOutliers = [hOutliers, hOutliersPoints];
    end

    % Legend
    legend_entries = {'Pelvis Defect','Reference pelvis','Main Cluster Vertices', 'Main Cluster Points'};
    legend_handles = [p1, p2, hMainVertices, hMainPoints];
    if ~isempty(hOtherClusters)
        legend_entries = [legend_entries, 'Other Clusters Vertices/Points'];
        legend_handles = [legend_handles, hOtherClusters];
    end
    if ~isempty(hOutliers)
        legend_entries = [legend_entries, 'Outliers Vertices/Points'];
        legend_handles = [legend_handles, hOutliers];
    end
    legend(legend_handles, legend_entries);

    % Reference acetabulum centre
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
        '.','Color','m', 'MarkerSize', 80); % pelvis(i).transform.trafo.(weightKabsch{w}).acentre
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
        'x', 'MarkerEdgeColor', 'm', 'MarkerSize', 35, 'LineWidth', 10); % pelvis(i).transform.trafo.(weightKabsch{w}).acentre

    % Display edge (acetabulum edge) (optional)
    % for k = 1:pelvisDefect(i).volume.patches.all.numComponents
    %     plot3(pelvisDefect(i).volume.edge.verticesLoops{k}(:,1),...
    %         pelvisDefect(i).volume.edge.verticesLoops{k}(:,2),...
    %         pelvisDefect(i).volume.edge.verticesLoops{k}(:,3),...
    %         '.','Color', 'r','MarkerSize',30) % Figure Paper 'r'
    % end

    % Edges of patches (optional)

    numLoops = size(pelvisDefect(i).volume.edge.edgeLoops,1);
    colors = jet(numLoops);  % Get distinct colors for each loop
    if size(pelvisDefect(i).volume.edge.edgeLoops,1) > 1
        for k = (pelvisDefect(i).volume.patches.all.numComponents+1):numLoops % Other patches
            currentColor = 'k';  % Get color from the colormap colors(k,:)
            plot3(pelvisDefect(i).volume.edge.verticesLoops{k}(:,1),pelvisDefect(i).volume.edge.verticesLoops{k}(:,2),...
                pelvisDefect(i).volume.edge.verticesLoops{k}(:,3),...
                '.', 'Color', currentColor, 'MarkerSize', 30); % Figure Paper 'k'
        end
    end

    % Format and display properties
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['DBSCAN Clustering Result: pelvis defect ' num2str(i)]);
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    view(80, -10) % Figure Paper
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')clusterPoints.fig'])
%end

clear p1 p2 hMainVertices hMainPoints hOtherClusters num_clusters colors hOtherClustersVertices hOtherClustersPoints hOutliers ...
    hOutliersVertices hOutliersPoints legend_entries legend_handles

%% Meshing / triangulation - generate volume 
% Combine methods (hybrid) - called inside, because of points inside

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% ShrinkWrap (alphaShape) for reference pelvis same as alpha
pelvisDefect(1).volume.shrink.refinedAlpha = pelvisDefect(1).volume.shrink.alpha; 
usedCount(1,1) = pelvisDefect(1).volume.shrink.alpha.usedDefectVerticesCount;
volumeAlpha(1,1) = pelvisDefect(1).volume.shrink.alpha.volume;
alphaRadius(1,1) = pelvisDefect(1).volume.shrink.alpha.radius;
timeRefinedAlpha(1,1) = NaN;

% Reduced vertices/points 
% Identify boundary points with point density (point density at the boundary lower)
boundaryThreshold = 33.33; 
parfor i = 1:dataCountDefect
    insideVerticesPoints = [pelvisDefect(i).volume.shrink.inside.inVertices; pelvisDefect(i).volume.shrink.inside.inPoints];
    thresholdLowDensity = prctile(pelvisDefect(i).volume.shrink.inside.pointCloudNormDensity,boundaryThreshold);
    lowDensityBoundary = pelvisDefect(i).volume.shrink.inside.pointCloudNormDensity <= thresholdLowDensity;
    % Remove filtered out vertices/points (density)
    lowDensityBoundaryReduced = lowDensityBoundary(pelvisDefect(i).volume.shrink.inside.inoutDensity,:);
    insideVerticesPointsReduced = insideVerticesPoints(pelvisDefect(i).volume.shrink.inside.inoutDensity,:);
    % Remove filtered out vertices/points (cluster)
    pelvisDefect(i).volume.shrink.inside.mainClusterMask = ...
        (pelvisDefect(i).volume.shrink.inside.clusterLabels == pelvisDefect(i).volume.shrink.inside.clusterMainLabel);
    boundaryMask = lowDensityBoundaryReduced(pelvisDefect(i).volume.shrink.inside.mainClusterMask,:);
    insideVerticesPointsReduced = insideVerticesPointsReduced(pelvisDefect(i).volume.shrink.inside.mainClusterMask,:);
    pelvisDefect(i).volume.shrink.refinedAlpha.boundaryVerticesPoints = insideVerticesPointsReduced(boundaryMask,:);
end
clear insideVerticesPoints thresholdLowDensity lowDensityBoundary lowDensityBoundaryReduced insideVerticesPointsReduced boundaryMask boundaryThreshold

% ShrinkWrap (alphaShape) without user interface (afterwards check with user input)
% Alpha radius: first time mesh is closed and manifold 
shrinkType = 'refinedAlpha'; 
aRadius = 1; % Start value for alpha Radius
parfor  i = 1:dataCountDefect
    tic
    numData = size(pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,1);
    addPoints = [pelvisDefect(i).volume.shrink.inside.mainClusterVertices; pelvisDefect(i).volume.shrink.inside.mainClusterPoints; ...
    pelvisDefect(i).volume.shrink.inside.verticesThicken]; % all points
    pelvisDefect(i).volume = pelvisDefect(i).volume.shrinkAlpha(i,shrinkType,aRadius,numData,...
        pelvisDefect(i).volume.patches.all, addPoints); 
    usedCount(i,1) = pelvisDefect(i).volume.shrink.refinedAlpha.usedDefectVerticesCount;
    volumeAlpha(i,1) = pelvisDefect(i).volume.shrink.refinedAlpha.volume;
    alphaRadius(i,1) = pelvisDefect(i).volume.shrink.refinedAlpha.radius;
    timeRefinedAlpha(i,1) = toc;
end 
allDefect.time.refinedAlpha = timeRefinedAlpha;
allDefect.time.refinedAlphaMean = mean(timeRefinedAlpha,'omitnan');
allDefect.time.refinedAlphaStd = std(timeRefinedAlpha,'omitnan');

% ShrinkWrap (alphaShape) with user interface (post-correction) 
shrinkType = 'refinedAlpha'; 
timeRefinedAlphaUI(1:dataCountDefect,1) = NaN;
for  i = 2:dataCountDefect % for-Loop because of plot and user input  
    tic
    aRadius = pelvisDefect(i).volume.shrink.refinedAlpha.radius; % Start value alpha radius
    numData = size(pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,1);
    addPoints = [pelvisDefect(i).volume.shrink.inside.mainClusterVertices; pelvisDefect(i).volume.shrink.inside.mainClusterPoints; ...
    pelvisDefect(i).volume.shrink.inside.verticesThicken]; % all points
    pelvisDefect(i).volume = pelvisDefect(i).volume.shrinkAlphaUI(i,shrinkType,aRadius,numData,...
        pelvisDefect(i).volume.patches.all, addPoints); 
    usedCount(i,1) = pelvisDefect(i).volume.shrink.refinedAlpha.usedDefectVerticesCount;
    volumeAlpha(i,1) = pelvisDefect(i).volume.shrink.refinedAlpha.volume;
    alphaRadius(i,1) = pelvisDefect(i).volume.shrink.refinedAlpha.radius;
    timeRefinedAlphaUI(i,1) = toc;
end 
allDefect.time.refinedAlphaUI = timeRefinedAlphaUI;
allDefect.time.refinedAlphaUIMean = mean(allDefect.time.refinedAlphaUI,'omitnan');
allDefect.time.refinedAlphaUIStd = std(allDefect.time.refinedAlphaUI,'omitnan');

% AlphaRadius, used vertices for alphaShape and generated volume
allDefect.volume.refinedAlpha.radius = alphaRadius;
allDefect.volume.refinedAlpha.radiusMean = mean(alphaRadius);
allDefect.volume.refinedAlpha.radiusStd = std(alphaRadius);
allDefect.volume.refinedAlpha.radiusMax = max(alphaRadius);
allDefect.volume.refinedAlpha.radiusMin = min(alphaRadius);
allDefect.volume.refinedAlpha.usedDefectVerticesCount = usedCount;
allDefect.volume.refinedAlpha.usedDefectVerticesCountMean = mean(usedCount);
allDefect.volume.refinedAlpha.usedDefectVerticesCountStd = std(usedCount);
allDefect.volume.refinedAlpha.usedDefectVerticesCountMax = max(usedCount);
allDefect.volume.refinedAlpha.usedDefectVerticesCountMin = min(usedCount);
allDefect.volume.refinedAlpha.volume = volumeAlpha;
allDefect.volume.refinedAlpha.volumeMean = mean(volumeAlpha);
allDefect.volume.refinedAlpha.volumeStd = std(volumeAlpha);
allDefect.volume.refinedAlpha.volumeMax = max(volumeAlpha);
allDefect.volume.refinedAlpha.volumeMin = min(volumeAlpha);

% Define the different Paprosky types and their corresponding variable names
paproskyTypes = {'2a', '2b', '2c', '3a', '3b', 'NaN'};
paproskyNames = {'IIa', 'IIb', 'IIc', 'IIIa', 'IIIb', 'NaN'};
% Initialize a structure to store the volumes for each Paprosky type
allDefect.volume.refinedAlpha.paproskyVol = struct();
for i = 1:length(paproskyNames)
    allDefect.volume.refinedAlpha.paproskyVol.(paproskyNames{i}) = [];
end

% Collect the volumes for each Paprosky type
for i = 1:dataCountDefect
    paproskyType = pelvisDefect(i).import.processed.paprosky;
    idx = find(strcmp(paproskyTypes, paproskyType));
    if ~isempty(idx)
        varName = paproskyNames{idx};
        allDefect.volume.refinedAlpha.paproskyVol.(varName) = [allDefect.volume.refinedAlpha.paproskyVol.(varName); ...
            pelvisDefect(i).volume.shrink.refinedAlpha.volume];
    end
end
% Calculate the mean volume for each Paprosky type
for i = 1:length(paproskyNames)
    varName = paproskyNames{i};
    allDefect.volume.refinedAlpha.paproskyVolMean.(varName) = mean(allDefect.volume.refinedAlpha.paproskyVol.(varName), 'omitnan');
    allDefect.volume.refinedAlpha.paproskyVolStd.(varName) = std(allDefect.volume.refinedAlpha.paproskyVol.(varName), 'omitnan');
    allDefect.volume.refinedAlpha.paproskyVolMax.(varName) = max(allDefect.volume.refinedAlpha.paproskyVol.(varName));
    allDefect.volume.refinedAlpha.paproskyVolMin.(varName) = min(allDefect.volume.refinedAlpha.paproskyVol.(varName));
end

%% Comparison alpha - refined alpha

% Comparison refined alpha shape to alpha shape with vertex to nearest neighbour 
compAlphas(1) = NaN;
parfor i = 2:dataCountDefect
    % Find the nearest neighbour in alpha (ref) for each vertex in refined alpha 
    % (alpha vertices partial amount from refined alpha)
    pelvisDefect(i).volume = pelvisDefect(i).volume.shapeVertexNeighbour(i, ...
                pelvisDefect(i).volume.shrink.alpha.usedVerticesPoints, ... % alpha (ref)
                pelvisDefect(i).volume.shrink.refinedAlpha.usedVerticesPoints,... % refined alpha 
                pelvisDefect(i).volume.shrink.refinedAlpha.usedFaces);
    compAlphas(i,1) = pelvisDefect(i).volume.shrink.comp.nearVertexMean;
end
allDefect.volume.comp.nearVertex = compAlphas;
allDefect.volume.comp.nearVertexMean = mean(compAlphas,'omitnan'); 
allDefect.volume.comp.nearVertexStd = std(compAlphas,'omitnan'); 

clear usedCount volumeAlpha alphaRadius timeRefinedAlpha numData addPoints aRadius shrinkType compAlphas paproskyType ...
    paproskyTypes paproskyNames varName timeRefinedAlphaUI


% Data storage of alpha data
alphaBytes = 0;
for i = 1:numel(pelvisDefect)
    tmp = pelvisDefect(i).volume.shrink.alpha;   
    alphaBytes = alphaBytes + whos('tmp').bytes; 
end
allDefect.volume.comp.alphaSize = alphaBytes;              
% Data storage of refined alpha data
refBytes = 0;
for i = 1:numel(pelvisDefect)
    tmp = pelvisDefect(i).volume.edge;             refBytes = refBytes + whos('tmp').bytes;
    tmp = pelvisDefect(i).volume.patches;          refBytes = refBytes + whos('tmp').bytes;
    tmp = pelvisDefect(i).volume.shrink.alphaCupBack; refBytes = refBytes + whos('tmp').bytes;
    tmp = pelvisDefect(i).volume.shrink.inside;    refBytes = refBytes + whos('tmp').bytes;
    tmp = pelvisDefect(i).volume.shrink.refinedAlpha; refBytes = refBytes + whos('tmp').bytes;
end
tmp = allPelvis.refPoints.gridPoints;             refBytes = refBytes + whos('tmp').bytes;
tmp = allDefect.boundaries.cuboid.hull;           refBytes = refBytes + whos('tmp').bytes;
allDefect.volume.comp.refinedAlphaSize = refBytes;           

clear alphaBytes refBytes tmp

%% Fill pelvis defect with points (refined alpha)
% For comparison: same point grid for all defects - use the point grid from filling the reference pelvis + defects with points
% Pre-filter with bounding box (time-saving)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Boundary box / cuboid hull of pelvis defects (refined alphaShape)
% For class intersect
% Pre-filter with bounding box (time-saving: less points for inpolyhedron)
shrinkType = 'refinedAlpha'; 
level = 1;
parfor i = 1:dataCountDefect 
    tic
    % Boundary box / cuboid: hull
    pelvisDefect(i).boundaries = pelvisDefect(i).boundaries.boundBoxHull(i,pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,level);
    timeBoundDefect(i,1) = toc;
end
allDefect.time.boundDefect = timeBoundDefect;
allDefect.time.boundDefectMean = mean(timeBoundDefect);
allDefect.time.boundDefectStd = std(timeBoundDefect);

% Fill pelvis defects (refined alpha shape) with points (incl. pre-filter)
parfor i = 1:dataCountDefect
    tic
    pelvisDefect(i).volume = pelvisDefect(i).volume.fillDefect(i,allPelvis.refPoints.gridPoints.pointsBox,... % Points of box around pelvis + all defects
        pelvisDefect(i).boundaries.cuboid.hull,... % Bounding Box
        pelvisDefect(i).volume.shrink.refinedAlpha.faces,pelvisDefect(i).volume.shrink.refinedAlpha.allVerticesPoints); % Mesh: refined alpha
    timeDefectPointsInside(i,1) = toc;
end 
allDefect.time.defectPointsInside = timeDefectPointsInside;
allDefect.time.defectPointsInsideMean = mean(timeDefectPointsInside);
allDefect.time.defectPointsInsideStd = std(timeDefectPointsInside);

clear shrinkType level timeBoundDefect timeDefectPointsInside

%% Save and load properties of Class Volume (Defect) 

% Save class volume data (properties)  
savePelvisDefectDataVolume = struct(); 
% Meta information of volume
metaVolumeDefect = metaclass(pelvisDefect(1).volume);
propertiesVolumeDefect = {metaVolumeDefect.PropertyList.Name};
for i = 1:dataCountDefect
    for j = 1:length(propertiesVolumeDefect)
        propertyName = propertiesVolumeDefect{j};
        % Cache
        savePelvisDefectDataVolume(i).(propertyName) = pelvisDefect(i).volume.(propertyName);
    end
end
% Save data
save('.\Workspace\pelvisDefectDataVolume.mat', 'savePelvisDefectDataVolume', '-v7.3'); % Adapt storage location
clear savePelvisDefectDataVolume metaVolumeDefect propertiesVolumeDefect propertyName

% Clear class volume (selected)
propertiesVolumeDefectToClear = {'edge'}; % Adapt
for i = 1:dataCountDefect
    % Empty fields in pelvisDefect(i).transform
    for j = 1:length(propertiesVolumeDefectToClear)
        propertyName = propertiesVolumeDefectToClear{j};
        fieldType = class(pelvisDefect(i).volume.(propertyName));
        switch fieldType
            case 'double'
                pelvisDefect(i).volume.(propertyName) = []; % Empty double
            case 'struct'
                pelvisDefect(i).volume.(propertyName) = struct(); % Empty struct
            case 'cell'
                pelvisDefect(i).volume.(propertyName) = {}; % Empty cell
            otherwise
                pelvisDefect(i).volume.(propertyName) = []; % Default to empty
        end
    end
end

% Clear class volume (selected)
propertiesVolumeDefectToClear  = {'alphaCupBack','comp'}; % Adapt
for i = 1:dataCount
    % Empty fields in pelvis(i).volume.shrink
    for j = 1:length(propertiesVolumeDefectToClear)
        propertyName = propertiesVolumeDefectToClear{j};
        if isfield(pelvisDefect(i).volume.shrink, propertyName) % Check if property exists
            fieldType = class(pelvisDefect(i).volume.shrink.(propertyName));
            switch fieldType
                case 'double'
                    pelvisDefect(i).volume.shrink.(propertyName) = []; % Empty double
                case 'struct'
                    pelvisDefect(i).volume.shrink.(propertyName) = struct(); % Empty struct
                case 'cell'
                    pelvisDefect(i).volume.shrink.(propertyName) = {}; % Empty cell
                otherwise
                    pelvisDefect(i).volume.shrink.(propertyName) = []; % Default to empty
            end
        end
    end
end

% Clear class volume (selected)
propertiesVolumeDefectToClear  = {'refVerticesInside','refVerticesInsideMask','refVerticesInsideIdx','refFacesInside','refPointsInside',...
    'refPointsInsideMask','refPointsInsideIdx','inoutVerticesThicken','verticesThicken','inoutVertices','outVertices','inVertices','inVerticesMask',...
    'inVerticesIdx','inFaces','inoutPoints','outPoints','inPoints','inPointsMask','inPointsIdx','executionTimeDensity','pointCloudDensity',...
    'pointCloudNormDensity','inoutDensity','filteredOutAll','filteredInAll','filteredInVertices','filteredOutVertices','filteredInPoints',...
    'filteredOutPoints','filteredInVerticesMask','filteredInVerticesIdx','filteredInFaces','filteredInPointsMask','filteredInPointsIdx','densityRGB',...
    'executionTimeCluster','clusterLabels','clusterMainLabel','otherClustersVertices','otherClustersPoints','outliersVertices','outliersPoints',...
    'mainClusterMask'}; % Adapt
for i = 1:dataCount
    % Empty fields in pelvis(i).volume.shrink
    for j = 1:length(propertiesVolumeDefectToClear)
        propertyName = propertiesVolumeDefectToClear {j};
        fieldType = class(pelvisDefect(i).volume.shrink.inside.(propertyName));
        switch fieldType
            case 'double'
                pelvisDefect(i).volume.shrink.inside.(propertyName) = []; % Empty double
            case 'struct'
                pelvisDefect(i).volume.shrink.inside.(propertyName) = struct(); % Empty struct
            case 'cell'
                pelvisDefect(i).volume.shrink.inside.(propertyName) = {}; % Empty cell
            otherwise
                pelvisDefect(i).volume.shrink.inside.(propertyName) = []; % Default to empty
        end
    end
end
clear propertiesVolumeDefectToClear propertyName fieldType

% Load class volume data (properties) %%%
loadPelvisDefectDataVolume = load('.\pelvisDefectDataVolumeMin.mat', 'savePelvisDefectDataVolume'); % Adapt storage location
%pelvisDefect(dataCountDefect) = Defect; % Initialisation
metaVolumeDefect = metaclass(pelvisDefect(1).volume);
propertiesVolumeDefect = {metaVolumeDefect.PropertyList.Name};
%propertiesTransform = {'edge', 'patches', 'shrink', 'gridPoints'}; % load selected properties
for i = 1:dataCountDefect
    for j = 1:length(propertiesVolumeDefect)
        propertyName = propertiesVolumeDefect{j};
        if isfield(loadPelvisDefectDataVolume.savePelvisDefectDataVolume(i), propertyName)
            pelvisDefect(i).volume.(propertyName) = ...
                loadPelvisDefectDataVolume.savePelvisDefectDataVolume(i).(propertyName);
        end
    end
end
clear loadPelvisDefectDataVolume metaVolumeDefect propertiesVolumeDefect propertyName


% Save data inside 
savePelvisDefectDataInside = struct();
for i = 1:dataCountDefect
    % Cache 
    savePelvisDefectDataInside(i).inside = pelvisDefect(i).volume.shrink.inside;
end
save('.\pelvisDefectDataInside.mat', 'savePelvisDefectDataInside', '-v7.3'); % Adapt storage location
clear savePelvisDefectDataInside

% Load data inside 
loadPelvisDefectDataInside = load('.\pelvisDefectDataInside.mat'); % Adapt storage location
savePelvisDefectDataInside = loadPelvisDefectDataInside.savePelvisDefectDataInside;
for i = 1:dataCountDefect
    %pelvisDefect(i).volume.shrink.inside = struct();
    pelvisDefect(i).volume.shrink.inside = savePelvisDefectDataInside(i).inside;
end
clear savePelvisDefectDataInside

% Save data isPoint/VertexInside
savePelvisDefectDataIsInside = struct();
for i = 1:dataCountDefect
    % Cache 
    savePelvisDefectDataIsInside(i).isVertexInside = pelvisDefect(i).volume.shrink.inside.isVertexInside;
    savePelvisDefectDataIsInside(i).isPointInside = pelvisDefect(i).volume.shrink.inside.isPointInside;
end
save('.\PelvisDefectDataIsInside.mat', 'savePelvisDefectDataIsInside', '-v7.3'); % Adapt storage location
clear savePelvisDefectDataIsInside

% Load data isPoint/VertexInside
loadPelvisDefectDataIsInside = load('.\PelvisDefectDataIsInside.mat', 'savePelvisDefectDataIsInside'); % Adapt storage location
savePelvisDefectDataIsInside = loadPelvisDefectDataIsInside.savePelvisDefectDataIsInside;
for i = 1:dataCountDefect
    %pelvisDefect(i).volume.shrink.inside = struct();
    pelvisDefect(i).volume.shrink.inside.isVertexInside = savePelvisDefectDataIsInside(i).isVertexInside;
    pelvisDefect(i).volume.shrink.inside.isPointInside = savePelvisDefectDataIsInside(i).isPointInside;
end
clear savePelvisDefectDataIsInside


% Boundary of refined alpha shape
% Save class boundary defect data (properties)
savePelvisDefectDataBound = struct();
% Meta information of boundary
metaBoundDefect = metaclass(pelvisDefect(1).boundaries);
propertiesBoundDefect = {metaBoundDefect.PropertyList.Name};
for i = 1:dataCountDefect
    for j = 1:length(propertiesBoundDefect)
        propertyName = propertiesBoundDefect{j};
        % Cache
        savePelvisDefectDataBound(i).(propertyName) = pelvisDefect(i).boundaries.(propertyName);
    end
end
% Save data
save('.\Workspace\pelvisDefectDataBound.mat', 'savePelvisDefectDataBound', '-v7.3'); % Adapt storage location
clear savePelvisDefectDataBound metaBoundDefect propertiesBoundDefect propertyName

% Load class boundary data (properties) %%%
loadPelvisDefectDataBound = load('.\pelvisDefectDataBound.mat', 'savePelvisDefectDataBound'); % Adapt storage location
%pelvisDefect(dataCount) = Defect; % Initialisation
metaBoundDefect = metaclass(pelvisDefect(1).boundaries);
%propertiesBoundDefect = {metaBoundDefect.PropertyList.Name}; % load all properties
propertiesBoundDefect = {'cuboid'}; % load selected properties
for i = 1:dataCountDefect
    for j = 1:length(propertiesBoundDefect)
        propertyName = propertiesBoundDefect{j};
        if isfield(loadPelvisDefectDataBound.savePelvisDefectDataBound(i), propertyName)
            pelvisDefect(i).boundaries.(propertyName) = loadPelvisDefectDataBound.savePelvisDefectDataBound(i).(propertyName);
        end
    end
end
clear loadPelvisDefectDataBound metaBoundDefect propertiesBoundDefect propertyName

%% Display scaled, transformed defect with boundary points of reference pelvis points/vertices inside (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%parfor i = 2:dataCountDefect
    i = 1; % Pelvis number
    %figure('Visible','off')   
    figure
    hold on

    % Pelvis defect remeshed (scaled and transformed)
    patch('Faces',pelvisDefect(i).volume.patches.all.comFaces,...
        'Vertices',pelvisDefect(i).volume.patches.all.comVertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...    % Edge color TUMcolors.grey50
        'EdgeAlpha',0.25);                  % Transparency of the edges
    light('Position', [1 1 5], 'Style', 'infinite');

    % Reference pelvis mesh inside
    patch('Faces',pelvisDefect(i).volume.shrink.inside.mainClusterFaces,...
        'Vertices',pelvis(1).transform.trafo.(weightKabsch{w}).vertices,... % pelvis(1).import.processed
        'FaceColor',[0.8 0.8 0.8], ...      % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',TUMcolors.grey50,...    % Edge color
        'EdgeAlpha',0.25);                  % Transparency of the edges

    % Main cluster: vertices
    plot3(pelvisDefect(i).volume.shrink.inside.mainClusterVertices(:,1), ...
        pelvisDefect(i).volume.shrink.inside.mainClusterVertices(:,2), ...
        pelvisDefect(i).volume.shrink.inside.mainClusterVertices(:,3), '.','Color', TUMcolors.blue300, 'MarkerSize', 2);
    % Main cluster: points
    plot3(pelvisDefect(i).volume.shrink.inside.mainClusterPoints(:,1), ...
        pelvisDefect(i).volume.shrink.inside.mainClusterPoints(:,2), ...
        pelvisDefect(i).volume.shrink.inside.mainClusterPoints(:,3), '.','Color', TUMcolors.blue300, 'MarkerSize', 2);
    % Boundary vertices/points
    plot3(pelvisDefect(i).volume.shrink.refinedAlpha.boundaryVerticesPoints(:,1), ...
    pelvisDefect(i).volume.shrink.refinedAlpha.boundaryVerticesPoints(:,2), ...
        pelvisDefect(i).volume.shrink.refinedAlpha.boundaryVerticesPoints(:,3), '.','Color', 'r', 'MarkerSize', 5);

    % Format and display properties
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['Boundary Point Cloud: Pelvis Defect ' num2str(i)]);
    legend('Pelvis Defect','Reference pelvis','Reference vertices inside', 'Reference points inside','Boundary vertices/points of reference');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')boundaryRefPointsInside.fig'])
%end

%% Display scaled, transformed and refined alphaShaped defect (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Loop with save figure to save the data/figures
%for i = 1:dataCountDefect
    i = 1; % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    
    % Pelvis defect alphaShape
    shrinkType = 'refinedAlpha'; % type: 'alpha, 'refinedAlpha', ...
    patch('Faces',pelvisDefect(i).volume.shrink.(shrinkType).faces,...
        'Vertices',pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,...
        'FaceColor',[0.8 0.8 0.8], ...    % Face color
        'FaceAlpha',0.7,...                  % Transparency of the faces % 0.7 to check whether mesh is inside % Figure Paper 0.5
        'EdgeColor',TUMcolors.grey80 ,...    % Edge color % Figure Paper TUMcolors.grey50
        'EdgeAlpha',0.25 ,...                % Transparency of the edges % Figure Paper 0.25
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Pelvis defect
    patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color [0.502 0.502 0.502]
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Display pelvis defect with all vertices/points (optional)
    % plot3(pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints(:,1), ...
    %     pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints(:,2), ...
    %     pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints(:,3),...
    %     '.','Color',[0.2 0.2 0.2],'MarkerSize', 1)

    % Display edge (acetabulum edge) (optional)
    % for k = 1:pelvisDefect(i).volume.patches.all.numComponents
    %     plot3(pelvisDefect(i).volume.edge.verticesLoops{k}(:,1),...
    %         pelvisDefect(i).volume.edge.verticesLoops{k}(:,2),...
    %         pelvisDefect(i).volume.edge.verticesLoops{k}(:,3),...
    %         '.','Color', 'r','MarkerSize',30) %TUMcolors.orange
    % end

    % Reference acetabulum centre
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
        '.','Color','m', 'MarkerSize', 80); % pelvis(i).transform.trafo.(weightKabsch{w}).acentre
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
        'x', 'MarkerEdgeColor', 'm', 'MarkerSize', 35, 'LineWidth', 10); % pelvis(i).transform.trafo.(weightKabsch{w}).acentre

    % Used vertices of original defect for alphaShape (optional)
    % plot3(pelvisDefect(i).volume.shrink.(shrinkType).usedDefectVertices(:,1), ...
    %     pelvisDefect(i).volume.shrink.(shrinkType).usedDefectVertices(:,2), ....
    %     pelvisDefect(i).volume.shrink.(shrinkType).usedDefectVertices(:,3), '.','MarkerSize',20, 'Color', 'g');

    % Format and display properties
    title(['Refined alpha shape: pelvis defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    view(80, -10) % Figure Paper
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')refinedAlpha.fig'])
%end

%% Display scaled, transformed and (refined) alpha Shaped defect with reference pelvis (for control)

% Loop with save figure to save the data/figures
%parfor i = 1:dataCountDefect
    i = 1;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on
    
    % Pelvis defect alphaShape
    shrinkType = 'alpha'; % type: 'alpha, 'refinedAlpha' ... %%%
    % Pelvis defect scaled and transformed
    % patch('Faces',pelvisDefect(i).volume.shrink.(shrinkType).faces,...
    %     'Vertices',pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,...
    %     'FaceColor',[0.9 0.75 0.68], ...    % Face color
    %     'FaceAlpha',1,...                   % Transparency of the faces
    %     'EdgeColor',TUMcolors.grey50,...    % Edge color
    %     'EdgeAlpha',0.25);                  % Transparency of the edges
    % Pelvis defect alphaShape (lightning)
    patch('Faces',pelvisDefect(i).volume.shrink.(shrinkType).faces,...
        'Vertices',pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,...
        'FaceColor', [0.9 0.75 0.68], ...    % Face color [0.9 0.75 0.68]
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',[0.4 0.4 0.4],...              % Edge color
        'EdgeAlpha',0.1,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);

    % Reference pelvis
    patch('Faces',pelvis(1).import.processed.faces,...
        'Vertices',pelvis(1).import.processed.vertices,...
        'FaceColor',TUMcolors.grey20, ...    % Face color
        'FaceAlpha',0.25,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...              % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);

    % Display pelvis defect with all vertices/points (optional)
    % plot3(pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints(:,1), ...
    %     pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints(:,2), ...
    %     pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints(:,3),...
    %     '.','Color',[0.2 0.2 0.2],'MarkerSize', 1)
    
    % Reference acetabulum centre
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
        '.','Color','m', 'MarkerSize', 32); % pelvis(i).transform.trafo.(weightKabsch{w}).acentre
    plot3(pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(1), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(2), ...
        pelvisDefect(i).transform.trafo.(weightKabsch{w}).acentre(3), ...
        'x', 'MarkerEdgeColor', 'm', 'MarkerSize', 14, 'LineWidth', 4); % pelvis(i).transform.trafo.(weightKabsch{w}).acentre
    
    % Format and display properties
    title(['Alpha Shape: Pelvis Defect ' num2str(i)]);
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(3);
    view(80, -10) % Figure Paper   
    %view(260, -10) % Figure Paper
    camlight; lighting gouraud;
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')RefAlphaShape.fig'])
%end
clear shrinkType

%% Display scaled, transformed defect with different alphaShape methods (for control)

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

%parfor i = 1:dataCountDefect
    i = 1;  % Pelvis number
    %figure('Visible','off')  
    figure; % Create a new figure

    % 1st subplot: Pelvis defect
    subplot(2, 2, 1);
    hold on;
    patch('Faces',pelvisDefect(i).transform.trafo.(weightKabsch{w}).faces,...
        'Vertices',pelvisDefect(i).transform.trafo.(weightKabsch{w}).vertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',... % Edge color [0.502 0.502 0.502]
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    % Reference acetabulum centre 
    plot3(pelvisDefect(1).import.processed.acentre(1), pelvisDefect(1).import.processed.acentre(2), ...
        pelvisDefect(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10);
    light('Position', [1 1 5], 'Style', 'infinite');
    title('Pelvis Defect');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]);
    view(3);
    hold off;

    % 2nd subplot: patch of pelvis defect + used vertices/points
    subplot(2, 2, 2);
    hold on;
    patch('Faces',pelvisDefect(i).volume.patches.all.comFaces,...
        'Vertices',pelvisDefect(i).volume.patches.all.comVertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',... % Edge color [0.502 0.502 0.502]
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');
    % Reference pelvis mesh inside (optional)
    patch('Faces',pelvisDefect(i).volume.shrink.inside.mainClusterFaces,...
        'Vertices',pelvis(1).transform.trafo.(weightKabsch{w}).vertices,... % pelvis(1).import.processed
        'FaceColor',[0.8 0.8 0.8], ...      % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',TUMcolors.grey50,...    % Edge color
        'EdgeAlpha',0.25);                  % Transparency of the edges
    % Reference points/vertices inside
    verticesPointsInside = [pelvisDefect(i).volume.shrink.inside.mainClusterVertices; pelvisDefect(i).volume.shrink.inside.mainClusterPoints];
    plot3(verticesPointsInside(:,1), ...
        verticesPointsInside(:,2), ...
        verticesPointsInside(:,3), '.','Color', TUMcolors.blue300, 'MarkerSize', 1);
    title('Pelvis Defect + Reference inside');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]);
    view(3);
    hold off;

    % 3rd subplot: alphaShape of pelvis defect 
    shrinkType = 'alpha'; 
    subplot(2, 2, 3);
    hold on;
    patch('Faces',pelvisDefect(i).volume.shrink.(shrinkType).faces,...
        'Vertices',pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,...
        'FaceColor',[0.8 0.8 0.8], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');
    % Alpha shape: used vertices/points
    plot3(pelvisDefect(i).volume.shrink.(shrinkType).usedVerticesPoints(:,1), ...
        pelvisDefect(i).volume.shrink.(shrinkType).usedVerticesPoints(:,2), ...
        pelvisDefect(i).volume.shrink.(shrinkType).usedVerticesPoints(:,3), '.','Color',[0.2 0.2 0.2], 'MarkerSize', 1);
    title('AlphaShape of Pelvis Defect');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]);
    view(3);
    hold off;

    % 4th subplot: refined alphaShape of pelvis defect
    shrinkType = 'refinedAlpha'; 
    subplot(2, 2, 4);
    hold on;
    patch('Faces',pelvisDefect(i).volume.shrink.(shrinkType).faces,...
        'Vertices',pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,...
        'FaceColor',[0.8 0.8 0.8], ...    % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');
    % Reference inside: used vertices/points
    plot3(pelvisDefect(i).volume.shrink.(shrinkType).usedVerticesPoints(:,1), ...
        pelvisDefect(i).volume.shrink.(shrinkType).usedVerticesPoints(:,2), ...
        pelvisDefect(i).volume.shrink.(shrinkType).usedVerticesPoints(:,3), '.','Color',[0.2 0.2 0.2], 'MarkerSize', 1);
    title('AlphaShape with Vertices');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    daspect([1, 1, 1]);
    view(3);
    hold off;

    % Set overall figure properties
    sgtitle(['Volume - Pelvis defect ' num2str(i)]);

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')alphaShapes.fig'])
%end

clear verticesPointsInside

%% Display alpha shape to refined alphashape with vertex-to-nearest-neighbour color (for control)

% Loop with save figure to save the data/figures
%parfor i = 2:dataCount
    i = 2;  % Pelvis number
    %figure('Visible','off')   
    figure
    hold on

    % alpha shape
    shrinkType = 'alpha'; 
    patch('Faces',pelvisDefect(i).volume.shrink.(shrinkType).faces,...
        'Vertices',pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,...
        'FaceColor',[0.8 0.8 0.8], ...    % Face color
        'FaceAlpha',0.5,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');
    % Used vertices and points
    plot3(pelvisDefect(i).volume.shrink.(shrinkType).usedVerticesPoints(:,1), ...
        pelvisDefect(i).volume.shrink.(shrinkType).usedVerticesPoints(:,2), ...
        pelvisDefect(i).volume.shrink.(shrinkType).usedVerticesPoints(:,3),...
        '.','Color',[0.2 0.2 0.2],'MarkerSize', 1);

    % Refined alpha shape
    shrinkType = 'refinedAlpha';
    patch('Faces',pelvisDefect(i).volume.shrink.(shrinkType).faces,...
        'Vertices',pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,...
        'FaceVertexCData', pelvisDefect(i).volume.shrink.comp.nearVertexFaceColour, ... % Face color
        'FaceColor','flat', ...    
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');
    % Used vertices and points
    plot3(pelvisDefect(i).volume.shrink.(shrinkType).usedVerticesPoints(:,1), ...
        pelvisDefect(i).volume.shrink.(shrinkType).usedVerticesPoints(:,2), ...
        pelvisDefect(i).volume.shrink.(shrinkType).usedVerticesPoints(:,3),...
        '.','Color',[0.2 0.2 0.2],'MarkerSize', 1);

    % Reference acetabulum centre 
    plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
        pelvis(1).import.processed.acentre(3), '.','Color','r', 'MarkerSize', 30); 
    
    % Format and display properties
    title(['Alpha vs. refinedAlpha (Vertex-to-Nearest-Neighbour): Pelvis Defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    legend('Alpha Shape', 'Used Vertices/Points','Refined Alpha Shape','Used Vertices/Points', 'Acentre');
    colormap('viridis');
    colorbar
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;
    
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisDefect(',num2str(i),')alphaNearVertex.fig'])
%end

%% Display scaled, transformed defect (refined alpha) with boundary box (convex hull) and its points inside (for control)

% Loop with save figure to save the data/figures
for i = 1:dataCountDefect
    %i = 2;  % Pelvis number
    figure('Visible','off')
    %figure; 
    hold on

    % Pelvis defect alphaShape
    shrinkType = 'refinedAlpha'; % type: 'alpha, 'refinedAlpha', ...
    h_defect = patch('Faces',pelvisDefect(i).volume.shrink.(shrinkType).faces,...
        'Vertices',pelvisDefect(i).volume.shrink.(shrinkType).allVerticesPoints,...
        'FaceColor',[0.8 0.8 0.8], ...    % Face color
        'FaceAlpha',0.5,...                   % Transparency of the faces
        'EdgeColor','none',...              % Edge color
        'EdgeAlpha',0.25,...                % Transparency of the edges
        ... % Ligthing for 3d effect
        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
        'AmbientStrength', 0.5);
    light('Position', [1 1 5], 'Style', 'infinite');

    % Reference acetabulum centre
   h_acentre = plot3(pelvisDefect(1).import.processed.acentre(1), pelvisDefect(1).import.processed.acentre(2), ...
        pelvisDefect(1).import.processed.acentre(3), '*','Color',TUMcolors.green, 'MarkerSize', 10);

    % Display bounding box
    % Define the 8 corners of the bounding box
    h_box = trisurf(pelvisDefect(i).boundaries.cuboid.hull.tri,...
        pelvisDefect(i).boundaries.cuboid.hull.cornerpoints(:,1),...
        pelvisDefect(i).boundaries.cuboid.hull.cornerpoints(:,2),...
        pelvisDefect(i).boundaries.cuboid.hull.cornerpoints(:,3),...
        'FaceColor',TUMcolors.blue300,'EdgeColor',TUMcolors.blue300,'FaceAlpha',0.25);

    % Display cuboid edges
    for k=1:3
        h_edge = quiver3(pelvisDefect(i).boundaries.cuboid.hull.cornerpoints(1,1),...
            pelvisDefect(i).boundaries.cuboid.hull.cornerpoints(1,2),...
            pelvisDefect(i).boundaries.cuboid.hull.cornerpoints(1,3),... % start point
            pelvisDefect(i).boundaries.cuboid.hull.edgeVector(k,1),...
            pelvisDefect(i).boundaries.cuboid.hull.edgeVector(k,2),...
            pelvisDefect(i).boundaries.cuboid.hull.edgeVector(k,3),0,'LineWidth',2);
    end

    % Plot points inside bounding box (optional)
    % plot3(pelvisDefect(i).volume.gridPoints.boxInside(:,1), ...
    %     pelvisDefect(i).volume.gridPoints.boxInside(:,2), ...
    %     pelvisDefect(i).volume.gridPoints.boxInside(:,3), '*','Color',TUMcolors.green, 'MarkerSize', 0.25);

    % Plot points inside defect (max bounding box)
    h_insideDefect = plot3(pelvisDefect(i).volume.gridPoints.inside(:,1), ...
        pelvisDefect(i).volume.gridPoints.inside(:,2), ...
        pelvisDefect(i).volume.gridPoints.inside(:,3), '*','Color',TUMcolors.orange, 'MarkerSize', 0.75);

    % Format and display properties
    title(['Filled with Points: Pelvis Defect ' num2str(i)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    legend([h_defect, h_acentre, h_box, h_edge, h_insideDefect], {'Refined Alpha Shape', 'Acentre', 'Bounding Box', 'Cuboid Edge', 'Points Inside Defect'});
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
    savefig(['./Figures/PelvisDefect(',num2str(i),')boundBoxPoints.fig'])
end

clear  h_defect h_acentre h_box h_insideDefect h_edge

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Class DefectAnalysis %%%%%%%%%%%%

% Pelvis defect categories (all or selected) are base for intersection
baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
numBaseIntersect = length(baseIntersect);

% Initalization: Generate object array
pelvisDefectAnalysis(numBaseIntersect,1) = DefectAnalysis;

%% %%%%%%%%%%%% Class Intersect %%%%%%%%%%%%

% Pelvis defect categories (all or selected) are base for intersection
baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
numBaseIntersect = length(baseIntersect);

% Intersection types
intersectType = {'alpha','refinedAlpha','inside','grid'};
numIntersectTypes = length(intersectType);

for j = 1:numIntersectTypes
    for i = 1:numBaseIntersect
        % Base of intersections
        pelvisDefectAnalysis(i).intersect.(intersectType{j}).base = baseIntersect{i};
    end
        % Pelvis IDs for the categories
        pelvisDefectAnalysis(1).intersect.(intersectType{j}).IDs = (2:dataCountDefect)'; % all
        pelvisDefectAnalysis(2).intersect.(intersectType{j}).IDs = [2;3;(5:9)';11;12;14;15;(17:21)';23;24;(26:37)';39;40;(42:47)']; % BTE
        pelvisDefectAnalysis(3).intersect.(intersectType{j}).IDs = [4;10;13;16;22;25;38;41]; % Other
        pelvisDefectAnalysis(4).intersect.(intersectType{j}).IDs = [3;4;5]; % 2A
        pelvisDefectAnalysis(5).intersect.(intersectType{j}).IDs = [2;13;14;17;19;25;38;42;45]; % 2B
        pelvisDefectAnalysis(6).intersect.(intersectType{j}).IDs = [11;40;44]; % 2C
        pelvisDefectAnalysis(7).intersect.(intersectType{j}).IDs = [6;12;20;21;22;26;28;29;30;31;32;34;35;36;39;47;48]; % 3A
        pelvisDefectAnalysis(8).intersect.(intersectType{j}).IDs = [7;8;9;10;15;16;18;23;24;27;33;37;41;43;46]; % 3B
        pelvisDefectAnalysis(9).intersect.(intersectType{j}).IDs = [3;5]; % 2A & BTE
        pelvisDefectAnalysis(10).intersect.(intersectType{j}).IDs = [2;14;17;19;42;45]; % 2B & BTE
        pelvisDefectAnalysis(11).intersect.(intersectType{j}).IDs = [11;40;44]; % 2C & BTE
        pelvisDefectAnalysis(12).intersect.(intersectType{j}).IDs = [6;12;20;21;26;28;29;30;31;32;34;35;36;39;47;48]; % 3A & BTE
        pelvisDefectAnalysis(13).intersect.(intersectType{j}).IDs = [7;8;9;15;18;23;24;27;33;37;43;46]; % 3B & BTE
end

% Intersections: Pairs (calculated), 5%, 10%, 25%, 33.33%, 50%, 66.67%, 75%, 100%
intersectPercent = [NaN; 5; 10; 25; 33.33; 50; 66.67; 75; 90; 95; 100]; % percent
intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
numIntersectPct = length(intersectPercent);
for i = 1:numBaseIntersect
    for j = 1:numIntersectTypes
        for k = 1:numIntersectPct
            pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}) = struct();
            if k == 1
                pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{1}).percent = ...
                    (2 / size(pelvisDefectAnalysis(i).intersect.(intersectType{j}).IDs,1)) * 100;
            else
                pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).percent = intersectPercent(k);
            end
            pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).count = ...
                round(pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).percent / 100 * ...
                size(pelvisDefectAnalysis(i).intersect.(intersectType{j}).IDs,1));
        end
    end
end

% Store data size
for j = 1:numIntersectTypes
    sizeIntersect = pelvisDefectAnalysis(1).intersect.(intersectType{j});
    sizeInfo = whos('sizeIntersect');
    allIntersect.(intersectType{j}).sizeBaseInfos(1) = sizeInfo.bytes;
end
clear sizeIntersect sizeInfo 

%% Intersection: pairwise intersection of boundary box of defect and defect vertices/points (pre fitler) 
% Pairwise intersection for all defects - pelvisIntersect(1)
% Intersection base: 'all'

weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Defect pairs of pelvis defect categorie 'all' 
for j = 1:numIntersectTypes
    typeName = intersectType{j};
    pairs = pelvisDefectAnalysis(1).intersect.(typeName).pairs;
    IDsFor = nchoosek(2:dataCountDefect, 2); % Forward combinations (excluding reference 1)
    IDsRev = [IDsFor(:, 2), IDsFor(:, 1)]; % Reverse combinations
    pairs.IDs = sortrows([IDsFor; IDsRev], [1,2]); % Combined and sorted
    % Update the pelvisIntersect structure
    pelvisDefectAnalysis(1).intersect.(typeName).pairs = pairs;
end
clear IDsFor IDsRev pairs

% Pairwise intersection of boundary box of defect and defect points (temporary object)
for j = 1:numIntersectTypes % intersectType
    typeName = intersectType{j};
    [pairRow, ~] = size(pelvisDefectAnalysis(1).intersect.(typeName).pairs.IDs); % same for all intersect types
    timeIntersectBoxDefect = NaN(pairRow, 1);

    for i = 1:pairRow
        tic
        % Extract defect points (u) and box (v)
        u = pelvisDefectAnalysis(1).intersect.(typeName).pairs.IDs(i, 1); % Defect points
        v = pelvisDefectAnalysis(1).intersect.(typeName).pairs.IDs(i, 2); % Box

        % Input Points and Masks Setup
        switch typeName
            case {'alpha', 'refinedAlpha'}
                if strcmp(typeName, 'alpha')
                    % Intersect type 'alpha':
                    % points = defect vertices; volume = alpha shape
                    inputPoints = pelvisDefect(u).transform.trafo.(weightKabsch{w}).vertices;   % alpha vertices (original defect vertices)
                else
                    % Intersect type 'refinedAlpha':
                    % points = defect vertices + refined patches; volume = refined alpha shape
                    inputPoints = pelvisDefect(u).volume.patches.all.comVertices;               % refined alpha vertices (comVertices with patches)
                end
                inputPointNum = size(inputPoints, 1);
                inputMasksIdx = (1:inputPointNum)';
                inputBaseNum = inputPointNum;
            case 'inside'
                % Intersect type 'inside':
                % points = defect vertices + refined patches + reference vertices&points inside ; volume = refined alpha shape
                inputPoints = [pelvisDefect(u).volume.patches.all.comVertices; ...              % refined alpha vertices (comVertices with patches)
                                pelvisDefect(u).volume.shrink.inside.mainClusterVertices; ...   % reference pelvis vertices inside defect
                                pelvisDefect(u).volume.shrink.inside.mainClusterPoints];        % reference pelvis points inside defect
                inputPointNum = [size(pelvisDefect(u).volume.patches.all.comVertices, 1); ...   % point num = filtered vertices/point of the base point set
                                size(pelvisDefect(u).volume.shrink.inside.mainClusterVertices, 1); ...
                                size(pelvisDefect(u).volume.shrink.inside.mainClusterPoints, 1)];
                inputMasksIdx = [(1:inputPointNum(1))'; ...
                                pelvisDefect(u).volume.shrink.inside.mainClusterVerticesIdx; ...
                                pelvisDefect(u).volume.shrink.inside.mainClusterPointsIdx];
                inputBaseNum = [size(pelvisDefect(u).volume.patches.all.comVertices,1); ...                 % mask num (logical mask) = all vertices/points of the base point set
                                size(pelvisDefect(u).volume.shrink.inside.mainClusterVerticesMask,1); ...   % vertices of reference pelvis
                                size(pelvisDefect(u).volume.shrink.inside.mainClusterPointsMask,1)];        % grid points
            case 'grid'
                % Intersect type 'grid':
                % points = defect vertices + refined patches + defect (refined alpha) filled with grid points; volume = refined alpha shape
                inputPoints = [pelvisDefect(u).volume.patches.all.comVertices; ...              % refined alpha vertices (comVertices with patches)
                                pelvisDefect(u).volume.shrink.inside.mainClusterVertices; ...   % reference pelvis vertices inside defect
                                pelvisDefect(u).volume.gridPoints.inside];                      % grid points inside refined alpha (includes mainClusterPoints)
                inputPointNum = [size(pelvisDefect(u).volume.patches.all.comVertices, 1); ...   % point num = filtered vertices/point base point set
                                size(pelvisDefect(u).volume.shrink.inside.mainClusterVertices, 1); ...
                                size(pelvisDefect(u).volume.gridPoints.inside, 1)];
                inputMasksIdx = [(1:inputPointNum(1))'; ...
                                pelvisDefect(u).volume.shrink.inside.mainClusterVerticesIdx; ...
                                pelvisDefect(u).volume.gridPoints.insideMaskIdx];
                inputBaseNum = [size(pelvisDefect(u).volume.patches.all.comVertices,1); ...                 % mask num (logical mask) = all vertices/points of the base point set
                                size(pelvisDefect(u).volume.shrink.inside.mainClusterVerticesMask,1); ...   % vertices of reference pelvis
                                size(pelvisDefect(u).volume.gridPoints.insideMask,1)];                      % grid points
        end

        % Perform intersection
        tempIntersect = Intersection();
        tempIntersect = tempIntersect.intersectBox(typeName, u, v, pelvisDefect(v).boundaries.cuboid.hull, ... % Boundary box of refinedAlpha volume
            inputPoints, inputPointNum, inputMasksIdx, inputBaseNum);

        % Assign the collected results to the main objects
        % All points/vertices inside
        %pelvisIntersect(1).(typeName).pairs.boxInside{i,1} = ...
        % tempIntersect.(typeName).pairs.boxInside;                         %%% to save memory
        % Defect mesh
        pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxDefectInsideIdx{i,1} = ...
            tempIntersect.(typeName).pairs.boxDefectInsideMaskIdx;
        %pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxDefectMask{i,1} = ...
        % tempIntersect.(typeName).pairs.boxDefectInsideMask;               %%% to save memory
        allIntersect.(typeName).numDefectBoxDefect(i,1) = ...
            length(tempIntersect.(typeName).pairs.boxDefectInsideMaskIdx);
        if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
            % Reference points inside defect (grid points)
            % inside: filtered alphaCupBack shape with reference points inside;
            % grid: filtered alphaCupBack shape + refinedAlpha filled with grid points
            % Vertices
            pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxVerticesInsideIdx{i,1} = ...
                tempIntersect.(typeName).pairs.boxVerticesInsideMaskIdx;
            %pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxVerticesMask{i,1} = ...
            % tempIntersect.(typeName).pairs.boxVerticesInsideMask;         %%% to save memory
            allIntersect.(typeName).numVerticesBoxDefect(i,1) = ...
                length(tempIntersect.(typeName).pairs.boxVerticesInsideMaskIdx);
            % Points
            pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxPointsInsideIdx{i,1} = ...
                tempIntersect.(typeName).pairs.boxPointsInsideMaskIdx;
            %pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxPointsMask{i,1} = ...
            % tempIntersect.(typeName).pairs.boxPointsInsideMask;           %%% to save memory
            allIntersect.(typeName).numPointsBoxDefect(i,1) = ...
                length(tempIntersect.(typeName).pairs.boxPointsInsideMaskIdx);
        end
        timeIntersectBoxDefect(i, 1) = toc;
        % Store data size
        sizeIntersectDefect = pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxDefectInsideIdx(i);
        sizeInfoDefect = whos('sizeIntersectDefect');
        allIntersect.(typeName).sizePairBoxDefect(i,1) = sizeInfoDefect.bytes;
        if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
            sizeIntersectVertices = pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxVerticesInsideIdx(i);
            sizeInfoVertices = whos('sizeIntersectVertices');
            allIntersect.(typeName).sizePairBoxVertices(i,1) = sizeInfoVertices.bytes;
            sizeIntersectPoints = pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxPointsInsideIdx(i);
            sizeInfoPoints = whos('sizeIntersectPoints');
            allIntersect.(typeName).sizePairBoxPoints(i,1) = sizeInfoPoints.bytes;
        end
    end

    % Store computation time in allIntersect
    allIntersect.time.(['intersect' typeName 'BoxDefect']) = timeIntersectBoxDefect;
    allIntersect.time.(['intersect' typeName 'BoxDefectMean']) = mean(timeIntersectBoxDefect, 'omitnan');
    allIntersect.time.(['intersect' typeName 'BoxDefectStd']) = std(timeIntersectBoxDefect, 'omitnan');
    % Store data size
    sizeIntersect = pelvisDefectAnalysis(1).intersect.(typeName).pairs;
    sizeInfo = whos('sizeIntersect');
    allIntersect.(typeName).('sizeBoxDefect')(1) = sizeInfo.bytes;
   
    % Clear temporary variables
    clear inputPoints inputPointNum inputMasksIdx inputBaseNum tempIntersect sizeIntersect sizeInfo timeIntersectBoxDefect u v 
    clear sizeIntersect sizeIntersectDefect sizeIntersectVertices sizeIntersectPoints sizeInfo sizeInfoDefect ...
        sizeInfoVertices sizeInfoPoints
end

%% Intersection: pairwise intersection of defects (inpolyhedron) 
% Intersection base: 'all'

% weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting

% Pairwise intersection of defect volume and defect points
for j = 1:numIntersectTypes % intersectType
    typeName = intersectType{j};
    [pairRow, ~] = size(pelvisDefectAnalysis(1).intersect.(intersectType{j}).pairs.IDs); % same for all intersect types
    timeIntersectDefect(1:pairRow,1) = NaN;

    for i = 1:pairRow
        tic
        u = pelvisDefectAnalysis(1).intersect.(typeName).pairs.IDs(i,1); % Defect point
        v = pelvisDefectAnalysis(1).intersect.(typeName).pairs.IDs(i,2); % Defect volume

        % Input points and masks setup
        switch typeName
            case {'alpha', 'refinedAlpha'}
                if strcmp(typeName, 'alpha')
                    % Intersect type 'alpha': points = defect vertices; volume = alpha shape
                    pointBase = pelvisDefect(u).transform.trafo.(weightKabsch{w}).vertices;        % alpha vertices (original defect vertices)
                else % 'refinedAlpha'
                    % Intersect type 'refined Alpha': points = defect vertices + refined patches; volume = refined alpha shape
                    pointBase = pelvisDefect(u).volume.patches.all.comVertices;                    % refined alpha vertices (comVertices with patches)
                end
                inputFiltered = pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxDefectInsideIdx{i};         % filtered input points (defect-box)
                inputFilteredNums = size(inputFiltered,1);
                inputBaseNum = size(pointBase,1);
                inputVolume.faces = pelvisDefect(v).volume.shrink.(typeName).faces;
                inputVolume.vertices = pelvisDefect(v).volume.shrink.(typeName).allVerticesPoints;
            case {'inside', 'grid'}
                % Intersect type 'inside': points = defect vertices + refined patches + reference vertices&points inside ; volume = refined alpha shape
                % Intersect type 'grid': points = defect vertices + refined patches + defect (refined alpha) filled with grid points; volume = refined alpha shape
                % Filtered input points (defect-box)
                inputFiltered = [pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxDefectInsideIdx{i};...      % Mask: defect vertices + refined patches
                    pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxVerticesInsideIdx{i};...                 % Mask: vertices of reference pelvis
                    pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxPointsInsideIdx{i}];                     % Mask: grid points
                inputFilteredNums = [size(pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxDefectInsideIdx{i},1);...
                    size(pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxVerticesInsideIdx{i},1);...
                    size(pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxPointsInsideIdx{i},1)];
                % Point base
                pointBase = [pelvisDefect(u).volume.patches.all.comVertices; ...                    % refined alpha vertices (comVertices with patches)
                    pelvis(1).import.processed.vertices; ...                                        % vertices of reference pelvis
                    allPelvis.refPoints.gridPoints.pointsBox];                                      % reference pelvis points (base gridPoints)
                % Num of all vertices/points of the base point set
                inputBaseNum = [size(pelvisDefect(u).volume.patches.all.comVertices,1); ...
                    size(pelvis(1).import.processed.vertices,1); ...
                    size(allPelvis.refPoints.gridPoints.pointsBox,1)];
                % Volume
                inputVolume.faces = pelvisDefect(v).volume.shrink.refinedAlpha.faces;
                inputVolume.vertices = pelvisDefect(v).volume.shrink.refinedAlpha.allVerticesPoints;
        end

        % Perform intersection
        tempIntersect = Intersection();
        % Pairwise intersection
        tempIntersect = tempIntersect.intersectPairs(typeName, u, v, ...
            inputFiltered, inputFilteredNums, pointBase, inputBaseNum, inputVolume);

        % Assign the collected results to the main objects
        % All points/vertices inside
        %pelvisDefectAnalysis(1).intersect.(typeName).pairs.inside{i,1} = ...
        % tempIntersect.(typeName).pairs.inside;                         %%% to save memory
        % Defect mesh
        pelvisDefectAnalysis(1).intersect.(typeName).pairs.defectInsideIdx{i,1} = ...
            tempIntersect.(typeName).pairs.defectInsideMaskIdx;
        %pelvisDefectAnalysis(1).intersect.(typeName).pairs.defectMask{i,1} = ...
        % tempIntersect.(typeName).pairs.defectInsideMask;              %%% to save memory
        allIntersect.(typeName).numDefectPairs(i,1) = ...
            length(tempIntersect.(typeName).pairs.defectInsideMaskIdx);
        if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
            % Reference points inside defect (grid points)
            % inside: filtered alphaCupBack shape with reference points inside;
            % grid: filtered alphaCupBack shape + refinedAlpha filled with grid points
            pelvisDefectAnalysis(1).intersect.(typeName).pairs.verticesInsideIdx{i,1} = ...
                tempIntersect.(typeName).pairs.verticesInsideMaskIdx;
            %pelvisDefectAnalysis(1).intersect.(typeName).pairs.verticesMask{i,1} = ...
            % tempIntersect.(typeName).pairs.verticesInsideMask;        %%% to save memory
            allIntersect.(typeName).numVerticesPairs(i,1) = ...
                length(tempIntersect.(typeName).pairs.verticesInsideMaskIdx);
            pelvisDefectAnalysis(1).intersect.(typeName).pairs.pointsInsideIdx{i,1} = ...
                tempIntersect.(typeName).pairs.pointsInsideMaskIdx;
            %pelvisDefectAnalysis(1).intersect.(typeName).pairs.pointsMask{i,1} = ...
            % tempIntersect.(typeName).pairs.pointsInsideMask;          %%% to save memory
            allIntersect.(typeName).numPointsPairs(i,1) = ...
                length(tempIntersect.(typeName).pairs.pointsInsideMaskIdx);
        end
        timeIntersectDefect(i,1) = toc;
        % Store data size
        sizeIntersectDefect = pelvisDefectAnalysis(1).intersect.(typeName).pairs.defectInsideIdx(i);
        sizeInfoDefect = whos('sizeIntersectDefect');
        allIntersect.(typeName).sizePairDefect(i,1) = sizeInfoDefect.bytes;
        if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
            sizeIntersectVertices = pelvisDefectAnalysis(1).intersect.(typeName).pairs.verticesInsideIdx(i);
            sizeInfoVertices = whos('sizeIntersectVertices');
            allIntersect.(typeName).sizePairVertices(i,1) = sizeInfoVertices.bytes;
            sizeIntersectPoints = pelvisDefectAnalysis(1).intersect.(typeName).pairs.pointsInsideIdx(i);
            sizeInfoPoints = whos('sizeIntersectPoints');
            allIntersect.(typeName).sizePairPoints(i,1) = sizeInfoPoints.bytes;
        end
    end
    % Computation time
    allIntersect.time.(['intersect' (typeName) 'DefectPairs']) = timeIntersectDefect;
    allIntersect.time.(['intersect' (typeName) 'DefectPairsMean']) = mean(timeIntersectDefect);
    allIntersect.time.(['intersect' (typeName) 'DefectPairsStd']) = std(timeIntersectDefect);
    % Data storage of intersection data
    % Size of pelvisIntersect(1).(typeName).pairs variable
    sizeIntersect = pelvisDefectAnalysis(1).intersect.(typeName).pairs; % cache
    sizeInfo = whos('sizeIntersect');
    allIntersect.(typeName).sizePairDefectTotal(1) = sizeInfo.bytes;

    clear pairRow tempIntersect inputFiltered inputFilteredNums pointBase inputBaseNum inputVolume sizeIntersect ...
        sizeInfo u v timeIntersectDefect
    clear sizeIntersect sizeIntersectDefect sizeIntersectVertices sizeIntersectPoints sizeInfo sizeInfoDefect ...
        sizeInfoVertices sizeInfoPoints
end

%% Display pairwise intersection with boundary box and points inside/outside (for control)
% Intersection base: 'all'

% Intersection type
% intersectType = {'alpha','refinedAlpha','inside','grid'};
% parfor j = 1:length(intersectType)
    j = 1; %%%
    typeName = intersectType{j};

    % figure('Visible','off')   
    figure
    hold on
    h = [];  lab = {}; % for legend

    % Pair number -> rows
    i = 1; % pairRow %%%
    u = pelvisDefectAnalysis(1).intersect.(typeName).pairs.IDs(i,1); % Defect point
    v = pelvisDefectAnalysis(1).intersect.(typeName).pairs.IDs(i,2); % Defect Box

    switch typeName
        case {'alpha', 'refinedAlpha'}
            if strcmp(typeName, 'alpha')
                inputAll = pelvisDefect(u).transform.trafo.(weightKabsch{w}).vertices;
                defectPointsFaces = pelvisDefect(u).transform.trafo.(weightKabsch{w}).faces;
                defectPointsVertices = pelvisDefect(u).transform.trafo.(weightKabsch{w}).vertices;
                defectBoxFaces = pelvisDefect(v).transform.trafo.(weightKabsch{w}).faces;
                defectBoxVertices = pelvisDefect(v).transform.trafo.(weightKabsch{w}).vertices;
            else % 'refinedAlpha'
                inputAll = pelvisDefect(u).volume.patches.all.comVertices;
                defectPointsFaces = pelvisDefect(u).volume.patches.all.comFaces;
                defectPointsVertices = pelvisDefect(u).volume.patches.all.comVertices;
                defectBoxFaces = pelvisDefect(v).volume.patches.all.comFaces;
                defectBoxVertices = pelvisDefect(v).volume.patches.all.comVertices;
            end
            boxInside = inputAll(pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxDefectInsideIdx{i}, :);
        case {'inside', 'grid'}
            inputDefect = pelvisDefect(u).volume.patches.all.comVertices;
            inputVertices = pelvisDefect(u).volume.shrink.inside.mainClusterVertices;
            if strcmp(typeName, 'inside')
                inputPoints = pelvisDefect(u).volume.shrink.inside.mainClusterPoints;
            else % 'grid'
                inputPoints = pelvisDefect(u).volume.gridPoints.inside;
            end
            inputAll = [inputDefect; inputVertices; inputPoints];
            boxDefectInside = inputDefect(pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxDefectInsideIdx{i}, :);
            boxVerticesInside = pelvis(1).import.processed.vertices(pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxVerticesInsideIdx{i}, :);
            boxPointsInside = allPelvis.refPoints.gridPoints.pointsBox(pelvisDefectAnalysis(1).intersect.(typeName).pairs.boxPointsInsideIdx{i}, :);
            boxInside = [boxDefectInside; boxVerticesInside; boxPointsInside];
            defectPointsFaces = pelvisDefect(u).volume.patches.all.comFaces;
            defectPointsVertices = pelvisDefect(u).volume.patches.all.comVertices;
            defectBoxFaces = pelvisDefect(v).volume.patches.all.comFaces;
            defectBoxVertices = pelvisDefect(v).volume.patches.all.comVertices;
    end

    % Defect of points (optional)
    h(end+1) = patch('Faces',defectPointsFaces,...
            'Vertices',defectPointsVertices,...
            'FaceColor',TUMcolors.grey50, ...   % Face color
            'FaceAlpha',0.75,...                % Transparency of the faces
            'EdgeColor','none',...              % Edge color
            'EdgeAlpha',0.25);                  % Transparency of the edges
    lab{end+1} = 'Defect of Points';
    % Defect reference inside
    if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
        h(end+1) = patch('Faces',pelvisDefect(u).volume.shrink.inside.mainClusterFaces,...
            'Vertices',pelvis(1).import.processed.vertices,... % pelvis(1).transform.trafo.(weightKabsch{w}).vertices
            'FaceColor',[0.8 0.8 0.8], ...   % Face color
            'FaceAlpha',0.75,...             % Transparency of the faces % Figure Paper 0.25
            'EdgeColor','none',...           % Edge color
            'EdgeAlpha',0.25);               % Transparency of the edges
        lab{end+1} = 'Reference within Defect of Points';
    end
    % Input points (all)
    % h(end+1) = plot3(inputAll(:,1), inputAll(:,2), inputAll(:,3), '.','Color','r', 'MarkerSize', 2.5);  % TUMcolors.orange  
    % lab{end+1} = 'Points Outside';
    % Points outside
    [~, ia] = setdiff(round(inputAll,6), round(boxInside,6), 'rows'); 
    boxOutside = inputAll(ia, :);
    h(end+1) = plot3(boxOutside(:,1), boxOutside(:,2), boxOutside(:,3), '.', 'Color', 'r', 'MarkerSize', 2.5);
    lab{end+1} = 'Points Outside';
    % Points inside
    h(end+1) = plot3(boxInside(:,1), boxInside(:,2), boxInside(:,3), '.','Color','g', 'MarkerSize', 2.5); % TUMcolors.green
    lab{end+1} = 'Points Inside';
    % h(end+1) = scatter3(boxInside(:,1), boxInside(:,2), boxInside(:,3), 10, 'g', 'filled', 'MarkerFaceAlpha', 1);
    % lab{end+1} = 'Points Inside'; % scatter with transparency

    % Defect of box (optional)
    h(end+1) = patch('Faces',defectBoxFaces,...
        'Vertices',defectBoxVertices,...
        'FaceColor',[0.9 0.75 0.68], ...   % Face color 
        'FaceAlpha',0.75,...               % Transparency of the faces
        'EdgeColor','none',...             % Edge color 
        'EdgeAlpha',0.1);                  % Transparency of the edges
    light('Position', [1 1 1], 'Style', 'infinite');
    lab{end+1} = 'Defect of Box';
    % Box
    h(end+1) = trisurf(pelvisDefect(v).boundaries.cuboid.hull.tri,...
        pelvisDefect(v).boundaries.cuboid.hull.cornerpoints(:,1),...
        pelvisDefect(v).boundaries.cuboid.hull.cornerpoints(:,2),...
        pelvisDefect(v).boundaries.cuboid.hull.cornerpoints(:,3),...
        'FaceColor',[0.9 0.75 0.68],'EdgeColor',[0.9 0.75 0.68],'FaceAlpha',0.1); 
    lab{end+1} = 'Boundary Box';

    % Reference acetabulum centre
    h(end+1) = plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
        pelvis(1).import.processed.acentre(3), '.','Color','m', 'MarkerSize', 50); % pelvisDefect(1).import.processed.acentre
    lab{end+1} = 'Acentre';

    % Format and display properties
    title(['Intersection "' intersectType{j} '": pelvis defect points ' num2str(u) ' - defect boundary box ' num2str(v)]);
    legend(h, lab, 'Location','bestoutside');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid off
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    view(80, -10) % Figure Paper
    hold off;

    % Save figure (figure unvisible, but saved visible)
    % set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    % savefig(['./Figures/PelvisIntersect(',intersectType{j},')BoxPair(',num2str(u) '-' num2str(v),').fig']) %%%
% end

clear inputDefect inputVertices inputPoints inputAll boxDefectInside boxVerticesInside boxPointsInside boxInside ...
    defectBoxFaces defectBoxVertices defectPointsFaces defectPointsVertices u v lab  h boxOutside

%% Display pairwise intersection with defect volume and points inside/outside (for control)
% Intersection base: 'all'

% Intersection type
% intersectType = {'alpha','refinedAlpha','inside','grid'};
% parfor j = 1:length(intersectType)
    j = 1; %%%
    typeName = intersectType{j};

    % figure('Visible','off')   
    figure
    hold on
    h = [];  lab = {}; % for legend
    
    % Pair number -> rows
    i = 9; % pairRow %%%
    u = pelvisDefectAnalysis(1).intersect.(typeName).pairs.IDs(i,1); % Defect point
    v = pelvisDefectAnalysis(1).intersect.(typeName).pairs.IDs(i,2); % Defect volume

    switch typeName
        case {'alpha', 'refinedAlpha'}
            if strcmp(typeName, 'alpha')
                % alpha vertices (original defect vertices)
                inputBase = pelvisDefect(u).transform.trafo.(weightKabsch{w}).vertices;
                defectPointsFaces = pelvisDefect(u).transform.trafo.(weightKabsch{w}).faces;
                defectPointsVertices = pelvisDefect(u).transform.trafo.(weightKabsch{w}).vertices;
            else % refinedAlpha
                inputBase = pelvisDefect(u).volume.patches.all.comVertices;
                defectPointsFaces = pelvisDefect(u).volume.patches.all.comFaces;
                defectPointsVertices = pelvisDefect(u).volume.patches.all.comVertices;
            end
            inside = inputBase(pelvisDefectAnalysis(1).intersect.(typeName).pairs.defectInsideIdx{i},:);
            % Alpha shape or refined alpha shape (depend on type)
            defectVolFaces = pelvisDefect(v).volume.shrink.(typeName).faces;
            defectVolVertices = pelvisDefect(v).volume.shrink.(typeName).allVerticesPoints;
        case {'inside', 'grid'}
            if strcmp(typeName, 'inside')
                inputPoints = pelvisDefect(u).volume.shrink.inside.mainClusterPoints;     % reference pelvis points inside defect
            else % grid
                inputPoints = pelvisDefect(u).volume.gridPoints.inside;                   % reference pelvis points inside defect
            end
            inputDefect = pelvisDefect(u).volume.patches.all.comVertices;                 % refined alpha vertices (comVertices with patches)
            inputVertices = pelvisDefect(u).volume.shrink.inside.mainClusterVertices;     % reference pelvis vertices inside defect
            inputBase = [inputDefect; inputVertices; inputPoints]; % All
            defectInside = pelvisDefect(u).volume.patches.all.comVertices(pelvisDefectAnalysis(1).intersect.(typeName).pairs.defectInsideIdx{i},:);
            verticesInside = pelvis(1).import.processed.vertices(pelvisDefectAnalysis(1).intersect.(typeName).pairs.verticesInsideIdx{i},:); 
            pointsInside = allPelvis.refPoints.gridPoints.pointsBox(pelvisDefectAnalysis(1).intersect.(typeName).pairs.pointsInsideIdx{i},:);
            inside = [defectInside; verticesInside; pointsInside]; % All
            % Refined alpha shape
            defectPointsFaces = pelvisDefect(u).volume.patches.all.comFaces;
            defectPointsVertices = pelvisDefect(u).volume.patches.all.comVertices;
            defectVolFaces = pelvisDefect(v).volume.shrink.refinedAlpha.faces;
            defectVolVertices = pelvisDefect(v).volume.shrink.refinedAlpha.allVerticesPoints;
    end

    % Defect of points (optional)
    h(end+1) = patch('Faces',defectPointsFaces,...
            'Vertices',defectPointsVertices,...
            'FaceColor',TUMcolors.grey50, ...   % Face color
            'FaceAlpha',0.75,...                % Transparency of the faces
            'EdgeColor','none',...              % Edge color
            'EdgeAlpha',0.25);                  % Transparency of the edges
    lab{end+1} = 'Defect of Points';
    % Defect reference inside
    if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
        h(end+1) = patch('Faces',pelvisDefect(u).volume.shrink.inside.mainClusterFaces,...
            'Vertices',pelvis(1).import.processed.vertices,... % pelvis(1).transform.trafo.(weightKabsch{w}).vertices
            'FaceColor',[0.8 0.8 0.8], ...   % Face color
            'FaceAlpha',0.75,...             % Transparency of the faces % Figure Paper 0.25
            'EdgeColor','none',...           % Edge color
            'EdgeAlpha',0.25);               % Transparency of the edges
        lab{end+1} = 'Reference within Defect of Points';
    end
    % Input points (all)
    h(end+1) = plot3(inputBase(:,1), inputBase(:,2), inputBase(:,3), '.','Color','r', 'MarkerSize', 5);   % TUMcolors.orange
    lab{end+1} = 'Points Outside';
    % Points inside
    h(end+1) = plot3(inside(:,1), inside(:,2), inside(:,3), '.','Color','g', 'MarkerSize', 5); % TUMcolors.green
    lab{end+1} = 'Points Inside';

    % Defect volume
    h(end+1) = patch('Faces',defectVolFaces,...
        'Vertices',defectVolVertices,...
        'FaceColor',[0.9 0.75 0.68], ...    % Face color
        'FaceAlpha',0.75,...                % Transparency of the faces
        'EdgeColor','none',...              % Edge color 
        'EdgeAlpha',0.75);                  % Transparency of the edges
    light('Position', [1 1 1], 'Style', 'infinite');
    lab{end+1} = 'Defect Volume';

    % Reference acetabulum centre
    h(end+1) = plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
        pelvis(1).import.processed.acentre(3), '.','Color','m', 'MarkerSize', 80); % pelvisDefect(1).import.processed.acentre
    lab{end+1} = 'Acentre';

    % Format and display properties
    title(['Intersection "' intersectType{j} '": pelvis defect points ' num2str(u) ' - defect volume ' num2str(v)]);
    legend(h, lab, 'Location','bestoutside');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid off
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    view(80, -10) % Figure Paper
    hold off;

    % Save figure (figure unvisible, but saved visible)
    % set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    % savefig(['./Figures/PelvisIntersect(',intersectType{j},')Pair(',num2str(u) '-' num2str(v),').fig']) %%%
% end

clear inputDefect inputVertices inputPoints inputBase defectInside verticesInside pointsInside inside defectVolFaces ...
    defectVolVertices defectPointsFaces defectPointsVertices h lab

%% Intersection: combinations of the pairwise intersections

% Defect pairs and indices for all pelvis defect categories
for i = 1:numBaseIntersect
    %for j = 1:numIntersectTypes
        % Extract pelvis IDs
        inputIDs = pelvisDefectAnalysis(1).intersect.(intersectType{j}).IDs;
        numIDs = length(inputIDs);
        pairs = pelvisDefectAnalysis(1).intersect.(intersectType{j}).pairs;

        if numIDs > 1
            IDsFor = nchoosek(inputIDs, 2);                  % Forward combinations
            IDsRev = [IDsFor(:, 2), IDsFor(:, 1)];           % Reverse combinations
            pairs.IDs = sortrows([IDsFor; IDsRev], [1, 2]);  % Combined and sorted
        end
        % Find identical pairs of pelvis IDs in 'all'
        [~, rowIdx] = ismember(pairs.IDs, ...
            pelvisDefectAnalysis(1).intersect.(intersectType{j}).pairs.IDs, 'rows');
        pairs.idxOfAll = rowIdx;
        % Update the pelvisIntersect structure
        pelvisDefectAnalysis(i).intersect.(intersectType{j}).pairs = pairs;
    %end
end
clear IDsFor IDsRev pairs rowIdx numIDs inputIDs

% Point frequency of pairs (function intersectCount)
%for j = 1:numIntersectTypes
timeIntersectFreq = NaN(numBaseIntersect, 1);
    typeName = intersectType{j};
    for i = 1:numBaseIntersect
        t1 = tic;
        base = baseIntersect{i};
        % Input of pairs
        referenceDefectIdx = pelvisDefectAnalysis(1).intersect.(typeName).pairs.defectInsideIdx;
        if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
            referenceVerticesIdx = pelvisDefectAnalysis(1).intersect.(typeName).pairs.verticesInsideIdx;
            referencePointsIdx = pelvisDefectAnalysis(1).intersect.(typeName).pairs.pointsInsideIdx;
        else
            referenceVerticesIdx = [];
            referencePointsIdx = [];
        end
        % Frequency table of vertices/points
        pelvisDefectAnalysis(i).intersect = pelvisDefectAnalysis(i).intersect.intersectCount(base, typeName, intersectPctName,...
            referenceDefectIdx, referenceVerticesIdx, referencePointsIdx);
        % Computation time
        timeIntersectFreq(i,1) = toc(t1); 
        allIntersect.time.(['intersect' (typeName) 'FreqPct'])(:,i) = pelvisDefectAnalysis(i).intersect.(typeName).timeIntersectFilter;
        % Data storage of intersection data
        sizeIntersect = pelvisDefectAnalysis(i).intersect.(typeName); % Cache
    end
    % Computation time
    allIntersect.time.(['intersect' (typeName) 'Freq']) = timeIntersectFreq;
    allIntersect.time.(['intersect' (typeName) 'FreqMean']) = mean(timeIntersectFreq, 'omitnan');
    allIntersect.time.(['intersect' (typeName) 'FreqStd']) = std(timeIntersectFreq, 'omitnan');
    % Data storage of intersection data
    % Size of pelvisDefectAnalysis(i).intersect.(typeName) variable; depend on typeName
    allIntersect.(typeName).('sizeFreqBytes') = numel(getByteStreamFromArray(sizeIntersect));
%end
clear sizeIntersect timeIntersectFreq referenceDefectIdx referenceVerticesIdx referencePointsIdx t1 base

% Evaluation: Data storage of intersection data (per pct/intersect)
%typeName = intersectType{j};
for i = 1:numBaseIntersect
    for k = 1:numIntersectPct
        pctCount = pelvisDefectAnalysis(i).intersect.(typeName).(intersectPctName{k}).count;
        if pctCount > 1
            sizeIntersectDefect = pelvisDefectAnalysis(i).intersect.(typeName).(intersectPctName{k}).filteredFreqDefect; % Cache
            allIntersect.(typeName).sizeFreqPctDefect(k,i) = numel(getByteStreamFromArray(sizeIntersectDefect));
            if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
                sizeIntersectVertices = pelvisDefectAnalysis(i).intersect.(typeName).(intersectPctName{k}).filteredFreqVertices; % Cache
                allIntersect.(typeName).sizeFreqPctVertices(k,i) = numel(getByteStreamFromArray(sizeIntersectVertices));
                sizeIntersectPoints  = pelvisDefectAnalysis(i).intersect.(typeName).(intersectPctName{k}).filteredFreqPoints; % Cache
                allIntersect.(typeName).sizeFreqPctPoints(k,i) = numel(getByteStreamFromArray(sizeIntersectPoints));
            end
        else
            allIntersect.(typeName).sizeFreqPctDefect(k,i) = NaN;
            if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
                allIntersect.(typeName).sizeFreqPctVertices(k,i) = NaN;
                allIntersect.(typeName).sizeFreqPctPoints(k,i) = NaN;
            end
        end
    end
end
clear sizeIntersectDefect sizeIntersectVertices sizeIntersectPoints pctCount

% Allocation of coordinates to defects/vertices/points inside and input/output statistics
%for j = 1:numIntersectTypes
timeIntersect = NaN(numIntersectPct, numBaseIntersect);
sizeIntersectBytes = NaN(numBaseIntersect, 1);
    typeName = intersectType{j};
    allDIntersect.(typeName).sizeIntersectAllo = NaN(numIntersectPct,numBaseIntersect); 
    for i = 1:numBaseIntersect
        base = baseIntersect{i};
        for k = 1:numIntersectPct
            currentPctName = intersectPctName{k};
            currentPct = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName);
            if currentPct.count > 1
                tic;
                inputIDs = pelvisDefectAnalysis(i).intersect.(typeName).IDs;
                pelvisDefectAnalysis(i).intersect = pelvisDefectAnalysis(i).intersect.allocateCoords(base, typeName, inputIDs, ...
                    pelvisDefect, pelvis(1).import.processed, allPelvis.refPoints.gridPoints.pointsBox, ...
                    weightKabsch{w}, currentPct, currentPctName);
                timeIntersect(k,i) = toc;
                allIntersect.(typeName).sizeIntersectAllo(k,i) = ...
                pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).sizeIntersectAllo;
            end
        end
        % Data storage of intersection data
        sizeIntersect = pelvisDefectAnalysis(i).intersect.(typeName); 
        sizeInfo = whos('sizeIntersect');
        sizeIntersectBytes(i,1) = sizeInfo.bytes;
    end
    % Computation time
    allIntersect.time.(['intersect' (typeName)]) = timeIntersect;
    allIntersect.time.(['intersect' (typeName) 'Mean']) = mean(timeIntersect, 'omitnan');
    allIntersect.time.(['intersect' (typeName) 'Std']) = std(timeIntersect, 'omitnan');
    % Data storage of intersection data
    % Size of pelvisIntersect(i).(typeName) variable; depend on typeName
    allIntersect.(typeName).sizeIntersect = sizeIntersectBytes;
%end
clear timeIntersect sizeIntersect currentPctName currentPct size Intersect sizeInfo sizeIntersectBytes

%% Evaluate Intersections (vertices/points)

% Input vertices/points
%for j = 1:numIntersectType
    typeName = intersectType{j};
    for i = 1:numBaseIntersect
        IDs = pelvisDefectAnalysis(i).intersect.(typeName).IDs;
        numIDs = length(IDs);
        numDefect = zeros(numIDs, 1);
        for n = 1:numIDs
            id = IDs(n);
            if strcmp(typeName, 'alpha')
                %colDefect(n) = allDefect.stl.vertices(id);  
                numDefect(n) = size(pelvisDefect(id).transform.trafo.(weightKabsch{w}).vertices, 1);  
            else
                numDefect(n) = size(pelvisDefect(id).volume.patches.all.comVertices, 1); 
            end
        end
        allIntersect.(typeName).sumInputDefect(i,1) = sum(numDefect);
        if ismember(typeName, {'inside', 'grid'})
            numVertices = zeros(numIDs, 1);
            numPoints   = zeros(numIDs, 1);
            for n = 1:numIDs
                id = IDs(n);
                numVertices(n) = size(pelvisDefect(id).volume.shrink.inside.mainClusterVertices, 1);
                if strcmp(typeName, 'inside')
                    numPoints(n) = size(pelvisDefect(id).volume.shrink.inside.mainClusterPoints, 1);
                else  % grid
                    numPoints(n) = size(pelvisDefect(id).volume.gridPoints.inside, 1);
                end
            end
            allIntersect.(typeName).sumInputVertices(i,1) = sum(numVertices);
            allIntersect.(typeName).sumInputPoints(i,1) = sum(numPoints); 
        end
        clear IDs numIDs numDefect numVertices numPoints
    end
%end

% Vertices/Points in intersection
%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    for i = 1:numBaseIntersect
        for k = 1:numIntersectPct
            currentPctName = intersectPctName{k};
            if pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).count > 1
                alIntersect.(typeName).numInsideDefect(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).sumInsideDefect;
                allIntersect.(typeName).ratioDefectMean(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).usedDefectMean;
                allIntersect.(typeName).ratioDefectStd(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).usedDefectStd;
                allIntersect.(typeName).ratioOrigDefectMean(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).usedOrigDefectMean;
                allIntersect.(typeName).ratioOrigDefectStd(k,i) =pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).usedOrigDefectStd;
                if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
                    allIntersect.(typeName).numInsideVertices(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).sumInsideVertices;
                    allIntersect.(typeName).numInsidePoints(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).sumInsidePoints;
                    allIntersect.(typeName).ratioVerticesMean(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).usedVerticesMean;
                    allIntersect.(typeName).ratioVerticesStd(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).usedVerticesStd;
                    allIntersect.(typeName).ratioPointsMean(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).usedPointsMean;
                    allIntersect.(typeName).ratioPointsStd(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).usedPointsStd;
                end
            else
                allIntersect.(typeName).numInsideDefect(k,i) = NaN;
                allIntersect.(typeName).ratioDefectMean(k,i) = NaN;
                allIntersect.(typeName).ratioDefectStd(k,i) = NaN;
                allIntersect.(typeName).ratioOrigDefectMean(k,i) = NaN;
                allIntersect.(typeName).ratioOrigDefectStd(k,i) = NaN;
                if strcmp(typeName, 'inside') || strcmp(typeName, 'grid')
                    allIntersect.(typeName).numInsideVertices(k,i) = NaN;
                    allIntersect.(typeName).numInsidePoints(k,i) = NaN;
                    allIntersect.(typeName).ratioVerticesMean(k,i) = NaN;
                    allIntersect.(typeName).ratioVerticesStd(k,i) = NaN;
                    allIntersect.(typeName).ratioPointsMean(k,i) = NaN;
                    allIntersect.(typeName).ratioPointsStd(k,i) = NaN;
                end
            end
        end
    end
%end

%% Generate volumes of intersections - alpha shape

% ShrinkWrap (alphaShape) without user interface (afterwards check with user input) 
% Alpha radius: first time mesh is closed and manifold 

% Loop over intersect types and baseIntersect
%for j = 1:numIntersectTypes
% Preallocation
alphaRadius = NaN(numIntersectPct, numBaseIntersect);
usedDefectCount = NaN(numIntersectPct, numBaseIntersect);
usedCount = NaN(numIntersectPct, numBaseIntersect);
aVolume = NaN(numIntersectPct, numBaseIntersect);
timeIntersectAlpha = NaN(numBaseIntersect,1); 
timeIntersectPctAlpha = NaN(numIntersectPct,numBaseIntersect);
sizeIntersectBytes = NaN(numBaseIntersect,1);
aRadius = 1; % start value
    typeName = intersectType{j};
    for i = 1:numBaseIntersect 
        tOuter = tic;
        for k = 1:numIntersectPct
            currentPctName = intersectPctName{k};
            if pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).count > 1
                tInner = tic;
                % Extract relevant data
                inputDefect = cell2mat(pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).insideDefect);
                numInputDefect = size(inputDefect,1); % for usedCount (only defect for comparison) 
                pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).numInputDefect = numInputDefect; 
                if ismember(typeName, {'inside', 'grid'})
                    %inputVertices = cell2mat(pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).insideVertices);
                    %inputPoints = cell2mat(pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).insidePoints);
                    % for calculation time unique
                    inputVertices = unique(cell2mat(pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).insideVertices), 'rows', 'stable');
                    inputPoints   = unique(cell2mat(ppelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).insidePoints), 'rows', 'stable');
                    addPoints = [inputVertices; inputPoints];
                    numInputVertices = size(inputVertices,1);
                    pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).numInputVertices = numInputVertices; % unique
                    numInputPoints = size(inputPoints,1);
                    pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).numInputPoints = numInputPoints; % unique
                else
                    numInputVertices = 0;
                    numInputPoints = 0;
                    addPoints = [];
                end
                % Call shrinkAlpha function
                pelvisDefectAnalysis(i).intersect = pelvisDefectAnalysis(i).intersect.shrinkAlpha(baseIntersect{i}, typeName, currentPctName, ...
                    aRadius, numInputDefect, numInputVertices, inputDefect, addPoints);
                % Save results
                usedDefectCount(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).usedDefectCount;
                usedCount(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).usedCount;
                aVolume(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).volume;
                alphaRadius(k,i) = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).radius;
                timeIntersectPctAlpha(k,i) = toc(tInner);
            end
        end
        timeIntersectAlpha(i,1) = toc(tOuter);
        % Data storage
        sizeIntersect = pelvisDefectAnalysis(i).intersect.(typeName); % Cache
        sizeInfo = whos('sizeIntersect');
        sizeIntersectBytes(i,1) = sizeInfo.bytes;
    end
    % Computation time
    allIntersect.time.(['intersect' (typeName) 'VolPct']) = timeIntersectPctAlpha;
    allIntersect.time.(['intersect' (typeName) 'Vol']) = timeIntersectAlpha;
    allIntersect.time.(['intersect' (typeName) 'VolMean']) = mean(timeIntersectAlpha, 'omitnan');
    allIntersect.time.(['intersect' (typeName) 'VolStd']) = std(timeIntersectAlpha, 'omitnan');
    % Data storage of intersection data
    allIntersect.intersect.(typeName).('sizeVolBytes') = sizeIntersectBytes;
    % Store
    allIntersect.(typeName).alphaUsedDefectCount = usedDefectCount;
    allIntersect.(typeName).alphaUsedCount = usedCount;
    allIntersect.(typeName).alphaVolume = aVolume;
    allIntersect.(typeName).alphaRadius = alphaRadius;
%end
clear inputDefect numInputDefect inputVertices inputPoints addPoints usedDefectCount usedCount aVolume aRadius alphaRadius ...
    sizeIntersect sizeInfo sizeIntersectBytes tOuter tInner timeIntersectAlpha timeIntersectPctAlpha

% Data storage of intersection data
%for j = 1:numIntersectTypes
sizeSelectedBytes = NaN(numIntersectPct, numBaseIntersect);
    typeName = intersectType{j}; 
    for i = 1:numBaseIntersect
        for k = 1:numIntersectPct
            currentPctName = intersectPctName{k};
            if isfield(pelvisDefectAnalysis(i).intersect.(typeName), currentPctName) && ...
               isfield(pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName), 'radius')
                selected = struct();
                selected.radius             = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).radius;
                selected.faces              = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).faces;
                selected.allVerticesPoints = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).allVerticesPoints;
                selected.usedDefectCount   = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).usedDefectCount;
                selected.volume             = pelvisDefectAnalysis(i).intersect.(typeName).(currentPctName).volume;
                sizeInfo = whos('selected');
                sizeSelectedBytes(k,i) = sizeInfo.bytes;
            else
                sizeSelectedBytes(k,i) = NaN;
            end
        end
    end
    allIntersect.(typeName).sizeVolSelectedBytes = sizeSelectedBytes;
%end
clear sizeSelectedBytes selected sizeInfo

%% Display intersections (for control)

% Select type 
% baseIntersect = {'all','BTE','Other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
i = 7; %%%
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 3; %%%

% Colors
inputIDs = pelvisDefectAnalysis(i).intersect.(intersectType{j}).IDs;
numIDs = length(inputIDs);
colorList = lines(numIDs);

% intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
%for k = 1:numIntersectPct
    k = 4; %%%

    if pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).count > 1
        figure
        %figure('Visible','off')
        hold on

        inputIDs = pelvisDefectAnalysis(i).intersect.(intersectType{j}).IDs;
        numIDs = size(inputIDs,1);

        for l = 1:numIDs
            % Colors
            color = colorList(mod(l - 1, size(colorList, 1)) + 1, :);
            %color = TUMcolors.blue300;

            % Input
            u = pelvisDefectAnalysis(i).intersect.(intersectType{j}).IDs(l); % pelvisDefect ID
            switch intersectType{j}
                case 'alpha'
                    pointBaseDefect = pelvisDefect(u).transform.trafo.(weightKabsch{w}).vertices;
                    % Volume
                    defectFaces = pelvisDefect(u).volume.shrink.alpha.faces;
                    defectVertices = pelvisDefect(u).volume.shrink.alpha.allVerticesPoints;
                case 'refinedAlpha'
                    pointBaseDefect = pelvisDefect(u).volume.patches.all.comVertices;
                    % Volume
                    defectFaces = pelvisDefect(u).volume.shrink.refinedAlpha.faces;
                    defectVertices = pelvisDefect(u).volume.shrink.refinedAlpha.allVerticesPoints;
                case {'inside', 'grid'}
                    pointBaseDefect = pelvisDefect(u).volume.patches.all.comVertices;
                    pointBaseVertices = pelvis(1).import.processed.vertices;
                    pointBasePoints = allPelvis.refPoints.gridPoints.pointsBox;
                    % Volume
                    defectFaces = pelvisDefect(u).volume.shrink.refinedAlpha.faces;
                    defectVertices = pelvisDefect(u).volume.shrink.refinedAlpha.allVerticesPoints;
            end

            % Defect plot
            pointBaseDefectIdx = ...
                pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).filteredFreqDefect{l}.PointBaseIndex;
            insideDefect = pointBaseDefect(pointBaseDefectIdx, :);
            plot3(insideDefect(:, 1), insideDefect(:, 2), insideDefect(:, 3), '*', 'Color', color, 'MarkerSize', 5);

            % Vertices and Points plot
            if ismember(intersectType{j}, {'inside', 'grid'})
                % Vertices
                pointBaseVerticesIdx = ...
                    pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).filteredFreqVertices{l}.PointBaseIndex;
                insideVertices = pointBaseVertices(pointBaseVerticesIdx, :);
                plot3(insideVertices(:, 1), insideVertices(:, 2), insideVertices(:, 3), 'o', ...
                    'Color', color, 'MarkerSize', 3);
                % Reference pelvis inside
                patch('Faces',pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).facesVerticesInside{l},...
                    'Vertices',pointBaseVertices,...
                    'FaceColor',color, ...              % Face color
                    'FaceAlpha',0.5,...                 % Transparency of the faces
                    'EdgeColor','none',...              % Edge color
                    'EdgeAlpha',0.25);                  % Transparency of the edges

                % Points
                pointBasePointsIdx = ...
                    pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).filteredFreqPoints{l}.PointBaseIndex;
                insidePoints = pointBasePoints(pointBasePointsIdx, :);
                plot3(insidePoints(:, 1), insidePoints(:, 2), insidePoints(:, 3), 'x', ...
                    'Color', color, 'MarkerSize', 3);
            end

            % Pelvis defect volumes (optional)
            patch('Faces',defectFaces,...
                'Vertices',defectVertices,...
                'FaceColor',TUMcolors.grey20, ...   % Face color
                'FaceAlpha',0.05,...                % Transparency of the faces
                'EdgeColor','none',...              % Edge color
                'EdgeAlpha',0.25);                  % Transparency of the edges
            % Pelvis defect inside
            patch('Faces',pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).facesDefectInside{l},...
                'Vertices',pointBaseDefect,...
                'FaceColor',color, ...              % Face color
                'FaceAlpha',0.5,...                 % Transparency of the faces
                'EdgeColor','none',...              % Edge color
                'EdgeAlpha',0.25);                  % Transparency of the edges

            % Reference acetabulum centre
            plot3(pelvis(1).import.processed.acentre(1), ...
                pelvis(1).import.processed.acentre(2), ...
                pelvis(1).import.processed.acentre(3), ...
                '.','Color','m', 'MarkerSize', 62); % 58 62
            plot3(pelvis(1).import.processed.acentre(1), ...
               pelvis(1).import.processed.acentre(2), ...
                pelvis(1).import.processed.acentre(3), ...
                'x', 'MarkerEdgeColor', 'm', 'MarkerSize', 27, 'LineWidth', 8); % 25/7 27/8

        end
        light('Position', [1 1 5], 'Style', 'infinite');

        % Format and display properties
        title(['Intersection of type "' intersectType{j} '" for defect category "' baseIntersect{i} '" and percentage ' intersectPctName{k}]);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        grid off
        daspect([1, 1, 1]); % Equal aspect ratio for the axes
        view(3);
        view(80, -10) % Figure Paper
        hold off;

    else
        disp('Intersect count too low')
    end

    % Save figure (figure unvisible, but saved visible)
    % set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
    % savefig(['./Figures/PelvisIntersect(',intersectType{j},')Base(',baseIntersect{i},')Pct(',intersectPctName{k},').fig'])
%end

clear inputIDs numIDs colorList  color ...
    pointBaseDefect defectFaces defectVertices pointBaseVertices pointBasePoints ...
    pointBaseDefectIdx insideDefect pointBaseVerticesIdx insideVertices pointBasePointsIdx insidePoints

%% Display alphaShape of intersections (for control)

% Select type 
% baseIntersect = {'all','BTE','Other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
i = 1; %%%
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 3; %%%

% intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
%for k = 1:numIntersectPct
    k = 1; %%%
    
    if pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).count > 1
        %figure('Visible','off')
        figure
        hold on
    
        % Alpha shape of intersection (lightning)
        patch('Faces',pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).faces,...
            'Vertices',pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).allVerticesPoints,...
            'FaceColor',[0.8 0.8 0.8], ...      % Face color
            'FaceAlpha',0.25,...                % Transparency of the faces 
            'EdgeColor','none',...              % Edge color
            'EdgeAlpha',0.25,...                % Transparency of the edges
            ... % Ligthing for 3d effect
            'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
            'AmbientStrength', 0.5);
        light('Position', [1 1 5], 'Style', 'infinite');
    
        % Display all vertices/points
        plot3(pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).allVerticesPoints(:,1), ...
            pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).allVerticesPoints(:,2), ...
            pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).allVerticesPoints(:,3),...
            '.','Color',TUMcolors.blue300,'MarkerSize', 10)

        % Display points used for alphaShape
        usedIdx = unique(pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).faces(:));
        usedPoints = pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).allVerticesPoints(usedIdx, :);
        plot3(usedPoints(:,1), usedPoints(:,2), usedPoints(:,3), '.','Color',TUMcolors.orange,'MarkerSize', 7)

        % Display vertices from original defects used for alphaShape
        inputIDs = pelvisDefectAnalysis(i).intersect.(intersectType{j}).IDs;
        numIDs = length(inputIDs);
        validOrig = [];
        for l = 1:numIDs % pointBase / Defect
            u = inputIDs(l);
            pointBaseOrigDefect = pelvisDefect(u).transform.trafo.(weightKabsch{w}).vertices;
            pointBaseDefectIdx = pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).filteredFreqDefect{l}.PointBaseIndex;
            validOrigIdx = pointBaseDefectIdx(pointBaseDefectIdx <= size(pointBaseOrigDefect,1));
            origPoints = pointBaseOrigDefect(validOrigIdx, :);
            [~, ia] = intersect(origPoints, usedPoints, 'rows');
            validOrig = [validOrig; origPoints(ia,:)];
        end
        plot3(validOrig(:,1), validOrig(:,2), validOrig(:,3), '.','Color','g','MarkerSize', 10)

        % Reference acetabulum centre
        plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
                pelvis(1).import.processed.acentre(3), '.','Color','r', 'MarkerSize', 50);
    
        % Format and display properties
        title(['Alpha shape of type "' intersectType{j} '" for defect category "' baseIntersect{i} '" and percentage ' intersectPctName{k}]);
        legend('Alpha Shape', 'Points inside', 'Points on Alpha Shape', 'Original Defect Points on Alpha Shape', 'Acetabulum Centre');
        daspect([1, 1, 1]); % Equal aspect ratio for the axes
        xlabel('X'); ylabel('Y'); zlabel('Z');
        view(3);
        view(80, -10) % Figure Paper
        hold off;
    
    else
        disp('Intersect count too low')
    end

    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    %savefig(['./Figures/PelvisIntersect(',intersectType{j},')Base(',baseIntersect{i},')Pct(',intersectPctName{k},')AlphaShape.fig'])
%end
clear usedIdx usedPoints inputIDs numIDs pointBaseOrigDefect pointBaseDefectIdx validOrigIdx origPoints validOrig ia

%% Intersection: comparisons 

% Preallocate storage for comparison results
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
    'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
% Initialize temporary object
% Iterate through all combinations of intersect types
%for j = 1:numIntersectTypes % type1
    %for k = j+1:numIntersectTypes % type2
    tempIntersect = Intersection();
        type1 = intersectType{j};
        type2 = intersectType{k};
        combi = [type1 '_' type2];

        for i = 1:numBaseIntersect
            baseName = groupNames{i};
            inputIDs = pelvisDefectAnalysis(i).intersect.(type1).IDs;
            numIDs = length(inputIDs);

            for p = 1:numIntersectPct
                pctName = intersectPctName{p};
                currentPct1 = pelvisDefectAnalysis(i).intersect.(type1).(pctName);
                currentPct2 = pelvisDefectAnalysis(i).intersect.(type2).(pctName);

                if currentPct1.count > 1 && currentPct2.count > 1
                    % Used count
                    if isfield(pelvisDefectAnalysis(i).intersect.(type1).(pctName), 'usedDefect') % alpha, refinedAlpha
                        tempIntersect.comp.(combi).(baseName).(pctName).usedOrigDefectType1 = pelvisDefectAnalysis(i).intersect.(type1).(pctName).usedOrigDefect;
                        tempIntersect.comp.(combi).(baseName).(pctName).usedOrigDefectType2 = pelvisDefectAnalysis(i).intersect.(type2).(pctName).usedOrigDefect;
                    end
                    if isfield(pelvisDefectAnalysis(i).intersect.(type1).(pctName), 'usedVertices') % inside, grid
                        tempIntersect.comp.(combi).(baseName).(pctName).usedVerticesType1 = pelvisDefectAnalysis(i).intersect.(type1).(pctName).usedVertices;
                        tempIntersect.comp.(combi).(baseName).(pctName).usedVerticesType2 = pelvisDefectAnalysis(i).intersect.(type2).(pctName).usedVertices;
                        tempIntersect.comp.(combi).(baseName).(pctName).usedPointsType1 = pelvisDefectAnalysis(i).intersect.(type1).(pctName).usedPoints;
                        tempIntersect.comp.(combi).(baseName).(pctName).usedPointsType2 = pelvisDefectAnalysis(i).intersect.(type2).(pctName).usedPoints;
                    end

                    % Perform comparison
                    tempIntersect = tempIntersect.compareIntersect(type1, type2, baseName, pctName, numIDs, ...
                        pelvisDefectAnalysis(i).intersect.(type1).(pctName), pelvisDefectAnalysis(i).intersect.(type2).(pctName));
                    % Allocation of coordinates
                    tempIntersect = tempIntersect.allocateCoordsComp(combi, baseName, pctName, inputIDs, ...
                        pelvisDefect, pelvis, allPelvis, type1, type2, weightKabsch{w});
                    % Save results
                    allIntersect.(combi) = tempIntersect.comp.(combi);
                end
            end
        end
    %end
%end
clear type1 type2 currentPct1 currentPct2 combi numIDs baseName pctName tempIntersect

%for j = 1:numIntersectTypes % type1
    %for k = j+1:numIntersectTypes % type2
        type1 = intersectType{j};
        type2 = intersectType{k};
        combi = [type1 '_' type2];
        for i = 1:numBaseIntersect
            baseName = groupNames{i};
            for p = 1:numIntersectPct
                pctName = intersectPctName{p};
                currentPct1 = pelvisDefectAnalysis(i).intersect.(type1).(pctName);
                currentPct2 = pelvisDefectAnalysis(i).intersect.(type2).(pctName);
                if currentPct1.count > 1 && currentPct2.count > 1
                    % Defect
                    sumUniqueDefectType1(p,i) = allIntersect.(combi).(baseName).(pctName).sumUniqueDefectType1;
                    sumUniqueDefectType2(p,i) = allIntersect.(combi).(baseName).(pctName).sumUniqueDefectType2;
                    sumSharedDefect(p,i) = allIntersect.(combi).(baseName).(pctName).sumSharedDefect;
                    if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
                        % Vertices
                        sumUniqueVerticesType1(p,i) = allIntersect.(combi).(baseName).(pctName).sumUniqueVerticesType1;
                        sumUniqueVerticesType2(p,i) = allIntersect.(combi).(baseName).(pctName).sumUniqueVerticesType2;
                        sumSharedVertices(p,i) = allIntersect.(combi).(baseName).(pctName).sumSharedVertices;
                        % Points
                        sumUniquePointsType1(p,i) = allIntersect.(combi).(baseName).(pctName).sumUniquePointsType1;
                        sumUniquePointsType2(p,i) = allIntersect.(combi).(baseName).(pctName).sumUniquePointsType2;
                        sumSharedPoints(p,i) = allIntersect.(combi).(baseName).(pctName).sumSharedPoints;
                    end
                end
            end
        end
    %end
%end
clear sumUniqueDefectType1 sumUniqueDefectType2 sumSharedDefect sumUniqueVerticesType1 sumUniqueVerticesType2 sumSharedVertices...
    sumUniquePointsType1 sumUniquePointsType2 sumSharedPoints

%% Display differences of intersections (for control)

% Inputs
% groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
%   'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 7; % Base intersect index %%%
baseName = groupNames{i}; % Base intersect name
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 1; % type1 %%%
k = 3; % type2 %%%
type1 = intersectType{j};
type2 = intersectType{k};
combi = [type1 '_' type2];

% intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
%for p = 1:numIntersectPct
    p = 3; %%%
    pctName = intersectPctName{p}; % Intersect percentage name
    % Extract data from allIntersect
    data = allIntersect.(combi).(baseName).(pctName);
    numIDs = length(pelvisDefectAnalysis(i).intersect.(type1).IDs);

    %figure('Visible','off')
    figure
    hold on
    
    % Define colors and markers
    colors = {[0.9, 0.1, 0.1], [0.1, 0.1, 0.9], [0.1, 0.8, 0.1], [1, 1, 0]}; % Red, Blue, Green
    labels = {'Unique to Type1', 'Unique to Type2', 'Shared', 'Acetabulum Centre'};
    markerSize = 30;
    markerStyle = 'o';

    % Define field groups
    fieldsDefect   = {'insideDefectType1', 'insideDefectType2', 'insideDefectShared'};
    fieldsVertices = {'insideVerticesType1', 'insideVerticesType2', 'insideVerticesShared'};
    fieldsPoints   = {'insidePointsType1', 'insidePointsType2', 'insidePointsShared'};

    % Create dummy points for legend
    for f = 1:length(fieldsDefect)
    scatter3(NaN, NaN, NaN, markerSize, ...
        'Marker', markerStyle, ...
        'MarkerFaceColor', colors{f}, ...
        'MarkerEdgeColor', 'k', ...
        'DisplayName', labels{f});
    end
    % Dummy for Acetabulum Centre
    scatter3(NaN, NaN, NaN, markerSize, ...
        'Marker', markerStyle, ...
        'MarkerFaceColor', colors{4}, ...
        'MarkerEdgeColor', 'k', ...
        'DisplayName', labels{4});

    % Plot Vertices and Points if applicable
    if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
        for f = [1, 2, 3]  % Shared is polottet last
            % Vertices
            for l = 1:numIDs
                points = data.(fieldsVertices{f}){l};
                if ~isempty(points)
                    scatter3(points(:, 1), points(:, 2), points(:, 3), ...
                        markerSize, ...
                        'Marker', markerStyle, ...
                        'MarkerFaceColor', colors{f}, ...
                        'MarkerEdgeColor', 'k', ...
                        'LineWidth', 0.5, ...
                        'HandleVisibility', 'off');
                end
            end
            % Points
            for l = 1:numIDs
                points = data.(fieldsPoints{f}){l};
                if ~isempty(points)
                    scatter3(points(:, 1), points(:, 2), points(:, 3), ...
                        markerSize, ...
                        'Marker', markerStyle, ...
                        'MarkerFaceColor', colors{f}, ...
                        'MarkerEdgeColor', 'k', ...
                        'LineWidth', 0.5, ...
                        'HandleVisibility', 'off');
                end
            end
        end
    end

    % Plot each type of defect vertices
    for f = [1, 2, 3]  % Shared is polottet last
        for l = 1:numIDs
            points = data.(fieldsDefect{f}){l};
            if ~isempty(points)
                scatter3(points(:, 1), points(:, 2), points(:, 3), ...
                    markerSize, ...
                    'Marker', markerStyle, ...
                    'MarkerFaceColor', colors{f}, ...
                    'MarkerEdgeColor', 'k', ...
                    'LineWidth', 0.5, ...
                    'HandleVisibility', 'off');
            end
        end
    end

    % Reference acetabulum centre
    acentre = pelvis(1).import.processed.acentre;
    scatter3(acentre(1), acentre(2), acentre(3), ...
        200, ...                     
        'o', ...                     
        'MarkerFaceColor', 'm', ... 
        'MarkerEdgeColor', 'k', ... 
        'LineWidth', 0.8, ...       
        'HandleVisibility', 'off'); 

    % Plot title and labels
    title({['Unique/Shared Points: ', type1, ' vs ', type2], ...
           ['Percentage: ', pctName, ' | Base: ', baseName]});
    xlabel('X'); ylabel('Y'); zlabel('Z');
    legend('Location', 'best', 'Interpreter', 'none');
    grid off;
    axis equal;
    hold off;
    view(3);
    view(80, -10) % Figure Paper
        
    % Save figure (figure unvisible, but saved visible)
    %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
    %savefig(['./Figures/PelvisIntersect(',combi,')Base(',baseIntersect{i},')Pct(',intersectPctName{k},').fig'])
%end

clear data  numIDs colors labels markerSize markerStyle fieldsDefect fieldsVertices fieldsPoints points acentre

%% Display differences of intersections (3D) with slice plane (for control)

% Inputs
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
  'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 7; % Base intersect index %%%
baseName = groupNames{i}; % Base intersect name
intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
p = 4; % Intersect percentage index %%%
pctName = intersectPctName{p}; % Intersect percentage name
intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 1; % type1 %%%
k = 3; % type2 %%%
type1 = intersectType{j};
type2 = intersectType{k};
combi = [type1 '_' type2];

% Extract data from allIntersect
data = allIntersect.(combi).(baseName).(pctName);

% Marker settings
markerSize = 30;
markerStyle = 'o';
colors = {[1, 0, 0], [0, 0, 1],[0, 1, 0],[1, 0, 1], [0, 1, 1]}; 
labels = {'Unique to Type1', 'Unique to Type2', 'Shared', 'Acetabulum Centre', 'Slice Plane'};

% Fields
fieldsDefect   = {'insideDefectType1', 'insideDefectType2', 'insideDefectShared'};
fieldsVertices = {'insideVerticesType1', 'insideVerticesType2', 'insideVerticesShared'};
fieldsPoints   = {'insidePointsType1', 'insidePointsType2', 'insidePointsShared'};

% Collect vertices/points: first type1, then type2, then shared (in the foreground)
allPoints = [];
allColors = [];

% inside/grid: Add vertices and points
if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
    for f_idx = 1:2
        fieldDataV = data.(fieldsVertices{f_idx});
        fieldPointsV = cat(1, fieldDataV{:});
        if ~isempty(fieldPointsV)
            allPoints = [allPoints; fieldPointsV];
            allColors = [allColors; repmat(colors{f_idx}, size(fieldPointsV,1), 1)];
        end
        fieldDataP = data.(fieldsPoints{f_idx});
        fieldPointsP = cat(1, fieldDataP{:});
        if ~isempty(fieldPointsP)
            allPoints = [allPoints; fieldPointsP];
            allColors = [allColors; repmat(colors{f_idx}, size(fieldPointsP,1), 1)];
        end
    end
    % Shared Vertices/Points
    fieldPointsV = cat(1, data.(fieldsVertices{3}){:});
    fieldPointsP = cat(1, data.(fieldsPoints{3}){:});
    if ~isempty(fieldPointsV)
        allPoints = [allPoints; fieldPointsV];
        allColors = [allColors; repmat(colors{3}, size(fieldPointsV,1), 1)];
    end
    if ~isempty(fieldPointsP)
        allPoints = [allPoints; fieldPointsP];
        allColors = [allColors; repmat(colors{3}, size(fieldPointsP,1), 1)];
    end
end

% Defect 
% Plot defects at the end so that they are in the foreground
for f_idx = 1:2
    fieldData = data.(fieldsDefect{f_idx});
    fieldPoints = cat(1, fieldData{:});
    if ~isempty(fieldPoints)
        allPoints = [allPoints; fieldPoints];
        allColors = [allColors; repmat(colors{f_idx}, size(fieldPoints,1), 1)];
    end
end
% Shared
fieldData = data.(fieldsDefect{3});
fieldPoints = cat(1, fieldData{:});
if ~isempty(fieldPoints)
    allPoints = [allPoints; fieldPoints];
    allColors = [allColors; repmat(colors{3}, size(fieldPoints,1), 1)];
end

% Figure
figure
hold on

% For all types refinedAlpha volume
% %if pelvisIntersect(i).(intersectType{j}).(intersectPctName{k}).count > 1
% inputIDs = pelvisDefectAnalysis(i).intersect.refinedAlpha.IDs;
% numIDs = size(inputIDs,1);
% for l = 1:numIDs
%     u = pelvisDefectAnalysis(i).intersect.refinedAlpha.IDs(l); % pelvisDefect ID
%     defectFaces    = pelvisDefect(u).volume.shrink.refinedAlpha.faces;
%     defectVertices = pelvisDefect(u).volume.shrink.refinedAlpha.allVerticesPoints;
%     patch('Faces', defectFaces, ...
%         'Vertices', defectVertices, ...
%         'FaceColor', TUMcolors.grey20, ...
%         'FaceAlpha', 0.25, ...
%         'EdgeColor', 'none', ...
%         'EdgeAlpha', 0.2); 
% end

% Main scatter points
scatter3( ...
    allPoints(:,1), allPoints(:,2), allPoints(:,3), ...
    markerSize, ...
    'Marker', markerStyle, ...
    'MarkerFaceColor', 'flat', ...
    'MarkerEdgeColor', [0.2, 0.2, 0.2], ...
    'LineWidth', 0.25, ...
    'CData', allColors); 

% Reference acetabulum centre
acentre = pelvis(1).import.processed.acentre;
% scatter3(acentre(1), acentre(2), acentre(3), ...
%     1800, ...
%     'Marker', markerStyle, ...
%     'MarkerFaceColor', 'm', ...
%     'MarkerEdgeColor', 'none', ...
%     'LineWidth', 0.5, ...
%     'DisplayName', 'Acetabulum Centre');
plot3(acentre(1), acentre(2), acentre(3), ...
    '.','Color','m', 'MarkerSize', 130); % 120, 100, 140 ,130
plot3(acentre(1), acentre(2), acentre(3), ...
    'x', 'MarkerEdgeColor', 'm', 'MarkerSize', 52, 'LineWidth', 16); % 50/15 42/13 55/17 52/16

% Slice plane at currentZ
%currentZ = (min(allPoints(:,3))+ max(allPoints(:,3)))/2; % plane in centre of plot
currentZ = pelvis(1).import.processed.acentre(3); % plane through acentre 
[X, Y] = meshgrid(linspace(min(allPoints(:,1)), max(allPoints(:,1)), 10), ...
    linspace(min(allPoints(:,2)), max(allPoints(:,2)), 10));
Z = ones(size(X)) * currentZ;
surf(X, Y, Z, ...
    'FaceAlpha', 0.75, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'c', ...
    'DisplayName', 'Slice Plane'); 

% Dummy-Scatter for legend
legendHandles = gobjects(1, 5);
for f = 1:4
    legendHandles(f) = scatter3(nan, nan, nan, markerSize, ...
        'Marker', markerStyle, ...
        'MarkerFaceColor', colors{f}, ...
        'MarkerEdgeColor', [0.2, 0.2, 0.2], ...
        'LineWidth', 0.25, ...
        'DisplayName', labels{f});
end
legendHandles(5) = patch(nan, nan, nan, 'c', ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'Slice Plane', ...
    'Tag', 'LegendDummy');
legend(legendHandles, labels, 'Location', 'best', 'Interpreter', 'none');

% Format
xlabel('X'); ylabel('Y'); zlabel('Z');
title( {['3D View – Slice Z = ' num2str(currentZ, '%.2f')], ...
    ['Unique/Shared Points: ', type1, ' vs ', type2], ...
    ['Base: ', baseName, ' | Percentage: ', pctName]});
axis('equal');
grid( 'off');
view( 80, -10);  % Paper View

clear X Y Z fieldsDefect  fieldsVertices fieldsPoints allPoints allColors fieldData fieldDataV fieldPointsV fieldDataP...
    fieldPointsP fieldData fieldPoints inputIDs numIds u defectFaces defectVertices acentre currentZ legendHandles 

%% Display differences of intersections in slices (for control)

% Inputs
% groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
%   'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 7; % Base intersect index %%%
baseName = groupNames{i}; % Base intersect name
% intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
p = 4; % Intersect percentage index %%%
pctName = intersectPctName{p}; % Intersect percentage name
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 1; % type1 %%%
k = 3; % type2 %%%
type1 = intersectType{j};
type2 = intersectType{k};
combi = [type1 '_' type2];

% Extract data from allIntersect
data = allIntersect.(combi).(baseName).(pctName);

% Marker settings
markerSize = 30;
markerStyle = 'o';
colors = {[1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 1], [0, 1, 1]}; 
labels = {'Unique to Type1', 'Unique to Type2', 'Shared', 'Acetabulum Centre', 'Slice Plane'};

% Fields
fieldsDefect   = {'insideDefectType1', 'insideDefectType2', 'insideDefectShared'};
fieldsVertices = {'insideVerticesType1', 'insideVerticesType2', 'insideVerticesShared'};
fieldsPoints   = {'insidePointsType1', 'insidePointsType2', 'insidePointsShared'};

% Collect vertices/points: first type1, then type2, then shared (in the foreground)
allPoints = [];
allColors = [];

% inside/grid: Add vertices and points
if ismember(type1, {'inside', 'grid'}) || ismember(type2, {'inside', 'grid'})
    for f_idx = 1:2
        fieldDataV = data.(fieldsVertices{f_idx});
        fieldPointsV = cat(1, fieldDataV{:});
        if ~isempty(fieldPointsV)
            allPoints = [allPoints; fieldPointsV];
            allColors = [allColors; repmat(colors{f_idx}, size(fieldPointsV,1), 1)];
        end
        fieldDataP = data.(fieldsPoints{f_idx});
        fieldPointsP = cat(1, fieldDataP{:});
        if ~isempty(fieldPointsP)
            allPoints = [allPoints; fieldPointsP];
            allColors = [allColors; repmat(colors{f_idx}, size(fieldPointsP,1), 1)];
        end
    end
    % Shared Vertices/Points
    fieldPointsV = cat(1, data.(fieldsVertices{3}){:});
    fieldPointsP = cat(1, data.(fieldsPoints{3}){:});
    if ~isempty(fieldPointsV)
        allPoints = [allPoints; fieldPointsV];
        allColors = [allColors; repmat(colors{3}, size(fieldPointsV,1), 1)];
    end
    if ~isempty(fieldPointsP)
        allPoints = [allPoints; fieldPointsP];
        allColors = [allColors; repmat(colors{3}, size(fieldPointsP,1), 1)];
    end
end

% Defect 
% Plot defects at the end so that they are in the foreground
for f_idx = 1:2
    fieldData = data.(fieldsDefect{f_idx});
    fieldPoints = cat(1, fieldData{:});
    if ~isempty(fieldPoints)
        allPoints = [allPoints; fieldPoints];
        allColors = [allColors; repmat(colors{f_idx}, size(fieldPoints,1), 1)];
    end
end
% Shared
fieldData = data.(fieldsDefect{3});
fieldPoints = cat(1, fieldData{:});
if ~isempty(fieldPoints)
    allPoints = [allPoints; fieldPoints];
    allColors = [allColors; repmat(colors{3}, size(fieldPoints,1), 1)];
end

% Position of start plane in 3D
z_min = min(allPoints(:,3));
z_max = max(allPoints(:,3));
tolerance = 0.5;
%currentZ = (z_min + z_max)/2; % plane in centre of plot
currentZ = pelvis(1).import.processed.acentre(3); % plane through acentre 

% Figure + slider
figure('Name','2D- and 3D-View','NumberTitle','off');
ax1 = subplot(1,2,1); 
ax2 = subplot(1,2,2); 
hold(ax2, 'on');

% Create dummy points for legend (including Acetabulum Centre)
legendHandles = gobjects(1,5);
for f = 1:4
    legendHandles(f) = scatter3(NaN, NaN, NaN, markerSize, ...
        'Marker', markerStyle, ...
        'MarkerFaceColor', colors{f}, ...
        'MarkerEdgeColor', [0.2, 0.2, 0.2], ...
        'LineWidth', 0.25, ...
        'DisplayName', labels{f}, ...
        'Tag', 'LegendDummy');
end
% Add dummy surface for slice plane
legendHandles(5) = patch(nan, nan, nan, 'c', ...
    'FaceAlpha', 0.75, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'Slice Plane', ...
    'Tag', 'LegendDummy');
legend(ax2, legendHandles, labels, 'Location', 'best', 'Interpreter', 'none');

% Draw 3D 
draw3D = 'true';
% Slider
slider = uicontrol('Style', 'slider', ...
    'Min', z_min, 'Max', z_max, 'Value', currentZ, ...
    'SliderStep', [0.01, 0.1], ...
    'Units', 'normalized', ...
    'Position', [0.2 0.01 0.6 0.05], ...
    'Callback', @(src,evt) updateZSlice(ax1, ax2, allPoints, allColors, src.Value, tolerance, ...
                                        type1, type2, pctName, baseName, markerSize, markerStyle, pelvis, draw3D));
% Slice Update
updateZSlice(ax1, ax2, allPoints, allColors, currentZ, tolerance, type1, type2, pctName, baseName, markerSize, markerStyle, ...
    pelvis, draw3D);

%clear data numIDs colors labels markerSize markerStyle fieldsDefect fieldsVertices fieldsPoints allPoints allColors acentre ...
%    fieldDataV fieldPointsV fieldDataP fieldPointsP fieldData fieldPoints z_min z_max tolerance ax1 ax2 legendHandles slider ...
%    xlim2D ylim2D zlim3D draw3D p f_idx

% Function for slices
function updateZSlice(ax1, ax2, points, colors, currentZ, tolerance, type1, type2, pctName, baseName, markerSize, ...
    markerStyle, pelvis, draw3D)
    idx = abs(points(:,3) - currentZ) < tolerance;
    filteredPoints = points(idx,:);
    filteredColors = colors(idx,:);

    % 2D Plot
    axes(ax1); cla(ax1); hold(ax1, 'on');
    if ~isempty(filteredPoints)
        scatter(ax1, ...
            filteredPoints(:,1), filteredPoints(:,2), ...
            markerSize, ...
            'Marker', markerStyle, ...
            'MarkerFaceColor', 'flat', ...
            'MarkerEdgeColor', [0.2, 0.2, 0.2], ...
            'LineWidth', 0.25, ...
            'CData', filteredColors);
    end
    % Acetabulum Centre in 2D View
    acentre = pelvis(1).import.processed.acentre;
    if abs(acentre(3) - currentZ) < tolerance
        % scatter(ax1, acentre(1), acentre(2), ...
        %     600, ...
        %     'Marker', markerStyle, ...
        %     'MarkerFaceColor', 'm', ...
        %     'MarkerEdgeColor', 'none', ...
        %     'LineWidth', 0.25);
        plot3(acentre(1), acentre(2), acentre(3), ...
            '.','Color','m', 'MarkerSize', 95); % 80 90 85 95
        plot3(acentre(1), acentre(2), acentre(3), ...
            'x', 'MarkerEdgeColor', 'm', 'MarkerSize', 40, 'LineWidth', 13); % 32/11 37/12 35/12 40/13
    end
    xlabel(ax1,'X'); ylabel(ax1,'Y');
    title(ax1, {['2D View Z = ' num2str(currentZ, '%.2f') ' ± ' num2str(tolerance, '%.2f')], ...
        ['Unique/Shared Points: ', type1, ' vs ', type2], ...
        ['Percentage: ', pctName, ' | Base: ', baseName]});
    % Format
    xlim2D = [min(points(:,1))-5, max(points(:,1))+5];
    ylim2D = [min(points(:,2))-5, max(points(:,2))+5];
    axis(ax1, 'equal');
    xlim(ax1, xlim2D);
    ylim(ax1, ylim2D);
    grid(ax1, 'off');

    % 3D Plot
    if draw3D
        delete(findall(ax2, 'Type', 'Scatter', '-not', 'Tag', 'LegendDummy'));
        delete(findall(ax2, 'Type', 'Surface'));
        if ~isempty(points)
            scatter3(ax2, ...
                points(:,1), points(:,2), points(:,3), ...
                markerSize, ...
                'Marker', markerStyle, ...
                'MarkerFaceColor', 'flat', ...
                'MarkerEdgeColor', [0.2, 0.2, 0.2], ...
                'LineWidth', 0.25, ...
                'CData', colors, ...
                'HandleVisibility', 'off');
        end
        % Plot Acetabulum Centre
        acentre = pelvis(1).import.processed.acentre;
        scatter3(ax2, acentre(1), acentre(2), acentre(3), ...
            150, ...
            'Marker', markerStyle, ...
            'MarkerFaceColor', [1, 1, 0], ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.5, ...
            'DisplayName', 'Acetabulum Centre', ...
            'Tag', 'LegendDummy', ...
            'HandleVisibility', 'off');
        % Plane in 3D
        [X, Y] = meshgrid(linspace(xlim2D(1), xlim2D(2), 10), ...
            linspace(ylim2D(1), ylim2D(2), 10));
        Z = ones(size(X)) * currentZ;
        surf(ax2, X, Y, Z, 'FaceAlpha', 0.75, 'EdgeColor', 'none', 'FaceColor', 'c','HandleVisibility', 'off');
        % Format
        xlabel(ax2,'X'); ylabel(ax2,'Y'); zlabel(ax2,'Z');
        title(ax2,'3D View');
        axis(ax2,'equal');
        xlim(ax2, xlim2D);
        ylim(ax2, ylim2D);
        zlim3D = [min(points(:,3)), max(points(:,3))];
        zlim(ax2, zlim3D);
        grid(ax2,'off');
        view(ax2, 37.5, 30);
        %view(ax2,80, -10) % Figure Paper
    end

    drawnow;
end

%% Video: Display differences of intersections in slices 
% Dependency to Display differences of intersections in slices 

% Inputs
% groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
%   'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 7; % Base intersect index %%%
baseName = groupNames{i}; % Base intersect name
% intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
p = 3; % Intersect percentage index %%%
pctName = intersectPctName{p}; % Intersect percentage name
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 1; % type1 %%%
k = 2; % type2 %%%
type1 = intersectType{j};
type2 = intersectType{k};
combi = [type1 '_' type2];

% Video Setup
videoName = ['IntersectComp_' combi '_Base_' baseName '_' pctName '.mp4'];
myVideo = VideoWriter(videoName, 'MPEG-4');
myVideo.FrameRate = 10;
myVideo.Quality = 100;
open(myVideo);

% Figure 
f = figure('Visible', 'off', 'Position', [100, 100, 1920, 1080]);  % Full-HD (Position)
ax1 = subplot(1,2,1, 'Parent', f);
ax2 = subplot(1,2,2, 'Parent', f);
hold(ax2, 'on');

% Z-Loop: Frames
zSteps = z_min:0.5:z_max;
%zSteps = z_min:0.5:-168; 
for idxModel = 1:length(zSteps)
    currentZ = zSteps(idxModel);

    % Draw 3D in first frame or not inside/grid is included
    draw3D = (idxModel == 1) || ~(ismember(type1, {'inside','grid'}) || ismember(type2, {'inside','grid'}));

    % Slice Update
    updateZSlice(ax1, ax2, allPoints, allColors, currentZ, tolerance, ...
        type1, type2, pctName, baseName, markerSize, markerStyle, pelvis, draw3D);

    % Legend
    legendHandles = gobjects(1,5);
    for fIdx = 1:4
        legendHandles(fIdx) = scatter3(ax2, 1e10, 1e10, 1e10, markerSize, ...
            'Marker', markerStyle, ...
            'MarkerFaceColor', colors{fIdx}, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.5, ...
            'DisplayName', labels{fIdx});
    end
    legendHandles(5) = patch(ax2, nan, nan, nan, 'm', ...
        'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', ...
        'DisplayName', 'Slice Plane');
    legend(ax2, legendHandles, ...
        [labels, {'Slice Plane'}], ...
        'Location', 'best', ...
        'Interpreter', 'none');

    drawnow;
    frame = getframe(f);  
    writeVideo(myVideo, frame);
end

close(myVideo);
close(f);  
disp(['Video saved: ', videoName]);

% clear data numIDs colors labels markerSize markerStyle fieldsDefect fieldsVertices fieldsPoints allPoints allColors acentre ...
%     fieldDataV fieldPointsV fieldDataP fieldPointsP fieldData fieldPoints z_min z_max tolerance ax1 ax2 legendHandles slider ...
%    xlim2D ylim2D zlim3D draw3D p f_idx
% clear videoName myVideo f zSteps idx currentZ frame
% clear i j k p pctName baseName combi baseIntersect currentPct1 currentPct2 intersectPctName intersectPercent intersectType ...
%     numBaseIntersect numIntersectPct numIntersectTypes type1 type2 typeName % final clear of intersect

%% Intersection: comparison (points inside) - inside vs grid

% Initialisation
defectInsideMask = false(size(allPelvis.refPoints.gridPoints.pointsBox, 1), numel(pelvisDefect));  
defectInsideCount = zeros(size(allPelvis.refPoints.gridPoints.pointsBox, 1), 1);          
% Grid points inside (all defects)
for i = 1:dataCountDefect
    thisMaskIdx = pelvisDefect(i).volume.gridPoints.insideMaskIdx;
    if isempty(thisMaskIdx)
        continue;
    end
    defectInsideMask(thisMaskIdx, i) = true;
end
allIntersect.inside_grid.gridPoints.defectPresence = defectInsideMask;

% Defect count
defectInsideCount = sum(defectInsideMask, 2);
allIntersect.inside_grid.gridPoints.defectCount = defectInsideCount;
% Points inside reference pelvis and additional points in pelvis defects
combinedMask = allPelvis.refPoints.gridPoints.insideMask | any(defectInsideMask, 2);
combinedMaskIdx = find(combinedMask);
allIntersect.inside_grid.gridPoints.combinedMask = combinedMask;
allIntersect.inside_grid.gridPoints.combinedMaskIdx = combinedMaskIdx;

% Additional points in pelvis defects
addMask = ~allPelvis.refPoints.gridPoints.insideMask & any(defectInsideMask, 2);
addMaskIdx = find(addMask);
addDefectCount = defectInsideCount(addMaskIdx);
allIntersect.inside_grid.gridPoints.addMask = addMask;
allIntersect.inside_grid.gridPoints.addMaskIdx = addMaskIdx;
allIntersect.inside_grid.gridPoints.addDefectCount = addDefectCount;
% Additional points threshold
threshold = 5; % all: 1
validIdx = addDefectCount >= threshold;
addDefectCountThresholdIdx = addMaskIdx(validIdx);
addDefectCountThresholdCount = defectInsideCount(addDefectCountThresholdIdx);
allIntersect.inside_grid.gridPoints.addDefectCountThresholdIdx = addDefectCountThresholdIdx;
allIntersect.inside_grid.gridPoints.addDefectCountThresholdCount = addDefectCountThresholdCount;


% Display (for control)
figure('Name', 'Combined Grid Points');
hold on;

% Reference pelvis
patch('Faces', pelvis(1).import.processed.faces, ...
      'Vertices', pelvis(1).import.processed.vertices, ...
      'FaceColor', TUMcolors.grey20, ...
      'FaceAlpha', 0.25, ...
      'EdgeColor', 'none', ...
      'FaceLighting', 'gouraud', ...
      'AmbientStrength', 0.5);
light('Position', [1 1 5], 'Style', 'infinite');

if isempty(addDefectCountThresholdCount)
    warning('No points with defectCount >= %d', threshold);
else
    nColors = max(addDefectCountThresholdCount)+1;
    viridisMap = viridis(nColors);
    colours = viridisMap(addDefectCountThresholdCount+1, :);  % +1 for 0-based

    % Additional defect points
    scatter3(allPelvis.refPoints.gridPoints.pointsBox(addDefectCountThresholdIdx,1), ...
             allPelvis.refPoints.gridPoints.pointsBox(addDefectCountThresholdIdx,2), ...
             allPelvis.refPoints.gridPoints.pointsBox(addDefectCountThresholdIdx,3), ...
             10, colours, 'filled');

    % Color bar
    colormap(viridisMap);
    cb = colorbar;
    cb.Ticks = linspace(0, 1, nColors);
    cb.TickLabels = 0:max(addDefectCountThresholdCount);
    clim([0 1]);
end

% Format
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Additional Grid Points (threshold ≥ %d)', threshold));
daspect([1 1 1]);
view(3);
grid off;
hold off;

clear thisMaskIdx defectInsideMask defectInsideCount combinedMask combinedMaskIdx addMask addMaskIdx threshold validIdx...
    addDefectCount addDefectCountThresholdIdx addDefectCountThresholdCount nColors viridisMap colours cb

%% Save and load properties of Class Intersection - per intersection type (Defect) 

%%% For intersect type and pairs ('all')
% Save pelvisIntersect 'all' pairs data for each intersect type
savePelvisIntersectPairs = struct();
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 4; %%%
% Get meta-information about fields in pairs
metaPairs = fieldnames(pelvisDefectAnalysis(1).intersect.(intersectType{j}).pairs);
% Iterate over each field in pairs
for i = 1:length(metaPairs)
    fieldName = metaPairs{i};
    % Save the field data to the struct
    savePelvisIntersectPairs.(intersectType{j}).(fieldName) = pelvisDefectAnalysis(1).intersect.(intersectType{j}).pairs.(fieldName); % 'all'
end
% Save data to file
save(['.\pelvisIntersect(' intersectType{j} ')PairsData.mat'], 'savePelvisIntersectPairs', '-v7.3'); % Adapt storage location
% Clear the temporary save structure
clear savePelvisIntersectPairs metaPairs fieldName

% Load pelvisIntersect 'all' pairs data for each intersect type
% intersectType = {'alpha', 'refinedAlpha', 'inside', 'grid'};
j = 1; %%%
load(['.\pelvisIntersect(' intersectType{j} ')PairsData.mat'], 'savePelvisIntersectPairs'); % Adapt storage location
metaPairs = fieldnames(savePelvisIntersectPairs.(intersectType{j}));
% Restore fields to the structure
for i = 1:length(metaPairs)
    fieldName = metaPairs{i};
    pelvisDefectAnalysis(1).intersect.(intersectType{j}).pairs.(fieldName) = savePelvisIntersectPairs.(intersectType{j}).(fieldName);
end
% Clear temporary loaded data
clear savePelvisIntersectPairs metaPairs fieldName


%%% For intersect type
% Save pelvisIntersect data for each intersect type
savePelvisIntersectType = struct();
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 1; %%%
% Meta infos
metaIntersectType = fieldnames(pelvisDefectAnalysis(1).intersect.(intersectType{j}));
% Iterate over each base intersect
for i = 1:numBaseIntersect
    for l = 1:length(metaIntersectType)
        fieldName = metaIntersectType{l};
        % Save the field data to the struct
        savePelvisIntersectType.(intersectType{j})(i).(fieldName) = pelvisDefectAnalysis(i).intersect.(intersectType{j}).(fieldName);
    end
end
% Save data to file
save(['.\pelvisIntersect(' intersectType{j} ').mat'], 'savePelvisIntersectType', '-v7.3'); % Adjust path if needed
% Clear the temporary save structure
clear savePelvisIntersectType metaIntersectType fieldName

% Load pelvisIntersect data for each intersect type
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 3; %%%
load(['C:\Users\ga56man\#Projects\IntersectionPelvis\Workspace\pelvisIntersect(' intersectType{j} ')Reduced.mat']);
metaIntersectType = fieldnames(savePelvisIntersectType.(intersectType{j})(1));
% Iterate over each base intersect
for i = 1:numBaseIntersect
    for l = 1:length(metaIntersectType)
        fieldName = metaIntersectType{l};
        pelvisDefectAnalysis(i).intersect.(intersectType{j}).(fieldName) = ...
            savePelvisIntersectType.(intersectType{j})(i).(fieldName);
    end
end
% Clear the temporary save structure
clear savePelvisIntersectType metaIntersectType fieldName


%%% For intersect type per base
% Save pelvisIntersect data for each intersect type
savePelvisIntersectType = struct();
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 1; %%%
% Meta infos
metaIntersectType = fieldnames(pelvisDefectAnalysis(1).intersect.(intersectType{j}));
% Iterate over each base intersect
i = 1;
for l = 1:length(metaIntersectType)
    fieldName = metaIntersectType{l};
    % Save the field data to the struct
    savePelvisIntersectType.(intersectType{j})(i).(fieldName) = pelvisDefectAnalysis(i).intersect.(intersectType{j}).(fieldName);
end
% Save data to file
save(['.\pelvisIntersect(' intersectType{j} ')' baseIntersect{i} '.mat'], 'savePelvisIntersectType', '-v7.3'); % Adjust path if needed
% Clear the temporary save structure
clear savePelvisIntersectType metaIntersectType fieldName


% Clear for intersect type / pct 
% Keep fields
keepFields = {'filteredFreqDefect', 'filteredFreqVertices', 'filteredFreqPoints', ...
              'usedOrigDefect', 'usedDefect', 'usedVertices', 'usedPoints', ...
              'percent', 'count', ...
              'faces', 'allVerticesPoints'};
%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    for i = 1:numBaseIntersect
        for k = 1:numIntersectPct
            pctName = intersectPctName{k};
            currentStruct = pelvisDefectAnalysis(i).intersect.(typeName).(pctName);
            fieldNames = fieldnames(currentStruct);
            % Delete all other fields
            for f = 1:numel(fieldNames)
                fname = fieldNames{f};
                if ~ismember(fname, keepFields)
                    pelvisDefectAnalysis(i).intersect.(typeName).(pctName) = rmfield(pelvisDefectAnalysis(i).intersect.(typeName).(pctName), fname);
                end
            end
        end
    end
%end
clear keepFields currentStruct fieldNames fieldNames fname f

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Defect Distribution %%%%%%%%%%%%%%
% Defect distribution for intersect type inside

% Pelvis defect categories (all or selected) are base for intersection
baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
numBaseIntersect = length(baseIntersect);

% Intersection type only inside
intersectType = {'inside'};
numIntersectTypes = length(intersectType);
j=1;
for i = 1:numBaseIntersect
    % Base of intersections
    pelvisDefectAnalysis(i).distribut.(intersectType{j}).base = baseIntersect{i};
end
% Pelvis IDs for the categories
pelvisDefectAnalysis(1).distribut.(intersectType{j}).IDs = (2:dataCountDefect)'; % all
pelvisDefectAnalysis(2).distribut.(intersectType{j}).IDs = [2;3;(5:9)';11;12;14;15;(17:21)';23;24;(26:37)';39;40;(42:47)']; % BTE
pelvisDefectAnalysis(3).distribut.(intersectType{j}).IDs = [4;10;13;16;22;25;38;41]; % Other
pelvisDefectAnalysis(4).distribut.(intersectType{j}).IDs = [3;4;5]; % 2A
pelvisDefectAnalysis(5).distribut.(intersectType{j}).IDs = [2;13;14;17;19;25;38;42;45]; % 2B
pelvisDefectAnalysis(6).distribut.(intersectType{j}).IDs = [11;40;44]; % 2C
pelvisDefectAnalysis(7).distribut.(intersectType{j}).IDs = [6;12;20;21;22;26;28;29;30;31;32;34;35;36;39;47;48]; % 3A
pelvisDefectAnalysis(8).distribut.(intersectType{j}).IDs = [7;8;9;10;15;16;18;23;24;27;33;37;41;43;46]; % 3B
pelvisDefectAnalysis(9).distribut.(intersectType{j}).IDs = [3;5]; % 2A & BTE
pelvisDefectAnalysis(10).distribut.(intersectType{j}).IDs = [2;14;17;19;42;45]; % 2B & BTE
pelvisDefectAnalysis(11).distribut.(intersectType{j}).IDs = [11;40;44]; % 2C & BTE
pelvisDefectAnalysis(12).distribut.(intersectType{j}).IDs = [6;12;20;21;26;28;29;30;31;32;34;35;36;39;47;48]; % 3A & BTE
pelvisDefectAnalysis(13).distribut.(intersectType{j}).IDs = [7;8;9;15;18;23;24;27;33;37;43;46]; % 3B & BTE

%% Save and load relevant properties of Class Volume (Defect) for defect distribution

% Relevant Data from class volume
% Save data isPoint/VertexInside
savePelvisDefectDataIsInside = struct();
for i = 1:dataCountDefect
    % Cache 
    savePelvisDefectDataIsInside(i).isVertexInside = pelvisDefect(i).volume.shrink.inside.isVertexInside;
    savePelvisDefectDataIsInside(i).isPointInside = pelvisDefect(i).volume.shrink.inside.isPointInside;
end
save('.\PelvisDefectDataIsInside.mat', 'savePelvisDefectDataIsInside', '-v7.3'); % Adapt storage location
clear savePelvisDefectDataIsInside

% Load data isPoint/VertexInside
loadPelvisDefectDataIsInside = load('.\PelvisDefectDataIsInside.mat', 'savePelvisDefectDataIsInside'); % Adapt storage location
savePelvisDefectDataIsInside = loadPelvisDefectDataIsInside.savePelvisDefectDataIsInside;
for i = 1:dataCountDefect
    %pelvisDefect(i).volume.shrink.inside = struct();
    pelvisDefect(i).volume.shrink.inside.isVertexInside = savePelvisDefectDataIsInside(i).isVertexInside;
    pelvisDefect(i).volume.shrink.inside.isPointInside = savePelvisDefectDataIsInside(i).isPointInside;
end
clear savePelvisDefectDataIsInside loadPelvisDefectDataIsInside

%% Pelvis defect distribution: counting and normalisation

% Pelvis defect categories (all or selected) are base for intersection
baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
numBaseIntersect = length(baseIntersect);
% % Initialisation
% pelvisDistribution = struct();

% ID groups
groups      = cell(numBaseIntersect,1);
for i = 1:numBaseIntersect
    groups{i}     = pelvisDefectAnalysis(i).distribut.(intersectType{j}).IDs;                                     
end
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
    'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};

% Counting (points inside defects)
for i = 1:numBaseIntersect
    groupName = groupNames{i};
    currentGroup = groups{i};
    % Initialisation
    pointsInsideBoxCount = zeros(size(allPelvis.refPoints.gridPoints.pointsBox, 1), 1);
    % Point count
    for k = 1:length(currentGroup)
        id = currentGroup(k);
        pointsInsideBoxCount = pointsInsideBoxCount + pelvisDefect(id).volume.shrink.inside.isPointInside;
    end
    pelvisDefectAnalysis(i).distribut.inside.pointsInsideBoxCount = pointsInsideBoxCount;
    pointsInsideCount = pointsInsideBoxCount(allPelvis.refPoints.gridPoints.insideMask,:);
    pelvisDefectAnalysis(i).distribut.inside.pointsInsideCount = pointsInsideCount;
    % Apply colormap (recommended: viridis)
    colourNum = 32768;
    colourMap = viridis(colourNum);
    maxVal = colourNum;
    viridisCountData = pointsInsideCount / maxVal;
    colourIdx = round(viridisCountData * (colourNum - 1)) + 1;
    rgbColour = colourMap(colourIdx, :);    
    pelvisDefectAnalysis(i).distribut.inside.rgbColour = rgbColour;

    % Linear Normalisation
    maxVal = length(currentGroup); % minVal = 0;
    linNormalData = pointsInsideCount / maxVal;
    pelvisDefectAnalysis(i).distribut.inside.linNormalData = linNormalData;
    % Apply colormap
    colourIdxLin = round(linNormalData * (colourNum - 1)) + 1;
    rgbColourLin = colourMap(colourIdxLin, :);    
    pelvisDefectAnalysis(i).distribut.inside.rgbLinNormal = rgbColourLin;

    % Normalisation for visualisation (percentile)
    lowerPrctile = 70;   %%%
    upperPrctile = 97.5; %%%
    lowerBound = prctile(pointsInsideCount, lowerPrctile);
    upperBound = prctile(pointsInsideCount, upperPrctile);
    prctNormalData = (pointsInsideCount - lowerBound) / (upperBound - lowerBound);
    prctNormalData(prctNormalData < 0) = 0;  % Clamp values below the range to 0
    prctNormalData(prctNormalData > 1) = 1;  % Clamp values above the range to 1
    pelvisDefectAnalysis(i).distribut.inside.lowerPrctile = lowerPrctile;
    pelvisDefectAnalysis(i).distribut.inside.upperPrctile = upperPrctile;
    pelvisDefectAnalysis(i).distribut.inside.lowerBound = lowerBound;
    pelvisDefectAnalysis(i).distribut.inside.upperBound = upperBound;
    pelvisDefectAnalysis(i).distribut.inside.prctNormalData = prctNormalData;
    % Apply colormap
    colourIdxPrct = round(prctNormalData * (colourNum - 1)) + 1;
    rgbColourPrct = colourMap(colourIdxPrct, :);    
    pelvisDefectAnalysis(i).distribut.inside.rgbPrctNormal = rgbColourPrct;

    % Normalisation for visualisation (logarithmic)
    logNormalData = log1p(pointsInsideCount);
    logNormalData  = logNormalData / max(logNormalData);
    pelvisDefectAnalysis(i).distribut.inside.logNormalData = logNormalData;
    % Apply colormap
    colourIdxLog = round(logNormalData * (colourNum - 1)) + 1;
    rgbColourLog = colourMap(colourIdxLog, :);    
    pelvisDefectAnalysis(i).distribut.inside.rgbLogNormal = rgbColourLog;

end
clear id pointsInsideBoxCount pointsInsideCount maxVal linNormalData colourNum colourMap viridisCountData colourIdx rgbColour ...
    colourIdxLin rgbColourLin lowerPrctile upperPrctile lowerBound upperBound prctNormalData colourIdxPrct rgbColourPrct ...
    logNormalData colourIdxLog rgbColourLog

%% Pelvis defect distribution: histogram
% for all (whole data set)

% Histogram 
binEdgesOrig = -0.5:1:(dataCountDefect+0.5); % for original data
binEdgesNorm = (-0.5:1:dataCountDefect+0.5) / dataCountDefect;  

% Plot Histograms
figure('Name', 'Histogram Comparison');
% Original data
subplot(2,2,1);
histogram(pelvisDefectAnalysis(1).distribut.inside.pointsInsideCount, 'BinEdges', binEdgesOrig, 'FaceColor', [0.2 0.6 0.8]);
title('Original Counts');
xlabel('Points Inside'); ylabel('Frequency');
grid on;
% Linear Normalisation
subplot(2,2,2);
histogram(pelvisDefectAnalysis(1).distribut.inside.linNormalData, 'BinEdges', binEdgesNorm, 'FaceColor', [0.8 0.4 0.4]);%
title('Linear Normalized');
xlabel('Normalized Value'); ylabel('Frequency');
grid on;
% Percentil Normalisation
subplot(2,2,3);
histogram(pelvisDefectAnalysis(1).distribut.inside.prctNormalData, 'BinEdges', binEdgesNorm, 'FaceColor', [0.4 0.8 0.4]);
title(['Percentile Normalized (' num2str(pelvisDefectAnalysis(1).distribut.inside.lowerPrctile) '%–' num2str(pelvisDefectAnalysis(1).distribut.inside.upperPrctile) '%)']);
xlabel('Normalized Value'); ylabel('Frequency');
grid on;
% Logarithmic Normalisation
subplot(2,2,4);
histogram(pelvisDefectAnalysis(1).distribut.inside.logNormalData, 'BinEdges', binEdgesNorm, 'FaceColor', [0.6 0.4 0.8]);
title('Logarithmic Normalized');
xlabel('Normalized Value'); ylabel('Frequency');
grid on;

clear binEdgesOrig binEdgesNorm countsOrig countsLin countsPrct countsLog 

%% Display pelvis defect distribution (for control)
% Comparison of data visualisartion/normalisation 

% Groups
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
    'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 1; % group number

% Norm type (options: 'org', 'lin', 'prct', 'log')
normType = 'lin'; %%%
groupData = pelvisDefectAnalysis(i).distribut.inside;
colourMap = viridis(32768);  % Color map

switch normType
    case 'org'
        colourData = groupData.pointsInsideCount;
        cbTicks = 0:5:max(groupData.pointsInsideCount);
        if cbTicks(end) < max(colourData)
            cbTicks(end+1) = max(colourData);
        end
        cbLabels = arrayfun(@(x) sprintf('%d', x), cbTicks, 'UniformOutput', false);
        normName = 'original';
    case 'lin'
        colourData = groupData.linNormalData;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        minVal = min(groupData.pointsInsideCount);
        maxVal = max(groupData.pointsInsideCount);
        tickCounts = minVal + cbTicks * (maxVal - minVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'lin.normalized';
    case 'prct'
        colourData = groupData.prctNormalData;
        ticks = 0:0.25:1;
        lb = groupData.lowerBound;
        ub = groupData.upperBound;
        cbTicks = ticks;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        tickCounts = lb + ticks * (ub - lb);
        cbLabelsRight = arrayfun(@(x) sprintf('%.0f', x), tickCounts, 'UniformOutput', false);
        cbLabelsRight{1} = sprintf('0–%.0f', tickCounts(1));
        cbLabelsRight{end} = sprintf('%.0f–%d', tickCounts(end), dataCountDefect-1);
        normName = 'prct.normalized';
    case 'log'
        colourData = groupData.logNormalData;
        ticks = linspace(0,1,5);
        maxVal = log1p(max(groupData.pointsInsideCount));
        cbTicks = ticks;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        tickCounts = expm1(ticks * maxVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'log.normalized';
    otherwise
        error('Invalid normalisation type.');
end

% Plot
%figure('Visible','off')
figure('Name', [groupNames{i} ' - ' normName], 'NumberTitle', 'off');
hold on

% Reference pelvis
patch('Faces', pelvis(1).import.processed.faces,...
    'Vertices', pelvis(1).import.processed.vertices,...
    'FaceColor', TUMcolors.grey20, ...
    'FaceAlpha', 0.5,...
    'EdgeColor', 'none', ...
    'EdgeAlpha', 0.25, ...
    'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.5);
light('Position', [1 1 5], 'Style', 'infinite');

% Acetabulum Centre
plot3(pelvis(1).import.processed.acentre(1), ...
      pelvis(1).import.processed.acentre(2), ...
      pelvis(1).import.processed.acentre(3), '.', ...
      'Color', 'm', 'MarkerSize', 50);

% Defect Distribution Points
scatter3(allPelvis.refPoints.gridPoints.inside(:,1), ...
         allPelvis.refPoints.gridPoints.inside(:,2), ...
         allPelvis.refPoints.gridPoints.inside(:,3), ...
         10, colourData, 'filled');

% Colorbar
colormap(colourMap);
cb = colorbar;
cb.Ticks = cbTicks;

if ismember(normType, {'org'})
    cb.TickLabels = cbLabels;
elseif ismember(normType, {'lin', 'prct', 'log'})
    cb.TickLabels = cbLabelsLeft;
    cb.TickDirection = 'out';
    % Right-side annotation for real count values
    cbPos = cb.Position;
    for t = 1:length(cbTicks)
        y = cbPos(2) + cbTicks(t) * cbPos(4);
        annotation('textbox', [cbPos(1)+cbPos(3)+0.06, y-0.01, 0.05, 0.03], ...
            'String', cbLabelsRight{t}, 'EdgeColor', 'none', 'FontSize', 9);
    end
    annotation('textbox', [cbPos(1)+cbPos(3)+0.06, cbPos(2)+cbPos(4)+0.015, 0.05, 0.03], ...
        'String', 'Counts', 'EdgeColor', 'none', 'FontWeight', 'bold');
end


% Format
title(['Point Inside Count for ' groupNames{i} ' (' normName ' data)']);
xlabel('X'); ylabel('Y'); zlabel('Z');
daspect([1 1 1]);
view(3);
%view(ax2,80, -10) % Figure Paper

grid off;
hold off;

% Save figure (figure unvisible, but saved visible)
%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%savefig(['./Figures/PelvisDefectDistribution(', groupNames{i}, ')-', normName, '.fig'])

clear normType groupData colourMap colourData cbTicks cbLabels normName ticks lb ub cbLabelsLeft tickCounts cbLabelsRight ...
    maxVal cb cbPos y

%% Display pelvis defect distribution (3D) in slices-axial view (for control)

% Parameters
normType = 'lin';  % Options: 'org', 'lin', 'prct', 'log'
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
    'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 1;
groupData = pelvisDefectAnalysis(i).distribut.inside;
% Highlight intersection points (pct)
useHighlight = false;  %%%
highlightIndices = [];
if useHighlight
    intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
    k = 4; %%%
    numL = numel(pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints);
    for l = 1:numL
        thisCell = pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints{l}.PointBaseIndex;
        if ~isempty(thisCell)
            highlightIndices = [highlightIndices; thisCell(:)];
        end
    end
    highlightIndices = unique(highlightIndices);
end


% Select data based on normalisation
colourMap = viridis(32768);
switch normType
    case 'org'
        colourData = groupData.pointsInsideCount;
        cbTicks = 0:5:max(groupData.pointsInsideCount);
        if cbTicks(end) < max(colourData)
            cbTicks(end+1) = max(colourData);
        end
        cbLabels = arrayfun(@(x) sprintf('%d', x), cbTicks, 'UniformOutput', false);
        normName = 'original';
    case 'lin'
        colourData = groupData.linNormalData;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        minVal = min(groupData.pointsInsideCount);
        maxVal = max(groupData.pointsInsideCount);
        tickCounts = minVal + cbTicks * (maxVal - minVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'lin.normalized';
    case 'prct'
        colourData = groupData.prctNormalData;
        ticks = 0:0.25:1;
        lb = groupData.lowerBound;
        ub = groupData.upperBound;
        cbTicks = ticks;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        tickCounts = lb + ticks * (ub - lb);
        cbLabelsRight = arrayfun(@(x) sprintf('%.0f', x), tickCounts, 'UniformOutput', false);
        cbLabelsRight{1} = sprintf('0–%.0f', tickCounts(1));
        cbLabelsRight{end} = sprintf('%.0f–%d', tickCounts(end), dataCountDefect-1);
        normName = 'prct.normalized';
    case 'log'
        colourData = groupData.logNormalData;
        ticks = linspace(0,1,5);
        maxVal = log1p(max(groupData.pointsInsideCount));
        cbTicks = ticks;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        tickCounts = expm1(ticks * maxVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'log.normalized';
    otherwise
        error('Invalid normalisation type.');
end


% Create figure
figure('Name','Z-Slice View','NumberTitle','off');
hold on;

% Plot reference pelvis
patch('Faces', pelvis(1).import.processed.faces, ...
      'Vertices', pelvis(1).import.processed.vertices, ...
      'FaceColor', TUMcolors.grey20, ...
      'FaceAlpha', 0.25, ...
      'EdgeColor', 'none', ...
      'FaceLighting', 'gouraud', ...
      'AmbientStrength', 0.5);
light('Position', [1 1 5], 'Style', 'infinite');

% Plot acetabulum center
plot3(pelvis(1).import.processed.acentre(1), ...
      pelvis(1).import.processed.acentre(2), ...
      pelvis(1).import.processed.acentre(3), ...
      '*', 'Color', 'm', 'MarkerSize', 10); % pelvisDefect(1).import.processed.acentre

% Plot defect points
scatter3(allPelvis.refPoints.gridPoints.inside(:,1), ...
         allPelvis.refPoints.gridPoints.inside(:,2), ...
         allPelvis.refPoints.gridPoints.inside(:,3), ...
         5, colourData, 'filled');

% Colorbar
if ismember(normType, {'lin', 'prct', 'log'})
    clim([0 1]);
else
    clim([0 max(colourData)]);
end 
colourMap = viridis(32768);  % Color map
colormap(colourMap);
cb = colorbar;
cb.Ticks = cbTicks;
if ismember(normType, {'org'})
    cb.TickLabels = cbLabels;
elseif ismember(normType, {'lin', 'prct', 'log'})
    cb.TickLabels = cbLabelsLeft;
    cb.TickDirection = 'out';
    % Right-side annotation for real count values
    cbPos = cb.Position;
    for t = 1:length(cbTicks)
        y = cbPos(2) + cbTicks(t) * cbPos(4);
        annotation('textbox', [cbPos(1)+cbPos(3)+0.07, y-0.01, 0.05, 0.03], ...
            'String', cbLabelsRight{t}, 'EdgeColor', 'none', 'FontSize', 9);
    end
    annotation('textbox', [cbPos(1)+cbPos(3)+0.07, cbPos(2)+cbPos(4)+0.015, 0.05, 0.03], ...
        'String', 'Counts', 'EdgeColor', 'none', 'FontWeight', 'bold');
end

% Format
xlabel('X'); ylabel('Y'); zlabel('Z');
if useHighlight
    title(['Point Inside Count for ' groupNames{i} ' and ' intersectPctName{k} ' (' normName ' data)']);
else
    title(['Point Inside Count for ' groupNames{i} ' (' normName ' data)']);
end
view(3);
%view(ax2,80, -10) % Figure Paper
daspect([1 1 1]);
grid off;
hold off;

% Add slider for Z-slice navigation
slider = uicontrol('Style', 'slider', ...
                   'Min', min(allPelvis.refPoints.gridPoints.inside(:,3)), ...
                   'Max', max(allPelvis.refPoints.gridPoints.inside(:,3)), ...
                   'Value', mean(allPelvis.refPoints.gridPoints.inside(:,3)), ...
                   'Units', 'normalized', ...
                   'Position', [0.2 0.01 0.6 0.05], ...
                   'Callback', @(src, ~) pelvisDefectAnalysis(i).distribut.updateSlice(src, ...
                   allPelvis.refPoints.gridPoints.inside, colourData, pelvis, ...
                   highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx));

% Save figure (figure unvisible, but saved visible)
%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%savefig(['./Figures/PelvisDefectDistribution(', groupNames{i}, ')-', normName, '_3Dslice.fig'])

clear normType groupData colourMap colourData cbTicks cbLabels normName ticks lb ub cbLabelsLeft tickCounts cbLabelsRight ...
    maxVal cb cbPos y slider thisCell useHighlight highlightIndices numL minVal

%% Video: Display pelvis defect distribution (3D) in slices-axial view (for control) 

% Pelvis group
normType = 'log';  % Options: 'org', 'lin', 'prct', 'log' %%%
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
    'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 1;
groupName = groupNames{i};
groupData = pelvisDefectAnalysis(i).distribut.inside; 

% Highlight intersection points
useHighlight = false; %%%
highlightIndices = [];
if useHighlight
    intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
    k = 4; %%%
    numL = numel(pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints);
    for l = 1:numL
        thisCell = pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints{l}.PointBaseIndex;
        if ~isempty(thisCell)
            highlightIndices = [highlightIndices; thisCell(:)];
        end
    end
    highlightIndices = unique(highlightIndices);
end

% Select data based on normalisation
colourMap = viridis(32768);
switch normType
    case 'org'
        colourData = groupData.pointsInsideCount;
        cbTicks = 0:5:max(groupData.pointsInsideCount);
        if cbTicks(end) < max(colourData)
            cbTicks(end+1) = max(colourData);
        end
        cbLabels = arrayfun(@(x) sprintf('%d', x), cbTicks, 'UniformOutput', false);
        normName = 'original';
    case 'lin'
        colourData = groupData.linNormalData;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        minVal = min(groupData.pointsInsideCount);
        maxVal = max(groupData.pointsInsideCount);
        tickCounts = minVal + cbTicks * (maxVal - minVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'lin.normalized';
    case 'prct'
        colourData = groupData.prctNormalData;
        lb = groupData.lowerBound;
        ub = groupData.upperBound;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        tickCounts = lb + cbTicks * (ub - lb);
        cbLabelsRight = arrayfun(@(x) sprintf('%.0f', x), tickCounts, 'UniformOutput', false);
        cbLabelsRight{1} = sprintf('0–%.0f', tickCounts(1));
        cbLabelsRight{end} = sprintf('%.0f–%d', tickCounts(end), dataCountDefect-1);
        normName = 'prct.normalized';
    case 'log'
        colourData = groupData.logNormalData;
        maxVal = log1p(max(groupData.pointsInsideCount));
        cbTicks = linspace(0,1,5);
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        tickCounts = expm1(cbTicks * maxVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'log.normalized';
    otherwise
        error('Invalid normalisation type.');
end

% Global limits
x_min = min(allPelvis.refPoints.gridPoints.inside(:,1));
x_max = max(allPelvis.refPoints.gridPoints.inside(:,1));
y_min = min(allPelvis.refPoints.gridPoints.inside(:,2));
y_max = max(allPelvis.refPoints.gridPoints.inside(:,2));
z_min = min(allPelvis.refPoints.gridPoints.inside(:,3));
z_max = max(allPelvis.refPoints.gridPoints.inside(:,3));

% Video Setup
videoName = ['PelvisDistribution3D_' groupName '_' normName '.mp4'];
myVideo   = VideoWriter(videoName, 'MPEG-4');
myVideo.FrameRate = 10;
myVideo.Quality   = 100; 
open(myVideo);

% Figure
if useHighlight
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off', 'Position', [100 100 1920 1080]);  % Full-HD
    sgtitle(['Pelvis Group: ' groupName ' and ' intersectPctName{k} ' (' normName ')']);
else
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off', 'Position', [100 100 1920 1080]);  % Full-HD
    sgtitle(['Pelvis Group: ' groupName ' (' normName ')']);
end

% Add slider for Z-slice navigation
slider = uicontrol('Style', 'slider', ...
                   'Min', min(allPelvis.refPoints.gridPoints.inside(:,3)), ...
                   'Max', max(allPelvis.refPoints.gridPoints.inside(:,3)), ...
                   'Value', mean(allPelvis.refPoints.gridPoints.inside(:,3)), ...
                   'Units', 'normalized', ...
                   'Position', [0.2 0.01 0.6 0.05], ...
                   'Callback', @(src, ~) pelvisDefectAnalysis(i).distribut.updateSlice(src, ...
                   allPelvis.refPoints.gridPoints.inside, colourData, pelvis, ...
                   highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx));

% Save figure handle
guidata(f, struct('slider', slider, ...
    'D', allPelvis.refPoints.gridPoints.inside, ...
    'rgbColour', colourData, ...
    'x_min', x_min, 'x_max', x_max, ...
    'y_min', y_min, 'y_max', y_max, ...
    'z_min', z_min, 'z_max', z_max));
% First update (for initial view)
pelvisDefectAnalysis(i).distribut.updateSlice(slider, allPelvis.refPoints.gridPoints.inside, colourData, pelvis, highlightIndices, ...
    allPelvis.refPoints.gridPoints.insideMaskIdx);

% Colorbar
if ismember(normType, {'lin', 'prct', 'log'})
    clim([0 1]);
else
    clim([0 max(colourData)]);
end 
colormap(colourMap);
cb = colorbar;
cb.Ticks = cbTicks;
if ismember(normType, {'org'})
    cb.TickLabels = cbLabels;
elseif ismember(normType, {'lin', 'prct', 'log'})
    cb.TickLabels = cbLabelsLeft;
    cb.TickDirection = 'out';
    % Right-side annotation for real count values
    cbPos = cb.Position;
    for t = 1:length(cbTicks)
        y = cbPos(2) + cbTicks(t) * cbPos(4);
        annotation('textbox', [cbPos(1)+cbPos(3)+0.07, y-0.01, 0.05, 0.03], ...
            'String', cbLabelsRight{t}, 'EdgeColor', 'none', 'FontSize', 9);
    end
    annotation('textbox', [cbPos(1)+cbPos(3)+0.07, cbPos(2)+cbPos(4)+0.015, 0.05, 0.03], ...
        'String', 'Counts', 'EdgeColor', 'none', 'FontWeight', 'bold');
end

% Video
zVals = z_min:0.5:z_max; 
%zVals = z_min:0.5:-212; % for tests
for zPos = zVals
    slider.Value = zPos;
    pelvisDefectAnalysis(i).distribut.updateSlice(slider, allPelvis.refPoints.gridPoints.inside, colourData, pelvis, highlightIndices, ...
        allPelvis.refPoints.gridPoints.insideMaskIdx);

    drawnow;
    frame = getframe(f);
    writeVideo(myVideo, frame);
end

close(myVideo);
close(f);
disp(['Video saved: ', videoName]);

clear normType groupData colourMap colourData cbTicks cbLabels normName ticks lb ub cbLabelsLeft tickCounts cbLabelsRight ...
    maxVal cb cbPos y slider thisCell useHighlight highlightIndices numL minVal zVals zPos videoName myVideo frame

%% Display pelvis defect distribution in slices: all views 

% Pelvis group
normType = 'log';  % Options: 'org', 'lin', 'prct', 'log'
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
    'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 7;
groupName = groupNames{i};
groupData = pelvisDefectAnalysis(i).distribut.inside;

% Highlight intersection points (pct)
useHighlight = false;  %%%
highlightIndices = [];
if useHighlight
    intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
    k = 5; %%%
    numL = numel(pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints);
    for l = 1:numL
        thisCell = pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints{l}.PointBaseIndex;
        if ~isempty(thisCell)
            highlightIndices = [highlightIndices; thisCell(:)];
        end
    end
    highlightIndices = unique(highlightIndices);
end

% Select data based on normalisation
colourMap = viridis(32768);
switch normType
    case 'org'
        colourData = groupData.pointsInsideCount;
        cbTicks = 0:5:max(groupData.pointsInsideCount);
        if cbTicks(end) < max(colourData)
            cbTicks(end+1) = max(colourData);
        end
        cbLabels = arrayfun(@(x) sprintf('%d', x), cbTicks, 'UniformOutput', false);
        normName = 'original';
        rgbColour = groupData.rgbColour;
    case 'lin'
        colourData = groupData.linNormalData;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        minVal = min(groupData.pointsInsideCount);
        maxVal = max(groupData.pointsInsideCount);
        tickCounts = minVal + cbTicks * (maxVal - minVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'lin.normalized';
        rgbColour = groupData.rgbLinNormal;
    case 'prct'
        colourData = groupData.prctNormalData;
        ticks = 0:0.25:1;
        lb = groupData.lowerBound;
        ub = groupData.upperBound;
        cbTicks = ticks;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        tickCounts = lb + ticks * (ub - lb);
        cbLabelsRight = arrayfun(@(x) sprintf('%.0f', x), tickCounts, 'UniformOutput', false);
        cbLabelsRight{1} = sprintf('0–%.0f', tickCounts(1));
        cbLabelsRight{end} = sprintf('%.0f–%d', tickCounts(end), dataCountDefect-1);
        normName = 'prct.normalized';
        rgbColour = groupData.rgbPrctNormal;
    case 'log'
        colourData = groupData.logNormalData;
        ticks = linspace(0,1,5);
        maxVal = log1p(max(groupData.pointsInsideCount));
        cbTicks = ticks;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        tickCounts = expm1(ticks * maxVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'log.normalized';
        rgbColour = groupData.rgbLogNormal;
    otherwise
        error('Invalid normalisation type.');
end

% Global limits
x_min = min(allPelvis.refPoints.gridPoints.inside(:,1));
x_max = max(allPelvis.refPoints.gridPoints.inside(:,1));
y_min = min(allPelvis.refPoints.gridPoints.inside(:,2));
y_max = max(allPelvis.refPoints.gridPoints.inside(:,2));
z_min = min(allPelvis.refPoints.gridPoints.inside(:,3));
z_max = max(allPelvis.refPoints.gridPoints.inside(:,3));

% Create figure and subplots
if useHighlight
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off');
    sgtitle(['Pelvis Group: ' groupName ' and ' intersectPctName{k} ' (' normName ')']);
else
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off');
    sgtitle(['Pelvis Group: ' groupName ' (' normName ')']);
end
% Subplots
ax1 = subplot(2, 2, 1); % Axial 
ax2 = subplot(2, 2, 2); % Sagittal
ax3 = subplot(2, 2, 3); % Coronal 
ax4 = subplot(2, 2, 4); % 3D 

% Subplot 1: Axial (X-Y-Ebene)
axialSlider = uicontrol('Style', 'slider', ...
    'Min', z_min, 'Max', z_max, ...
    'Value', pelvis(1).import.processed.acentre(3), ... % pelvis(1).import.processed.acentre(3) % mean([z_min z_max])
    'Units', 'normalized', ...
    'Position', [0.1 0.01 0.2 0.05], ...
    'Callback', @(src, event) pelvisDefectAnalysis(i).distribut.updateViews(src, ax1, ax4, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'axial',...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx));
uicontrol('Style', 'text', ...
    'String', 'Axial', ...
    'Units', 'normalized', ...
    'Position', [0.1 0.065 0.2 0.02], ...
    'FontWeight', 'bold', ...
    'BackgroundColor', f.Color);

% Subplot 2: Sagittal (Y-Z-Ebene)
sagittalSlider = uicontrol('Style', 'slider', ...
    'Min', x_min, 'Max', x_max, ...
    'Value', pelvis(1).import.processed.acentre(1), ... % pelvis(1).import.processed.acentre(1) % mean([x_min x_max])
    'Units', 'normalized', ...
    'Position', [0.4 0.01 0.2 0.05], ...
    'Callback', @(src, event) pelvisDefectAnalysis(i).distribut.updateViews(src, ax2, ax4, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'sagittal',...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx));
uicontrol('Style', 'text', ...
    'String', 'Sagittal', ...
    'Units', 'normalized', ...
    'Position', [0.4 0.065 0.2 0.02], ...
    'FontWeight', 'bold', ...
    'BackgroundColor', f.Color);

% Subplot 3: Coronal (X-Z-Ebene)
coronalSlider = uicontrol('Style', 'slider', ...
    'Min', y_min, 'Max', y_max, ...
    'Value', pelvis(1).import.processed.acentre(2), ... % pelvis(1).import.processed.acentre(2) % mean([y_min y_max])
    'Units', 'normalized', ...
    'Position', [0.7 0.01 0.2 0.05], ...
    'Callback', @(src, event) pelvisDefectAnalysis(i).distribut.updateViews(src, ax3, ax4, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'coronal', ...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx));
uicontrol('Style', 'text', ...
    'String', 'Coronal', ...
    'Units', 'normalized', ...
    'Position', [0.7 0.065 0.2 0.02], ...
    'FontWeight', 'bold', ...
    'BackgroundColor', f.Color);

% Save handles
guidata(f, struct('sagittalSlider', sagittalSlider, ...
                  'coronalSlider', coronalSlider, ...
                  'axialSlider', axialSlider, ...
                  'D', allPelvis.refPoints.gridPoints.inside, ...
                  'rgbColour', rgbColour));

% Subplot 4: 3D View
axes(ax4); % subplot 4 active
patch('Faces', pelvis(1).import.processed.faces, ...
      'Vertices', pelvis(1).import.processed.vertices, ...
      'FaceColor', TUMcolors.grey20, ...
      'FaceAlpha', 0.5, ...
      'EdgeColor', 'none', ...
      'FaceLighting', 'gouraud', ...
      'AmbientStrength', 0.5);
light('Position', [1 1 5], 'Style', 'infinite');
scatter3(ax4, NaN, NaN, NaN, 1, colourData(1), 'filled', 'Visible', 'off'); % Dummy points for colorbar
% Colorbar
colormap(ax4, viridis(32768));
if ismember(normType, {'org'})
    clim(ax4, [0 max(colourData)]);
else
    clim(ax4, [0 1]);
end
cb = colorbar(ax4);
cb.Ticks = cbTicks;
if ismember(normType, {'org'})
    cb.TickLabels = cbLabels;
elseif ismember(normType, {'lin', 'prct', 'log'})
    cb.TickLabels = cbLabelsLeft;
    cb.TickDirection = 'out';
    % Right-side annotation for real count values
    cbPos = cb.Position;
    for t = 1:length(cbTicks)
        y = cbPos(2) + cbTicks(t) * cbPos(4);
        annotation('textbox', [cbPos(1)+cbPos(3)+0.07, y-0.01, 0.05, 0.03], ...
            'String', cbLabelsRight{t}, 'EdgeColor', 'none', 'FontSize', 9);
    end
    annotation('textbox', [cbPos(1)+cbPos(3)+0.07, cbPos(2)+cbPos(4)+0.015, 0.05, 0.03], ...
        'String', 'Counts', 'EdgeColor', 'none', 'FontWeight', 'bold');
end

% Format
xlabel('X'); ylabel('Y'); zlabel('Z');
if useHighlight
    title(['Point Inside Count for ' groupNames{i} ' and ' intersectPctName{k} ' (' normName ' data)']);
else
    title(['Point Inside Count for ' groupNames{i} ' (' normName ' data)']);
end
view(3);
%view(ax2,80, -10) % Figure Paper
daspect([1 1 1]);
grid off;
hold off;

% Update all views
pelvisDefectAnalysis(i).distribut.updateViews(axialSlider, ax1, ax4, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'axial',...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx);
pelvisDefectAnalysis(i).distribut.updateViews(sagittalSlider, ax2, ax4, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'sagittal',...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx);
pelvisDefectAnalysis(i).distribut.updateViews(coronalSlider, ax3, ax4, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'coronal',...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx);

% Save figure (figure unvisible, but saved visible)
%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%savefig(['./Figures/PelvisDefectDistribution(', groupNames{i}, ')-', normName, 'sliceAll.fig'])

clear normType groupData colourMap colourData cbTicks cbLabels normName ticks lb ub cbLabelsLeft tickCounts cbLabelsRight ...
    maxVal cb cbPos y sagittalSlider coronalSlider axialSlider ax1 ax2 ax3 ax4 rgbColour x_min x_max y_min y_max z_min z_max...
    thisCell useHighlight highlightIndices numL minVal

%% Display pelvis defect distribution in sagittal view / Y-Z-plane (for control)

% Pelvis group
normType = 'log';  % Options: 'org', 'lin', 'prct', 'log'
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
    'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 7;
groupName = groupNames{i};
groupData = pelvisDefectAnalysis(i).distribut.inside;

% Highlight intersection points
useHighlight = false; %%%
highlightIndices = [];
if useHighlight
    intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
    k = 4; %%%
    numL = numel(pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints);
    for l = 1:numL
        thisCell = pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints{l}.PointBaseIndex;
        if ~isempty(thisCell)
            highlightIndices = [highlightIndices; thisCell(:)];
        end
    end
    highlightIndices = unique(highlightIndices);
end

% Select data based on normalisation
colourMap = viridis(32768);
switch normType
    case 'org'
        colourData = groupData.pointsInsideCount;
        cbTicks = 0:5:max(groupData.pointsInsideCount);
        if cbTicks(end) < max(colourData)
            cbTicks(end+1) = max(colourData);
        end
        cbLabels = arrayfun(@(x) sprintf('%d', x), cbTicks, 'UniformOutput', false);
        normName = 'original';
        rgbColour = groupData.rgbColour;
    case 'lin'
        colourData = groupData.linNormalData;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        minVal = min(groupData.pointsInsideCount);
        maxVal = max(groupData.pointsInsideCount);
        tickCounts = minVal + cbTicks * (maxVal - minVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'lin.normalized';
        rgbColour = groupData.rgbLinNormal;
    case 'prct'
        colourData = groupData.prctNormalData;
        ticks = 0:0.25:1;
        lb = groupData.lowerBound;
        ub = groupData.upperBound;
        cbTicks = ticks;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        tickCounts = lb + ticks * (ub - lb);
        cbLabelsRight = arrayfun(@(x) sprintf('%.0f', x), tickCounts, 'UniformOutput', false);
        cbLabelsRight{1} = sprintf('0–%.0f', tickCounts(1));
        cbLabelsRight{end} = sprintf('%.0f–%d', tickCounts(end), dataCountDefect-1);
        normName = 'prct.normalized';
        rgbColour = groupData.rgbPrctNormal;
    case 'log'
        colourData = groupData.logNormalData;
        ticks = linspace(0,1,5);
        maxVal = log1p(max(groupData.pointsInsideCount));
        cbTicks = ticks;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        tickCounts = expm1(ticks * maxVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'log.normalized';
        rgbColour = groupData.rgbLogNormal;
    otherwise
        error('Invalid normalisation type.');
end

% Global limits
x_min = min(allPelvis.refPoints.gridPoints.inside(:,1));
x_max = max(allPelvis.refPoints.gridPoints.inside(:,1));
y_min = min(allPelvis.refPoints.gridPoints.inside(:,2));
y_max = max(allPelvis.refPoints.gridPoints.inside(:,2));
z_min = min(allPelvis.refPoints.gridPoints.inside(:,3));
z_max = max(allPelvis.refPoints.gridPoints.inside(:,3));

% Create figure and subplots
if useHighlight
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off');
    sgtitle(['Pelvis Group: ' groupName ' and ' intersectPctName{k} ' (' normName ')']);
else
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off');
    sgtitle(['Pelvis Group: ' groupName ' (' normName ')']);
end
ax1 = subplot(1, 2, 1); % Sagittal view
ax2 = subplot(1, 2, 2); % 3D view

% Sagittal view with slider (Y-Z-plane)
stepSize = 0.5;
range = x_max - x_min;
sagittalSlider = uicontrol('Style', 'slider', ...
    'Min', x_min, 'Max', x_max, ...
    'Value', pelvis(1).import.processed.acentre(1), ... % pelvis(1).import.processed.acentre(1) % mean([x_min x_max])
    'SliderStep', [stepSize / range, 0.1], ... 
    'Units', 'normalized', ... 
    'Position', [0.2 0.01 0.2 0.05], ... 
    'Callback', @(src, event) pelvisDefectAnalysis(i).distribut.updateViews(src, ax1, ax2, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'sagittal', ...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx));

% Saving the sliders in the figure's handles
guidata(f, struct('sagittalSlider', sagittalSlider, ...
    'D', allPelvis.refPoints.gridPoints.inside, 'rgbColour', rgbColour, ...
    'x_min', x_min, 'x_max', x_max, 'y_min', y_min, 'y_max', y_max, 'z_min', z_min, 'z_max', z_max));
% Update view
pelvisDefectAnalysis(i).distribut.updateViews(sagittalSlider, ax1, ax2, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'sagittal',...
   allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx);

% Colorbar
hold on
scatter3(ax2, 0, 0, 0, 1, 0, 'Visible', 'off'); % Dummy points for colorbar
colormap(ax2, viridis(32768));
if ismember(normType, {'org'})
    clim(ax2, [0 max(colourData)]);
else
    clim(ax2, [0 1]);
end
cb = colorbar(ax2); 
cb.Ticks = cbTicks;
if ismember(normType, {'org'})
    cb.TickLabels = cbLabels;
elseif ismember(normType, {'lin', 'prct', 'log'})
    cb.TickLabels = cbLabelsLeft;
    cb.TickDirection = 'out';
    % Right-side annotation for real count values
    cbPos = cb.Position;
    for t = 1:length(cbTicks)
        y = cbPos(2) + cbTicks(t) * cbPos(4);
        annotation('textbox', [cbPos(1)+cbPos(3)+0.07, y-0.01, 0.05, 0.03], ...
            'String', cbLabelsRight{t}, 'EdgeColor', 'none', 'FontSize', 9);
    end
    annotation('textbox', [cbPos(1)+cbPos(3)+0.07, cbPos(2)+cbPos(4)+0.015, 0.05, 0.03], ...
        'String', 'Counts', 'EdgeColor', 'none', 'FontWeight', 'bold');
end

% Save figure (figure unvisible, but saved visible)
%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%savefig(['./Figures/PelvisDefectDistribution(', groupNames{i}, ')-', normName, 'sliceSagittal.fig'])

clear normType groupData colourMap colourData cbTicks cbLabels normName ticks lb ub cbLabelsLeft tickCounts cbLabelsRight ...
    maxVal cb cbPos y sagittalSlider ax1 ax2 ax3 ax4 rgbColour x_min x_max y_min y_max z_min z_max...
    stepSize range thisCell useHighlight highlightIndices numL minVal

%% Video: Display pelvis defect distribution in sagittal view / Y-Z-plane (for control)

% Pelvis group
normType = 'lin';  % Options: 'org', 'lin', 'prct', 'log' %%%
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
    'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 7;
groupName = groupNames{i};
groupData = pelvisDefectAnalysis(i).distribut.inside;

% Highlight intersection points
useHighlight = true; %%%
highlightIndices = [];
if useHighlight
    intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
    k = 4; %%%
    numL = numel(pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints);
    for l = 1:numL
        thisCell = pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints{l}.PointBaseIndex;
        if ~isempty(thisCell)
            highlightIndices = [highlightIndices; thisCell(:)];
        end
    end
    highlightIndices = unique(highlightIndices);
end

% Select data based on normalisation
colourMap = viridis(32768);
switch normType
    case 'org'
        colourData = groupData.pointsInsideCount;
        cbTicks = 0:5:max(groupData.pointsInsideCount);
        if cbTicks(end) < max(colourData)
            cbTicks(end+1) = max(colourData);
        end
        cbLabels = arrayfun(@(x) sprintf('%d', x), cbTicks, 'UniformOutput', false);
        normName = 'original';
        rgbColour = groupData.rgbColour; 
    case 'lin'
        colourData = groupData.linNormalData;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        minVal = min(groupData.pointsInsideCount);
        maxVal = max(groupData.pointsInsideCount);
        tickCounts = minVal + cbTicks * (maxVal - minVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'lin.normalized';
        rgbColour = groupData.rgbLinNormal;
    case 'prct'
        colourData = groupData.prctNormalData;
        lb = groupData.lowerBound;
        ub = groupData.upperBound;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        tickCounts = lb + cbTicks * (ub - lb);
        cbLabelsRight = arrayfun(@(x) sprintf('%.0f', x), tickCounts, 'UniformOutput', false);
        cbLabelsRight{1} = sprintf('0–%.0f', tickCounts(1));
        cbLabelsRight{end} = sprintf('%.0f–%d', tickCounts(end), dataCountDefect-1);
        normName = 'prct.normalized';
        rgbColour = groupData.rgbPrctNormal;
    case 'log'
        colourData = groupData.logNormalData;
        maxVal = log1p(max(groupData.pointsInsideCount));
        cbTicks = linspace(0,1,5);
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        tickCounts = expm1(cbTicks * maxVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'log.normalized';
        rgbColour = groupData.rgbLogNormal;
    otherwise
        error('Invalid normalisation type.');
end

% Global limits
x_min = min(allPelvis.refPoints.gridPoints.inside(:,1));
x_max = max(allPelvis.refPoints.gridPoints.inside(:,1));
y_min = min(allPelvis.refPoints.gridPoints.inside(:,2));
y_max = max(allPelvis.refPoints.gridPoints.inside(:,2));
z_min = min(allPelvis.refPoints.gridPoints.inside(:,3));
z_max = max(allPelvis.refPoints.gridPoints.inside(:,3));

% Video Setup
videoName = ['PelvisDistributionSagittal_' groupName '_' normName '.mp4'];
myVideo   = VideoWriter(videoName, 'MPEG-4');
myVideo.FrameRate = 10;
myVideo.Quality   = 100; 
open(myVideo);

% Figure
if useHighlight
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off', 'Position', [100 100 1920 1080]);  % Full-HD
    sgtitle(['Pelvis Group: ' groupName ' and ' intersectPctName{k} ' (' normName ')']);
else
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off', 'Position', [100 100 1920 1080]);  % Full-HD
    sgtitle(['Pelvis Group: ' groupName ' (' normName ')']);
end
ax1 = subplot(1, 2, 1); % Sagittal view
ax2 = subplot(1, 2, 2); % 3D view

% Slider 
stepSize = 0.5;
range = x_max - x_min;
sagittalSlider = uicontrol('Style', 'slider', ...
    'Min', x_min, ...
    'Max', x_max, ...
    'Value', mean([x_min x_max]), ... % pelvis(1).import.processed.acentre(1) % mean([x_min x_max])
    'SliderStep', [stepSize / range, 0.1], ... 
    'Units', 'normalized', ...
    'Position', [0.2 0.01 0.2 0.05], ...
    'Callback', @(src, event) pelvisDefectAnalysis(i).distribut.updateViews(src, ax1, ax2, pelvis, ...
                                         x_min, x_max, y_min, y_max, ...
                                         z_min, z_max, 'sagittal',...
                                         allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, ...
                                         highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx));
% Save figure handle
guidata(f, struct('sagittalSlider', sagittalSlider, ...
    'D', allPelvis.refPoints.gridPoints.inside, ...
    'rgbColour', rgbColour, ...
    'x_min', x_min, 'x_max', x_max, ...
    'y_min', y_min, 'y_max', y_max, ...
    'z_min', z_min, 'z_max', z_max));
% First update (for initial view)
pelvisDefectAnalysis(i).distribut.updateViews(sagittalSlider, ax1, ax2, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'sagittal',...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx);

% Colorbar
hold on
scatter3(ax2, 0, 0, 0, 1, 0, 'Visible', 'off'); % Dummy points for colorbar
colormap(ax2, viridis(32768));
if ismember(normType, {'org'})
    clim(ax2, [0 max(colourData)]);
else
    clim(ax2, [0 1]);
end
cb = colorbar(ax2); 
cb.Ticks = cbTicks;
if ismember(normType, {'org'})
    cb.TickLabels = cbLabels;
elseif ismember(normType, {'lin', 'prct', 'log'})
    cb.TickLabels = cbLabelsLeft;
    cb.TickDirection = 'out';
    % Right-side annotation for real count values
    cbPos = cb.Position;
    for t = 1:length(cbTicks)
        y = cbPos(2) + cbTicks(t) * cbPos(4);
        annotation('textbox', [cbPos(1)+cbPos(3)+0.07, y-0.01, 0.05, 0.03], ...
            'String', cbLabelsRight{t}, 'EdgeColor', 'none', 'FontSize', 9);
    end
    annotation('textbox', [cbPos(1)+cbPos(3)+0.07, cbPos(2)+cbPos(4)+0.015, 0.05, 0.03], ...
        'String', 'Counts', 'EdgeColor', 'none', 'FontWeight', 'bold');
end

% Video
xVals = x_min:0.5:x_max;
%xVals = x_min:0.5:5; % for tests
for xPos = xVals
    sagittalSlider.Value = xPos;
    pelvisDefectAnalysis(i).distribut.updateViews(sagittalSlider, ax1, ax2, pelvis, ...
                x_min, x_max, y_min, y_max, z_min, z_max, 'sagittal',...
                allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx);
    
    drawnow;
    frame = getframe(f);
    writeVideo(myVideo, frame);
end

close(myVideo);
close(f);  
disp(['Video saved: ', videoName]);

clear normType groupData colourMap colourData cbTicks cbLabels normName ticks lb ub cbLabelsLeft tickCounts cbLabelsRight ...
    maxVal cb cbPos y sagittalSlider ax1 ax2 rgbColour x_min x_max y_min y_max z_min z_max stepSize range f ...
    xVals xPos videoName myVideo frame thisCell useHighlight highlightIndices numL minVal

%% Display pelvis defect distribution in coronal view / X-Z-plane (for control)

% Pelvis group
normType = 'log';  % Options: 'org', 'lin', 'prct', 'log'
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
    'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 7;
groupName = groupNames{i};
groupData = pelvisDefectAnalysis(i).distribut.inside;

% Highlight intersection points
useHighlight = false; %%%
highlightIndices = [];
if useHighlight
    intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
    k = 4; %%%
    numL = numel(pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints);
    for l = 1:numL
        thisCell = pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints{l}.PointBaseIndex;
        if ~isempty(thisCell)
            highlightIndices = [highlightIndices; thisCell(:)];
        end
    end
    highlightIndices = unique(highlightIndices);
end

% Select data based on normalisation
colourMap = viridis(32768);
switch normType
    case 'org'
        colourData = groupData.pointsInsideCount;
        cbTicks = 0:5:max(groupData.pointsInsideCount);
        if cbTicks(end) < max(colourData)
            cbTicks(end+1) = max(colourData);
        end
        cbLabels = arrayfun(@(x) sprintf('%d', x), cbTicks, 'UniformOutput', false);
        normName = 'original';
        rgbColour = groupData.rgbColour;
    case 'lin'
        colourData = groupData.linNormalData;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        minVal = min(groupData.pointsInsideCount);
        maxVal = max(groupData.pointsInsideCount);
        tickCounts = minVal + cbTicks * (maxVal - minVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'lin.normalized';
        rgbColour = groupData.rgbLinNormal;
    case 'prct'
        colourData = groupData.prctNormalData;
        ticks = 0:0.25:1;
        lb = groupData.lowerBound;
        ub = groupData.upperBound;
        cbTicks = ticks;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        tickCounts = lb + ticks * (ub - lb);
        cbLabelsRight = arrayfun(@(x) sprintf('%.0f', x), tickCounts, 'UniformOutput', false);
        cbLabelsRight{1} = sprintf('0–%.0f', tickCounts(1));
        cbLabelsRight{end} = sprintf('%.0f–%d', tickCounts(end), dataCountDefect-1);
        normName = 'prct.normalized';
        rgbColour = groupData.rgbPrctNormal;
    case 'log'
        colourData = groupData.logNormalData;
        ticks = linspace(0,1,5);
        maxVal = log1p(max(groupData.pointsInsideCount));
        cbTicks = ticks;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        tickCounts = expm1(ticks * maxVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'log.normalized';
        rgbColour = groupData.rgbLogNormal;
    otherwise
        error('Invalid normalisation type.');
end

% Global limits
x_min = min(allPelvis.refPoints.gridPoints.inside(:,1));
x_max = max(allPelvis.refPoints.gridPoints.inside(:,1));
y_min = min(allPelvis.refPoints.gridPoints.inside(:,2));
y_max = max(allPelvis.refPoints.gridPoints.inside(:,2));
z_min = min(allPelvis.refPoints.gridPoints.inside(:,3));
z_max = max(allPelvis.refPoints.gridPoints.inside(:,3));

% Create figure and subplots
if useHighlight
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off');
    sgtitle(['Pelvis Group: ' groupName ' and ' intersectPctName{k} ' (' normName ')']);
else
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off');
    sgtitle(['Pelvis Group: ' groupName ' (' normName ')']);
end
ax1 = subplot(1, 2, 1); % Coronal view
ax2 = subplot(1, 2, 2); % 3D view

% Coronal view with slider (Y-Z-plane)
stepSize = 0.5;
range = y_max - y_min;
coronalSlider = uicontrol('Style', 'slider', ...
    'Min', y_min, 'Max', y_max, ...
    'Value', pelvis(1).import.processed.acentre(2), ... % pelvis(1).import.processed.acentre(2) % mean([y_min y_max])
    'SliderStep', [stepSize / range, 0.1], ...
    'Units', 'normalized', ...
    'Position', [0.2 0.01 0.2 0.05], ...
    'Callback', @(src, event) pelvisDefectAnalysis(i).distribut.updateViews(src, ax1, ax2, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'coronal',...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx));

% Saving the sliders in the figure's handles
guidata(f, struct('coronalSlider', coronalSlider, ...
    'D', allPelvis.refPoints.gridPoints.inside, 'rgbColour', rgbColour, ...
    'x_min', x_min, 'x_max', x_max, 'y_min', y_min, 'y_max', y_max, 'z_min', z_min, 'z_max', z_max));
% Update view
pelvisDefectAnalysis(i).distribut.updateViews(coronalSlider, ax1, ax2, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'coronal',...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx);

% Colorbar
hold on
scatter3(ax2, 0, 0, 0, 1, 0, 'Visible', 'off'); % Dummy points for colorbar
colormap(ax2, viridis(32768));
if ismember(normType, {'org'})
    clim(ax2, [0 max(colourData)]);
else
    clim(ax2, [0 1]);
end
cb = colorbar(ax2); 
cb.Ticks = cbTicks;
if ismember(normType, {'org'})
    cb.TickLabels = cbLabels;
elseif ismember(normType, {'lin', 'prct', 'log'})
    cb.TickLabels = cbLabelsLeft;
    cb.TickDirection = 'out';
    % Right-side annotation for real count values
    cbPos = cb.Position;
    for t = 1:length(cbTicks)
        y = cbPos(2) + cbTicks(t) * cbPos(4);
        annotation('textbox', [cbPos(1)+cbPos(3)+0.07, y-0.01, 0.05, 0.03], ...
            'String', cbLabelsRight{t}, 'EdgeColor', 'none', 'FontSize', 9);
    end
    annotation('textbox', [cbPos(1)+cbPos(3)+0.07, cbPos(2)+cbPos(4)+0.015, 0.05, 0.03], ...
        'String', 'Counts', 'EdgeColor', 'none', 'FontWeight', 'bold');
end

% Save figure (figure unvisible, but saved visible)
%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%savefig(['./Figures/PelvisDefectDistribution(', groupNames{i}, ')-', normName, 'sliceCoronal.fig'])

clear normType groupData colourMap colourData cbTicks cbLabels normName ticks lb ub cbLabelsLeft tickCounts cbLabelsRight ...
    maxVal cb cbPos y sagittalSlider coronalSlider axialSlider ax1 ax2 ax3 ax4 rgbColour x_min x_max y_min y_max z_min z_max...
    stepSize range thisCell useHighlight highlightIndices numL minVal

%% Video: Display pelvis defect distribution in coronal view / X-Z-plane (for control)

% Pelvis group
normType = 'lin';  % Options: 'org', 'lin', 'prct', 'log' %%%
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
    'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 7;
groupName = groupNames{i};
groupData = pelvisDefectAnalysis(i).distribut.inside;

% Highlight intersection points
useHighlight = true; %%%
highlightIndices = [];
if useHighlight
    intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
    k = 4; %%%
    numL = numel(pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints);
    for l = 1:numL
        thisCell = pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints{l}.PointBaseIndex;
        if ~isempty(thisCell)
            highlightIndices = [highlightIndices; thisCell(:)];
        end
    end
    highlightIndices = unique(highlightIndices);
end

% Select data based on normalisation
colourMap = viridis(32768);
switch normType
    case 'org'
        colourData = groupData.pointsInsideCount;
        cbTicks = 0:5:max(groupData.pointsInsideCount);
        if cbTicks(end) < max(colourData)
            cbTicks(end+1) = max(colourData);
        end
        cbLabels = arrayfun(@(x) sprintf('%d', x), cbTicks, 'UniformOutput', false);
        normName = 'original';
        rgbColour = groupData.rgbColour;
    case 'lin'
        colourData = groupData.linNormalData;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        minVal = min(groupData.pointsInsideCount);
        maxVal = max(groupData.pointsInsideCount);
        tickCounts = minVal + cbTicks * (maxVal - minVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'lin.normalized';
        rgbColour = groupData.rgbLinNormal;
    case 'prct'
        colourData = groupData.prctNormalData;
        lb = groupData.lowerBound;
        ub = groupData.upperBound;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        tickCounts = lb + cbTicks * (ub - lb);
        cbLabelsRight = arrayfun(@(x) sprintf('%.0f', x), tickCounts, 'UniformOutput', false);
        cbLabelsRight{1} = sprintf('0–%.0f', tickCounts(1));
        cbLabelsRight{end} = sprintf('%.0f–%d', tickCounts(end), dataCountDefect-1);
        normName = 'prct.normalized';
        rgbColour = groupData.rgbPrctNormal;
    case 'log'
        colourData = groupData.logNormalData;
        maxVal = log1p(max(groupData.pointsInsideCount));
        cbTicks = linspace(0,1,5);
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        tickCounts = expm1(cbTicks * maxVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'log.normalized';
        rgbColour = groupData.rgbLogNormal;
    otherwise
        error('Invalid normalisation type.');
end

% Global limits
x_min = min(allPelvis.refPoints.gridPoints.inside(:,1));
x_max = max(allPelvis.refPoints.gridPoints.inside(:,1));
y_min = min(allPelvis.refPoints.gridPoints.inside(:,2));
y_max = max(allPelvis.refPoints.gridPoints.inside(:,2));
z_min = min(allPelvis.refPoints.gridPoints.inside(:,3));
z_max = max(allPelvis.refPoints.gridPoints.inside(:,3));

% Video Setup
videoName = ['PelvisDistributionCoronal_' groupName '_' normName '.mp4'];
myVideo   = VideoWriter(videoName, 'MPEG-4');
myVideo.FrameRate = 10;
myVideo.Quality   = 100; 
open(myVideo);

% Figure
if useHighlight
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off', 'Position', [100 100 1920 1080]);  % Full-HD
    sgtitle(['Pelvis Group: ' groupName ' and ' intersectPctName{k} ' (' normName ')']);
else
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off', 'Position', [100 100 1920 1080]);  % Full-HD
    sgtitle(['Pelvis Group: ' groupName ' (' normName ')']);
end
ax1 = subplot(1, 2, 1); % Coronal view
ax2 = subplot(1, 2, 2); % 3D view

% Slider 
stepSize = 0.5;
range = y_max - y_min;
coronalSlider = uicontrol('Style', 'slider', ...
    'Min', y_min, ...
    'Max', y_max, ...
    'Value', mean([y_min y_max]), ... % pelvis(1).import.processed.acentre(2) % mean([y_min y_max])
    'SliderStep', [stepSize / range, 0.1], ... 
    'Units', 'normalized', ...
    'Position', [0.2 0.01 0.2 0.05], ...
    'Callback', @(src, event) pelvisDefectAnalysis(i).distribut.updateViews(src, ax1, ax2, pelvis, ...
                                         x_min, x_max, y_min, y_max, ...
                                         z_min, z_max, 'coronal',...
                                         allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, ...
                                         highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx));
% Save figure handle
guidata(f, struct('coronalSlider', coronalSlider, ...
    'D', allPelvis.refPoints.gridPoints.inside, ...
    'rgbColour', rgbColour, ...
    'x_min', x_min, 'x_max', x_max, ...
    'y_min', y_min, 'y_max', y_max, ...
    'z_min', z_min, 'z_max', z_max));
% First update (for initial view)
pelvisDefectAnalysis(i).distribut.updateViews(coronalSlider, ax1, ax2, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'coronal',...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx);

% Colorbar
hold on
scatter3(ax2, 0, 0, 0, 1, 0, 'Visible', 'off'); % Dummy points for colorbar
colormap(ax2, viridis(32768));
if ismember(normType, {'org'})
    clim(ax2, [0 max(colourData)]);
else
    clim(ax2, [0 1]);
end
cb = colorbar(ax2); 
cb.Ticks = cbTicks;
if ismember(normType, {'org'})
    cb.TickLabels = cbLabels;
elseif ismember(normType, {'lin', 'prct', 'log'})
    cb.TickLabels = cbLabelsLeft;
    cb.TickDirection = 'out';
    % Right-side annotation for real count values
    cbPos = cb.Position;
    for t = 1:length(cbTicks)
        y = cbPos(2) + cbTicks(t) * cbPos(4);
        annotation('textbox', [cbPos(1)+cbPos(3)+0.07, y-0.01, 0.05, 0.03], ...
            'String', cbLabelsRight{t}, 'EdgeColor', 'none', 'FontSize', 9);
    end
    annotation('textbox', [cbPos(1)+cbPos(3)+0.07, cbPos(2)+cbPos(4)+0.015, 0.05, 0.03], ...
        'String', 'Counts', 'EdgeColor', 'none', 'FontWeight', 'bold');
end

% Video
yVals = y_min:0.5:y_max;
%yVals = y_min:0.5:5; % for tests
for yPos = yVals
    coronalSlider.Value = yPos;
    pelvisDefectAnalysis(i).distribut.updateViews(coronalSlider, ax1, ax2, pelvis, ...
                x_min, x_max, y_min, y_max, z_min, z_max, 'coronal',...
                allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx);
    
    drawnow;
    frame = getframe(f);
    writeVideo(myVideo, frame);
end

close(myVideo);
close(f);  
disp(['Video saved: ', videoName]);

clear normType groupData colourMap colourData cbTicks cbLabels normName ticks lb ub cbLabelsLeft tickCounts cbLabelsRight ...
    maxVal cb cbPos y sagittalSlider ax1 ax2 rgbColour x_min x_max y_min y_max z_min z_max stepSize range f ...
    yVals yPos videoName myVideo frame thisCell useHighlight highlightIndices numL minVal

%% Display pelvis defect distribution in axial view / X-Y-plane (for control)

% Pelvis group
normType = 'log';  % Options: 'org', 'lin', 'prct', 'log'
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
    'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 7;
groupName = groupNames{i};
groupData = pelvisDefectAnalysis(i).distribut.inside;

% Highlight intersection points
useHighlight = false; %%%
highlightIndices = [];
if useHighlight
    intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
    k = 4; %%%
    numL = numel(pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints);
    for l = 1:numL
        thisCell = pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints{l}.PointBaseIndex;
        if ~isempty(thisCell)
            highlightIndices = [highlightIndices; thisCell(:)];
        end
    end
    highlightIndices = unique(highlightIndices);
end

% Select data based on normalisation
colourMap = viridis(32768);
switch normType
    case 'org'
        colourData = groupData.pointsInsideCount;
        cbTicks = 0:5:max(groupData.pointsInsideCount);
        if cbTicks(end) < max(colourData)
            cbTicks(end+1) = max(colourData);
        end
        cbLabels = arrayfun(@(x) sprintf('%d', x), cbTicks, 'UniformOutput', false);
        normName = 'original';
        rgbColour = groupData.rgbColour; 
    case 'lin'
        colourData = groupData.linNormalData;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        minVal = min(groupData.pointsInsideCount);
        maxVal = max(groupData.pointsInsideCount);
        tickCounts = minVal + cbTicks * (maxVal - minVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'lin.normalized';
        rgbColour = groupData.rgbLinNormal;
    case 'prct'
        colourData = groupData.prctNormalData;
        ticks = 0:0.25:1;
        lb = groupData.lowerBound;
        ub = groupData.upperBound;
        cbTicks = ticks;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        tickCounts = lb + ticks * (ub - lb);
        cbLabelsRight = arrayfun(@(x) sprintf('%.0f', x), tickCounts, 'UniformOutput', false);
        cbLabelsRight{1} = sprintf('0–%.0f', tickCounts(1));
        cbLabelsRight{end} = sprintf('%.0f–%d', tickCounts(end), dataCountDefect-1);
        normName = 'prct.normalized';
        rgbColour = groupData.rgbPrctNormal;
    case 'log'
        colourData = groupData.logNormalData;
        ticks = linspace(0,1,5);
        maxVal = log1p(max(groupData.pointsInsideCount));
        cbTicks = ticks;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        tickCounts = expm1(ticks * maxVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'log.normalized';
        rgbColour = groupData.rgbLogNormal;
    otherwise
        error('Invalid normalisation type.');
end

% Global limits
x_min = min(allPelvis.refPoints.gridPoints.inside(:,1));
x_max = max(allPelvis.refPoints.gridPoints.inside(:,1));
y_min = min(allPelvis.refPoints.gridPoints.inside(:,2));
y_max = max(allPelvis.refPoints.gridPoints.inside(:,2));
z_min = min(allPelvis.refPoints.gridPoints.inside(:,3));
z_max = max(allPelvis.refPoints.gridPoints.inside(:,3));

% Create figure and subplots
if useHighlight
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off');
    sgtitle(['Pelvis Group: ' groupName ' and ' intersectPctName{k} ' (' normName ')']);
else
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off');
    sgtitle(['Pelvis Group: ' groupName ' (' normName ')']);
end
ax1 = subplot(1, 2, 1); % Axial view
ax2 = subplot(1, 2, 2); % 3D view

% Axial view with slider (X-Y-plane)
stepSize = 0.5;
range = z_max - z_min;
axialSlider = uicontrol('Style', 'slider', ...
    'Min', z_min, ...
    'Max', z_max, ...
    'Value', pelvis(1).import.processed.acentre(3), ... % pelvis(1).import.processed.acentre(3) % mean([z_min z_max])
    'SliderStep', [stepSize / range, 0.1], ...
    'Units', 'normalized', ...
    'Position', [0.2 0.01 0.2 0.05], ...
    'Callback', @(src, event) pelvisDefectAnalysis(i).distribut.updateViews(src, ax1, ax2, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'axial',...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx));

% Saving the sliders in the figure's handles
guidata(f, struct('axialSlider', axialSlider, ...
    'D', allPelvis.refPoints.gridPoints.inside, 'rgbColour', rgbColour, ...
    'x_min', x_min, 'x_max', x_max, 'y_min', y_min, 'y_max', y_max, 'z_min', z_min, 'z_max', z_max));
% Update view
pelvisDefectAnalysis(i).distribut.updateViews(axialSlider, ax1, ax2, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'axial',...
allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx);

% Colorbar
hold on
scatter3(ax2, 0, 0, 0, 1, 0, 'Visible', 'off'); % Dummy points for colorbar
colormap(ax2, viridis(32768));
if ismember(normType, {'org'})
    clim(ax2, [0 max(colourData)]);
else
    clim(ax2, [0 1]);
end
cb = colorbar(ax2); 
cb.Ticks = cbTicks;
if ismember(normType, {'org'})
    cb.TickLabels = cbLabels;
elseif ismember(normType, {'lin', 'prct', 'log'})
    cb.TickLabels = cbLabelsLeft;
    cb.TickDirection = 'out';
    % Right-side annotation for real count values
    cbPos = cb.Position;
    for t = 1:length(cbTicks)
        y = cbPos(2) + cbTicks(t) * cbPos(4);
        annotation('textbox', [cbPos(1)+cbPos(3)+0.07, y-0.01, 0.05, 0.03], ...
            'String', cbLabelsRight{t}, 'EdgeColor', 'none', 'FontSize', 9);
    end
    annotation('textbox', [cbPos(1)+cbPos(3)+0.07, cbPos(2)+cbPos(4)+0.015, 0.05, 0.03], ...
        'String', 'Counts', 'EdgeColor', 'none', 'FontWeight', 'bold');
end

% Save figure (figure unvisible, but saved visible)
%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%savefig(['./Figures/PelvisDefectDistribution(', groupNames{i}, ')-', normName, 'sliceAxial.fig'])

clear normType groupData colourMap colourData cbTicks cbLabels normName ticks lb ub cbLabelsLeft tickCounts cbLabelsRight ...
    maxVal cb cbPos y sagittalSlider coronalSlider axialSlider ax1 ax2 ax3 ax4 rgbColour x_min x_max y_min y_max z_min z_max...
    stepSize range thisCell useHighlight highlightIndices numL minVal

%% Video: Display pelvis defect distribution in axial view / X-Y-plane (for control)

% Pelvis group
normType = 'lin';  % Options: 'org', 'lin', 'prct', 'log' %%%
groupNames = {'all', 'BTE', 'other', 'Pap2a', 'Pap2b', 'Pap2c', 'Pap3a', 'Pap3b', ...
    'Pap2aBTE', 'Pap2bBTE', 'Pap2cBTE', 'Pap3aBTE', 'Pap3bBTE'};
i = 7;
groupName = groupNames{i};
groupData = pelvisDefectAnalysis(i).distribut.inside;

% Highlight intersection points
useHighlight = true; %%%
highlightIndices = [];
if useHighlight
    intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
    k = 4; %%%
    numL = numel(pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints);
    for l = 1:numL
        thisCell = pelvisDefectAnalysis(i).intersect.inside.(intersectPctName{k}).filteredFreqPoints{l}.PointBaseIndex;
        if ~isempty(thisCell)
            highlightIndices = [highlightIndices; thisCell(:)];
        end
    end
    highlightIndices = unique(highlightIndices);
end

% Select data based on normalisation
colourMap = viridis(32768);
switch normType
    case 'org'
        colourData = groupData.pointsInsideCount;
        cbTicks = 0:5:max(groupData.pointsInsideCount);
        if cbTicks(end) < max(colourData)
            cbTicks(end+1) = max(colourData);
        end
        cbLabels = arrayfun(@(x) sprintf('%d', x), cbTicks, 'UniformOutput', false);
        normName = 'original';
        rgbColour = groupData.rgbColour; 
    case 'lin'
        colourData = groupData.linNormalData;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        minVal = min(groupData.pointsInsideCount);
        maxVal = max(groupData.pointsInsideCount);
        tickCounts = minVal + cbTicks * (maxVal - minVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'lin.normalized';
        rgbColour = groupData.rgbLinNormal;
    case 'prct'
        colourData = groupData.prctNormalData;
        lb = groupData.lowerBound;
        ub = groupData.upperBound;
        cbTicks = 0:0.25:1;
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        tickCounts = lb + cbTicks * (ub - lb);
        cbLabelsRight = arrayfun(@(x) sprintf('%.0f', x), tickCounts, 'UniformOutput', false);
        cbLabelsRight{1} = sprintf('0–%.0f', tickCounts(1));
        cbLabelsRight{end} = sprintf('%.0f–%d', tickCounts(end), dataCountDefect-1);
        normName = 'prct.normalized';
        rgbColour = groupData.rgbPrctNormal;
    case 'log'
        colourData = groupData.logNormalData;
        maxVal = log1p(max(groupData.pointsInsideCount));
        cbTicks = linspace(0,1,5);
        cbLabelsLeft = arrayfun(@(x) sprintf('%.2f', x), cbTicks, 'UniformOutput', false);
        tickCounts = expm1(cbTicks * maxVal);
        cbLabelsRight = arrayfun(@(x) sprintf('%d', round(x)), tickCounts, 'UniformOutput', false);
        normName = 'log.normalized';
        rgbColour = groupData.rgbLogNormal;
    otherwise
        error('Invalid normalisation type.');
end

% Global limits
x_min = min(allPelvis.refPoints.gridPoints.inside(:,1));
x_max = max(allPelvis.refPoints.gridPoints.inside(:,1));
y_min = min(allPelvis.refPoints.gridPoints.inside(:,2));
y_max = max(allPelvis.refPoints.gridPoints.inside(:,2));
z_min = min(allPelvis.refPoints.gridPoints.inside(:,3));
z_max = max(allPelvis.refPoints.gridPoints.inside(:,3));

% Video Setup
videoName = ['PelvisDistributionAxial_' groupName '_' normName '.mp4'];
myVideo   = VideoWriter(videoName, 'MPEG-4');
myVideo.FrameRate = 10;
myVideo.Quality   = 100; 
open(myVideo);

% Figure
if useHighlight
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off', 'Position', [100 100 1920 1080]);  % Full-HD
    sgtitle(['Pelvis Group: ' groupName ' and ' intersectPctName{k} ' (' normName ')']);
else
    f = figure('Name', [groupName ' - ' normName], 'NumberTitle', 'off', 'Position', [100 100 1920 1080]);  % Full-HD
    sgtitle(['Pelvis Group: ' groupName ' (' normName ')']);
end
ax1 = subplot(1, 2, 1); % Axial view
ax2 = subplot(1, 2, 2); % 3D view

% Slider 
stepSize = 0.5;
range = z_max - z_min;
axialSlider = uicontrol('Style', 'slider', ...
    'Min', z_min, ...
    'Max', z_max, ...
    'Value', mean([z_min z_max]), ... % pelvis(1).import.processed.acentre(2) % mean([y_min y_max])
    'SliderStep', [stepSize / range, 0.1], ... 
    'Units', 'normalized', ...
    'Position', [0.2 0.01 0.2 0.05], ...
    'Callback', @(src, event) pelvisDefectAnalysis(i).distribut.updateViews(src, ax1, ax2, pelvis, ...
                                         x_min, x_max, y_min, y_max, ...
                                         z_min, z_max, 'axial',...
                                         allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, ...
                                         highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx));
% Save figure handle
guidata(f, struct('axialSlider', axialSlider, ...
    'D', allPelvis.refPoints.gridPoints.inside, ...
    'rgbColour', rgbColour, ...
    'x_min', x_min, 'x_max', x_max, ...
    'y_min', y_min, 'y_max', y_max, ...
    'z_min', z_min, 'z_max', z_max));
% First update (for initial view)
pelvisDefectAnalysis(i).distribut.updateViews(axialSlider, ax1, ax2, pelvis, x_min, x_max, y_min, y_max, z_min, z_max, 'axial',...
    allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx);

% Colorbar
hold on
scatter3(ax2, 0, 0, 0, 1, 0, 'Visible', 'off'); % Dummy points for colorbar
colormap(ax2, viridis(32768));
if ismember(normType, {'org'})
    clim(ax2, [0 max(colourData)]);
else
    clim(ax2, [0 1]);
end
cb = colorbar(ax2); 
cb.Ticks = cbTicks;
if ismember(normType, {'org'})
    cb.TickLabels = cbLabels;
elseif ismember(normType, {'lin', 'prct', 'log'})
    cb.TickLabels = cbLabelsLeft;
    cb.TickDirection = 'out';
    % Right-side annotation for real count values
    cbPos = cb.Position;
    for t = 1:length(cbTicks)
        y = cbPos(2) + cbTicks(t) * cbPos(4);
        annotation('textbox', [cbPos(1)+cbPos(3)+0.07, y-0.01, 0.05, 0.03], ...
            'String', cbLabelsRight{t}, 'EdgeColor', 'none', 'FontSize', 9);
    end
    annotation('textbox', [cbPos(1)+cbPos(3)+0.07, cbPos(2)+cbPos(4)+0.015, 0.05, 0.03], ...
        'String', 'Counts', 'EdgeColor', 'none', 'FontWeight', 'bold');
end

% Video
zVals = z_min:0.5:z_max;
%zVals = z_min:0.5:-210; % for tests
for zPos = zVals
    axialSlider.Value = zPos;
    pelvisDefectAnalysis(i).distribut.updateViews(axialSlider, ax1, ax2, pelvis, ...
                x_min, x_max, y_min, y_max, z_min, z_max, 'axial',...
                allPelvis.refPoints.gridPoints.inside, rgbColour, useHighlight, highlightIndices, allPelvis.refPoints.gridPoints.insideMaskIdx);
    
    drawnow;
    frame = getframe(f);
    writeVideo(myVideo, frame);
end

close(myVideo);
close(f);  
disp(['Video saved: ', videoName]);

clear normType groupData colourMap colourData cbTicks cbLabels normName ticks lb ub cbLabelsLeft tickCounts cbLabelsRight ...
    maxVal cb cbPos y sagittalSlider ax1 ax2 rgbColour x_min x_max y_min y_max z_min z_max stepSize range f ...
    zVals zPos videoName myVideo frame thisCell useHighlight highlightIndices numL minVal

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Class APPROXIMATE %%%%%%%%%%

% Pelvis defect categories (all or selected) are base for intersection
baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
numBaseIntersect = length(baseIntersect);

% Intersection types
intersectType = {'alpha','refinedAlpha','inside','grid'};
numIntersectTypes = length(intersectType);

for j = 1:numIntersectTypes
    for i = 1:numBaseIntersect
        % Base of intersections
        pelvisDefectAnalysis(i).approx.(intersectType{j}).base = baseIntersect{i};
    end
    % Pelvis IDs for the categories
    pelvisDefectAnalysis(1).approx.(intersectType{j}).IDs = (2:dataCountDefect)'; % all
    pelvisDefectAnalysis(2).approx.(intersectType{j}).IDs = [2;3;(5:9)';11;12;14;15;(17:21)';23;24;(26:37)';39;40;(42:47)']; % BTE
    pelvisDefectAnalysis(3).approx.(intersectType{j}).IDs = [4;10;13;16;22;25;38;41]; % Other
    pelvisDefectAnalysis(4).approx.(intersectType{j}).IDs = [3;4;5]; % 2A
    pelvisDefectAnalysis(5).approx.(intersectType{j}).IDs = [2;13;14;17;19;25;38;42;45]; % 2B
    pelvisDefectAnalysis(6).approx.(intersectType{j}).IDs = [11;40;44]; % 2C
    pelvisDefectAnalysis(7).approx.(intersectType{j}).IDs = [6;12;20;21;22;26;28;29;30;31;32;34;35;36;39;47;48]; % 3A
    pelvisDefectAnalysis(8).approx.(intersectType{j}).IDs = [7;8;9;10;15;16;18;23;24;27;33;37;41;43;46]; % 3B
    pelvisDefectAnalysis(9).approx.(intersectType{j}).IDs = [3;5]; % 2A & BTE
    pelvisDefectAnalysis(10).approx.(intersectType{j}).IDs = [2;14;17;19;42;45]; % 2B & BTE
    pelvisDefectAnalysis(11).approx.(intersectType{j}).IDs = [11;40;44]; % 2C & BTE
    pelvisDefectAnalysis(12).approx.(intersectType{j}).IDs = [6;12;20;21;26;28;29;30;31;32;34;35;36;39;47;48]; % 3A & BTE
    pelvisDefectAnalysis(13).approx.(intersectType{j}).IDs = [7;8;9;15;18;23;24;27;33;37;43;46]; % 3B & BTE
end

%% Approximation: Preparation of the stl

% Pelvis defect categories (all or selected) are base for intersection
baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 

% Intersection types
intersectType = {'alpha','refinedAlpha','inside','grid'};
numIntersectTypes = length(intersectType);

% Save stl files (intersection alpha shape)
%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    for i = 1:numBaseIntersect 
        baseName = baseIntersect{i};
        for k = 1:numIntersectPct
            pctName = intersectPctName{k};
            if pelvisDefectAnalysis(i).intersect.(typeName).(pctName).count > 1
                % Extract relevant data (mesh)
                intersectFaces = pelvisDefectAnalysis(i).intersect.(typeName).(pctName).faces;
                intersectVertices = pelvisDefectAnalysis(i).intersect.(typeName).(pctName).allVerticesPoints;
                % Input fpr clumped spheres: stl files
                % create stl files from intersections
                intersectTri = triangulation(intersectFaces,intersectVertices);
                fileName = ['./GeometriesIntersection/intersect' typeName '_' baseName '_' pctName '.stl']; % base,typeName,pctName
                stlwrite(intersectTri,fileName); 
            end
        end
    end
%end
clear intersectFaces intersectVertices intersectTri fileName

%% Approximation: Clumped sphere (sihaeri/DEM-ClumpedSphere) for alphaShape from Class Intersection 
% with parameter variation

% Clumped spheres
scaleFactor = 1;
smoothFactorStart = 3;
smoothFactorMax = 5;
nNodesStart = 1;
nNodesMax = 10;
maxSpheres = 10;
%for j = 1:numIntersectTypes
typeName = intersectType{j};
for i = 1:numBaseIntersect
    baseName = baseIntersect{i};
    for k = 1:numIntersectPct
        pctName = intersectPctName{k};
        
        t2 = tic; 
        pelvisDefectAnalysis(i).approx.(typeName).(pctName) = struct(); 
        if pelvisDefectAnalysis(i).intersect.(typeName).(pctName).count > 1
            fileName = ['./GeometriesIntersection/intersect' typeName '_' baseName '_' pctName '.stl'];
            validModels = []; % Storage for all sphere models ≤10 
            fallbackModel = []; % Storage for smallest model >10
            fallbackNodes = NaN;
            
            for smoothFactor = smoothFactorStart:smoothFactorMax
                for nNodes = nNodesStart:nNodesMax               
                    try
                        t1 = tic; % time to generate sphere model
                        approx = pelvisDefectAnalysis(i).approx.clumpedSphere(... % temp
                            typeName, pctName, fileName, nNodes, scaleFactor, smoothFactor);
                        tocVal = toc(t1); 
                        numSpheres = size(approx.(typeName).(pctName).spheres, 2);
                        % Warning
                        if approx.(typeName).(pctName).volSpheres == 0
                            warning(['Zero volume with ' fileName ', nNodes=' num2str(nNodes) ', smooth=' num2str(smoothFactor)]);
                        end

                        % Case 1: sphere count ≤ 10 -> save
                        if numSpheres <= maxSpheres && approx.(typeName).(pctName).volSpheres > 0
                            validModels(end+1).smooth = smoothFactor;
                            validModels(end).nNodes = nNodes;
                            validModels(end).scaleFactor = scaleFactor;
                            validModels(end).spheres = approx.(typeName).(pctName).spheres;
                            validModels(end).volume = approx.(typeName).(pctName).volSpheres;
                            validModels(end).memoryBytes = whos('approx').bytes;
                            validModels(end).time = tocVal; % time only for valid models
                        else
                            % Case 2: Only save if no fallback has been set yet
                            if isempty(fallbackModel) && approx.(typeName).(pctName).volSpheres > 0
                                fallbackModel.smooth = smoothFactor;
                                fallbackModel.nNodes = nNodes;
                                fallbackModel.scaleFactor = scaleFactor;
                                fallbackModel.spheres = approx.(typeName).(pctName).spheres;
                                fallbackModel.volume = approx.(typeName).(pctName).volSpheres;
                                fallbackModel.memoryBytes = whos('approx').bytes;
                            end
                        end
                        disp(['Sphere model with ' fileName ', nNodes=' num2str(nNodes) ', smooth=' num2str(smoothFactor)]);

                    catch ME
                        % Ignore error, but display
                        disp(['Error with ' fileName ', nNodes=' num2str(nNodes) ', smooth=' num2str(smoothFactor)]);
                        disp(ME.message);
                        continue
                    end
                end

                % If at least one valid model (≤10 spheres) is found -> cancel
                if ~isempty(validModels)
                    pelvisDefectAnalysis(i).approx.(typeName).(pctName).models = validModels;
                    break;  % No further smooth values necessary
                end
            end
            % If no model with ≤10 balls found -> save fallback
            if isempty(validModels) && ~isempty(fallbackModel)
                pelvisDefectAnalysis(i).approx.(typeName).(pctName).fallbackModel = fallbackModel;
                disp(['fallbackModel for i=' num2str(i) ', typeName=' typeName ', pctName=' pctName]);
            end
        end
        tocPct(k,i) = toc(t2); 
    end
end
%end
clear ME approxi scaleFactor smoothFactorStart smoothFactorMax smoothFactor nNodes nNodesStart nNodesMax maxSpheres...
fileName validModels fallbackModel fallbackNodes approx numSpheres m t1 t2 tocVal tocPct %tocPct


% Store all possible models to different intersections and pct amount
allIntersect.(typeName).approxModels(numBaseIntersect,1) = struct('models', [], 'fallbackModels', []);
%for j = 1:numIntersectTypes
typeName = intersectType{j};
for i = 1:numBaseIntersect
    baseName = baseIntersect{i};
    for k = 1:numIntersectPct
        pctName = intersectPctName{k};
        % Save models
        if pelvisDefectAnalysis(i).intersect.(typeName).(pctName).count > 1 && ...
                isfield(pelvisDefectAnalysis(i).approx.(typeName), pctName) && ...
                isfield(pelvisDefectAnalysis(i).approx.(typeName).(pctName), 'models')

            models = pelvisDefectAnalysis(i).approx.(typeName).(pctName).models;
            for m = 1:length(models)
                modelEntry = models(m);
                modelEntry.sphereCount = size(pelvisDefectAnalysis(i).approx.(typeName).(pctName).models(m).spheres,2);
                modelEntry.pctName = pctName;
                if isempty(allIntersect.(typeName).approxModels(i).models)
                    allIntersect.(typeName).approxModels(i).models = modelEntry;
                else
                    allIntersect.(typeName).approxModels(i).models(end+1) = modelEntry;
                end
            end
        end

        % Save fallback models
        if pelvisDefectAnalysis(i).intersect.(typeName).(pctName).count > 1 && ...
                isfield(pelvisDefectAnalysis(i).approx.(typeName), pctName) && ...
                isfield(pelvisDefectAnalysis(i).approx.(typeName).(pctName), 'fallbackModel')

            fallbackEntry = pelvisDefectAnalysis(i).approx.(typeName).(pctName).fallbackModel;
            fallbackEntry.sphereCount = size(pelvisDefectAnalysis(i).approx.(typeName).(pctName).fallbackModel.spheres,2);
            fallbackEntry.pctName = pctName;
            if isempty(allIntersect.(typeName).approxModels(i).fallbackModels)
                allIntersect.(typeName).approxModels(i).fallbackModels = fallbackEntry;
            else
                allIntersect.(typeName).approxModels(i).fallbackModels(end+1) = fallbackEntry;
            end
        end
    end
end
%end
clear models m modelEntry fallbackEntry

%% Display generated approximated sphere model (for control)

% Inputs
% baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
i = 7; % Base intersect index %%%
baseName = baseIntersect{i}; % Base intersect name
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 3; %%%
typeName = intersectType{j};

% intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
%for k = 1:numIntersectPct
    k = 4; %%%
    pctName = intersectPctName{k};

    if isfield(pelvisDefectAnalysis(i).approx.(typeName), pctName) && ...
       isfield(pelvisDefectAnalysis(i).approx.(typeName).(pctName), 'models')

        models = pelvisDefectAnalysis(i).approx.(typeName).(pctName).models;
        numModels = numel(models);
    
         for modelNr = 1:numModels
            model = models(modelNr);
            spheres = model.spheres;
    
            %figure('Visible','off')
            figure;
            hold on;
            legendEntries = [];
            legendLabels = {};
            
            % Sphere model
            for s = 1:size(spheres, 2) 
                [x, y, z] = sphere;  
                r = spheres(4, s);       % Radius
                cx = spheres(1, s);      % centre x
                cy = spheres(2, s);      % centre y
                cz = spheres(3, s);      % centre z
                hSphere = surf(r * x + cx, r * y + cy, r * z + cz, ...
                    'FaceAlpha', 0.15, 'EdgeAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0.2 0.6 0.9]);
                if s == 1  % for legend
                    legendEntries(end+1) = hSphere;
                    legendLabels{end+1} = 'Approximated Spheres';
                end
            end
            
            % Alpha shape of intersection
            hPatch = patch('Faces', pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).faces, ...
                  'Vertices', pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).allVerticesPoints, ...
                  'FaceColor', [0.8 0.8 0.8], ...
                  'FaceAlpha', 1, ...
                  'EdgeColor', [0.7 0.7 0.7], ...
                  'EdgeAlpha', 0.25);
            
            % Format
            legendEntries(end+1) = hPatch;
            legendLabels{end+1} = 'Intersection Shape';
            legend(legendEntries, legendLabels, 'Location', 'best');
            daspect([1 1 1]);
            view(3);
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title({['Clumped Sphere Model #' num2str(modelNr)], ...
                   ['Type: ' typeName ', Base: ' baseName ', Pct: ' pctName], ...
                   ['nNodes = ' num2str(model.nNodes) ...
                   ', smooth = ' num2str(model.smooth)...
                   ', sphere count = ' num2str(size(spheres,2))]});  
            grid off;
            hold off;
            
            % Save figure (figure unvisible, but saved visible)
            %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
            %savefig(['./Figures/PelvisDefectApprox_' baseName '_' typeName '_' pctName '_model' num2str(modelNr) '.fig']);

        end
    end
%end

clear modelNr spheres s models numModels model x y z cx cy cz r hSphere hPatch legendEntries legendLabels

%% Volume approximation: intersection of sphere model and reference pelvis

% Parameter
baseIntersect     = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'};
intersectType     = {'alpha','refinedAlpha','inside','grid'};
intersectPctName  = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
reference = pelvis(1).import.processed;
voxelSpacing = 0.5;

%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    for i = 1:numBaseIntersect
        pelvisDefectAnalysis(i).approx = pelvisDefectAnalysis(i).approx.volSphereIntersect(i, reference, ...
            typeName, intersectPctName, voxelSpacing); % with function inpolyhedron (long runtime)
    end
%end
clear reference voxelSpacing

% Store all possible volume models to different intersections and pct amount
%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    allIntersect.(typeName).approxVolModels(numBaseIntersect,1) = struct('modelVol', []);
    for i = 1:numBaseIntersect
        baseName = baseIntersect{i};
        for k = 1:numIntersectPct
            pctName = intersectPctName{k};
            % Save modelVol 
            if pelvisDefectAnalysis(i).intersect.(typeName).(pctName).count > 1 && ...
               isfield(pelvisDefectAnalysis(i).approx.(typeName), pctName) && ...
               isfield(pelvisDefectAnalysis(i).approx.(typeName).(pctName), 'modelVol')
    
                vols = pelvisDefectAnalysis(i).approx.(typeName).(pctName).modelVol;
                for mv = 1:numel(vols)
                    volEntry = vols(mv);
                    volEntry.pctName = pctName;
                    if isempty(allIntersect.(typeName).approxVolModels(i).modelVol)
                        allIntersect.(typeName).approxVolModels(i).modelVol = volEntry;
                    else
                        allIntersect.(typeName).approxVolModels(i).modelVol(end+1) = volEntry;
                    end
                end
            end
        end
    end
%end
clear vols m mv volEntry 

%% Volume approximation with scaling: intersection of sphere model and reference pelvis 

% Parameter
tolerancePct = 0.01;     % Tolerance
maxIter = 10;            
initialScale = 0.75;     % Start value
scalingExp = 0.3;        % Damping, optional
refVolume = [70001.44662]; % ref all  %%% refVolume(1:numBaseIntersect) %%%
refVolume(7) = [63927.09906]; % 3A
voxelSpacing = 0.5;
voxelVol = voxelSpacing^3;

%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    for i = 1:numBaseIntersect
        targetVol = refVolume(i);
        gridPointsGlobal = pelvisDefectAnalysis(i).approx.(typeName).gridSpheres;
        for k = 1:numIntersectPct
            pctName = intersectPctName{k};
            pelvisDefectAnalysis(i).approx = pelvisDefectAnalysis(i).approx.volSphereScale(i, gridPointsGlobal, ...
                typeName, pctName, voxelVol, targetVol, tolerancePct, maxIter, initialScale, scalingExp);
        end
    end
%end
clear tolerancePct maxIter initialScale scalingExp refVolume voxelSpacing voxelVol gridPointsGlobal targetVol

% Store all possible scaled volume models to different intersections and pct amount
% for j = 1:numIntersectTypes
    typeName = intersectType{j};
    allIntersect.(typeName).approxScaleModels(numBaseIntersect,1) = struct('modelScaleVol', []);
    for i = 1:numBaseIntersect
        baseName = baseIntersect{i};
        for k = 1:numIntersectPct
            pctName = intersectPctName{k};
            % Save scaled modelVol 
            if pelvisDefectAnalysis(i).intersect.(typeName).(pctName).count > 1 && ...
               isfield(pelvisDefectAnalysis(i).approx.(typeName), pctName) && ...
               isfield(pelvisDefectAnalysis(i).approx.(typeName).(pctName), 'modelScale')
    
                scales = pelvisDefectAnalysis(i).approx.(typeName).(pctName).modelScale;
                for mv = 1:numel(scales)
                    scaleEntry = scales(mv);
                    scaleEntry.pctName = pctName;
                    if isempty(allIntersect.(typeName).approxScaleModels(i).modelScaleVol)
                        allIntersect.(typeName).approxScaleModels(i).modelScaleVol = scaleEntry;
                    else
                        allIntersect.(typeName).approxScaleModels(i).modelScaleVol(end+1) = scaleEntry;
                    end
                end
            end
        end
    end
% end
clear scales mv scaleEntry

%% Display grid points inside sphere model (for control)

% Inputs
% baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
i = 1; % Base intersect index %%%
baseName = baseIntersect{i}; % Base intersect name
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 1; %%%
typeName = intersectType{j};

% intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
%for k = 1:numIntersectPct
    k = 11; %%%
    pctName = intersectPctName{k};

    if isfield(pelvisDefectAnalysis(i).approx.(typeName), pctName) && ...
       isfield(pelvisDefectAnalysis(i).approx.(typeName).(pctName), 'modelVol')

        models = pelvisDefectAnalysis(i).approx.(typeName).(pctName).models;
        numModels = numel(models);

         for modelNr = 1:numModels
            model = models(modelNr);
            spheres = model.spheres;
    
            % Figure points inside spheres
            figure('Name','Voxel inside Sphere Model');
            hold on;
            legendEntries = [];
            legendLabels = {};
            
            % Sphere model
            for s = 1:size(spheres,2)
                [x, y, z] = sphere(20);
                x = x * spheres(4,s) + spheres(1,s);
                y = y * spheres(4,s) + spheres(2,s);
                z = z * spheres(4,s) + spheres(3,s);
                hSphere = surf(x, y, z, ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.1, ...
                    'FaceColor', TUMcolors.blue300);
                if s == 1
                    legendEntries(end+1) = hSphere;
                    legendLabels{end+1} = 'Approximated Spheres';
                end
            end
            % Points inside sphere mode
            pointsInSpheres = pelvisDefectAnalysis(i).approx.(typeName).(pctName).modelVol(modelNr).pointsInSpheres;
            hScatter = scatter3(pointsInSpheres(:,1), ...
                     pointsInSpheres(:,2), ...
                     pointsInSpheres(:,3), ...
                     1, 'r', 'filled');
            
            % Format
            legendEntries(end+1) = hScatter;
            legendLabels{end+1} = 'Points inside Spheres';
            legend(legendEntries, legendLabels, 'Location', 'best');
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title({['Voxel/Points inside sphere model: Sphere Model #' num2str(modelNr)], ...
                   ['Type: ' typeName ', Base: ' baseName ', Pct: ' pctName]});  
            axis equal;
            view(3);
            grid off;
            camlight; lighting gouraud;
            hold off;
            
            % Save figure (figure unvisible, but saved visible)
            %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
            %savefig(['./Figures/PelvisDefectApproxVoxel_' baseName '_' typeName '_' pctName '_model' num2str(modelNr) '.fig']);

         end
    end
%end

clear modelNr spheres s models numModels model x y z legendEntries legendLabels hSphere hScatter pointsInSpheres

%% Display grid points (volume) inside reference pelvis (for control)

% Inputs
% baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
i = 1; % Base intersect index %%%
baseName = baseIntersect{i}; % Base intersect name
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 1; %%%
typeName = intersectType{j};

% intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
%for k = 1:numIntersectPct
    k = 3; %%%
    pctName = intersectPctName{k};

    if isfield(pelvisDefectAnalysis(i).approx.(typeName), pctName) && ...
       isfield(pelvisDefectAnalysis(i).approx.(typeName).(pctName), 'modelVol')

        models = pelvisDefectAnalysis(i).approx.(typeName).(pctName).models;
        numModels = numel(models);
         for modelNr = 1:numModels
            model = models(modelNr);
            spheres = model.spheres;

            figure('Name','Intersection: Sphere model - reference pelvis');
            hold on;
            
            % Ref pelvis
            hPelvis = patch('Faces', pelvis(1).import.processed.faces, ...
                  'Vertices', pelvis(1).import.processed.vertices, ...
                  'FaceColor', [0.8 0.8 0.8], ...
                  'FaceAlpha', 0.2, ...
                  'EdgeColor', 'none');
            
            % Sphere model
            for s = 1:size(spheres,2)
                [x, y, z] = sphere(20);
                x = x * spheres(4,s) + spheres(1,s);
                y = y * spheres(4,s) + spheres(2,s);
                z = z * spheres(4,s) + spheres(3,s);
                hSphere = surf(x, y, z, ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.05, ...
                    'FaceColor', TUMcolors.blue300);
                if s == 1        % nur ein Eintrag für die Legende
                    hSphereLegend = hSphere;
                end
            end
            
            % Points inside (intersection) 
            intersectionPoints = pelvisDefectAnalysis(i).approx.(typeName).gridSpheres(...
                pelvisDefectAnalysis(i).approx.(typeName).(pctName).modelVol(modelNr).intersectIdxGlobal, :);
            hPts  = scatter3(intersectionPoints(:,1), ...
                     intersectionPoints(:,2), ...
                     intersectionPoints(:,3), ...
                     10, 'r', 'filled');
            % Alpha shape of intersection
            hAlpha = patch('Faces', pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).faces, ...
                  'Vertices', pelvisDefectAnalysis(i).intersect.(intersectType{j}).(intersectPctName{k}).allVerticesPoints, ...
                  'FaceColor', [0.8 0.8 0.8], ...
                  'FaceAlpha', 1, ...
                  'EdgeColor', [0.7 0.7 0.7], ...
                  'EdgeAlpha', 0.25);
            
            % Format
            legend([hPelvis, hSphereLegend, hPts, hAlpha], ...
           {'Reference pelvis', 'Sphere model', 'Intersection voxels', 'Alpha shape'}, 'Location', 'best');
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title({['Voxel/Points inside sphere model: Sphere Model #' num2str(modelNr)], ...
                   ['Type: ' typeName ', Base: ' baseName ', Pct: ' pctName]});  
            axis equal;
            view(3);
            grid off;
            camlight; lighting gouraud;
            hold off;

            % Save figure (figure unvisible, but saved visible)
            %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
            %savefig(['./Figures/PelvisDefectApproxIntersectVol_' baseName '_' typeName '_' pctName '_model' num2str(modelNr) '.fig']);
         end
    end
%end

clear modelNr spheres s models numModels model x y z intersectionPoints hPelvis hSphere hSphereLegend hPts hAlpha

%% Display grid points (scaled volume) inside reference pelvis (for control)

% Inputs
% baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
i = 1; % Base intersect index %%%
baseName = baseIntersect{i}; % Base intersect name
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 3; %%%
typeName = intersectType{j};

% intersectPctName = {'pairs','pct5','pct10','pct25','pct33','pct50','pct66','pct75','pct90','pct95','pct100'};
%for k = 1:numIntersectPct
    k = 9; %%%
    pctName = intersectPctName{k};

    if isfield(pelvisDefectAnalysis(i).approx.(typeName), pctName) && ...
       isfield(pelvisDefectAnalysis(i).approx.(typeName).(pctName), 'modelScale')
    
        models = pelvisDefectAnalysis(i).approx.(typeName).(pctName).models;
        modelScale = pelvisDefectAnalysis(i).approx.(typeName).(pctName).modelScale;
        numModels = numel(models);
        for modelNr = 1:numModels
            model = models(modelNr); 
            modelScaled = modelScale(modelNr); 
            origSpheres = model.spheres;
            scaledSpheres = modelScaled.scaledSpheres;
            intersectionPointsScaled = pelvisDefectAnalysis(i).approx.(typeName).gridSpheres(...
                pelvisDefectAnalysis(i).approx.(typeName).(pctName).modelScale(modelNr).intersectIdxGlobal, :);
    
            figure('Name', 'Scaled vs. Original Sphere Model');
            hold on;
    
            % Ref pelvis
            hPelvis = patch('Faces', pelvis(1).import.processed.faces, ...
                  'Vertices', pelvis(1).import.processed.vertices, ...
                  'FaceColor', [0.8 0.8 0.8], ...
                  'FaceAlpha', 0.2, ...
                  'EdgeColor', 'none');
    
            % Original sphere model (background) 
            for s = 1:size(origSpheres, 2)
                [x, y, z] = sphere(20);
                x = x * origSpheres(4,s) + origSpheres(1,s);
                y = y * origSpheres(4,s) + origSpheres(2,s);
                z = z * origSpheres(4,s) + origSpheres(3,s);
                hO = surf(x, y, z, ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.05, ...
                    'FaceColor', [0.7 0.7 0.7]);
                if s==1
                    hOrig = hO;
                end
            end
    
            % Scaled Sphere Model
            for s = 1:size(scaledSpheres, 2)
                [x, y, z] = sphere(20);
                x = x * scaledSpheres(4,s) + scaledSpheres(1,s);
                y = y * scaledSpheres(4,s) + scaledSpheres(2,s);
                z = z * scaledSpheres(4,s) + scaledSpheres(3,s);
                hS = surf(x, y, z, ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.25, ...
                    'FaceColor', TUMcolors.blue300);
                if s==1
                    hScaled = hS;
                end
            end

            % Intersection volume (points)
            hPts = scatter3(intersectionPointsScaled(:,1), ...
                     intersectionPointsScaled(:,2), ...
                     intersectionPointsScaled(:,3), ...
                     10, 'r', 'filled');

            % Format
            legend([hPelvis,hOrig,hScaled,hPts], {'Reference pelvis','Original spheres','Scaled spheres','Intersection voxels'}, ...
           'Location','best');
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title({['Skalierung: Modell #' num2str(modelNr)], ...
                   ['Typ: ' typeName ', Base: ' baseName ', Pct: ' pctName], ...
                   ['Scale Factor: ' num2str(modelScaled.scaleFactorToRef, '%.3f')]});
            axis equal;
            view(3);
            grid off;
            camlight; lighting gouraud;
            hold off;
    
            % Save figure (figure unvisible, but saved visible)
            %set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
            %savefig(['./Figures/PelvisDefectScaledApprox_' baseName '_' typeName '_' pctName '_model' num2str(modelNr) '.fig']);
        end
    end
%end

clear modelNr origSpheres scaledSpheres intersectionPointsScaled model models modelScale modelScaled numModels x y z s ....
    hPelvis hO hOrig hS hScaled hPts

%% Comparison sphere models - best model for base (point frequency, voxel-based) 
% Best fit sphere model with point frequency (with scaling)

% Sphere models frequency (per base)
%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    for i = 1:numBaseIntersect
        tic
        pelvisDefectAnalysis(i).approx = pelvisDefectAnalysis(i).approx.compSpheres(i, typeName, intersectPctName);

        % Time and memory
        pelvisDefectAnalysis(i).approx.(typeName).modelFitTime = toc;
        tmp = pelvisDefectAnalysis(i).approx.(typeName).modelFit;
        memInfo = whos('tmp');
        pelvisDefectAnalysis(i).approx.(typeName).modelFitMemory = memInfo.bytes;
    end
%end
clear tmp memInfo

% Store all frequency-based model scores per base and type
% for j = 1:numIntersectTypes
    typeName = intersectType{j};
    allIntersect.(typeName).approxFreqModels(numBaseIntersect,1) = struct('modelFreq', []);
    for i = 1:numBaseIntersect
        baseName = baseIntersect{i};
    
        if isfield(pelvisDefectAnalysis(i).approx.(typeName), 'modelFit')
            models = pelvisDefectAnalysis(i).approx.(typeName).modelFit;
            for m = 1:numel(models)
                freqEntry = models(m);  
                if isempty(allIntersect.(typeName).approxFreqModels(i).modelFreq)
                    allIntersect.(typeName).approxFreqModels(i).modelFreq = freqEntry;
                else
                    allIntersect.(typeName).approxFreqModels(i).modelFreq(end+1) = freqEntry;
                end
            end
        end
    end
% end
clear models m freqEntry


% Best fit sphere model
%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    for i = 1:numBaseIntersect
        modelList = pelvisDefectAnalysis(i).approx.(typeName).modelFit;
        [~, bestIdxSum]    = max([modelList.scoreSum]);
        [~, bestIdxMean]   = max([modelList.scoreMean]);
        [~, bestIdxGlobal] = max([modelList.scoreGlobal]);
    
        pelvisDefectAnalysis(i).approx.(typeName).bestModel.scoreSumModelNr = modelList(bestIdxSum).modelNr;
        pelvisDefectAnalysis(i).approx.(typeName).bestModel.scoreSumPct = modelList(bestIdxSum).pctName;
        pelvisDefectAnalysis(i).approx.(typeName).bestModel.scoreMeanModelNr = modelList(bestIdxMean).modelNr;
        pelvisDefectAnalysis(i).approx.(typeName).bestModel.scoreMeanPct = modelList(bestIdxMean).pctName;
        pelvisDefectAnalysis(i).approx.(typeName).bestModel.scoreGlobalModelNr = modelList(bestIdxGlobal).modelNr;
        pelvisDefectAnalysis(i).approx.(typeName).bestModel.scoreGlobalPct = modelList(bestIdxGlobal).pctName;

    end
%end
clear modelList bestIdxSum bestIdxMean bestIdxGlobal tmp memInfo

% Best fit model 
bestModelType = 'Mean';  % 'Sum', 'Mean', 'Global'
%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    allIntersect.(typeName).approxBestModel(numBaseIntersect,1) = struct('bestModel', []);
    for i = 1:numBaseIntersect
        bestModel = pelvisDefectAnalysis(i).approx.(typeName).bestModel;
        bestPct = bestModel.(['score' bestModelType 'Pct']);
        bestNr = bestModel.(['score' bestModelType 'ModelNr']);
        model = pelvisDefectAnalysis(i).approx.(typeName).(bestPct);
        bestMdl = struct();
        bestMdl.pctName = bestPct;
        bestMdl.modelNr = bestNr;
        bestMdl.modelType= bestModelType;
        bestMdl.smooth = model.models(bestNr).smooth;
        bestMdl.nNodes = model.models(bestNr).nNodes;
        bestMdl.idxGlobal = model.modelScale(bestNr).intersectIdxGlobal;
        bestMdl.scaledSpheres = model.modelScale(bestNr).scaledSpheres;
        bestMdl.nSpheres = size(model.modelScale(bestNr).scaledSpheres,2);
        bestMdl.volumeScaled = model.modelScale(bestNr).intersectVolumeScaled;
        bestMdl.scaleFactorToRef = model.modelScale(bestNr).scaleFactorToRef;
        % Save
        allIntersect.(typeName).approxBestModel(i).bestModel = bestMdl;
    end
% end
clear bestPct bestNr bestModel model bestMdl

%% Comparison sphere model to acetabular defect
% Comparison intersect types

% Create mesh for sphere model
spacing =  0.5;
bestModelType = 'Mean'; % 'Sum', 'Mean', 'Global'
%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    for i = 1:numBaseIntersect
        pelvisDefectAnalysis(i).approx.(typeName).bestModelDist.bestModelType = bestModelType;
        bestModelNr = pelvisDefectAnalysis(i).approx.(typeName).bestModel.(['score' bestModelType 'ModelNr']);
        bestModelPct = pelvisDefectAnalysis(i).approx.(typeName).bestModel.(['score' bestModelType 'Pct']);
        bestModelScaledSpheres = pelvisDefectAnalysis(i).approx.(typeName).(bestModelPct).modelScale(bestModelNr).scaledSpheres;
        pelvisDefectAnalysis(i).approx = pelvisDefectAnalysis(i).approx.sphereMesh(i, typeName, bestModelType, bestModelScaledSpheres, spacing);
    end
%end
clear spacing bestModelScaledSpheres bestModelNr bestModelPct bestModelType

% Vertex-to-nearest neighbour: sphere model - defect mesh %%%
weightKabsch = {'w1', 'w2', 'w3', 'w4', 'w5'};
w = 5; % Select weighting
bestModelType = 'Mean'; % 'Sum', 'Mean', 'Global'
%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    for i = 1:numBaseIntersect
        relevantIDs = pelvisDefectAnalysis(i).approx.(typeName).IDs;  % Defect indices
        for d = relevantIDs(:)'
            vertices = pelvisDefect(d).transform.trafo.(weightKabsch{w}).vertices;
            faces = pelvisDefect(d).transform.trafo.(weightKabsch{w}).faces;
            pelvisDefectAnalysis(i).approx = pelvisDefectAnalysis(i).approx.sphereDistance(i, d, ...
                vertices, faces, typeName, bestModelType);
        end
    end
%end
% Mean Distance
bestModelType = 'Mean'; % 'Sum', 'Mean', 'Global'
%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    for i = 1:numBaseIntersect
        relevantIDs = pelvisDefectAnalysis(i).approx.(typeName).IDs;
        distVals = pelvisDefectAnalysis(i).approx.(typeName).bestModelDist.(['nearDist' bestModelType 'Mean']);
        distFiltered = distVals(relevantIDs);
        pelvisDefectAnalysis(i).approx.(typeName).bestModelDist.(['meanDistance' bestModelType]) = mean(distFiltered); (['meanDistance' bestModelType])
    end
%end
% Store in allIntresect
bestModelType = 'Mean'; % 'Sum', 'Mean', 'Global'
%for j = 1:numIntersectTypes
    typeName = intersectType{j};
    allIntersect.(typeName).approxV2NNmodels(numBaseIntersect,1) = ...
    struct(['nearDist' bestModelType],[],['nearDist' bestModelType 'Mean'],[],['meanDistance' bestModelType],[]); 
    for i = 1:numBaseIntersect  
        allIntersect.(typeName).approxV2NNmodels(i).(['nearDist' bestModelType]) = ...
            pelvisDefectAnalysis(i).approx.(typeName).bestModelDist.(['nearDist' bestModelType]);   % per defect and per vertex
        allIntersect.(typeName).approxV2NNmodels(i).(['nearDist' bestModelType 'Mean']) = ...
            pelvisDefectAnalysis(i).approx.(typeName).bestModelDist.(['nearDist' bestModelType 'Mean']);  % mean per defect
        allIntersect.(typeName).approxV2NNmodels(i).(['meanDistance' bestModelType]) = ...
            pelvisDefectAnalysis(i).approx.(typeName).bestModelDist.(['meanDistance' bestModelType]); % mean relevant defects
    end
%end

clear vertices faces bestModelType  relevantIDs d distVals distFiltered

%% Display best fit sphere model (for control)

% Inputs
% baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
i = 1; % Base intersect index %%%
baseName = baseIntersect{i}; % Base intersect name
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 3; %%%
typeName = intersectType{j};     
bestModelType = 'Mean'; % 'Sum', 'Mean', 'Global' %%%

% Laod model
pctName = pelvisDefectAnalysis(i).approx.(typeName).bestModel.scoreMeanPct;
modelNr = pelvisDefectAnalysis(i).approx.(typeName).bestModel.scoreMeanModelNr;
models = pelvisDefectAnalysis(i).approx.(typeName).(pctName).models;
modelScale = pelvisDefectAnalysis(i).approx.(typeName).(pctName).modelScale;
origSpheres = models(modelNr).spheres;
scaledSpheres = modelScale(modelNr).scaledSpheres;
idxModel = modelScale(modelNr).intersectIdxGlobal;
ptsModel = pelvisDefectAnalysis(i).approx.(typeName).gridSpheres(idxModel,:);

figure('Name','Best-Fit scaled sphere model');
hold on;

% Reference
hPelvis = patch('Faces', pelvis(1).import.processed.faces, ...
                'Vertices', pelvis(1).import.processed.vertices, ...
                'FaceColor',[.8 .8 .8],'FaceAlpha',.20,'EdgeColor','none');

% Original sphere model (optional)
% hOrig = gobjects(size(origSpheres,2),1);  
% for s = 1:size(origSpheres,2)
%     [xs,ys,zs] = sphere(20);
%     hOrig(s) = surf(xs*origSpheres(4,s)+origSpheres(1,s), ...
%          ys*origSpheres(4,s)+origSpheres(2,s), ...
%          zs*origSpheres(4,s)+origSpheres(3,s), ...
%          'EdgeColor','none','FaceColor',[.7 .7 .7],'FaceAlpha',.05);
% end

% Scaled sphere model
% for s = 1:size(scaledSpheres,2)
%     [xs,ys,zs] = sphere(20);
%     surf(xs*scaledSpheres(4,s)+scaledSpheres(1,s), ...
%          ys*scaledSpheres(4,s)+scaledSpheres(2,s), ...
%          zs*scaledSpheres(4,s)+scaledSpheres(3,s), ...
%          'EdgeColor','none','FaceColor',TUMcolors.blue300,'FaceAlpha',.25);
% end
% hScaled = get(gca,'Children'); hScaled = hScaled(1);

% Intersection points/voxel
hPts = scatter3(ptsModel(:,1), ptsModel(:,2), ptsModel(:,3), ...
    20, 'o', 'MarkerFaceColor', 'g', ... % 20
    'MarkerEdgeColor', TUMcolors.grey50, 'LineWidth', 0.5); % 0.5

% Reference acetabulum centre
plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
    pelvis(1).import.processed.acentre(3), '.','Color','m', 'MarkerSize', 35);
plot3(pelvis(1).import.processed.acentre(1), pelvis(1).import.processed.acentre(2), ...
    pelvis(1).import.processed.acentre(3), 'x','Color','m', 'MarkerSize', 15, 'Linewidth', 4.5);

% Format
axis equal;
view(3);
view(80, -10) % Figure Paper
%view(260, -10) % Figure Paper
xlabel('X'); ylabel('Y'); zlabel('Z');
title({['Best-Fit Model (score' bestModelType ') – #' num2str(modelNr)], ... 
       ['Typ: ' typeName ', Base: ' baseName ', Pct: ' pctName], ...
       ['Scale Factor = ' num2str(modelScale(modelNr).scaleFactorToRef,'%.3f')]});
camlight; lighting gouraud;
% legend([hPelvis,hOrig(1),hScaled,hPts], ...
%        {'Reference pelvis','Original spheres','Scaled spheres','Intersection voxels'}, ...
%        'Location','best');
legend([hPelvis,hPts], {'Reference pelvis','Intersection voxels'}, 'Location','best');
grid off;
hold off;

clear pctName modelNr models modelScale origSpheres scaledSpheres idxModel ptsModel hPelvis hOrig hScaled hPts s xs ys zs

%% Display sphere model mesh (for control)

% Inputs
% baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
i = 1; % Base intersect index %%%
baseName = baseIntersect{i}; % Base intersect name
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 1; %%%
typeName = intersectType{j};     
% Best model
bestModelType = 'Mean'; % 'Sum', 'Mean', 'Global' %%%
bestModelNr = pelvisDefectAnalysis(i).approx.(typeName).bestModel.(['score' bestModelType 'ModelNr']);

figure
hold on;

% Sphere model mesh
pts =  pelvisDefectAnalysis(i).approx.(typeName).bestModelDist.(['sphere' bestModelType 'Mesh']);
scatter3(pts(:,1), pts(:,2), pts(:,3), 5, 'filled');

% Format
axis equal; 
view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title({['Fibonacci point cloud of sphere model (best fit model: ' bestModelType ') – '], ...
       [ 'model #' num2str(bestModelNr) ', Typ: ' typeName ', Base: ' baseName]});
camlight; lighting gouraud;
%legend([hPelvis,hPts], {'Reference pelvis','Intersection voxels'}, 'Location','best');
grid off;
hold off;

clear bestModelType bestModelNr pts

%% Display all defects with sphere model (for control)

% Inputs
% baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
i = 1; % Base intersect index %%%
baseName = baseIntersect{i}; % Base intersect name
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 3; %%%
typeName = intersectType{j};     
% Best model
bestModelType = 'Mean'; % 'Sum', 'Mean', 'Global' %%%
bestModelNr = pelvisDefectAnalysis(i).approx.(typeName).bestModel.(['score' bestModelType 'ModelNr']);
bestModelPct = pelvisDefectAnalysis(i).approx.(typeName).bestModel.(['score' bestModelType 'Pct']);
bestModelScaledSpheres = pelvisDefectAnalysis(i).approx.(typeName).(bestModelPct).modelScale(bestModelNr).scaledSpheres;
%spherePts = pelvisDefectAnalysis(i).approx.(typeName).bestModel.(['sphere' bestModelType 'Mesh']); %%%

relevantIDs = pelvisDefectAnalysis(i).approx.(typeName).IDs(:).';   % Zeilen-Vektor
colors = lines(numel(relevantIDs));    

figure('Name','Defects + Sphere Model');
hold on;

% Defects
for k = 1:numel(relevantIDs)
    d   = relevantIDs(k);
    vtx = pelvisDefect(d).transform.trafo.(weightKabsch{w}).vertices;
    fcs = pelvisDefect(d).transform.trafo.(weightKabsch{w}).faces;
    patch('Vertices', vtx, 'Faces', fcs, ...
          'FaceColor', colors(k,:), ...
          'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
          'DisplayName', ['Defect #' num2str(d)]);
end

% Sphere model 
%scatter3(spherePts(:,1), spherePts(:,2), spherePts(:,3), 1, 'k', 'filled', 'DisplayName', 'Sphere model');
for s = 1:size(bestModelScaledSpheres,2)
    [x, y, z] = sphere(50);
    x = x * bestModelScaledSpheres(4,s) + bestModelScaledSpheres(1,s);
    y = y * bestModelScaledSpheres(4,s) + bestModelScaledSpheres(2,s);
    z = z * bestModelScaledSpheres(4,s) + bestModelScaledSpheres(3,s);
    hS = surf(x, y, z, ...
              'EdgeColor','none', ...
              'FaceColor', 'k', ...
              'FaceAlpha', 0.75);
    if s == 1
        hSphereLegend = hS;  
    end
end

% Format
axis equal; view(3); grid off;
xlabel('X'); ylabel('Y'); zlabel('Z');
title({'Best-Fit Sphere Model vs. Defects', ...
       ['Type: ' typeName ', Base: ' baseName]});
%legend('show');
camlight; lighting gouraud;
hold off;

clear bestModelType bestModelNr bestModelPct bestModelScaledSpheres spherePts relevantIDs colors d vtx fcs s x y z hS hSphereLegend

%% Display sphere models of different intersect types (for control)

% Inputs
% baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
i = 7; % Base intersect index %%%
% Best model
bestModelType = 'Mean'; % 'Sum', 'Mean', 'Global' %%%
% Colours
cols = lines(numel(intersectType)); 

figure('Name','Best-Fit Sphere Models');
hold on;

% Reference pelvis
patch('Faces', pelvis(1).import.processed.faces, ...
      'Vertices', pelvis(1).import.processed.vertices, ...
      'FaceColor',[.85 .85 .85],...
      'FaceAlpha',.15,...
      'EdgeColor','none');

hLegend = [];  
lblLegend = {};

% Sphere models for intersect types
% intersectType = {'alpha','refinedAlpha','inside','grid'};   
for j = [1,3] %1:numel(intersectType)
    typeName = intersectType{j};
        % Sphere model
        bestModelNr = pelvisDefectAnalysis(i).approx.(typeName).bestModel.(['score' bestModelType 'ModelNr']);
        bestModelPct = pelvisDefectAnalysis(i).approx.(typeName).bestModel.(['score' bestModelType 'Pct']);
        spheres = pelvisDefectAnalysis(i).approx.(typeName).(bestModelPct).modelScale(bestModelNr).scaledSpheres;
        [xs,ys,zs] = sphere(50);
        for s = 1:size(spheres,2)
            h = surf(xs*spheres(4,s)+spheres(1,s), ...
                ys*spheres(4,s)+spheres(2,s), ...
                zs*spheres(4,s)+spheres(3,s), ...
                'EdgeColor','none', ...
                'FaceAlpha',0.25, ...
                'FaceColor',cols(j,:), ...
                'DisplayName',typeName);
            if s==1               
                hLegend(end+1) = h;
                lblLegend{end+1} = typeName;
            end
        end
end

axis equal; 
view(3); 
grid off;
xlabel('X'); ylabel('Y'); zlabel('Z');
title({'Best-Fit Sphere Models (all intersect types)', ...
       ['Base: ' baseIntersect{i} ', Score: ' bestModelType]});
legend(hLegend, lblLegend,'Location','bestoutside');
camlight; 
lighting gouraud;
hold off;

clear bestModelType cols bestModelNr bestModelPct spheres s xs ys zs hLegend lblLegend h

%% Display vertex-to-nearest neighbour: sphere model - defet mesh

% Inputs
% baseIntersect = {'all','BTE','other','2A','2B','2C','3A','3B','BTE2A','BTE2B','BTE2C','BTE3A','BTE3B'}; 
i = 1; % Base intersect index %%%
baseName = baseIntersect{i}; % Base intersect name
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 1; %%%
typeName = intersectType{j};     
% Defect
defectNr = 47; %%%
% Model
bestModelType = 'Mean'; % 'Sum', 'Mean', 'Global'
bestModel = pelvisDefectAnalysis(i).approx.(typeName).bestModel;
bestModelNr = bestModel.(['score' bestModelType 'ModelNr']);
bestModelPct = bestModel.(['score' bestModelType 'Pct']);
bestModelScaledSpheres = pelvisDefectAnalysis(i).approx.(typeName).(bestModelPct).modelScale(bestModelNr).scaledSpheres;
faceColour  = pelvisDefectAnalysis(i).approx.(typeName).bestModelDist.nearDistMeanFaceColour{defectNr};
% Defect mesh
facesDef = pelvisDefect(defectNr).transform.trafo.(weightKabsch{w}).faces;
vertsDef = pelvisDefect(defectNr).transform.trafo.(weightKabsch{w}).vertices;

figure('Name','NN-Distance: Sphere ↔ Defect');
hold on;

% Defect
patch('Faces',facesDef, 'Vertices',vertsDef, ...
      'FaceVertexCData',faceColour, ...
      'FaceColor','flat', ...
      'EdgeColor','none');

% Sphere model
for s = 1:size(bestModelScaledSpheres,2)
    [x, y, z] = sphere(50);
    x = x * bestModelScaledSpheres(4,s) + bestModelScaledSpheres(1,s);
    y = y * bestModelScaledSpheres(4,s) + bestModelScaledSpheres(2,s);
    z = z * bestModelScaledSpheres(4,s) + bestModelScaledSpheres(3,s);
    surf(x, y, z, ...
              'EdgeColor','none', ...
              'FaceColor', 'k', ...
              'FaceAlpha', 0.15);
end

% Format
colormap(viridis);
clim([0 1]); 
cb = colorbar;
cb.Ticks = 0:0.2:1; 
axis equal; 
view(3); 
grid off;
title({['V2NN – Base ' baseName ', Type ' typeName], ...
       ['Defect #' num2str(defectNr) ', Score: ' bestModelType]});
xlabel('X'); ylabel('Y'); zlabel('Z');
camlight; 
lighting gouraud;
hold off;

clear bestModelType bestModel bestModelNr bestModelPct bestModelScaledSpheres faceColour facesDef vertsDef s x y z cb

%% Save and load properties of Class Approxi (Defect) - per intersection type 

%%% For intersect type
% Save pelvisIntersect data for each intersect type
savePelvisApproxType = struct();
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 1; %%%
% Meta infos
metaApproxType = fieldnames(pelvisDefectAnalysis(1).approx.(intersectType{j}));
% Iterate over each base intersect
for i = 1:numBaseIntersect
    for l = 1:length(metaApproxType)
        fieldName = metaApproxType{l};
        % Save the field data to the struct
        if isfield(pelvisDefectAnalysis(i).approx.(intersectType{j}), fieldName)
            savePelvisApproxType.(intersectType{j})(i).(fieldName) = pelvisDefectAnalysis(i).approx.(intersectType{j}).(fieldName);
        end
    end
end
% Save data to file
save(['.\pelvisApprox(' intersectType{j} ').mat'], 'savePelvisApproxType', '-v7.3'); % Adjust path if needed
% Clear the temporary save structure
clear savePelvisApproxType metaApproxType fieldName

% Load pelvisIntersect data for each intersect type
% intersectType = {'alpha','refinedAlpha','inside','grid'};
j = 1; %%%
load(['C:\Users\ga56man\#Projects\IntersectionPelvis\Workspace\pelvisApprox(' intersectType{j} ').mat']);
metaApproxType = fieldnames(savePelvisApproxType.(intersectType{j})(1));
% Iterate over each base intersect
for i = 1:numBaseIntersect
    for l = 1:length(metaApproxType)
        fieldName = metaApproxType{l};
        if isfield(savePelvisApproxType.(intersectType{j})(i), fieldName)
            pelvisDefectAnalysis(i).approx.(intersectType{j}).(fieldName) = savePelvisApproxType.(intersectType{j})(i).(fieldName);
        end
    end
end
% Clear the temporary save structure
clear savePelvisApproxType metaApproxType fieldName

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
