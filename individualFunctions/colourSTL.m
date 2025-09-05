% Coloured stl
% For open and closed meshes

% Input:    TR: triangulated Data with TR.Points (vertices) and TR.ConnectivityList (faces)
%           attribute: curve values for coloured stl
%           pelvisID: ID of the pelvis (for naming)
%           stlName: name for stlwrite

% Output:   normCurve: normalized curve data
%                       non-linear normalization for non-linear colour distribution, 
%                       otherwise colour differentiation of the clustered value range is hardly possible
%                       colour distribution adjusted to 95% of the data (percentiles 2.5% and 97.5%)
%           RGBnormCurve: normalized curve data as the corresponding RGB colour value
%           stl-file: coloured stl with adjusted curvature data / colour scale

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function [normCurve, RGBnormCurve] = colourSTL(TR, attribute, pelvisID, stlName)

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
    normCurve = normalData;
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
    RGBnormCurve = rgbColour;
    
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
    
    disp('write generated');

end