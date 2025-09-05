%% Approximation with the clumped sphere algorithmn of Sihaeri

% Input:    typeName: 'alpha' | 'refinedAlpha' | 'inside' | 'grid'
%           pctName: Subfield label (e.g., 'p50', 'p75').
%           fileName: Path to the input mesh 
%           nNodes: Number of seed nodes used to populate spheres
%           scaleFactor: Global scale factor forwarded to populateSpheres
%           smoothFactor: Surface smoothing factor forwarded to populateSpheres

% Output:   approx: approx.(typeName).(pctName).spheres    % clumped-sphere assembly
%                   approx.(typeName).(pctName).volSpheres % volume estimated from spheres

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich



function [approx] = clumpedSphere(typeName,pctName,fileName,nNodes,scaleFactor,smoothFactor)

% sihaeri/DEM-ClumpedSphere: https://de.mathworks.com/matlabcentral/fileexchange/67754-sihaeri-dem-clumpedsphere
% Paramter for pupulateSpheres:
% fileName, nNodes, scaleFactor, smoothFact, writeTec, writeLammpsTemp, lammpsType, writeStlScaled, rho, findSTLmI
[assembly, ~, volSpheres, ~] = populateSpheres(fileName, nNodes, scaleFactor, smoothFactor, 0, 0, 1, 0, 1, 0);
approx.(typeName).(pctName).spheres = assembly;
approx.(typeName).(pctName).volSpheres = volSpheres;

end