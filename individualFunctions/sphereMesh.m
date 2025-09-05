%% Create mesh für sphere model (Fibonacci)

% Input:    typeName: 'alpha' | 'refinedAlpha' | 'inside' | 'grid'
%           i: Numeric identifier used only for logging
%           bestModelType: Suffix used to name the output field (e.g., 'P50', 'P75')
%           spheres: 4xN array of sphere parameters [x; y; z; r] (centers and radii)
%           spacing: spacing on sphere surfaces (same units as coordinates)

% Output:   approx: updated struct

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [approx] = sphereMesh(i, typeName,  bestModelType, spheres, spacing)

nSpheres = size(spheres,2);
outerPts = [];

% Fibonacci points on unit sphere
goldAng = pi * (3 - sqrt(5));   % golden angle

for s = 1:nSpheres
    centre = spheres(1:3,s)';
    radius = spheres(4,s);
    % Required number of points for desired spacing
    approxN = max(round(4*pi*radius^2 / spacing^2), 50);   % min. 50 points
    k = (0:approxN-1)';       % Fibonacci Index
    z = 1 - 2*(k+0.5)/approxN;
    rad = sqrt(1 - z.^2);
    phi = k * goldAng;
    ptsLocal = [rad .* cos(phi), rad .* sin(phi), z];  % unit sphere
    ptsWorld = radius * ptsLocal + centre;  % scale and transform

    % Filter points inside
    inOther = false(size(ptsWorld,1),1);
    for j = 1:nSpheres
        if j == s
            continue;
        end
        cj = spheres(1:3,j)';
        rj = spheres(4,j);
        d  = vecnorm(ptsWorld - cj, 2, 2);            % Distance to sphere j
        inOther = inOther | (d < rj - 1e-6);          % Tolerance
    end
    outerPts = [outerPts; ptsWorld(~inOther,:)];
end

approx.(typeName).bestModelDist.(['sphere' bestModelType 'Mesh']) = outerPts;

disp(['Created sphere mesh: i=' num2str(i) ', typeName=' typeName]);

end