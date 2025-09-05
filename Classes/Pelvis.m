 classdef Pelvis
    % Pelvis class for the pelvis attributes
    
    % Developed by C.Micheler,
    % Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
    % Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

    
    properties 
        import = Import(0);         % Object array for class Import: Loaded stl data 
        landmarks = struct();       % Struct for anatomical landmarks (class Landmark)
        boundaries = Boundary;      % Object array for class Boundary: Boundary body around pelvis (scaling)
        curve = Curvature;          % Object array for class Curvature: Curvature of the whole pelvis
        transform = Transform;      % Object array for class Transform: Transformed pelvis
    end
    
    methods
        %% Constructer: generate object
        function obj = Pelvis()
            disp('class pelvis initialized')
        end  
    end
    
 end