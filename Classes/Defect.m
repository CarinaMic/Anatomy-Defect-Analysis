 classdef Defect
    % Defect class for the defect attributes
    
    % Developed by C.Micheler,
    % Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
    % Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich
    
    
    properties 
        import = Import(0);     % Import: loaded stl data                   
        curve = Curvature;      % Object array for class Curvature: Curvature of defect area
        transform = Transform;  % Object array for class Transform: Transform defect area 
        volume = Volume;        % Object array for class Volume: Create volume of defect area 
        boundaries  = Boundary; % Object array for class boundary: Boundary body around pelvis (scaling) 
    end
    
    methods
        %% Constructer: generate object
        function obj = Defect()
            disp('class defect initialized')
        end        
    end
    
 end