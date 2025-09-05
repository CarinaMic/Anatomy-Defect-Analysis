 classdef DefectAnalysis
    % DefectAnalysis class for the analysis of the different defects and methods
    
    % Developed by C.Micheler,
    % Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
    % Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich
    
    
    properties              
        intersect = Intersection;      % Object array for class Intersection
        distribut = Distribution;      % Object array for class Distribution
        approx = Approxi;              % Object array for class Approxi
    end
    
    methods
        %% Constructer: generate object
        function obj = DefectAnalysis()
            disp('class defectAnalysis initialized')
        end        
    end
    
 end