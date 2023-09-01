classdef LineProfile
    % LINEPROFILE storage class used by the PaperSegmentation class
    % author: Felix Binder
    % e-mail: felix.binder@fmt.fau.de
    % organisation: Chair of Manufacturing Metrology,
    %               Friedrich-Alexander University,
    %               91052 Erlangen, Germany
    % date (dd.mm.yyyy): 01.09.2023
    % version: 1.0
    % description: This is a storage class for the PaperSegmentation class.
    
    %% Public properties
    properties (Access = public)
        lineProfile = 0; % [G]
        linePeaks   = 0; % [G] 
        linePeakIdx = 0; % [-]
        foundPeaks  = 0; % [-]
        lineValleys = 0; % [G]
        lineValleyIdx = 0; % [-]
        peakDistances = 0; % [um]
        peakDIPs = 0; % [pct]
    end
    
    %% Public methods
    methods (Access = public)
        % Constructor
        function obj = LineProfile(lineProfile)
            obj.lineProfile = lineProfile;
        end
    end
end