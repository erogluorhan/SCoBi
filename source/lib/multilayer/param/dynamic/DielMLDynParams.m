classdef DielMLDynParams < handle
% class DielMLDynParams 
%
%   Maintains dynamic dielectric profile parameters for multi-layered 
%   ground, if any. It keeps the parameters that are specific to 
%   multi-layer ground with row structures and are updated in each 
%   simulation iteration. It can have only one instance throughout the 
%   entire simulation thanks to Singleton Pattern. Its properties are 
%   updated in every simulation iteration and are then used by other 
%   entities by using the get() functions provided. 
%
%   See also generateDielMLProfiles.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


    
    properties (SetAccess = private, GetAccess = public)
        
        % 2nd Order Polynomial Fit
        eps_diel_z2nd
        
        % 3rd Order Polynomial Fit
        eps_diel_z3rd
        
        % Logistic Function
        eps_diel_zL
        
        % Discrete Slab
        eps_diel_zS
        
    end
    
    methods (Access = private)
    
        function obj = DielMLDynParams
            % DIELMLDYNPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = DielMLDynParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods        
        
        function initialize(obj, eps_diel_z2nd, eps_diel_z3rd, eps_diel_zL, eps_diel_zS )
            % INITIALIZE - Initializes all the properties
                     
            obj.eps_diel_z2nd = eps_diel_z2nd;
            obj.eps_diel_z3rd = eps_diel_z3rd;
            obj.eps_diel_zL = eps_diel_zL;
            obj.eps_diel_zS = eps_diel_zS;
                      
        end 
         
        
        function out = get.eps_diel_z2nd(obj)
            out = obj.eps_diel_z2nd;        
        end
         
        
        function out = get.eps_diel_z3rd(obj)
            out = obj.eps_diel_z3rd;        
        end
         
        
        function out = get.eps_diel_zL(obj)
            out = obj.eps_diel_zL;        
        end
         
        
        function out = get.eps_diel_zS(obj)
            out = obj.eps_diel_zS;        
        end
        
    end 
    
end

