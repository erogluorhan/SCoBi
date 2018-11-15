classdef GndDynParams < handle
% class GndDynParams
%
%   Maintains ground dynamic parameters. It keeps the parameters that are 
%   specific to the dynamic surface and each simulation. It can have only 
%   one instance throughout the whole simulation thanks to Singleton 
%   pattern. Its properties are update in every simulation iteration and 
%   then used by other entities by using the get() functions provided by it. 
%
%   See also updateGndDynParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0



    properties (SetAccess = private, GetAccess = public)
                
        % Effective surface roughness parameter
        h
        
        % Dielectric constant
        eps_g
        
        
    end
    
    
    methods (Access = private)
    
        function obj = GndDynParams
            % GROUNDPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = GndDynParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods
        
        function update(obj, h, eps_g )
            % INITIALIZE - Initializes all the properties
            
            obj.h = h;
                        
            obj.eps_g = eps_g;
            
        end
        
        
        function out = get.eps_g(obj)
            out = obj.eps_g;
        end
        
        
        function out = get.h(obj)
            out = obj.h;
        end
        
    end
    
end

