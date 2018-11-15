classdef RotMatDynParams < handle
% class RotMatDynParams
%
%   Maintains rotation matrices dynamic params. It keeps the parameters 
%   that are specific to the bistatic rotations and updated in the 
%   configuration of each simulation. It can have only one instance 
%   throughout the whole simulation thanks to Singleton pattern. Its 
%   properties are update in every simulation iteration and then used by 
%   other entities by using the get() functions provided by it. 
%
%   See also updateRotMatDynParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0    
    

    
    properties (SetAccess = private, GetAccess = public)
        
        % Rotation Matrix (Specular to Receiver)
        u_sr
        
        % Rotation Matrix (Transmitter to Receiver)
        u_tr
        
        % Rotation Matrix (Transmitter to Specular)
        u_ts
        
        % Rotation Matrix (Transmitter Image to Specular)
        u_tIs
        
    end
    
    
    methods (Access = private)
    
        function obj = RotMatDynParams
            % ROTMATDYNPARAMS - Private constructor
        end
    
    end
    
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = RotMatDynParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods
        
        function update(obj, u_sr, u_tr, u_ts, u_tIs )
            % INITIALIZE - Initializes all the properties
            
            obj.u_sr = u_sr;
            obj.u_tr = u_tr;
            obj.u_ts = u_ts;
            obj.u_tIs = u_tIs;
            
        end
        
        function out = get.u_sr(obj)
            out = obj.u_sr;        
        end 
        
        function out = get.u_tr(obj)
            out = obj.u_tr;        
        end 
        
        function out = get.u_ts(obj)
            out = obj.u_ts;        
        end 
        
        function out = get.u_tIs(obj)
            out = obj.u_tIs;        
        end 
        
    end
    
end

