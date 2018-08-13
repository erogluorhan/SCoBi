classdef DielMLParams < handle
    % DIELMLPARAMS Maintains specific Dielectric Profile parameters for 
    % MultiLayer soil simulations
    % It keeps the parameters that are specific to multi-layer ground with
    % row structures and each simulation. It can have only one instance 
    % throughout the whole simulation thanks to Singleton Pattern. Its 
    % properties should be initialized once in the simulation and then used
    % by other entities by using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        
        % 2nd Order Polynomial Fit
        eps_diel_z2nd
        
        % 3rd Order Polynomial Fit
        eps_diel_z3rd
        
        % Logistic Function
        eps_diel_zL
        
        % Discrete Slab
        eps_diel_zS
        
        % Figure 1 to draw VSM and reflectivity data 
        fig1
        
        % Figure 1 to draw VSM and reflectivity data
        fig2
        
    end
    
    methods (Access = private)
    
        function obj = DielMLParams
            % DIELMLPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = DielMLParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods        
        
        function initialize(obj, fig1, fig2, eps_diel_z2nd, eps_diel_z3rd, eps_diel_zL, eps_diel_zS )
            % INITIALIZE - Initializes all the properties
                     
            obj.fig1 = fig1;
            obj.fig2 = fig2;
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
               
        
        function out = get.fig1(obj)
            out = obj.fig1;        
        end
               
        
        function out = get.fig2(obj)
            out = obj.fig2;        
        end
        
    end 
    
end

