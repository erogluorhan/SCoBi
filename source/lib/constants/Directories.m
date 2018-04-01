classdef Directories < handle
    %DIRECTORIES Class to keep track of all source code and input directories
    %   This class has one attribute for each source code or input
    %   directory. Every attribute can be reached by a static getter method.
    
    
    properties (SetAccess = private, GetAccess = public)
        main_dir
        lib
        analysis
        bistatic
        constants
        gui
        init
        input
        input_sys
        input_veg
        input_veg_hom
        input_veg_vir
        input_veg_vir_plg
        input_veg_vir_row
        input_veg_vir_rnd
        monte_carlo
        param
        plot
        products
        SCoBi
        util
    end
    
    
    methods (Access = private)
        
        function obj = Directories
        end
        
    end
    
    
    methods (Static)
        function singleObj = getInstance
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = Directories;
                localObj.initialize();
             end
             
             singleObj = localObj;
        end
    end
        
    
    methods
        
        function initialize(obj)
            % Initialize attribute for each source code or input directory

            obj.main_dir = fileparts(mfilename('fullpath'));
            obj.main_dir = strrep(obj.main_dir, '\lib\constants', '');

            obj.lib = fullfile(obj.main_dir, 'lib');

            obj.analysis = fullfile(obj.lib, 'analysis');

            obj.bistatic = fullfile(obj.lib, 'bistatic');

            obj.constants = fullfile( obj.lib, 'constants');

            obj.gui = fullfile( obj.lib, 'gui');
            
            obj.init = fullfile( obj.lib, 'init');
            
            obj.input = fullfile( obj.main_dir, 'input');
            obj.input_sys = fullfile( obj.input, 'sys');
            obj.input_veg = fullfile( obj.input, 'veg');
            obj.input_veg_hom = fullfile( obj.input_veg, Constants.veg_methods.HOMOGENOUS);
            obj.input_veg_vir = fullfile( obj.input_veg, Constants.veg_methods.VIRTUAL);
            obj.input_veg_vir_plg = fullfile( obj.input_veg_vir, 'plg');
            obj.input_veg_vir_row = fullfile( obj.input_veg_vir, Constants.veg_vir_types.ROW );
            obj.input_veg_vir_rnd = fullfile( obj.input_veg_vir, Constants.veg_vir_types.RANDOM );

            obj.monte_carlo = fullfile( obj.lib, 'monte_carlo');

            obj.param = fullfile( obj.lib, 'param');

            obj.plot = fullfile( obj.lib, 'plot');

            obj.products = fullfile( obj.lib, 'products');

            obj.SCoBi = fullfile( obj.lib, 'SCoBi');

            obj.util = fullfile( obj.lib, 'util');
            
        end
        
        function out = get.main_dir(obj)
            out = obj.main_dir;        
        end
        
        function out = get.lib(obj)
            out = obj.lib;        
        end
        
        function out = get.analysis(obj)
            out = obj.analysis;
        end
        
        function out = get.bistatic(obj)
            out = obj.bistatic;
        end
        
        function out = get.constants(obj)
            out = obj.constants;
        end
        
        function out = get.gui(obj)
            out = obj.gui;
        end
        
        function out = get.init(obj)
            out = obj.init;
        end
        
        function out = get.input(obj)
            out = obj.input;
        end
        
        function out = get.input_sys(obj)
            out = obj.input_sys;
        end
            
        function out = get.input_veg(obj)
            out = obj.input_veg;
        end
            
        function out = get.input_veg_hom(obj)
            out = obj.input_veg_hom;
        end
            
        function out = get.input_veg_vir(obj)
            out = obj.input_veg_vir;
        end
            
        function out = get.input_veg_vir_plg(obj)
            out = obj.input_veg_vir_plg;
        end
            
        function out = get.input_veg_vir_row(obj)
            out = obj.input_veg_vir_row;
        end
            
        function out = get.input_veg_vir_rnd(obj)
            out = obj.input_veg_vir_rnd;
        end
        
        function out = get.monte_carlo(obj)
            out = obj.monte_carlo;
        end
        
        function out = get.param(obj)
            out = obj.param;
        end
        
        function out = get.plot(obj)
            out = obj.plot;
        end
        
        function out = get.products(obj)
            out = obj.products;
        end
        
        function out = get.SCoBi(obj)
            out = obj.SCoBi;
        end
        
        function out = get.util(obj)
            out = obj.util;
        end
    end
    
end

