classdef Directories < handle
    %DIRECTORIES Class to keep track of all source code and input directories
    %   This class has one attribute for each source code or input
    %   directory. Every attribute can be reached by a static getter method.
    
    
    properties (SetAccess = private, GetAccess = public)
        main_dir
        lib
        
        common
        common_bistatic
        common_constants
        common_gui
        common_gui_ml
        common_gui_veg
        common_init
        common_param
        common_products
        common_SCoBi
        common_util
        
        input
        input_dyn
        input_sys
        input_ml
        input_veg
        input_veg_hom
        input_veg_vir
        input_veg_vir_plg
        input_veg_vir_row
        input_veg_vir_rnd
        
        scobi_veg
        scobi_veg_analysis
        scobi_veg_init
        scobi_veg_monte_carlo
        scobi_veg_param
        scobi_veg_plot
        scobi_veg_products
        
        scobi_ml
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
            obj.main_dir = strrep(obj.main_dir, '\lib\common\constants', '');

            obj.lib = fullfile(obj.main_dir, 'lib');
            
            
            %% INPUT DIRECTORIES
            obj.input = fullfile( obj.main_dir, 'input');
            obj.input_dyn = fullfile( obj.input, 'dyn');
            obj.input_sys = fullfile( obj.input, 'sys');
            obj.input_ml = fullfile( obj.input, 'multi_layer');
            obj.input_veg = fullfile( obj.input, 'veg');
            obj.input_veg_hom = fullfile( obj.input_veg, Constants.veg_methods{Constants.id_veg_hom} );
            obj.input_veg_vir = fullfile( obj.input_veg, Constants.veg_methods{Constants.id_veg_vir} );
            obj.input_veg_vir_plg = fullfile( obj.input_veg_vir, 'plg');
            obj.input_veg_vir_row = fullfile( obj.input_veg_vir, Constants.veg_vir_orientations{Constants.id_veg_vir_row_crop} );
            obj.input_veg_vir_rnd = fullfile( obj.input_veg_vir, Constants.veg_vir_orientations{Constants.id_veg_vir_row_crop} );
            
            
            %% COMMON DIRECTORIES
            obj.common = fullfile(obj.lib, 'common');
            
            obj.common_bistatic = fullfile(obj.common, 'bistatic');
            
            obj.common_constants = fullfile( obj.common, 'constants');

            obj.common_gui = fullfile( obj.common, 'gui');

            obj.common_gui_ml = fullfile( obj.common_gui, 'multilayer');

            obj.common_gui_veg = fullfile( obj.common_gui, 'veg');
            
            obj.common_init = fullfile( obj.common, 'init');

            obj.common_param = fullfile( obj.common, 'param');

            obj.common_products = fullfile( obj.common, 'products');

            obj.common_util = fullfile( obj.common, 'util');

            obj.common_SCoBi = fullfile( obj.common, 'SCoBi');
            

            %% SCoBi-Veg SOURCE-CODE DIRECTORIES
            obj.scobi_veg = fullfile(obj.lib, 'scobi_veg');

            obj.scobi_veg_analysis = fullfile(obj.scobi_veg, 'analysis');
            
            obj.scobi_veg_init = fullfile( obj.scobi_veg, 'init');

            obj.scobi_veg_monte_carlo = fullfile( obj.scobi_veg, 'monte_carlo');

            obj.scobi_veg_param = fullfile( obj.scobi_veg, 'param');

            obj.scobi_veg_plot = fullfile( obj.scobi_veg, 'plot');

            obj.scobi_veg_products = fullfile( obj.scobi_veg, 'products');
            

            %% SCoBi-ML SOURCE-CODE DIRECTORIES
            obj.scobi_ml = fullfile(obj.lib, 'scobi_ml');
            
        end
        
        function out = get.main_dir(obj)
            out = obj.main_dir;        
        end
        
        function out = get.lib(obj)
            out = obj.lib;        
        end
        
        function out = get.common_bistatic(obj)
            out = obj.common_bistatic;
        end
        
        function out = get.common_constants(obj)
            out = obj.common_constants;
        end
        
        function out = get.common_init(obj)
            out = obj.common_init;
        end
        
        function out = get.common_param(obj)
            out = obj.common_param;
        end
        
        function out = get.common_products(obj)
            out = obj.common_products;
        end
        
        function out = get.common_SCoBi(obj)
            out = obj.common_SCoBi;
        end
        
        function out = get.common_util(obj)
            out = obj.common_util;
        end
        
        function out = get.input(obj)
            out = obj.input;
        end
        
        function out = get.input_dyn(obj)
            out = obj.input_dyn;
        end
        
        function out = get.input_ml(obj)
            out = obj.input_ml;
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
        
        function out = get.scobi_ml(obj)
            out = obj.scobi_ml;
        end
        
        function out = get.scobi_veg(obj)
            out = obj.scobi_veg;
        end
        
        function out = get.scobi_veg_analysis(obj)
            out = obj.scobi_veg_analysis;
        end
        
        function out = get.common_gui(obj)
            out = obj.common_gui;
        end
        
        function out = get.common_gui_ml(obj)
            out = obj.common_gui_ml;
        end
        
        function out = get.common_gui_veg(obj)
            out = obj.common_gui_veg;
        end
        
        function out = get.scobi_veg_init(obj)
            out = obj.scobi_veg_init;
        end
        
        function out = get.scobi_veg_monte_carlo(obj)
            out = obj.scobi_veg_monte_carlo;
        end
        
        function out = get.scobi_veg_param(obj)
            out = obj.scobi_veg_param;
        end
        
        function out = get.scobi_veg_plot(obj)
            out = obj.scobi_veg_plot;
        end
        
        function out = get.scobi_veg_products(obj)
            out = obj.scobi_veg_products;
        end
    end
    
end

