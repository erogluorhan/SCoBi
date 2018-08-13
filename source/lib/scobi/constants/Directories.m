classdef Directories < handle
    %DIRECTORIES Class to keep track of all source code and input directories
    %   This class has one attribute for each source code or input
    %   directory. Every attribute can be reached by a static getter method.
    
    
    properties (SetAccess = private, GetAccess = public)
        main_dir
        
        scobi        
        scobi_gui_multilayer
        scobi_gui_veg        
        scobi_ml
        scobi_veg
        
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
            obj.main_dir = strrep(obj.main_dir, '\lib\scobi\constants', '');

                        
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
            
            
            %% SCOBI DIRECTORIES
            lib = fullfile(obj.main_dir, 'lib');
            
            obj.scobi = fullfile(lib, 'scobi');
            
            obj.scobi_ml = fullfile(lib, 'multilayer');
            
            obj.scobi_veg = fullfile(lib, 'vegetation');

            
            %% GUI DIRECTORIES
            scobi_gui = fullfile(obj.scobi, 'gui');
            
            obj.scobi_gui_multilayer = fullfile( scobi_gui, 'multilayer');

            obj.scobi_gui_veg = fullfile( scobi_gui, 'vegetation');
            
        end
        
        function out = get.main_dir(obj)
            out = obj.main_dir;        
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
        
        function out = get.scobi_gui_multilayer(obj)
            out = obj.scobi_gui_multilayer;
        end
        
        function out = get.scobi_gui_veg(obj)
            out = obj.scobi_gui_veg;
        end
    end
    
end

