classdef Directories < handle
    %DIRECTORIES Class to keep track of all source code and input directories
    %   This class has one attribute for each source code or input
    %   directory. Every attribute can be reached by a static getter method.
    
    
    properties (SetAccess = private, GetAccess = public)
        
        main_dir
        
        scobi        
        scobi_gui_last_input    
        scobi_gui_scobi        
        multi_layer
        vegetation
        
        input
        input_config
        input_sys
        input_ant_pat_Rx
        input_veg
        
        sims_main
            
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
            obj.input_config = fullfile( obj.input, 'configuration');
            obj.input_sys = fullfile( obj.input, 'system');
            obj.input_ant_pat_Rx = fullfile( obj.input, 'Rx_antenna_pattern');
            obj.input_veg = fullfile( obj.input, 'vegetation');
            
            
            %% SCOBI DIRECTORIES
            lib = fullfile(obj.main_dir, 'lib');
            
            obj.scobi = fullfile(lib, 'scobi');
            
            obj.multi_layer = fullfile(lib, 'multilayer');
            
            obj.vegetation = fullfile(lib, 'vegetation');

            
            %% GUI DIRECTORIES
            scobi_gui = fullfile(obj.scobi, 'gui');

            obj.scobi_gui_scobi = fullfile( scobi_gui, 'scobi');

            obj.scobi_gui_last_input = fullfile( scobi_gui, 'last_input');

            
            %% GUI DIRECTORIES
            obj.sims_main = fullfile( obj.main_dir, 'sims');
            
        end
        
        function out = get.main_dir(obj)
            out = obj.main_dir;        
        end
        
        function out = get.input(obj)
            out = obj.input;
        end
        
        function out = get.input_config(obj)
            out = obj.input_config;
        end
        
        function out = get.input_sys(obj)
            out = obj.input_sys;
        end
        
        function out = get.input_ant_pat_Rx(obj)
            out = obj.input_ant_pat_Rx;
        end
            
        function out = get.input_veg(obj)
            out = obj.input_veg;
        end
        
        function out = get.multi_layer(obj)
            out = obj.multi_layer;
        end
        
        function out = get.vegetation(obj)
            out = obj.vegetation;
        end
        
        function out = get.scobi_gui_scobi(obj)
            out = obj.scobi_gui_scobi;
        end
        
        function out = get.scobi_gui_last_input(obj)
            out = obj.scobi_gui_last_input;
        end
        
        function out = get.sims_main(obj)
            out = obj.sims_main;
        end
    end
    
end

