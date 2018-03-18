classdef VegParams < handle
    %% VEGPARAMS CLASS - Maintains vegetation parameters
    % It keeps the parameters that are specific to the vegetation and each 
    % simulation. It can have only one instance throughout the whole
    % simulation thanks to Singleton Pattern. Its properties should be 
    % initialized once in the simulation and then used by other entities by 
    % using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        
        
        % Vegetation stage to separate simulations for different stages of 
        % the same vegetation 
        vegetation_stage
        
        % Vegetation particle types struct
        % L: Leaf
        % B: Branch
        % T: Trunk
        % N: Needle
        TYPES = struct('L', 1, 'B', 2, 'T', 3, 'N', 4);
        
        % Layer dimensions
        dim_layers = 0;
        
        % Number of layers
        num_layers = 1
        
        
        TYPKND
        
       
        sTYPKND
        
        
        num_types;
        
        
        num_kinds;
        
        
        scat_cal_veg;
        
        
        LTK
        
        
        dsty
        
        
        dim1
        
        
        dim2
        
        
        dim3
        
        
        epsr
        
        
        parm1
        
        
        parm2
        
        
    end
    
    
    methods (Access = protected)
    
        function obj = VegParams
            % VEGPARAMS - Protected constructor
        end
    
    end
    
    methods (Static)
        
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = VegParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods
        
        function initialize(obj, vegetation_stage, dim_layers, particleIDs, particlesCell, layersCell )
            % INITIALIZE - Initializes all the properties from input file           
        
                        
            obj.vegetation_stage = vegetation_stage;
            
            obj.dim_layers = dim_layers;
            obj.num_layers = length( obj.dim_layers );         
            
            %Initialize TYPKND
            obj.TYPKND = zeros(obj.num_layers, length(fieldnames( obj.TYPES )) );            
            
            % Fill TYPKND
            for ii = 1 : obj.num_layers
                layerParticles = layersCell{ii, 1};
                
                for jj = 1 : length( layerParticles )
                    
                    part = layerParticles{1, jj};
                    
                    if strcmp( part(1,1), 'L' ) 
                        obj.TYPKND(ii,1) = obj.TYPKND(ii,1) + 1;                        
                    elseif strcmp( part(1,1), 'B' ) 
                        obj.TYPKND(ii,2) = obj.TYPKND(ii,2) + 1;                        
                    elseif strcmp( part(1,1), 'T' ) 
                        obj.TYPKND(ii,3) = obj.TYPKND(ii,3) + 1;                        
                    elseif strcmp( part(1,1), 'N' ) 
                        obj.TYPKND(ii,4) = obj.TYPKND(ii,4) + 1;                    
                    end 
                end
            end  
            
            obj.sTYPKND = sum(obj.TYPKND);
            obj.num_types = length(obj.sTYPKND(obj.sTYPKND ~= 0)) ; % L, B, T
            obj.num_kinds = max(max(obj.TYPKND)) ;
                        
            obj.dsty = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
            obj.dim1 = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
            obj.dim2 = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
            obj.dim3 = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
            obj.epsr = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
            obj.parm1 = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
            obj.parm2 = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
          
            obj.LTK = cell(obj.num_kinds, obj.num_types, obj.num_layers) ;   
            
            % to set if we calculate scattering
            obj.scat_cal_veg = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;           
            
            
            for ii = 1 : obj.num_layers
                layerParticles = layersCell{ii, 1};
                
                for jj = 1 : length( layerParticles )
                    
                    part = layerParticles{1, jj};
                    
                    type_index = strfind( layerParticles(1,:), part(1,1) );
                    type_index = find(not(cellfun('isempty', type_index)));
                    
                    part_layer_index = strfind( layerParticles(1,:), part );
                    part_layer_index = find(not(cellfun('isempty', part_layer_index)));
                    
                    part_kind = part_layer_index - type_index(1) + 1;
                    
                    part_index = strfind( particleIDs(1,:), part );
                    part_index = find(not(cellfun('isempty', part_index)));
                    
                    part_layers = [];
                    index_down = 2;
                    while ~isempty( particleIDs{index_down, part_index} )
                        part_layers = [part_layers; particleIDs{index_down, part_index}];
                        index_down = index_down + 1;
                    end
                    
                    partStruct = particlesCell{part_index};
                    
                    if strcmp( part(1,1), 'L' )                         
                        part_type = 1;
                        %part_kind = obj.TYPKND(ii,1);
                        
                    elseif strcmp( part(1,1), 'B' )                         
                        part_type = 2;
                        %part_kind = obj.TYPKND(ii,2);
                        
                    elseif strcmp( part(1,1), 'T' )                         
                        part_type = 3;
                        %part_kind = obj.TYPKND(ii,3);
                        
                    elseif strcmp( part(1,1), 'N' )                         
                        part_type = 4;
                        %part_kind = obj.TYPKND(ii,4);
                    end   
                    
                    obj.LTK{part_kind, part_type, ii} = part;
                    
                    obj.dsty(part_kind, part_type, ii) = partStruct.DENSITY * sum(obj.dim_layers) / sum(dim_layers(part_layers));
                    obj.dim1(part_kind, part_type, ii) = partStruct.DIM1;
                    obj.dim2(part_kind, part_type, ii) = partStruct.DIM2;
                    obj.dim3(part_kind, part_type, ii) = partStruct.DIM3;
                    obj.epsr(part_kind, part_type, ii) = partStruct.EPSILON;
                    obj.parm1(part_kind, part_type, ii) = partStruct.PARM1;
                    obj.parm2(part_kind, part_type, ii) = partStruct.PARM2;
                    
                    obj.scat_cal_veg(part_kind, part_type, ii) = partStruct.IS_SCATTERER;
                end
            end
            
        end 
        
        function initialize2(obj, dim_layers, scat_cal_veg,...
            TYPKND, dsty, dim1, dim2, dim3, epsr, parm1, parm2)
            % INITIALIZE - Initializes all the properties excluding LTK

            obj.scat_cal_veg = scat_cal_veg;
            obj.dim_layers = dim_layers;
            obj.TYPKND = TYPKND;
            obj.dsty = dsty;
            obj.dim1 = dim1;
            obj.dim2 = dim2;
            obj.dim3 = dim3;
            obj.epsr = epsr;
            obj.parm1 = parm1;
            obj.parm2 = parm2;
            
        end
        
        function initialize3(obj, dim_layers, scat_cal_veg,...
            TYPKND, LTK, dsty, dim1, dim2, dim3, epsr, parm1, parm2)
            % INITIALIZE - Initializes all the properties

            obj.scat_cal_veg = scat_cal_veg;
            obj.dim_layers = dim_layers;
            obj.TYPKND = TYPKND;
            obj.LTK = LTK;
            obj.dsty = dsty;
            obj.dim1 = dim1;
            obj.dim2 = dim2;
            obj.dim3 = dim3;
            obj.epsr = epsr;
            obj.parm1 = parm1;
            obj.parm2 = parm2;
            
        end
        
        function initializeStage(obj, stage)

            obj.vegetation_stage = stage;
            
        end
         
        
        function out = get.dim_layers(obj)
            out = obj.dim_layers;        
        end
        
        function out = get.num_layers(obj)
            out = obj.num_layers;
        end
        
        function out = get.TYPKND(obj)
            out = obj.TYPKND;
        end
        
        function out = get.num_types(obj)
            out = obj.num_types;
        end
        
        function out = get.num_kinds(obj)
            out = obj.num_kinds;
        end
        
        function out = get.scat_cal_veg(obj)
            out = obj.scat_cal_veg;
        end
        
        function out = get.LTK(obj)
            out = obj.LTK;        
        end
        
        function out = get.dsty(obj)
            out = obj.dsty;        
        end
        
        function out = get.dim1(obj)
            out = obj.dim1;        
        end
        
        function out = get.dim2(obj)
            out = obj.dim2;        
        end
        
        function out = get.dim3(obj)
            out = obj.dim3;        
        end
        
        function out = get.epsr(obj)
            out = obj.epsr;        
        end
        
        function out = get.parm1(obj)
            out = obj.parm1;        
        end     
        
        function out = get.parm2(obj)
            out = obj.parm2;        
        end      
        
        function write(obj)
            
            pathName = SimulationFolders.getInstance.veg;
            
            writeVar(pathName, 'dim_layers', obj.dim_layers) ;
            
            writeVar(pathName, 'TYPKND', obj.TYPKND) ;
            
            writeVar(pathName, 'scat_cal_veg', obj.scat_cal_veg) ;
            
            writeVar(pathName, 'dsty', obj.dsty) ;
            
            writeVar(pathName, 'dim1', obj.dim1) ;
            
            writeVar(pathName, 'dim2', obj.dim2) ;
            
            writeVar(pathName, 'dim3', obj.dim3) ;
            
            writeComplexVar(pathName, 'epsr', obj.epsr) ;
            
            writeVar(pathName, 'parm1', obj.parm1) ;
            
            writeVar(pathName, 'parm2', obj.parm2) ;
            
            currentDir = pwd;
            cd( SimulationFolders.getInstance.veg )
            tmp = obj.LTK;
            save LTK tmp
            cd(currentDir)
        end      
        
    end
    
end

