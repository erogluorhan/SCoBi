classdef VegParams < handle
    %% VEGPARAMS CLASS - Maintains vegetation parameters
    % It keeps the parameters that are specific to the vegetation and each 
    % simulation. It can have only one instance throughout the whole
    % simulation thanks to Singleton Pattern. Its properties should be 
    % initialized once in the simulation and then used by other entities by 
    % using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        
        % Vegetation particle types struct
        % L: Leaf
        % B: Branch
        % T: Trunk
        % N: Needle
        TYPES = struct('L', 1, 'B', 2, 'T', 3, 'N', 4);
        
        % Layer dimensions
        dim_layers_m = 0;
        
        % Number of layers
        num_layers = 0;
        
        
        TYPKND;
        
       
        sTYPKND;
        
        
        num_types;
        
        
        num_kinds;
        
        
        LTK;
        
        
        dsty;
        
        
        dim1_m
        
        
        dim2_m
        
        
        dim3_m
        
        
        epsr
        
        
        parm1_deg
        
        
        parm2_deg
        
        
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
        
        function initialize( obj, dim_layers_m, num_layers, TYPKND, ...
                       num_types, num_kinds, LTK, dsty, dim1_m, dim2_m, ...
                       dim3_m, epsr, parm1_deg, parm2_deg )
      
            obj.dim_layers_m = dim_layers_m;
            obj.num_layers = num_layers;
            obj.TYPKND = TYPKND;
            obj.num_types = num_types;
            obj.num_kinds = num_kinds;
            obj.LTK = LTK;
            obj.dsty = dsty;
            obj.dim1_m = dim1_m;
            obj.dim2_m = dim2_m;
            obj.dim3_m = dim3_m;
            obj.epsr = epsr;
            obj.parm1_deg = parm1_deg;
            obj.parm2_deg = parm2_deg;
                   
        end
        
        function setup(obj, dim_layers_m, particleIDs, particlesCell, layersCell )
            % INITIALIZE - Initializes all the properties from input file  
            
            obj.dim_layers_m = dim_layers_m;
            obj.num_layers = length( obj.dim_layers_m );         
            
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
            obj.dim1_m = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
            obj.dim2_m = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
            obj.dim3_m = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
            obj.epsr = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
            obj.parm1_deg = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
            obj.parm2_deg = zeros(obj.num_kinds, obj.num_types, obj.num_layers) ;
          
            obj.LTK = cell(obj.num_kinds, obj.num_types, obj.num_layers) ;            
            
            
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
                    
                    obj.dsty(part_kind, part_type, ii) = partStruct.DENSITY * sum(obj.dim_layers_m) / sum(dim_layers_m(part_layers));
                    obj.dim1_m(part_kind, part_type, ii) = partStruct.DIM1;
                    obj.dim2_m(part_kind, part_type, ii) = partStruct.DIM2;
                    obj.dim3_m(part_kind, part_type, ii) = partStruct.DIM3;
                    obj.epsr(part_kind, part_type, ii) = partStruct.EPSILON;
                    obj.parm1_deg(part_kind, part_type, ii) = partStruct.PARM1;
                    obj.parm2_deg(part_kind, part_type, ii) = partStruct.PARM2;
                end
            end
            
        end 
         
        
        function out = get.dim_layers_m(obj)
            out = obj.dim_layers_m;        
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
        
        function out = get.LTK(obj)
            out = obj.LTK;        
        end
        
        function out = get.dsty(obj)
            out = obj.dsty;        
        end
        
        function out = get.dim1_m(obj)
            out = obj.dim1_m;        
        end
        
        function out = get.dim2_m(obj)
            out = obj.dim2_m;        
        end
        
        function out = get.dim3_m(obj)
            out = obj.dim3_m;        
        end
        
        function out = get.epsr(obj)
            out = obj.epsr;        
        end
        
        function out = get.parm1_deg(obj)
            out = obj.parm1_deg;        
        end     
        
        function out = get.parm2_deg(obj)
            out = obj.parm2_deg;        
        end    
        
    end
    
end

