classdef BistaticDynParams < handle
    %% BISTATICDYNPARAMS CLASS - Maintains bistatic dynamic parameters
    % It keeps the parameters that are specific to the bistatic geometry 
    % and updated in the configuration of each simulation. It can have only
    % one instance throughout the whole simulation thanks to Singleton 
    % Pattern. 
    % Its properties should be initialized once in the simulation and then 
    % used by other entities by using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        
        % Slant range - distance between ground and Transmitter (meters)
        rd_m 
        
        % Transformation matrix for transforming a vector from the ground 
        % frame to transmitter system
        Tgt
        
        % Transformation matrix for transforming a vector from the ground 
        % frame to Image transmitter system
        TgtI
        
        % propagation vector (i_d^-)
        idn
        
        % propagation vector (i_s^-)
        isn
        
        % propagation vector (o_s^+)
        osp
        
        % propagation vector (o_s^-)
        osn
        
        % Transformation Gnd -> Specular
        Tgs
        
        % Transformation Gnd -> Rx
        Tgr
        
        % Transformation Gnd -> Rx Image
        TgrI
        
        % Receiver Rotation about z-axis (Azimuth rotation)
        AntRotZ_Rx
        
        % Receiver Rotation about y-axis (Elevation rotation)
        AntRotY_Rx
        
        % Receiver Rotation in both azimuth and elevation
        AntRot_Rx
        
        % Transmitter Rotation about z-axis (Azimuth rotation)
        AntRotZ_Tx
        
        % pos_Tx_m, pos_SP_m, pos_Rx_m, pos_Gnd_m, pos_B_Rx_m, 
        % pos_FP_Rx_m, pos_FZ_m in ground (refrence) frame (G)
        AllPoints_m
        
        % Incidence angle (T -> R) in receiver frame (R)
        AngT2R_rf
        
        % Incidence angle (S -> R) in receiver frame (R)
        AngS2R_rf
        
        % Incidence angle (T -> S) in specular frame (S)
        AngT2S_sf
        
    end
    
    
    methods (Access = private)
    
        function obj = BistaticDynParams
            % BISTATICDYNPARAMS - Private constructor
        end
    
    end
    
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = BistaticDynParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods
        
        function update(obj, rd_m, idn, isn, osp, osn, Tgt, TgtI, ...
                Tgs, Tgr, TgrI, AntRotZ_Rx, AntRotY_Rx, AntRot_Rx, ...
                AntRotZ_Tx, AllPoints_m, AngT2R_rf, ...
                AngS2R_rf, AngT2S_sf )
            % INITIALIZE - Initializes all the properties
            
            obj.rd_m = rd_m;
            obj.idn = idn;
            obj.isn = isn;
            obj.osp = osp;
            obj.osn = osn; 
            obj.Tgt = Tgt;
            obj.TgtI = TgtI;
            obj.Tgs = Tgs;
            obj.Tgr = Tgr;
            obj.TgrI = TgrI;
            obj.AntRotZ_Rx = AntRotZ_Rx;
            obj.AntRotY_Rx = AntRotY_Rx;
            obj.AntRot_Rx = AntRot_Rx;
            obj.AntRotZ_Tx = AntRotZ_Tx;
            obj.AllPoints_m = AllPoints_m;
            obj.AngT2R_rf = AngT2R_rf;
            obj.AngS2R_rf = AngS2R_rf;
            obj.AngT2S_sf = AngT2S_sf;
            
        end
        
        function out = get.rd_m(obj)
            out = obj.rd_m;        
        end 
        
        function out = get.idn(obj)
            out = obj.idn;        
        end 
        
        function out = get.isn(obj)
            out = obj.isn;        
        end 
        
        function out = get.osp(obj)
            out = obj.osp;        
        end 
        
        function out = get.osn(obj)
            out = obj.osn;        
        end 
        
        function out = get.Tgt(obj)
            out = obj.Tgt;        
        end 
        
        function out = get.TgtI(obj)
            out = obj.TgtI;        
        end 
        
        function out = get.Tgr(obj)
            out = obj.Tgr;        
        end 
        
        function out = get.TgrI(obj)
            out = obj.TgrI;        
        end 
        
        function out = get.Tgs(obj)
            out = obj.Tgs;        
        end 
        
        function out = get.AntRotZ_Rx(obj)
            out = obj.AntRotZ_Rx;        
        end 
        
        function out = get.AntRotY_Rx(obj)
            out = obj.AntRotY_Rx;        
        end 
        
        function out = get.AntRot_Rx(obj)
            out = obj.AntRot_Rx;        
        end 
        
        function out = get.AntRotZ_Tx(obj)
            out = obj.AntRotZ_Tx;        
        end 
                
        function out = get.AllPoints_m(obj)
            out = obj.AllPoints_m;        
        end 
        
        function out = get.AngT2R_rf(obj)
            out = obj.AngT2R_rf;        
        end 
        
        function out = get.AngS2R_rf(obj)
            out = obj.AngS2R_rf;        
        end 
        
        function out = get.AngT2S_sf(obj)
            out = obj.AngT2S_sf;        
        end 
        
    end
    
end

