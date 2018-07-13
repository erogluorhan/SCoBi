classdef TxParams < handle
    %% TXPARAMS CLASS - Maintains transmitter parameters
    % It keeps the parameters that are specific to the transmitter in the
    % bistatic configuration of each simulation. It can have only one 
    % instance throughout the whole simulation thanks to Singleton Pattern. 
    % Its properties should be initialized once in the simulation and then 
    % used by other entities by using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        
        % Operating frequency (MHz)
        f_MHz 
        
        % Transmitter range from Earth's center (m)
        r_Tx_m 
        
        % Pt * G0t - Equivalent Isotropic Radiated Power
        EIRP_dB 
        
        % Transmitter Antenna pattern
        g_t = [ 1, 0; 0, 1 ];
        
        % Transmitter Polarization State 1
        e_t1 = [ 1; 0 ];
        
        % Transmitter Polarization State 2
        e_t2 = [ 0; 1 ];
        
        % Transmitter antenna polarization
        pol_Tx
    end
    
    
    methods (Access = private)
    
        function obj = TxParams
            % TXPARAMS - Private constructor
        end
    
    end
    
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = TxParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods
        
        function initialize(obj, f_MHz, r_Tx_m, EIRP_dB, pol_Tx )
            % INITIALIZE - Initializes all the properties
            
            obj.f_MHz = f_MHz;
            obj.r_Tx_m = r_Tx_m;
            obj.EIRP_dB = EIRP_dB; 
            obj.pol_Tx = pol_Tx;
        end
        
        function out = get.f_MHz(obj)
            out = obj.f_MHz;        
        end
        
        function out = get.r_Tx_m(obj)
            out = obj.r_Tx_m;
        end
        
        function out = get.EIRP_dB(obj)
            out = obj.EIRP_dB;
        end
        
        function out = get.g_t(obj)
            out = obj.g_t;        
        end
        
        function out = get.e_t1(obj)
            out = obj.e_t1;        
        end
        
        function out = get.e_t2(obj)
            out = obj.e_t2;        
        end
        
        function out = get.pol_Tx(obj)
            out = obj.pol_Tx;        
        end
    end
    
end

