classdef SatParams < handle
    %% SATPARAMS CLASS - Maintains transmitter satellite parameters
    % It keeps the parameters that are specific to the transmitter in the
    % bistatic configuration of each simulation. It can have only one 
    % instance throughout the whole simulation thanks to Singleton Pattern. 
    % Its properties should be initialized once in the simulation and then 
    % used by other entities by using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        
        % Operating frequency (MHz)
        f_MHz 
        
        % Satellite radius from Earth's center (m)
        rsat_m 
        
        % Incidence angle list - measured between ground zenith and the 
        % ground-satellite direction (degrees) 
        th0_list_deg
        
        % Elveation angle list (degrees)
        EL0_list_deg
        
        % Azimuth angle list (degrees)- the horizontal angle measured at the 
        % transmitter's position on the earth to Northpole, i.e. measured 
        % clockwise from the North pole IMPORTANT: It is used in Geo 
        % Staellite comm. and different from standard spherical azimuth 
        % direction that is measured from local East in the 
        % counter-clockwise direction.
        PH0_list_deg
        
        % Pt * G0t - Equivalent Isotropic Radiated Power
        EIRP_dB 
        
        % Transmitter Antenna pattern
        g_t = [ 1, 0; 0, 1 ];
        
        % Transmitter Polarization State 1
        e_t1 = [ 1; 0 ];
        
        % Transmitter Polarization State 2
        e_t2 = [ 0; 1 ];
        
        % Satellite polarization
        polT
    end
    
    
    methods (Access = private)
    
        function obj = SatParams
            % SATPARAMS - Private constructor
        end
    
    end
    
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = SatParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods
        
        function initialize(obj, f_MHz, rsat_m, th0_list_deg, PH0_list_deg, EIRP_dB, polT )
            % INITIALIZE - Initializes all the properties
            
            obj.f_MHz = f_MHz;
            obj.rsat_m = rsat_m;
            obj.th0_list_deg = th0_list_deg;
            obj.EL0_list_deg = 90 - th0_list_deg;
            obj.PH0_list_deg = PH0_list_deg;
            obj.EIRP_dB = EIRP_dB; 
            obj.polT = polT;
        end
        
        function out = get.f_MHz(obj)
            out = obj.f_MHz;        
        end
        
        function out = get.rsat_m(obj)
            out = obj.rsat_m;
        end
        
        function out = get.th0_list_deg(obj)
            out = obj.th0_list_deg;
        end
        
        function out = get.EL0_list_deg(obj)
            out = obj.EL0_list_deg;
        end
        
        function out = get.PH0_list_deg(obj)
            out = obj.PH0_list_deg;
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
        
        function out = get.polT(obj)
            out = obj.polT;        
        end
        
        function rep_Th( obj, factor )
            obj.th0_list_deg = repmat(obj.th0_list_deg, 1, factor);
        end
        
        function rep_Ph( obj, factor )
            obj.PH0_list_deg = repmat(obj.PH0_list_deg, 1, factor);
        end
    end
    
end

