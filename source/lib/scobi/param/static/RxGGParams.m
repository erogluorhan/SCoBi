
classdef RxGGParams < handle
% class RxGGParams
%
%   Maintains Generalized-Gaussian receiver antenna parameters. It can have 
%   only one instance throughout the whole simulation thanks to Singleton 
%   Pattern. Its properties should be initialized once in the simulation 
%   and then used by other entities by using the get() functions provided 
%   by it. 
%
%   See also initRxGGParams.

%   Copyright � 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


    
    properties (SetAccess = private, GetAccess = public)
        
        % Beamwidth (degrees)
        hpbw_deg 
        
        % Sidelobe Level (dB)
        SLL_dB 
        
        % X-pol level (dB)
        XPL_dB 
    end
    
    
    methods (Access = private)
    
        function obj = RxGGParams
            % RXGGPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = RxGGParams;
             end
             
             singleObj = localObj;
        end
             
    end
    
    methods
        
        function initialize( obj, hpbw_deg, SLL_dB, XPL_dB )
            % INITIALIZE - Initializes all the properties
            
            obj.hpbw_deg = hpbw_deg;
            obj.SLL_dB = SLL_dB;
            obj.XPL_dB = XPL_dB;
            
        end
        
        function out = get.hpbw_deg(obj)
            out = obj.hpbw_deg;
        end
        
        function out = get.SLL_dB(obj)
            out = obj.SLL_dB;
        end
        
        function out = get.XPL_dB(obj)
            out = obj.XPL_dB;
        end
        
        function ant_pat_struct = calcRxAntPatMatrix( obj, ant_pat_res_deg )

            ant_pat_struct = [];
            
            
            %% CALCULATIONS            
            % Calculate parameter "a" based on SLL
            % Use the known values for SLL and "a"
            sll_levels = [15; 20; 30; 40];
            a_values = [0.18; 0.11; 0.04; 0.01];
            % Fit a 3rd order polynomial to these known values
            p = polyfit(sll_levels, a_values, 3);            
            % Calculate the "a" value for the given SLL value
            a = p(1,1) * obj.SLL_dB ^ 3 + p(1,2) * obj.SLL_dB ^ 2 + ...
                 p(1,3) * obj.SLL_dB + p(1,4);
            
            % Calculate parameter V based on XPL
            xpl_levels = [15; 25; 30; 40];
            V_values = [sqrt(0.0316); sqrt(0.0100); sqrt(0.0010); sqrt(0.0001)];
            % Fit a 3rd order polynomial to these known values
            pV = polyfit(xpl_levels, V_values, 3);            
            % Calculate the "a" value for the given SLL value
            V = pV(1,1) * obj.XPL_dB ^ 3 + pV(1,2) * obj.XPL_dB ^ 2 + ...
                 pV(1,3) * obj.XPL_dB + pV(1,4);

            % Calculate Generalized-Gaussian pattern
            [th, ph, gg] = obj.GGpattern( obj.hpbw_deg, a, ant_pat_res_deg ) ;
            gXX = gg ; gYY = gg ;

            [~, ~, gg] = obj.GGpattern( 2 * obj.hpbw_deg, a, ant_pat_res_deg ) ;
            gXY = V * gg ; gYX = V * gg ;


            %% complex (voltage) co- and x- patterns : X-PORT - v-pol
            magg = sqrt(abs(gXX) .^2 + abs(gXY) .^2) ;
            maxg = max(max(magg)) ;

            % complex normalized (voltage) pattern
            gnXX = gXX / maxg ;  
            gnXY = gXY / maxg ;

            %% complex (voltage) co- and x- patterns : Y-PORT - h-pol
            magg = sqrt(abs(gYX) .^2 + abs(gYY) .^2) ;
            maxg = max(max(magg)) ;

            % complex normalized (voltage) pattern
            gnYX = gYX / maxg ;  
            gnYY = gYY / maxg ;

            %% Antenna Pattern Matrix Elements

            g11 = abs(gnXX) .^ 2 ;
            g12 = abs(gnXY) .^ 2 ;
            g13 = real(gnXX .* conj(gnXY)) ;
            g14 = -imag(gnXX .* conj(gnXY)) ;

            g21 = abs(gnYX) .^ 2 ;
            g22 = abs(gnYY) .^ 2 ;
            g23 = real(gnYX .* conj(gnYY)) ;
            g24 = -imag(gnYX .* conj(gnYY)) ;

            g31 = 2 * real(gnXX .* conj(gnYX)) ; 
            g32 = 2 * real(gnXY .* conj(gnYY)) ; 
            g33 = real(gnXX .* conj(gnYY) + gnXY .* conj(gnYX)) ;
            g34 = -imag(gnXX .* conj(gnYY) - gnXY .* conj(gnYX)) ;

            g41 = 2 * imag(gnXX .* conj(gnYX)) ; 
            g42 = 2 * imag(gnXY .* conj(gnYY)) ; 
            g43 = imag(gnXX .* conj(gnYY) + gnXY .* conj(gnYX)) ;
            g44 = real(gnXX .* conj(gnYY) - gnXY .* conj(gnYX)) ;

            %% Antenna Pattern Matrix

            G = cell(4) ;

            G{1,1} = g11 ; G{1,2} = g12 ; G{1,3} = g13 ; G{1,4} = g14 ;
            G{2,1} = g21 ; G{2,2} = g22 ; G{2,3} = g23 ; G{2,4} = g24 ;
            G{3,1} = g31 ; G{3,2} = g32 ; G{3,3} = g33 ; G{3,4} = g34 ;
            G{4,1} = g41 ; G{4,2} = g42 ; G{4,3} = g43 ; G{4,4} = g44 ;

            %%

            g = cell(2) ;

            g{1,1} = gnXX ; g{1,2} = gnXY ;
            g{2,1} = gnYX ; g{2,2} = gnYY ; 


            %% Half-Power Beanwidth
            % X-port
            Gn_co = G{1, 1} ;
            Gn_x = G{1, 2} ;

            Gn = Gn_co + Gn_x ; % co- + x- pols

            indhpbw = (Gn >= 0.49) & (Gn <= 0.51) ;
            hpbwX = mean(mean(th(indhpbw))) ;

            del = 0 ;
            while isnan(hpbwX)
                indhpbw = (Gn >= (0.50 - del)) & (Gn <= (0.50 + del)) ;
                hpbwX = mean(mean(th(indhpbw))) ;
                del = del + 0.01 ;
            end

            % Y-port
            Gn_co = G{2, 2} ;
            Gn_x = G{2, 1} ;

            Gn = Gn_co + Gn_x ; % co- + x- pols

            indhpbw = (Gn >= 0.49) & (Gn <= 0.51) ;
            hpbwY = mean(mean(th(indhpbw))) ;

            del = 0 ;
            while isnan(hpbwY)
                indhpbw = (Gn >= (0.50 - del)) & (Gn <= (0.50 + del)) ;
                hpbwY = mean(mean(th(indhpbw))) ;
                del = del + 0.01 ;
            end

            hpbw2 = hpbwX + hpbwY ;  % Main beam

            hpbw2 = hpbw2 * 180 / pi ;
            
            
            ant_pat_struct.G = G;
            ant_pat_struct.g = g;
            ant_pat_struct.th = th;
            ant_pat_struct.ph = ph;
            
        end  
        
    end
    
    methods( Access = private )
        
        function [th, ph, gg] = GGpattern( obj, ths_deg, a, ant_pat_res_deg )
            
            
            ant_pat_th_range_deg = Constants.ANT_PAT_TH_RANGE_DEG;
            ant_pat_th_range_rad = deg2rad( ant_pat_th_range_deg );
            
            ant_pat_ph_range_deg = Constants.ANT_PAT_PH_RANGE_DEG;
            ant_pat_ph_range_rad = deg2rad( ant_pat_ph_range_deg );

            % Beamwidth
            % ths_deg = 12, 6, 3

            % Sidelobe levels
            % a = 0.35 (15 dB), 0.25 (20 dB), 0.12 (30 dB), 0.055 (40 db)

            % TO-DO: Future work: Test for bad values (e.g. non-integer result)?
            Nth = floor( ant_pat_th_range_deg / ant_pat_res_deg ) + 1;
            Nph = floor( ant_pat_ph_range_deg / ant_pat_res_deg ) + 1;

            th_rad = linspace(0, ant_pat_th_range_rad, Nth) ;
            ph_rad = linspace(0, ant_pat_ph_range_rad, Nph) ;

            [th, ph] =  meshgrid(th_rad, ph_rad);

            % Angle span (Beamwidth)
            % ths_deg = 12 ;                     
            ths_rad = deg2rad(ths_deg) ; 

            % Generalized Gaussian Pattern parameters
            alpha = 0.2 ; % 0.5 ;  % sidelobe width                 
            % a = 0.35 ;  % sidelobe level

            % Generelized Gaussian
            gg = abs(1 / (1 - a) * exp( -(tan(th) / tan(ths_rad)) .^ 2)...
                    - a / (1 - a) * exp(-(alpha * tan(th) / tan(ths_rad)) .^ 2)) ;

            % To eliminate the angles that are over 0.4pi degrees from main beam
            XX = min( min( gg( :, (th_rad > 0.4*pi) & (th_rad < 0.45*pi) ) ) ) ;
            gg(:, th_rad > pi/2) = XX ;

            gg(gg < XX) = XX ;
        
        end
        
    end
    
end

