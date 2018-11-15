
function updateGndDynParams
% function updateGndDynParams 
%
%   Calculates the effective roughness parameter 
%   and calls the ground dielectric functions (Dobson, Mironow, Wang) to
%   calculate the dielectric constant, and updates GndDynParams class with 
%   those values in each simulation iteration.  
%
%   See also GndDynParams, dielDobson, dielMironov, dielWang.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0



%% GET GLOBAL PARAMETERS
sim_counter = ParamsManager.sim_counter;
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
f_Hz = f_MHz * Constants.MHZ_TO_HZ ;
% Configuration Parameters
VSM_list_cm3cm3 = ConfigParams.getInstance.VSM_list_cm3cm3;
VSM_cm3cm3 = VSM_list_cm3cm3( sim_counter, : )';
RMSH_list_cm = ConfigParams.getInstance.RMSH_list_cm;
RMSH_cm = RMSH_list_cm( sim_counter );
% Ground Parameters
% TO-DO: Future work: Implement for different ground layers
diel_model_id = GndParams.getInstance.diel_model_id;
sand_ratio = GndParams.getInstance.sand_ratio;
clay_ratio = GndParams.getInstance.clay_ratio;
rhob_gcm3 = GndParams.getInstance.rhob_gcm3;


%% CALCULATIONS
% Surface Roughness
lambda_cm = Constants.LIGHTSPEED / f_Hz * Constants.M_TO_CM ;        % in cm
ko = 2 * pi / lambda_cm ;
% effective roughness parameter
h = (2 * RMSH_cm * ko) ^ 2 ;

% Soil Dielectric Constant
if diel_model_id == Constants.ID_DIEL_DOBSON
    
    eps_g = dielDobson( f_Hz, VSM_cm3cm3, sand_ratio, clay_ratio, rhob_gcm3) ;
    
elseif diel_model_id == Constants.ID_DIEL_MIRONOV
    
    eps_g = dielMironov(f_Hz, VSM_cm3cm3, clay_ratio);
    
elseif diel_model_id == Constants.ID_DIEL_WANG
        
    eps_g = dielWang(VSM_cm3cm3, sand_ratio, clay_ratio, rhob_gcm3);
    
end

% Round the dielectric constant
eps_g = round(eps_g * 10) / 10 ;


% Initialize Ground Dynamic Parameters
GndDynParams.getInstance.update(h, eps_g);

end