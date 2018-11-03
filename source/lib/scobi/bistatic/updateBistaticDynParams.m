
function updateBistaticDynParams
% function updateBistaticDynParams 
%
%   Calls the functions that calculate the geometry and transmitter bistatic 
%   parameters and updates the bistatic dynamic parameters 
%   (BistaticDynParams) with those parameter values in each 
%   simulation iteration.
%   
%   - Calls transmitterGeometryManual function to calculate the transmitter 
%   geometry-related  transformation matrices.
%   - Calls bistaticGeometry function to calculate the bistatic geometry-
%   related transformation matrices, direction vectors, & rotation matrices.
%   - Updates BistaticDynParams class.  
%
%   See also bistaticGeometry, transmitterGeometryManuel, BistaticDynParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


%% TO-DO: Future work: Real transmitterGeometry
%% TRANSMITTER GEOMETRY
[Tgt, TgtI] = transmitterGeometryManuel();


%% BISTATIC GEOMETRY
[rd_m, idn, isn, osp, osn, Tgs, Tgr, TgrI, AntRotZ_Rx, AntRotY_Rx, ...
    AntRot_Rx, AntRotZ_Tx, AllPoints_m, AngT2R_rf, ...
    AngS2R_rf, AngT2S_sf] = bistaticGeometry();


% Initialize Bistatic Parameters
BistaticDynParams.getInstance.update(rd_m, idn, isn, osp, osn, Tgt, ...
    TgtI, Tgs, Tgr, TgrI, AntRotZ_Rx, AntRotY_Rx, AntRot_Rx, AntRotZ_Tx, ...
    AllPoints_m, AngT2R_rf, AngS2R_rf, AngT2S_sf);

end