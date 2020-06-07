function inputStruct = initBackwardsControl( inputStruct )
% function initBackwardsControl
%
%   Checks for version-specific additions to the SCoBi inputStruct in order
%   to maintain backwards compatibility. If an inputStruct field from a
%   newer version of SCoBi is not present, it is initialized with a null
%   value.
%
%   initBackwardsControl( inputStruct )
%
%   INPUTS:
%   inputStruct: Input structure that comes from GUI
%
%   See also initAllInputParams, initTxParams, initRxParams, initGndParams,
%   initConfigParams, initVegParams.

%   Copyright © 2017-2020 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.3

%% Version 1.03 Additions
if ~isfield(inputStruct, 'calculate_penetration_depth')
    inputStruct.calculate_penetration_depth = 0;
end

end