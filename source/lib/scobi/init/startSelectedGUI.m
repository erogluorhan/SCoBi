function inputStruct = startSelectedGUI( simulator_id )

inputStruct = [];

% TO-DO: Check this!
% If the OS is not UNIX OR it is MAC and Matlab version below 7.14
if (~isunix || (ismac && verLessThan('matlab', '7.14')))
        
    inputStruct = gui_SCoBi(simulator_id);
    
end

end