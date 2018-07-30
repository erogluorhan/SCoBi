
function mainSCoBiML


%% GET GLOBAL DIRECTORIES
dir_out_ml_ref = SimulationFolders.getInstance.out_ml_ref;


%% GET GLOBAL PARAMETERS
% Dynamic Parameters
DOYs = DynParams.getInstance.DOYs;
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;    

%% PLOT
[fig1, fig2] = plotSMdata();


%% MULTILAYER CALCULATIONS

% Initialize reflectivities for four different profiles
DoYlist = [];
Rp1_list = []; 
Rp2_list = [];
Rp1_L_list = [];
Rp2_L_list = [];
Rp1_2nd_list = [];
Rp2_2nd_list = [];
Rp1_3rd_list = [];
Rp2_3rd_list = [];

Ninterval = length(DOYs);

for ii = 1 : 10 : Ninterval
    
    ParamsManager.index_VSM( ii );

    [Rp1, Rp2, Rp1_L, Rp2_L, Rp1_2nd, Rp2_2nd, Rp1_3rd, Rp2_3rd] = multiLayerModel( fig1, fig2, VSM_list_cm3cm3(ii,:)', DOYs(ii) ) ;

    DoYlist = [DoYlist; DOYs(ii)];
    Rp1_list = [Rp1_list; Rp1]; 
    Rp2_list = [Rp2_list; Rp2];
    Rp1_L_list = [Rp1_L_list; Rp1_L];
    Rp2_L_list = [Rp2_L_list; Rp2_L];
    Rp1_2nd_list = [Rp1_2nd_list; Rp1_2nd];
    Rp2_2nd_list = [Rp2_2nd_list; Rp2_2nd];
    Rp1_3rd_list = [Rp1_3rd_list; Rp1_3rd];
    Rp2_3rd_list = [Rp2_3rd_list; Rp2_3rd];

    [M1(ii), M2(ii)] = plotAddSMpoint( fig1, fig2, DOYs(ii), VSM_list_cm3cm3(ii,1) );

end

% plotMovie( M1 );

close(fig1);
close(fig2);

%% SAVE
Rp1s = [Rp1_list, Rp1_L_list, Rp1_2nd_list, Rp1_3rd_list];
Rp2s = [Rp2_list, Rp2_L_list, Rp2_2nd_list, Rp2_3rd_list];

filename01 = strcat('ML_Ref01') ;
filename02 = strcat('ML_Ref02') ;
filenameDoY = strcat('DoY') ;

writeVar(dir_out_ml_ref, filename01, Rp1s);
writeVar(dir_out_ml_ref, filename02, Rp2s);
writeVar(dir_out_ml_ref, filenameDoY, DoYlist);

end
