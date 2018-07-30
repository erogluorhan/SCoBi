function plotDielProfiles( fig1, eps_diel_soil, eps_diel_z2, eps_diel_z3, eps_diel_zL, eps_diel_zS, z )

% 2nd Order Polyfit
plotDielProfile2ndOrder( fig1, eps_diel_z2, eps_diel_soil, z);

% 3rd Order Polyfit
plotDielProfile3rdOrder( fig1, eps_diel_z3, eps_diel_soil, z);

% Logistic
plotDielProfileLogisticFunc( fig1, eps_diel_zL, eps_diel_soil, z );

% Discrete Slab
plotDielProfileDiscreteSlab( fig1, eps_diel_zS, eps_diel_soil, z );

end

