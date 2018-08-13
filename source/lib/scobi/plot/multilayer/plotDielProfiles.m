function plotDielProfiles( fig1, eps_diel_z2nd, eps_diel_z3rd, eps_diel_zL, eps_diel_zS, z )

% 2nd Order Polyfit
plotDielProfile2ndOrder( fig1, eps_diel_z2nd, z);

% 3rd Order Polyfit
plotDielProfile3rdOrder( fig1, eps_diel_z3rd, z);

% Logistic
plotDielProfileLogisticFunc( fig1, eps_diel_zL, z );

% Discrete Slab
plotDielProfileDiscreteSlab( fig1, eps_diel_zS, z );

end

