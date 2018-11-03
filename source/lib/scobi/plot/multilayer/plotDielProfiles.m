function plotDielProfiles( fig1, selpath )


%% GET GLOBAL PARAMETER
% Ground ML Parameters
gnd_layer_depth_m = GndMLParams.getInstance.layer_depth_m;
calc_diel_profile_fit_functions = GndMLParams.getInstance.calc_diel_profile_fit_functions;
zA_m = GndMLParams.getInstance.zA_m;
z_m = GndMLParams.getInstance.z_m;    % Layer profile
% Ground Dynamic Params
eps_g = GndDynParams.getInstance.eps_g;

eps_g = ones (size(gnd_layer_depth_m));

eps_diel_profiles = [];

% 2nd Order Polyfit
eps_diel_z2nd = zeros;
if calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_2ND_ORDER, 1)
   
    % Generate the reflectivity folder
    dir_products_specular_reflectivity = strcat(selpath, '\products\specular\reflectivity\2nd-order' );

    % Read 2nd-order dielectric profile
    eps_diel_z2nd = ones (701, 1);
    
    eps_diel_profiles(:,1) = eps_diel_z2nd;
    
end

% 3rd Order Polyfit
eps_diel_z3rd = [];
if calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_3RD_ORDER, 1)
   
    % Generate the reflectivity folder
    dir_products_specular_reflectivity = strcat(selpath, '\products\specular\reflectivity\3rd-order' );

    % Read 3rd-order dielectric profile
    eps_diel_z3rd = 2 * ones (701, 1);
    
    eps_diel_profiles(:,2) = eps_diel_z3rd;

end

% Logistic
eps_diel_zL = [];
if calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_LOGISTIC, 1)
   
    % Generate the reflectivity folder
    dir_products_specular_reflectivity = strcat(selpath, '\products\specular\reflectivity\Logistic' );

    % Read Logistic-regression dielectric profile
    eps_diel_zL = 3 * ones (701, 1);
    
    eps_diel_profiles(:,3) = eps_diel_zL;
    
end

% Discrete Slab
eps_diel_zS = [];
if calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_SLAB, 1)
   
    % Generate the reflectivity folder
    dir_products_specular_reflectivity = strcat(selpath, '\products\specular\reflectivity\Discrete-slab' );

    % Read Discrete-slab dielectric profile
    eps_diel_zS = 4 * ones (701, 1);
    
    eps_diel_profiles(:,4) = eps_diel_zS;
    
end


%% PLOT
FigON = 1 ;
if FigON == 1
    
    figure(fig1)
    
    for ii = 1 : 4
    
        if ~ isempty( eps_diel_profiles(:, ii) )

            subplot(3,4,ii)
            plot( real(eps_diel_profiles(:, ii) ), (z_m - zA_m) * Constants.M_TO_CM, 'k', 'linewidth', 2)
            hold on
            plot(real(eps_g), gnd_layer_depth_m * Constants.M_TO_CM, 'o')
            title( Constants.DIEL_PROFILES{1,ii})

        end
    
%     subplot(3,4,2)
%     if ~ isempty( eps_diel_z3rd )
%         plot(real(eps_diel_z3rd), (z_m - zA_m) * Constants.M_TO_CM, 'k', 'linewidth', 2)
%         hold on
%         plot(real(eps_g), gnd_layer_depth_m * Constants.M_TO_CM, 'o')
%     end
%     title('3rd order polynomial')
%     
%     subplot(3,4,3)
%     if ~ isempty( eps_diel_zL )
%         plot(real(eps_diel_zL), (z_m - zA_m) * Constants.M_TO_CM, 'k', 'linewidth', 2)
%         hold on
%         plot(real(eps_g), gnd_layer_depth_m * Constants.M_TO_CM, 'o')
%     end
%     title('Logistic regression')
%     
%     subplot(3,4,4)
%     if ~ isempty( eps_diel_zS )
%         plot(real(eps_diel_zS), (z_m - zA_m) * Constants.M_TO_CM, 'k', 'linewidth', 2)
%         hold on
%         plot(real(eps_g), gnd_layer_depth_m * Constants.M_TO_CM, 'o')
%     end
%     title('Discrete-slab')    
    
        grid
        axis([0 30 0 z_m(end)])
        set(gca,'YDir','reverse')
        set(gca,'YLim',[0 - zA_m * Constants.M_TO_CM, z_m(end, :) * Constants.M_TO_CM - zA_m * Constants.M_TO_CM])
        % set(gca,'YTick',[0 zz z_m(end)] * Constants.M_TO_CM)
        % aa = sort([zz, z_m(end), zA_m + gnd_layer_depth_m]) - zA_m ;
        aa = [0; gnd_layer_depth_m] ;
        set(gca, 'YTick', aa * Constants.M_TO_CM)
        bb = real([Constants.EPS_DIEL_AIR; eps_g]) ;
        bb = sort(unique(bb)) ;
        set(gca,'XTick', bb)
        xlabel('\epsilon\prime - real part')
        ylabel('z [cm]')
        hold off
    
    end
    
end

FigON = 0 ;
if FigON == 1
    
    figure(10)
    if ~ isempty( eps_diel_z2nd )
        subplot(2,2,1)
        plot(real(eps_diel_z2nd), (z_m - zA_m) * Constants.M_TO_CM, 'k', 'linewidth', 2)
        title('2nd order polynomial')
    end
    hold on
    if ~ isempty( eps_diel_z3rd )
        subplot(2,2,2)
        plot(real(eps_diel_z3rd), (z_m - zA_m) * Constants.M_TO_CM, 'k', 'linewidth', 2)
        title('3rd order polynomial')
    end
    if ~ isempty( eps_diel_zL )
        subplot(2,2,3)
        plot(real(eps_diel_zL), (z_m - zA_m) * Constants.M_TO_CM, 'k', 'linewidth', 2)
        title('Logistic regression')
    end
    if ~ isempty( eps_diel_zS )
        subplot(2,2,4)
        plot(real(eps_diel_zS), (z_m - zA_m) * Constants.M_TO_CM, 'k', 'linewidth', 2)
        title('Discrete-slab')
    end
    
    plot(real(eps_g), gnd_layer_depth_m * Constants.M_TO_CM, 'o')
    
    grid
    axis([0 30 0 z_m(end)])
    set(gca,'YDir','reverse')
    set(gca,'YLim',[0 - zA_m * Constants.M_TO_CM, z_m(end, :) * Constants.M_TO_CM - zA_m * Constants.M_TO_CM])
    % set(gca,'YTick',[0 zz z_m(end)] * 1e2)
    % aa = sort([zz, z_m(end), zA_m + gnd_layer_depth_m]) - zA_m ;
    aa = [0, gnd_layer_depth_m] ;
    set(gca,'YTick',aa * Constants.M_TO_CM, 'FontSize', 6)
    bb = real([Constants.EPS_DIEL_AIR, eps_g]) ;
    bb = sort(unique(bb)) ;
    set(gca,'XTick',bb)
%     xlabel('\epsilon\prime - real part')
    ylabel('z [cm]')
    hold off
    
end


end

