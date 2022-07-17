
clearvars
clc
close all

%% INPUT

% Specify the set #
%Set_list = [1, 2, 3, 4, 5, 6];
Set_list = [1, 2, 3, 4];

% Scenario_list = [9, 10, 11, 13, 14, 15, 21, 22, 23, 25, 26];
% Kw_over_Ks_list = [2, 5, 10, 2, 5, 10, 2, 5, 10, 2, 5];
% Dip_list = [90, 90, 90, 90, 90, 90, -30, -30, -30, -30, -30];
Scenario_list = [21, 22, 23];
Kw_over_Ks_list = [2, 5, 10];
Dip_list = [-30, -30, -30];

Advection_rates_mpyr = [2.5E-04, 5E-04];

U_rate_mpyr = 1E-04;

% I'm using this to represent the horizontal motion of the contact due to
% uplift and the contact's dip (not including advection on top of that)
G = 0;%U_rate_mpyr / tan(deg2rad(30));

Advection_rates_mpyr = Advection_rates_mpyr + G;

Threshold_y_DD_vals = [1.02, 1.02];

plot_regression_option = 1;

use_stan_dev_option = 0;

Dip_val_title = -30;

Kw_over_Ks_vals = [2, 5, 10];

Kw_over_Ks_list_Extra_Scenarios = [2, 5, 10, 2, 5, 10];

%Dip_list = [90, 90, 90, -30, -30, -30];
Dip_vals = [90, -30];

K_colors = cbrewer('div', 'RdYlBu', 100, 'pchip');
K_colors = flipud(K_colors);

K_color_mid = [0.5, 0.5, 0.5];

markertype_list = {'o', '^', 's'};
markersize_list = [12, 12, 12];
linewidth_list = [1.5, 2, 1.5];
linestyle_list = {'-', ':'};

strong_unit_gone_marker = '*';
strong_unit_gone_color = 'k';
Sunit_gone_markersize_ref = 8;
Sunit_gone_linewidth_ref = 1;

Use_error_bars = 0;

saving_interval_ref_yrs = 5E+05;

if exist('Compiled_DIVIDEobj_Analysis_Sets1to4_Dip30_May29_95p','dir') == 0
    
    mkdir('Compiled_DIVIDEobj_Analysis_Sets1to4_Dip30_May29_95p');
    
end

Compiled_DIVIDEobj_Analysis_location = [pwd '/Compiled_DIVIDEobj_Analysis_Sets1to4_Dip30_May29_95p'];

%% COMPILE RESULTS

LowDiff_y_vals_fig = figure(1);
HighDiff_y_vals_fig = figure(2);

LowandHighDiff_dy_dt_vals_fig = figure(9);
HighDiff_dy_dt_vals_fig = figure(10);

LowandHighDiff_lag_time_fig = figure(11);
HighDiff_lag_time_fig = figure(12);

LowandHighDiff_y_vals_end_fig = figure(3);
HighDiff_y_vals_end_fig = figure(4);

LowandHighDiff_y_vals_max_fig = figure(13);
HighDiff_y_vals_max_fig = figure(14);

LowDiff_ksn_S_vals_fig = figure(5);
HighDiff_ksn_S_vals_fig = figure(6);

LowDiff_Relief_S_vals_fig = figure(7);
HighDiff_Relief_S_vals_fig = figure(8);

Original_location_ref = pwd;

for Set_ref = 1:1:length(Set_list)
    
    for Scenario_ref = 1:1:length(Scenario_list)
        
        Set = Set_list(1,Set_ref);
        Scenario = Scenario_list(1,Scenario_ref);
        
        %
        proceed_test = 1;
        
        %
        if proceed_test == 1
            
            if Kw_over_Ks_list(1,Scenario_ref) == Kw_over_Ks_vals(1,1)
                
                color_ref = K_colors(1,:);
                Threshold_y_DD_val = Threshold_y_DD_vals(1,1);
                
            elseif Kw_over_Ks_list(1,Scenario_ref) == Kw_over_Ks_vals(1,2)
                
                %color_ref = K_colors(ceil(length(K_colors) / 2),:);
                color_ref = K_color_mid;
                Threshold_y_DD_val = Threshold_y_DD_vals(1,2);
                
            elseif Kw_over_Ks_list(1,Scenario_ref) == Kw_over_Ks_vals(1,3)
                
                color_ref = K_colors(end,:);
                Threshold_y_DD_val = Threshold_y_DD_vals(1,2);
                
            end
            
            if Set == 1 || Set == 2
                
                marker_ref = markertype_list(1,1);
                markersize_ref = markersize_list(1,1);
                linewidth_ref = linewidth_list(1,1);
                linestyle_ref = linestyle_list(1,1);
                Advection_rate_m_per_yr_ref = Advection_rates_mpyr(1,2);
                
            elseif Set == 3 || Set == 4
                
                marker_ref = markertype_list(1,2);
                markersize_ref = markersize_list(1,2);
                linewidth_ref = linewidth_list(1,2);
                linestyle_ref = linestyle_list(1,2);
                Advection_rate_m_per_yr_ref = Advection_rates_mpyr(1,1);
                
            end
            
            %
            
            cd([Original_location_ref '/Set' num2str(Set,'%.0f') '\Output\Scenario_' num2str(Scenario,'%.0f') '/DIVIDEobj_Analysis'])
            
            clear ksn_selected_excprob_in_strong_N_over_time
            load('ksn_in_strong_Workspace_May25.mat')
            
            % Get time at which strong unit is gone
            time_strong_unit_gone = max(time_range(isnan(ksn_selected_excprob_in_strong_N_over_time) ~= 1)) + saving_interval_ref_yrs;
            
            clear mean_y_dist_km_DrnDiv_OrderAboveThresh std_y_dist_km_DrnDiv_OrderAboveThresh
            load('Drainage_Divide_Workspace.mat')
            
            close(figure(15))
            
            % Get dystar_dt vals
            
            Ystar_over_thresh = mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1);
            Y_over_thresh = mean_y_dist_km_DrnDiv_OrderAboveThresh;
            time_range_Ystar_over_thresh = time_range;
            %time_range_Ystar_over_thresh_copy = time_range_Ystar_over_thresh;
            
            time_max_Ystar = time_range_Ystar_over_thresh(Ystar_over_thresh == max(Ystar_over_thresh));
            
            Ystar_over_thresh = Ystar_over_thresh(time_range_Ystar_over_thresh <= time_max_Ystar);
            Y_over_thresh = Y_over_thresh(time_range_Ystar_over_thresh <= time_max_Ystar);
            time_range_Ystar_over_thresh = time_range_Ystar_over_thresh(time_range_Ystar_over_thresh <= time_max_Ystar);
            
            time_range_Ystar_over_thresh = time_range_Ystar_over_thresh(Ystar_over_thresh >= Threshold_y_DD_val);
            Y_over_thresh = Y_over_thresh(Ystar_over_thresh >= Threshold_y_DD_val);
            Ystar_over_thresh = Ystar_over_thresh(Ystar_over_thresh >= Threshold_y_DD_val);
            
            min_time_Ystar_over_thresh = min(time_range_Ystar_over_thresh);
            min_Y_over_thresh = min(Y_over_thresh);
            min_Ystar_over_thresh = min(Ystar_over_thresh);
            
            dYstar_dt_fitlm = fitlm(time_range_Ystar_over_thresh - min_time_Ystar_over_thresh, Ystar_over_thresh - min_Ystar_over_thresh, ...
                'intercept', false);
            Coeffs_temp = table2array(dYstar_dt_fitlm.Coefficients);
            Slope_temp = Coeffs_temp(1,1);
            R2_temp = dYstar_dt_fitlm.Rsquared.Adjusted;
            
            y_val_temp = mean_y_dist_km_DrnDiv_OrderAboveThresh(time_range == time_strong_unit_gone) ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1);
            
            %  DIVIDE Y VALS FIG
            if mod(Set,2) == 1
                
                figure(LowDiff_y_vals_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_mean_y_vals_LowDiff_HighAdvec = plot(time_range ./ 1E+06, mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    elseif Set == 3 || Set  == 4
                        
                        h_mean_y_vals_LowDiff_LowAdvec = plot(time_range ./ 1E+06, mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    if Set == 1 || Set == 2
                        
                        h_mean_y_vals_LowDiff_HighAdvec = errorbar(time_range ./ 1E+06, mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            std_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    elseif Set == 3 || Set  == 4
                        
                        h_mean_y_vals_LowDiff_LowAdvec = errorbar(time_range ./ 1E+06, mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            std_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    end
                    
                end
                
                h_S_gone_yvals_LD = plot(time_strong_unit_gone ./ 1E+06, y_val_temp, 'marker', char(strong_unit_gone_marker), ...
                    'color', char(strong_unit_gone_color), 'markersize', Sunit_gone_markersize_ref, 'linestyle', 'none', 'linewidth', Sunit_gone_linewidth_ref);
                hold on
                
                if plot_regression_option == 1
                    
                    plot([min(time_range_Ystar_over_thresh), max(time_range_Ystar_over_thresh)] ./ 1E+06, ...
                        ([min(time_range_Ystar_over_thresh) - min_time_Ystar_over_thresh, max(time_range_Ystar_over_thresh) - min_time_Ystar_over_thresh] .* Slope_temp) + min_Ystar_over_thresh, 'marker', 'none', ...
                        'color', char(strong_unit_gone_color), 'linestyle', char(linestyle_ref), 'linewidth', 3);
                    hold on
                    
                end
                
            elseif mod(Set,2) == 0
                
                figure(HighDiff_y_vals_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_mean_y_vals_HighDiff_HighAdvec = plot(time_range ./ 1E+06, mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_mean_y_vals_HighDiff_LowAdvec = plot(time_range ./ 1E+06, mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    if Set == 1 || Set == 2
                        
                        h_mean_y_vals_HighDiff_HighAdvec = errorbar(time_range ./ 1E+06, mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            std_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_mean_y_vals_HighDiff_LowAdvec = errorbar(time_range ./ 1E+06, mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            std_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    end
                    
                end
                
                h_S_gone_yvals_HD = plot(time_strong_unit_gone ./ 1E+06, y_val_temp, 'marker', char(strong_unit_gone_marker), ...
                    'color', char(strong_unit_gone_color), 'markersize', Sunit_gone_markersize_ref, 'linestyle', 'none', 'linewidth', Sunit_gone_linewidth_ref);
                hold on
                
                if plot_regression_option == 1
                    
                    plot([min(time_range_Ystar_over_thresh), max(time_range_Ystar_over_thresh)] ./ 1E+06, ...
                        ([min(time_range_Ystar_over_thresh) - min_time_Ystar_over_thresh, max(time_range_Ystar_over_thresh) - min_time_Ystar_over_thresh] .* Slope_temp) + min_Ystar_over_thresh, 'marker', 'none', ...
                        'color', char(strong_unit_gone_color), 'linestyle', char(linestyle_ref), 'linewidth', 3);
                    hold on
                    
                end
                
            end
            
            if use_stan_dev_option == 1
                
                x_vals_temp = [time_range, fliplr(time_range)] ./ 1E+06;
                y_vals_temp = [(mean_y_dist_km_DrnDiv_OrderAboveThresh + std_y_dist_km_DrnDiv_OrderAboveThresh) ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1); ...
                    flipud((mean_y_dist_km_DrnDiv_OrderAboveThresh - std_y_dist_km_DrnDiv_OrderAboveThresh) ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1))];
                y_vals_temp = y_vals_temp';
                
                pgon = polyshape(x_vals_temp, y_vals_temp, 'Simplify', false);
                
                plot(pgon, 'FaceColor', color_ref, 'EdgeColor', color_ref, ...
                    'FaceAlpha', 0.1, 'EdgeAlpha', 0.1)
                
            end
            
            
            
            % DY/DT VALS FIG
            if mod(Set,2) == 1
                
                figure(LowandHighDiff_dy_dt_vals_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_mean_dy_dt_vals_LowDiff_HighAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), Slope_temp, ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref .* 0.5, 'markerfacecolor', color_ref, 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_mean_dy_dt_vals_LowDiff_LowAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), Slope_temp, ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref .* 0.5, 'markerfacecolor', color_ref, 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    
                    
                end
                
                
            elseif mod(Set,2) == 0
                
                figure(LowandHighDiff_dy_dt_vals_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_mean_dy_dt_vals_HighDiff_HighAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), Slope_temp, ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref, 'markerfacecolor', 'none', 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_mean_dy_dt_vals_HighDiff_LowAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), Slope_temp, ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref, 'markerfacecolor', 'none', 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    
                    
                end
                
            end
            
            
            
            % LAG TIME FIG
            if mod(Set,2) == 1
                
                figure(LowandHighDiff_lag_time_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_mean_lag_time_LowDiff_HighAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), min_time_Ystar_over_thresh / 1E+06, ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref .* 0.5, 'markerfacecolor', color_ref, 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_mean_lag_time_LowDiff_LowAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), min_time_Ystar_over_thresh / 1E+06, ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref .* 0.5, 'markerfacecolor', color_ref, 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    
                    
                end
                
                
            elseif mod(Set,2) == 0
                
                figure(LowandHighDiff_lag_time_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_mean_lag_time_HighDiff_HighAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), min_time_Ystar_over_thresh / 1E+06, ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref, 'markerfacecolor', 'none', 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_mean_lag_time_HighDiff_LowAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), min_time_Ystar_over_thresh / 1E+06, ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref, 'markerfacecolor', 'none', 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    
                    
                end
                
            end
            
            
            
            % FINAL DIVIDE Y VALS FIG
            if mod(Set,2) == 1
                
                figure(LowandHighDiff_y_vals_end_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_y_vals_end_LowDiff_HighAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), mean_y_dist_km_DrnDiv_OrderAboveThresh(end,1) ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref .* 0.5, 'markerfacecolor', color_ref, 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_y_vals_end_LowDiff_LowAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), mean_y_dist_km_DrnDiv_OrderAboveThresh(end,1) ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref .* 0.5, 'markerfacecolor', color_ref, 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    
                    
                end
                
            elseif mod(Set,2) == 0
                
                figure(LowandHighDiff_y_vals_end_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_y_vals_end_HighDiff_HighAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), mean_y_dist_km_DrnDiv_OrderAboveThresh(end,1) ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref, 'markerfacecolor', 'none', 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_y_vals_end_HighDiff_LowAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), mean_y_dist_km_DrnDiv_OrderAboveThresh(end,1) ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1), ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref, 'markerfacecolor', 'none', 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    
                    
                end
                
            end
            
            
            
            
            % MAX DIVIDE Y VALS FIG
            if mod(Set,2) == 1
                
                figure(LowandHighDiff_y_vals_max_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_y_vals_max_LowDiff_HighAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), max(mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1)), ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref .* 0.5, 'markerfacecolor', color_ref, 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_y_vals_max_LowDiff_LowAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), max(mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1)), ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref .* 0.5, 'markerfacecolor', color_ref, 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    
                    
                end
                
            elseif mod(Set,2) == 0
                
                figure(LowandHighDiff_y_vals_max_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_y_vals_max_HighDiff_HighAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), max(mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1)), ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref, 'markerfacecolor', 'none', 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_y_vals_max_HighDiff_LowAdvec = plot(Kw_over_Ks_list(1,Scenario_ref), max(mean_y_dist_km_DrnDiv_OrderAboveThresh ./ mean_y_dist_km_DrnDiv_OrderAboveThresh(1,1)), ...
                            'marker', char(marker_ref), 'linewidth', 1.5, 'color', color_ref, 'markerfacecolor', 'none', 'linestyle', 'none', 'markersize', markersize_ref);
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    
                    
                end
                
            end
            
            
            
            clear ksn_selected_excprob_in_strong_S_over_time ksn_selected_excprob_in_strong_N_over_time ...
                std_ksn_in_strong_S_over_time std_ksn_in_strong_N_over_time
            load('ksn_in_strong_Workspace_May25.mat')
            
            %y_val_temp = ksn_selected_excprob_in_strong_N_over_time(time_range == time_strong_unit_gone) ./ ksn_selected_excprob_in_strong_N_over_time(1,1);
            
            % ksn in strong fig
            if mod(Set,2) == 1
                
                figure(LowDiff_ksn_S_vals_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_ksn_S_LowDiff_HighAdvec = plot(time_range ./ 1E+06, ksn_selected_excprob_in_strong_S_over_time ./ ksn_selected_excprob_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref * 2, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_ksn_N_LowDiff_HighAdvec = plot(time_range ./ 1E+06, ksn_selected_excprob_in_strong_N_over_time ./ ksn_selected_excprob_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_ksn_S_LowDiff_LowAdvec = plot(time_range ./ 1E+06, ksn_selected_excprob_in_strong_S_over_time ./ ksn_selected_excprob_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref * 2, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_ksn_N_LowDiff_LowAdvec = plot(time_range ./ 1E+06, ksn_selected_excprob_in_strong_N_over_time ./ ksn_selected_excprob_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    if Set == 1 || Set == 2
                        
                        h_ksn_S_LowDiff_HighAdvec = errorbar(time_range ./ 1E+06, ksn_selected_excprob_in_strong_S_over_time ./ ksn_selected_excprob_in_strong_S_over_time(1,1), ...
                            std_ksn_in_strong_S_over_time ./ ksn_selected_excprob_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 3, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_ksn_N_LowDiff_HighAdvec = errorbar(time_range ./ 1E+06, ksn_selected_excprob_in_strong_N_over_time ./ ksn_selected_excprob_in_strong_N_over_time(1,1), ...
                            std_ksn_in_strong_N_over_time ./ ksn_selected_excprob_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_ksn_S_LowDiff_LowAdvec = errorbar(time_range ./ 1E+06, ksn_selected_excprob_in_strong_S_over_time ./ ksn_selected_excprob_in_strong_S_over_time(1,1), ...
                            std_ksn_in_strong_S_over_time ./ ksn_selected_excprob_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 3, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_ksn_N_LowDiff_LowAdvec = errorbar(time_range ./ 1E+06, ksn_selected_excprob_in_strong_N_over_time ./ ksn_selected_excprob_in_strong_N_over_time(1,1), ...
                            std_ksn_in_strong_N_over_time ./ ksn_selected_excprob_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    end
                    
                end
                
%                 h_S_gone_ksn_LD = plot(time_strong_unit_gone ./ 1E+06, y_val_temp, 'marker', char(strong_unit_gone_marker), ...
%                     'color', char(strong_unit_gone_color), 'markersize', Sunit_gone_markersize_ref, 'linestyle', 'none', 'linewidth', Sunit_gone_linewidth_ref);
%                 hold on
                
            elseif mod(Set,2) == 0
                
                figure(HighDiff_ksn_S_vals_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_ksn_S_HighDiff_HighAdvec = plot(time_range ./ 1E+06, ksn_selected_excprob_in_strong_S_over_time ./ ksn_selected_excprob_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref * 2, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_ksn_N_HighDiff_HighAdvec = plot(time_range ./ 1E+06, ksn_selected_excprob_in_strong_N_over_time ./ ksn_selected_excprob_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_ksn_S_HighDiff_LowAdvec = plot(time_range ./ 1E+06, ksn_selected_excprob_in_strong_S_over_time ./ ksn_selected_excprob_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref * 2, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_ksn_N_HighDiff_LowAdvec = plot(time_range ./ 1E+06, ksn_selected_excprob_in_strong_N_over_time ./ ksn_selected_excprob_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    if Set == 1 || Set == 2
                        
                        h_ksn_S_HighDiff_HighAdvec = errorbar(time_range ./ 1E+06, ksn_selected_excprob_in_strong_S_over_time ./ ksn_selected_excprob_in_strong_S_over_time(1,1), ...
                            std_ksn_in_strong_S_over_time ./ ksn_selected_excprob_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 3, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_ksn_N_HighDiff_HighAdvec = errorbar(time_range ./ 1E+06, ksn_selected_excprob_in_strong_N_over_time ./ ksn_selected_excprob_in_strong_N_over_time(1,1), ...
                            std_ksn_in_strong_N_over_time ./ ksn_selected_excprob_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_ksn_S_HighDiff_LowAdvec = errorbar(time_range ./ 1E+06, ksn_selected_excprob_in_strong_S_over_time ./ ksn_selected_excprob_in_strong_S_over_time(1,1), ...
                            std_ksn_in_strong_S_over_time ./ ksn_selected_excprob_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 3, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_ksn_N_HighDiff_LowAdvec = errorbar(time_range ./ 1E+06, ksn_selected_excprob_in_strong_N_over_time ./ ksn_selected_excprob_in_strong_N_over_time(1,1), ...
                            std_ksn_in_strong_N_over_time ./ ksn_selected_excprob_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    end
                    
                end
                
%                 h_S_gone_ksn_HD = plot(time_strong_unit_gone ./ 1E+06, y_val_temp, 'marker', char(strong_unit_gone_marker), ...
%                     'color', char(strong_unit_gone_color), 'markersize', Sunit_gone_markersize_ref, 'linestyle', 'none', 'linewidth', Sunit_gone_linewidth_ref);
%                 hold on
                
            end
            
            if use_stan_dev_option == 1
                
                x_vals_temp = [time_range(isnan(ksn_selected_excprob_in_strong_N_over_time) ~= 1), fliplr(time_range(isnan(ksn_selected_excprob_in_strong_N_over_time) ~= 1))] ./ 1E+06;
                y_vals_temp = [(ksn_selected_excprob_in_strong_N_over_time(isnan(ksn_selected_excprob_in_strong_N_over_time) ~= 1) + std_ksn_in_strong_N_over_time(isnan(ksn_selected_excprob_in_strong_N_over_time) ~= 1)) ./ ksn_selected_excprob_in_strong_N_over_time(1,1); ...
                    flipud((ksn_selected_excprob_in_strong_N_over_time(isnan(ksn_selected_excprob_in_strong_N_over_time) ~= 1) - std_ksn_in_strong_N_over_time(isnan(ksn_selected_excprob_in_strong_N_over_time) ~= 1)) ./ ksn_selected_excprob_in_strong_N_over_time(1,1))];
                y_vals_temp = y_vals_temp';
                
                pgon = polyshape(x_vals_temp, y_vals_temp, 'Simplify', false);
                
                plot(pgon, 'FaceColor', color_ref, 'EdgeColor', color_ref, ...
                    'FaceAlpha', 0.1, 'EdgeAlpha', 0.1)
                
                
                x_vals_temp = [time_range(isnan(ksn_selected_excprob_in_strong_S_over_time) ~= 1), fliplr(time_range(isnan(ksn_selected_excprob_in_strong_S_over_time) ~= 1))] ./ 1E+06;
                y_vals_temp = [(ksn_selected_excprob_in_strong_S_over_time(isnan(ksn_selected_excprob_in_strong_S_over_time) ~= 1) + std_ksn_in_strong_S_over_time(isnan(ksn_selected_excprob_in_strong_S_over_time) ~= 1)) ./ ksn_selected_excprob_in_strong_S_over_time(1,1); ...
                    flipud((ksn_selected_excprob_in_strong_S_over_time(isnan(ksn_selected_excprob_in_strong_S_over_time) ~= 1) - std_ksn_in_strong_S_over_time(isnan(ksn_selected_excprob_in_strong_S_over_time) ~= 1)) ./ ksn_selected_excprob_in_strong_S_over_time(1,1))];
                y_vals_temp = y_vals_temp';
                
                pgon = polyshape(x_vals_temp, y_vals_temp, 'Simplify', false);
                
                plot(pgon, 'FaceColor', color_ref, 'EdgeColor', color_ref, ...
                    'FaceAlpha', 0.1, 'EdgeAlpha', 0.1)
                
            end
            
            
            
            clear relief_selected_excprob_window1_in_strong_N_over_time relief_selected_excprob_window1_in_strong_S_over_time ...
                std_relief_window1_in_strong_N_over_time std_relief_window1_in_strong_S_over_time
            load('Relief_in_strong_Workspace_May25.mat')
            
            %y_val_temp = relief_selected_excprob_window1_in_strong_N_over_time(time_range == time_strong_unit_gone)./ relief_selected_excprob_window1_in_strong_N_over_time(1,1);
            
            % Relief in strong fig
            if mod(Set,2) == 1
                
                figure(LowDiff_Relief_S_vals_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_Relief_S_LowDiff_HighAdvec = plot(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_S_over_time ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref * 2, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_Relief_N_LowDiff_HighAdvec = plot(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_N_over_time ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_Relief_S_LowDiff_LowAdvec = plot(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_S_over_time ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref * 2, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_Relief_N_LowDiff_LowAdvec = plot(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_N_over_time ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    if Set == 1 || Set == 2
                        
                        h_Relief_S_LowDiff_HighAdvec = errorbar(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_S_over_time ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1), ...
                            std_relief_window1_in_strong_S_over_time ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 3, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_Relief_N_LowDiff_HighAdvec = errorbar(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_N_over_time ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1), ...
                            std_relief_window1_in_strong_N_over_time ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_Relief_S_LowDiff_LowAdvec = errorbar(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_S_over_time ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1), ...
                            std_relief_window1_in_strong_S_over_time ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 3, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_Relief_N_LowDiff_LowAdvec = errorbar(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_N_over_time ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1), ...
                            std_relief_window1_in_strong_N_over_time ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    end
                    
                end
                
%                 h_S_gone_Relief_LD = plot(time_strong_unit_gone ./ 1E+06, y_val_temp, 'marker', char(strong_unit_gone_marker), ...
%                     'color', char(strong_unit_gone_color), 'markersize', Sunit_gone_markersize_ref, 'linestyle', 'none', 'linewidth', Sunit_gone_linewidth_ref);
%                 hold on
                
            elseif mod(Set,2) == 0
                
                figure(HighDiff_Relief_S_vals_fig)
                
                if Use_error_bars ~= 1
                    
                    if Set == 1 || Set == 2
                        
                        h_Relief_S_HighDiff_HighAdvec = plot(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_S_over_time ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref * 2, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_Relief_N_HighDiff_HighAdvec = plot(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_N_over_time ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_Relief_S_HighDiff_LowAdvec = plot(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_S_over_time ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref * 2, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_Relief_N_HighDiff_LowAdvec = plot(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_N_over_time ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', linewidth_ref, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    end
                    
                elseif Use_error_bars == 1
                    
                    if Set == 1 || Set == 2
                        
                        h_Relief_S_HighDiff_HighAdvec = errorbar(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_S_over_time ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1), ...
                            std_relief_window1_in_strong_S_over_time ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 3, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_Relief_N_HighDiff_HighAdvec = errorbar(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_N_over_time ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1), ...
                            std_relief_window1_in_strong_N_over_time ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    elseif Set == 3 || Set == 4
                        
                        h_Relief_S_HighDiff_LowAdvec = errorbar(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_S_over_time ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1), ...
                            std_relief_window1_in_strong_S_over_time ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                        h_Relief_N_HighDiff_LowAdvec = errorbar(time_range ./ 1E+06, relief_selected_excprob_window1_in_strong_N_over_time ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1), ...
                            std_relief_window1_in_strong_N_over_time ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1), ...
                            'marker', 'none', 'linewidth', 1.5, 'color', color_ref, 'linestyle', char(linestyle_ref));
                        hold on
                        
                    end
                    
                end
                
%                 h_S_gone_Relief_HD = plot(time_strong_unit_gone ./ 1E+06, y_val_temp, 'marker', char(strong_unit_gone_marker), ...
%                     'color', char(strong_unit_gone_color), 'markersize', Sunit_gone_markersize_ref, 'linestyle', 'none', 'linewidth', Sunit_gone_linewidth_ref);
%                 hold on
                
            end
            
            if use_stan_dev_option == 1
                
                x_vals_temp = [time_range(isnan(relief_selected_excprob_window1_in_strong_N_over_time) ~= 1), fliplr(time_range(isnan(relief_selected_excprob_window1_in_strong_N_over_time) ~= 1))] ./ 1E+06;
                y_vals_temp = [(relief_selected_excprob_window1_in_strong_N_over_time(isnan(relief_selected_excprob_window1_in_strong_N_over_time) ~= 1) + std_relief_window1_in_strong_N_over_time(isnan(relief_selected_excprob_window1_in_strong_N_over_time) ~= 1)) ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1); ...
                    flipud((relief_selected_excprob_window1_in_strong_N_over_time(isnan(relief_selected_excprob_window1_in_strong_N_over_time) ~= 1) - std_relief_window1_in_strong_N_over_time(isnan(relief_selected_excprob_window1_in_strong_N_over_time) ~= 1)) ./ relief_selected_excprob_window1_in_strong_N_over_time(1,1))];
                y_vals_temp = y_vals_temp';
                
                pgon = polyshape(x_vals_temp, y_vals_temp, 'Simplify', false);
                
                plot(pgon, 'FaceColor', color_ref, 'EdgeColor', color_ref, ...
                    'FaceAlpha', 0.1, 'EdgeAlpha', 0.1)
                
                x_vals_temp = [time_range(isnan(relief_selected_excprob_window1_in_strong_S_over_time) ~= 1), fliplr(time_range(isnan(relief_selected_excprob_window1_in_strong_S_over_time) ~= 1))] ./ 1E+06;
                y_vals_temp = [(relief_selected_excprob_window1_in_strong_S_over_time(isnan(relief_selected_excprob_window1_in_strong_S_over_time) ~= 1) + std_relief_window1_in_strong_S_over_time(isnan(relief_selected_excprob_window1_in_strong_S_over_time) ~= 1)) ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1); ...
                    flipud((relief_selected_excprob_window1_in_strong_S_over_time(isnan(relief_selected_excprob_window1_in_strong_S_over_time) ~= 1) - std_relief_window1_in_strong_S_over_time(isnan(relief_selected_excprob_window1_in_strong_S_over_time) ~= 1)) ./ relief_selected_excprob_window1_in_strong_S_over_time(1,1))];
                y_vals_temp = y_vals_temp';
                
                pgon = polyshape(x_vals_temp, y_vals_temp, 'Simplify', false);
                
                plot(pgon, 'FaceColor', color_ref, 'EdgeColor', color_ref, ...
                    'FaceAlpha', 0.1, 'EdgeAlpha', 0.1)
                
            end
            
            
            
        end
        
    end
    
end

%%

cd(Compiled_DIVIDEobj_Analysis_location)

%%

for Iter = 1:2
    %%
    
    figure(LowDiff_y_vals_fig)
    
    if Iter == 1
        
        set(gca, 'fontsize', 8)
        
        ylim_ref = ylim;
        
        ylim_min_y_vals = ylim_ref(1,1);
        ylim_max_y_vals = ylim_ref(1,2);
        
        title(['Drainage Divide Position, Lower D, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
        xlabel('t (Myr)', 'fontsize', 10)
        ylabel('Y_{DD avg}*(t)', 'fontsize', 10)
        
    elseif Iter == 2
        
        ylim([ylim_min_y_vals, ylim_max_y_vals])
        
        %         h = colorbar;
        %         colormap(K_colors)
        %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
        %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
        
        legend([h_mean_y_vals_LowDiff_HighAdvec, h_mean_y_vals_LowDiff_LowAdvec, h_S_gone_yvals_LD], ...
            'High Advec.', 'Low Advec.', 'S Unit Gone', 'location', 'best')
        
        set(gcf, 'renderer', 'Painters')
        
        if Use_error_bars ~= 1
            
            saveas(LowDiff_y_vals_fig,'DD_Position_over_time_Fig_LowDiff.fig')
            saveas(LowDiff_y_vals_fig,'DD_Position_over_time_Fig_LowDiff.png')
            saveas(LowDiff_y_vals_fig,'DD_Position_over_time_Fig_LowDiff','epsc')
            
        elseif Use_error_bars == 1
            
            saveas(LowDiff_y_vals_fig,'DD_Position_over_time_Fig_LowDiff_errorbars.fig')
            saveas(LowDiff_y_vals_fig,'DD_Position_over_time_Fig_LowDiff_errorbars.png')
            saveas(LowDiff_y_vals_fig,'DD_Position_over_time_Fig_LowDiff_errorbars','epsc')
            
        end
        
    end
    
    %%
    
    figure(LowandHighDiff_dy_dt_vals_fig)
    
    if Iter == 1
        
        set(gca, 'fontsize', 8)
        
        ylim_ref = ylim;
        
        ylim_min_dy_dt_vals = ylim_ref(1,1);
        ylim_max_dy_dt_vals = ylim_ref(1,2);
        
        title(['Drainage Divide Migration, Higher D: Hollow, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
        xlabel('K_W / K_S', 'fontsize', 10)
        ylabel('(dY_{DD avg}*/dt) (yr^{-1})', 'fontsize', 10)
        
    elseif Iter == 2
        
        ylim([ylim_min_dy_dt_vals, ylim_max_dy_dt_vals])
        
        %         h = colorbar;
        %         colormap(K_colors)
        %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
        %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
        
        legend([h_mean_y_vals_LowDiff_HighAdvec, h_mean_y_vals_LowDiff_LowAdvec], ...
            'High Advec.', 'Low Advec.', 'S Unit Gone', 'location', 'best')
        
        set(gcf, 'renderer', 'Painters')
        
        if Use_error_bars ~= 1
            
            saveas(LowandHighDiff_dy_dt_vals_fig,'DD_dydt_Fig_LowDiff.fig')
            saveas(LowandHighDiff_dy_dt_vals_fig,'DD_dydt_Fig_LowDiff.png')
            saveas(LowandHighDiff_dy_dt_vals_fig,'DD_dydt_Fig_LowDiff','epsc')
            
        elseif Use_error_bars == 1
            
            saveas(LowandHighDiff_dy_dt_vals_fig,'DD_dydt_Fig_LowDiff_errorbars.fig')
            saveas(LowandHighDiff_dy_dt_vals_fig,'DD_dydt_Fig_LowDiff_errorbars.png')
            saveas(LowandHighDiff_dy_dt_vals_fig,'DD_dydt_Fig_LowDiff_errorbars','epsc')
            
        end
        
    end
    
    %%
    
    figure(LowandHighDiff_lag_time_fig)
    
    if Iter == 1
        
        set(gca, 'fontsize', 8)
        
        ylim_ref = ylim;
        
        ylim_min_lag_time_vals = ylim_ref(1,1);
        ylim_max_lag_time_vals = ylim_ref(1,2);
        
        title(['Drainage Divide Migration, Higher D: Hollow, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
        xlabel('K_W / K_W', 'fontsize', 10)
        ylabel('Lag Time (Myr)', 'fontsize', 10)
        
    elseif Iter == 2
        
        ylim([ylim_min_lag_time_vals, ylim_max_lag_time_vals])
        
        %         h = colorbar;
        %         colormap(K_colors)
        %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
        %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
        
        legend([h_mean_lag_time_LowDiff_HighAdvec, h_mean_lag_time_LowDiff_LowAdvec], ...
            'High Advec.', 'Low Advec.', 'S Unit Gone', 'location', 'best')
        
        set(gcf, 'renderer', 'Painters')
        
        if Use_error_bars ~= 1
            
            saveas(LowandHighDiff_lag_time_fig,'Lag_Time_Fig_LowDiff.fig')
            saveas(LowandHighDiff_lag_time_fig,'Lag_Time_Fig_LowDiff.png')
            saveas(LowandHighDiff_lag_time_fig,'Lag_Time_Fig_LowDiff','epsc')
            
        elseif Use_error_bars == 1
            
            saveas(LowandHighDiff_lag_time_fig,'Lag_Time_Fig_LowDiff_errorbars.fig')
            saveas(LowandHighDiff_lag_time_fig,'Lag_Time_Fig_LowDiff_errorbars.png')
            saveas(LowandHighDiff_lag_time_fig,'Lag_Time_Fig_LowDiff_errorbars','epsc')
            
        end
        
    end
    
    %%
    
    figure(LowandHighDiff_y_vals_end_fig)
    
    if Iter == 1
        
        set(gca, 'fontsize', 8)
        
        ylim_ref = ylim;
        
        ylim_min_y_vals_end = ylim_ref(1,1);
        ylim_max_y_vals_end = ylim_ref(1,2);
        
        title(['Drainage Divide Final Position, Higher D: Hollow, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
        xlabel('K_W / K_S', 'fontsize', 10)
        ylabel('Y_{DD avg}*(t_f)', 'fontsize', 10)
        
    elseif Iter == 2
        
        ylim([ylim_min_y_vals_end, ylim_max_y_vals_end])
        
        %         h = colorbar;
        %         colormap(K_colors)
        %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
        %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
        
        legend([h_y_vals_end_LowDiff_HighAdvec, h_y_vals_end_LowDiff_LowAdvec], ...
            'High Advec.', 'Low Advec.', 'location', 'best')
        
        set(gcf, 'renderer', 'Painters')
        
        if Use_error_bars ~= 1
            
            saveas(LowandHighDiff_y_vals_end_fig,'DD_Position_Final_Fig_LowDiff.fig')
            saveas(LowandHighDiff_y_vals_end_fig,'DD_Position_Final_Fig_LowDiff.png')
            saveas(LowandHighDiff_y_vals_end_fig,'DD_Position_Final_Fig_LowDiff','epsc')
            
        elseif Use_error_bars == 1
            
            saveas(LowandHighDiff_y_vals_end_fig,'DD_Position_Final_Fig_LowDiff_errorbars.fig')
            saveas(LowandHighDiff_y_vals_end_fig,'DD_Position_Final_Fig_LowDiff_errorbars.png')
            saveas(LowandHighDiff_y_vals_end_fig,'DD_Position_Final_Fig_LowDiff_errorbars','epsc')
            
        end
        
    end
    
    %%
    
    figure(LowandHighDiff_y_vals_max_fig)
    
    if Iter == 1
        
        set(gca, 'fontsize', 8)
        
        ylim_ref = ylim;
        
        ylim_min_y_vals_max = ylim_ref(1,1);
        ylim_max_y_vals_max = ylim_ref(1,2);
        
        title(['Drainage Divide Max. Position, Higher D: Hollow, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
        xlabel('K_W / K_S', 'fontsize', 10)
        ylabel('max(Y_{DD avg}*)', 'fontsize', 10)
        
    elseif Iter == 2
        
        ylim([ylim_min_y_vals_max, ylim_max_y_vals_max])
        
        %         h = colorbar;
        %         colormap(K_colors)
        %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
        %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
        
        legend([h_y_vals_max_LowDiff_HighAdvec, h_y_vals_max_LowDiff_LowAdvec], ...
            'High Advec.', 'Low Advec.', 'location', 'best')
        
        set(gcf, 'renderer', 'Painters')
        
        if Use_error_bars ~= 1
            
            saveas(LowandHighDiff_y_vals_max_fig,'DD_Position_Max_Fig_LowDiff.fig')
            saveas(LowandHighDiff_y_vals_max_fig,'DD_Position_Max_Fig_LowDiff.png')
            saveas(LowandHighDiff_y_vals_max_fig,'DD_Position_Max_Fig_LowDiff','epsc')
            
        elseif Use_error_bars == 1
            
            saveas(LowandHighDiff_y_vals_max_fig,'DD_Position_Max_Fig_LowDiff_errorbars.fig')
            saveas(LowandHighDiff_y_vals_max_fig,'DD_Position_Max_Fig_LowDiff_errorbars.png')
            saveas(LowandHighDiff_y_vals_max_fig,'DD_Position_Max_Fig_LowDiff_errorbars','epsc')
            
        end
        
    end
    
    %%
    
    figure(LowDiff_ksn_S_vals_fig)
    
    if Iter == 1
        
        set(gca, 'fontsize', 8)
        
        ylim_ref = ylim;
        
        ksn_S_min_y_vals = ylim_ref(1,1);
        ksn_S_max_y_vals = ylim_ref(1,2);
        
        title(['k_{sn} in Strong Unit, Lower D, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
        xlabel('t (Myr)', 'fontsize', 10)
        ylabel('k_{sn S avg}*(t)', 'fontsize', 10)
        
    elseif Iter == 2
        
        ylim([ksn_S_min_y_vals, ksn_S_max_y_vals])
        
        %         h = colorbar;
        %         colormap(K_colors)
        %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
        %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
        
        legend([h_ksn_S_LowDiff_HighAdvec, h_ksn_N_LowDiff_HighAdvec, ...
            h_ksn_S_LowDiff_LowAdvec, h_ksn_N_LowDiff_LowAdvec], ...
            'High Advec., S Basins', 'High Advec., N Basins', ...
            'Low Advec., S Basins', 'Low Advec., N Basins', 'location', 'best')
        
        set(gcf, 'renderer', 'Painters')
        
        if Use_error_bars ~= 1
            
            saveas(LowDiff_ksn_S_vals_fig,'ksn_in_S_over_time_Fig_LowDiff.fig')
            saveas(LowDiff_ksn_S_vals_fig,'ksn_in_S_over_time_Fig_LowDiff.png')
            saveas(LowDiff_ksn_S_vals_fig,'ksn_in_S_over_time_Fig_LowDiff','epsc')
            
        elseif Use_error_bars == 1
            
            saveas(LowDiff_ksn_S_vals_fig,'ksn_in_S_over_time_Fig_LowDiff_errorbars.fig')
            saveas(LowDiff_ksn_S_vals_fig,'ksn_in_S_over_time_Fig_LowDiff_errorbars.png')
            saveas(LowDiff_ksn_S_vals_fig,'ksn_in_S_over_time_Fig_LowDiff_errorbars','epsc')
            
        end
        
    end
    
    %%
    
    figure(LowDiff_Relief_S_vals_fig)
    
    if Iter == 1
        
        set(gca, 'fontsize', 8)
        
        ylim_ref = ylim;
        
        Relief_S_min_y_vals = ylim_ref(1,1);
        Relief_S_max_y_vals = ylim_ref(1,2);
        
        title(['Relief in Strong Unit, Lower D, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
        xlabel('t (Myr)', 'fontsize', 10)
        ylabel('Relief_{S avg}*(t)', 'fontsize', 10)
        
    elseif Iter == 2
        
        ylim([Relief_S_min_y_vals, Relief_S_max_y_vals])
        
        %         h = colorbar;
        %         colormap(K_colors)
        %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
        %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
        
        legend([h_Relief_S_LowDiff_HighAdvec, h_Relief_N_LowDiff_HighAdvec, ...
            h_Relief_S_LowDiff_LowAdvec, h_Relief_N_LowDiff_LowAdvec], ...
            'High Advec., S Basins', 'High Advec., N Basins', 'Low Advec., S Basins', 'Low Advec., N Basins', 'location', 'best')
        
        set(gcf, 'renderer', 'Painters')
        
        if Use_error_bars ~= 1
            
            saveas(LowDiff_Relief_S_vals_fig,'Relief_in_S_over_time_Fig_LowDiff.fig')
            saveas(LowDiff_Relief_S_vals_fig,'Relief_in_S_over_time_Fig_LowDiff.png')
            saveas(LowDiff_Relief_S_vals_fig,'Relief_in_S_over_time_Fig_LowDiff','epsc')
            
        elseif Use_error_bars == 1
            
            saveas(LowDiff_Relief_S_vals_fig,'Relief_in_S_over_time_Fig_LowDiff_errorbars.fig')
            saveas(LowDiff_Relief_S_vals_fig,'Relief_in_S_over_time_Fig_LowDiff_errorbars.png')
            saveas(LowDiff_Relief_S_vals_fig,'Relief_in_S_over_time_Fig_LowDiff_errorbars','epsc')
            
        end
        
    end
    
    %%
    
    figure(HighDiff_y_vals_fig)
    
    if Iter == 1
        
        set(gca, 'fontsize', 8)
        
        ylim_ref = ylim;
        
        if ylim_ref(1,1) < ylim_min_y_vals
            
            ylim_min_y_vals = ylim_ref(1,1);
            
        end
        
        if ylim_ref(1,2) > ylim_max_y_vals
            
            ylim_max_y_vals = ylim_ref(1,2);
            
        end
        
        title(['Drainage Divide Position, Higher D, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
        xlabel('t (Myr)', 'fontsize', 10)
        ylabel('Y_{DD}*(t)', 'fontsize', 10)
        
    elseif Iter == 2
        
        ylim([ylim_min_y_vals, ylim_max_y_vals])
        
        %         h = colorbar;
        %         colormap(K_colors)
        %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
        %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
        
        legend([h_mean_y_vals_HighDiff_HighAdvec, h_mean_y_vals_HighDiff_LowAdvec, h_S_gone_yvals_HD], ...
            'High Advec.', 'Low Advec.', 'S Unit Gone', 'location', 'best')
        
        set(gcf, 'renderer', 'Painters')
        
        if Use_error_bars ~= 1
            
            saveas(HighDiff_y_vals_fig,'DD_Position_over_time_Fig_HighDiff.fig')
            saveas(HighDiff_y_vals_fig,'DD_Position_over_time_Fig_HighDiff.png')
            saveas(HighDiff_y_vals_fig,'DD_Position_over_time_Fig_HighDiff','epsc')
            
        elseif Use_error_bars == 1
            
            saveas(HighDiff_y_vals_fig,'DD_Position_over_time_Fig_HighDiff_errorbars.fig')
            saveas(HighDiff_y_vals_fig,'DD_Position_over_time_Fig_HighDiff_errorbars.png')
            saveas(HighDiff_y_vals_fig,'DD_Position_over_time_Fig_HighDiff_errorbars','epsc')
            
        end
        
    end
    
    %%
    
%     figure(HighDiff_dy_dt_vals_fig)
%     
%     if Iter == 1
%         
%         set(gca, 'fontsize', 8)
%         
%         ylim_ref = ylim;
%         
%         if ylim_ref(1,1) < ylim_min_dy_dt_vals
%             
%             ylim_min_dy_dt_vals = ylim_ref(1,1);
%             
%         end
%         
%         if ylim_ref(1,2) > ylim_max_dy_dt_vals
%             
%             ylim_max_dy_dt_vals = ylim_ref(1,2);
%             
%         end
%         
%         title(['Drainage Divide Migration, Higher D, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
%         xlabel('K_W / K_S', 'fontsize', 10)
%         ylabel('(dY_{DD avg}*/dt) (yr^{-1})', 'fontsize', 10)
%         
%     elseif Iter == 2
%         
%         ylim([ylim_min_dy_dt_vals, ylim_max_dy_dt_vals])
%         
%         %         h = colorbar;
%         %         colormap(K_colors)
%         %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
%         %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
%         
%         legend([h_mean_y_vals_HighDiff_HighAdvec, h_mean_y_vals_HighDiff_LowAdvec], ...
%             'High Advec.', 'Low Advec.', 'S Unit Gone', 'location', 'best')
%         
%         set(gcf, 'renderer', 'Painters')
%         
%         if Use_error_bars ~= 1
%             
%             saveas(HighDiff_dy_dt_vals_fig,'DD_dydt_Fig_HighDiff.fig')
%             saveas(HighDiff_dy_dt_vals_fig,'DD_dydt_Fig_HighDiff.png')
%             saveas(HighDiff_dy_dt_vals_fig,'DD_dydt_Fig_HighDiff','epsc')
%             
%         elseif Use_error_bars == 1
%             
%             saveas(HighDiff_dy_dt_vals_fig,'DD_dydt_Fig_HighDiff_errorbars.fig')
%             saveas(HighDiff_dy_dt_vals_fig,'DD_dydt_Fig_HighDiff_errorbars.png')
%             saveas(HighDiff_dy_dt_vals_fig,'DD_dydt_Fig_HighDiff_errorbars','epsc')
%             
%         end
%         
%     end
    
    %%
    
%     figure(HighDiff_lag_time_fig)
%     
%     if Iter == 1
%         
%         set(gca, 'fontsize', 8)
%         
%         ylim_ref = ylim;
%         
%         if ylim_ref(1,1) < ylim_min_lag_time_vals
%             
%             ylim_min_lag_time_vals = ylim_ref(1,1);
%             
%         end
%         
%         if ylim_ref(1,2) > ylim_max_lag_time_vals
%             
%             ylim_max_lag_time_vals = ylim_ref(1,2);
%             
%         end
%         
%         title(['Drainage Divide Migration, Higher D, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
%         xlabel('K_W / K_W', 'fontsize', 10)
%         ylabel('Lag Time (Myr)', 'fontsize', 10)
%         
%     elseif Iter == 2
%         
%         ylim([ylim_min_lag_time_vals, ylim_max_lag_time_vals])
%         
%         %         h = colorbar;
%         %         colormap(K_colors)
%         %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
%         %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
%         
%         legend([h_mean_lag_time_HighDiff_HighAdvec, h_mean_lag_time_HighDiff_LowAdvec], ...
%             'High Advec.', 'Low Advec.', 'S Unit Gone', 'location', 'best')
%         
%         set(gcf, 'renderer', 'Painters')
%         
%         if Use_error_bars ~= 1
%             
%             saveas(HighDiff_lag_time_fig,'Lag_Time_Fig_HighDiff.fig')
%             saveas(HighDiff_lag_time_fig,'Lag_Time_Fig_HighDiff.png')
%             saveas(HighDiff_lag_time_fig,'Lag_Time_Fig_HighDiff','epsc')
%             
%         elseif Use_error_bars == 1
%             
%             saveas(HighDiff_lag_time_fig,'Lag_Time_Fig_HighDiff_errorbars.fig')
%             saveas(HighDiff_lag_time_fig,'Lag_Time_Fig_HighDiff_errorbars.png')
%             saveas(HighDiff_lag_time_fig,'Lag_Time_Fig_HighDiff_errorbars','epsc')
%             
%         end
%         
%     end
    
    %%
    
%     figure(HighDiff_y_vals_end_fig)
%     
%     if Iter == 1
%         
%         set(gca, 'fontsize', 8)
%         
%         ylim_ref = ylim;
%         
%         if ylim_ref(1,1) < ylim_min_y_vals_end
%             
%             ylim_min_y_vals_end = ylim_ref(1,1);
%             
%         end
%         
%         if ylim_ref(1,2) > ylim_max_y_vals_end
%             
%             ylim_max_y_vals_end = ylim_ref(1,2);
%             
%         end
%         
%         title(['Drainage Divide Final Position, Higher D, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
%         xlabel('K_W / K_S', 'fontsize', 10)
%         ylabel('Y_{DD avg}*(t_f)', 'fontsize', 10)
%         
%     elseif Iter == 2
%         
%         ylim([ylim_min_y_vals_end, ylim_max_y_vals_end])
%         
%         %         h = colorbar;
%         %         colormap(K_colors)
%         %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
%         %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
%         
%         legend([h_y_vals_end_HighDiff_HighAdvec, h_y_vals_end_HighDiff_LowAdvec], ...
%             'High Advec.', 'Low Advec.', 'location', 'best')
%         
%         set(gcf, 'renderer', 'Painters')
%         
%         if Use_error_bars ~= 1
%             
%             saveas(HighDiff_y_vals_end_fig,'DD_Position_Final_Fig_HighDiff.fig')
%             saveas(HighDiff_y_vals_end_fig,'DD_Position_Final_Fig_HighDiff.png')
%             saveas(HighDiff_y_vals_fig,'DD_Position_Final_Fig_HighDiff','epsc')
%             
%         elseif Use_error_bars == 1
%             
%             saveas(HighDiff_y_vals_end_fig,'DD_Position_Final_Fig_HighDiff_errorbars.fig')
%             saveas(HighDiff_y_vals_end_fig,'DD_Position_Final_Fig_HighDiff_errorbars.png')
%             saveas(HighDiff_y_vals_fig,'DD_Position_Final_Fig_HighDiff_errorbars','epsc')
%             
%         end
%         
%     end
    
    %%
    
%     figure(HighDiff_y_vals_max_fig)
%     
%     if Iter == 1
%         
%         set(gca, 'fontsize', 8)
%         
%         ylim_ref = ylim;
%         
%         if ylim_ref(1,1) < ylim_min_y_vals_max
%             
%             ylim_min_y_vals_max = ylim_ref(1,1);
%             
%         end
%         
%         if ylim_ref(1,2) > ylim_max_y_vals_max
%             
%             ylim_max_y_vals_max = ylim_ref(1,2);
%             
%         end
%         
%         title(['Drainage Divide Max. Position, Higher D, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
%         xlabel('K_W / K_S', 'fontsize', 10)
%         ylabel('max(Y_{DD avg}*)', 'fontsize', 10)
%         
%     elseif Iter == 2
%         
%         ylim([ylim_min_y_vals_max, ylim_max_y_vals_max])
%         
%         %         h = colorbar;
%         %         colormap(K_colors)
%         %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
%         %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
%         
%         legend([h_y_vals_max_HighDiff_HighAdvec, h_y_vals_max_HighDiff_LowAdvec], ...
%             'High Advec.', 'Low Advec.', 'location', 'best')
%         
%         set(gcf, 'renderer', 'Painters')
%         
%         if Use_error_bars ~= 1
%             
%             saveas(HighDiff_y_vals_max_fig,'DD_Position_Max_Fig_HighDiff.fig')
%             saveas(HighDiff_y_vals_max_fig,'DD_Position_Max_Fig_HighDiff.png')
%             saveas(HighDiff_y_vals_fig,'DD_Position_Max_Fig_HighDiff','epsc')
%             
%         elseif Use_error_bars == 1
%             
%             saveas(HighDiff_y_vals_max_fig,'DD_Position_Max_Fig_HighDiff_errorbars.fig')
%             saveas(HighDiff_y_vals_max_fig,'DD_Position_Max_Fig_HighDiff_errorbars.png')
%             saveas(HighDiff_y_vals_fig,'DD_Position_Max_Fig_HighDiff_errorbars','epsc')
%             
%         end
%         
%     end
    
    %%
    
    figure(HighDiff_ksn_S_vals_fig)
    
    if Iter == 1
        
        set(gca, 'fontsize', 8)
        
        ylim_ref = ylim;
        
        if ylim_ref(1,1) < ksn_S_min_y_vals
            
            ksn_S_min_y_vals = ylim_ref(1,1);
            
        end
        
        if ylim_ref(1,2) > ksn_S_max_y_vals
            
            ksn_S_max_y_vals = ylim_ref(1,2);
            
        end
        
        title(['k_{sn} in Strong Unit, Higher D, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
        xlabel('t (Myr)', 'fontsize', 10)
        ylabel('k_{sn S avg}*(t)', 'fontsize', 10)
        
    elseif Iter == 2
        
        ylim([ksn_S_min_y_vals, ksn_S_max_y_vals])
        
        %         h = colorbar;
        %         colormap(K_colors)
        %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
        %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
        
        legend([h_ksn_S_HighDiff_HighAdvec, h_ksn_N_HighDiff_HighAdvec, ...
            h_ksn_S_HighDiff_LowAdvec, h_ksn_N_HighDiff_LowAdvec], ...
            'High Advec., S Basins', 'High Advec., N Basins', ...
            'Low Advec., S Basins', 'Low Advec., N Basins', 'location', 'best')
        
        set(gcf, 'renderer', 'Painters')
        
        if Use_error_bars ~= 1
            
            saveas(HighDiff_ksn_S_vals_fig,'ksn_in_S_over_time_Fig_HighDiff.fig')
            saveas(HighDiff_ksn_S_vals_fig,'ksn_in_S_over_time_Fig_HighDiff.png')
            saveas(HighDiff_ksn_S_vals_fig,'ksn_in_S_over_time_Fig_HighDiff','epsc')
            
        elseif Use_error_bars == 1
            
            saveas(HighDiff_ksn_S_vals_fig,'ksn_in_S_over_time_Fig_HighDiff_errorbars.fig')
            saveas(HighDiff_ksn_S_vals_fig,'ksn_in_S_over_time_Fig_HighDiff_errorbars.png')
            saveas(HighDiff_ksn_S_vals_fig,'ksn_in_S_over_time_Fig_HighDiff_errorbars','epsc')
            
        end
        
    end
    
    %%
    
    figure(HighDiff_Relief_S_vals_fig)
    
    if Iter == 1
        
        set(gca, 'fontsize', 8)
        
        ylim_ref = ylim;
        
        if ylim_ref(1,1) < Relief_S_min_y_vals
            
            Relief_S_min_y_vals = ylim_ref(1,1);
            
        end
        
        if ylim_ref(1,2) > Relief_S_max_y_vals
            
            Relief_S_max_y_vals = ylim_ref(1,2);
            
        end
        
        title(['Relief in Strong Unit, Higher D, Dip = ' num2str(Dip_val_title, '%.0f') char(176)], 'fontsize', 10)
        xlabel('t (Myr)', 'fontsize', 10)
        ylabel('Relief_{S avg}*(t)', 'fontsize', 10)
        
    elseif Iter == 2
        
        ylim([Relief_S_min_y_vals, Relief_S_max_y_vals])
        
        %         h = colorbar;
        %         colormap(K_colors)
        %         caxis([min(Kw_over_Ks_vals) max(Kw_over_Ks_vals)])
        %         ylabel(h, 'K_W / K_S', 'fontsize', 12)
        
        legend([h_Relief_S_HighDiff_HighAdvec, h_Relief_N_HighDiff_HighAdvec, ...
            h_Relief_S_HighDiff_LowAdvec, h_Relief_N_HighDiff_LowAdvec], ...
            'High Advec., S Basins', 'High Advec., N Basins', 'Low Advec., S Basins', 'Low Advec., N Basins', 'location', 'best')
        
        set(gcf, 'renderer', 'Painters')
        
        if Use_error_bars ~= 1
            
            saveas(HighDiff_Relief_S_vals_fig,'Relief_in_S_over_time_Fig_HighDiff.fig')
            saveas(HighDiff_Relief_S_vals_fig,'Relief_in_S_over_time_Fig_HighDiff.png')
            saveas(HighDiff_Relief_S_vals_fig,'Relief_in_S_over_time_Fig_HighDiff','epsc')
            
        elseif Use_error_bars == 1
            
            saveas(HighDiff_Relief_S_vals_fig,'Relief_in_S_over_time_Fig_HighDiff_errorbars.fig')
            saveas(HighDiff_Relief_S_vals_fig,'Relief_in_S_over_time_Fig_HighDiff_errorbars.png')
            saveas(HighDiff_Relief_S_vals_fig,'Relief_in_S_over_time_Fig_HighDiff_errorbars','epsc')
            
        end
        
    end
    
end
