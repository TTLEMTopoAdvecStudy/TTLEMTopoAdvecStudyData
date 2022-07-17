% Written by Nate Mitchell (Indiana University, Dept. of Earth and
% Atmospheric Sciences). Uses Topotoolbox v2 scripts to produce .txt and
% shapefile output of river steepness and concavity using both slope-area
% and slope-integral (Chi) methods. The user can select a number of reaches
% to evaluate - note that if reaches overlap, the shapefile output will
% show the values for higher number reach instead of the lower number
% reach. The values and boundaries for each reach will still be in the .txt
% output, but try to avoid having the reaches overlap. Areas not included
% in reaches (for steepness or concavity) will have zeros in the shapefile
% output. If you are displaying these reach-wide values in ArcMap, you may
% want to set zero values as a separate color, so they don't look like low
% values. You will also need to set the projection of the shapefile output
% in ArcMap ("Define Projection" tool under Data Management, Projections
% and Transformations).

% If the automatic placement of figures doesn't work well on your screen,
% just use find and replace to turn all instances of "set(figure(" into
% "% set(figure(."

% Last updated 2.27.2018 by Nate Mitchell
% Topotoolbox v2 made by Wolfgang Schwanghart
% https://topotoolbox.wordpress.com/

clear
close all
clc

%% INPUT TO CHANGE

% Name of the ASCII .txt file of the DEM. Must have a projection that uses
% meters (not decimal degrees)! Basin_name_figures is how the study area
% will be described in figure titles.
Set = 6;
Scenario = 22;
time_Myr = 12;
Basin_name_figures = ['Set ' num2str(Set, '%.0f') ...
    ' Scenario ' num2str(Scenario, '%.0f') ...
    ', t = ' num2str(time_Myr, '%.1f') ' Myr'];

enforce_longprofile_lims = 1;
longprofile_xmin_km = 0;
longprofile_xmax_km = 20;
longprofile_zmin_m = 100;
longprofile_zmax_m = 900;
longprofile_ksnmin_val = 50;
longprofile_ksnmax_val = 400;

Initialization_run = 0;
KS_ref = 2.5E-09;

smooth_ksn_data = 1;
ksn_smooth_window = 10;

A_cr_ref = 1E+06;
manual_Acr_option = 0;

make_shapefiles_option = 0;

% You can have the code move to a different folder to find the input.
Input_location = pwd;
% Input_location = 'C:\Users\Nate\Documents\UI\Research\Fall_2017\transient\fall';

% Location where the code makes a new folder for just the output created
% here. More conventient then a ton of files created at different times,
% all in the same location. You may want to consolidate your ouput folders
% to one location, if so change "Output_location"
Output_location = pwd;
% Output_location = 'C:\Users\Nate\Documents\UI\Research\Summer_2017\Transient_Profile_Fitting\Output_folder_GSA2017\Selected_files\Chi_Projections';

% Option to fill (== 1) or not fill (~=1) the DEM file. If extracted from
% an already filled DEM, filling it again isn't necessary.
fillsinks_option = 0;

% Option to enable (== 1) or diable (~= 1) the enforcement of a theta
% value, avoiding the measurement of concavity. The shapefile output will
% have all zeros for concavity, and no .txt output willbe created for the
% concavity measurements.
use_ref_theta = 1;
theta_ref = 0.5;

% This is the reference concavity for the gridded, "at-a-point" ksn values.
% It is separate b/c the normal reference theta can be changed.
theta_ref_ksn_pnt = 0.5;

% Enables (== 1) or disables (~= 1) the use of the DrEICH method for
% identifying critical drainage areas. The method was created by Clubb et
% al. (2014), "Objective extraction of channel heads from high-resolution
% topographic data." If you are working on a transient profile, you may
% want to disable DrEICH and pick the critical drainage area manually. If
% you use DrEICH or pick a value manually, you should extend the chosen
% reach as far upstream as possible. Note that the DrEICH method uses
% theta_ref, but the manual measurement of theta occurs later in the code.
% If you want to use an exact theta value for DrEICH, set it as theta_ref.
use_DrEICH = 0;

% The number of points over which the raw profile will be smoothed. This
% value should be changed with DEM resolution. Ensure that the smoothing
% performed does not obliterate the profile's shape. The value used should
% be odd. Smoothing is not performed on the resampled profile, since this
% would be somewhat redundant. Set to 1 for no smoothing.
smooth_window = 1;

% Option to smooth (== 1) or not smooth (~= 1) elevations before using the
% DrEICH method. Uses smooth_window. Smoothing is only performed on the
% separate Z_sort_DrEICH variable, not Z_sort itself (that is done later in
% the code).
smooth_for_DrEICH = 0;

% Length of each river segment in the shapefile output. Should be changed
% with DEM resolution (might want a larger value for coarser resolutions).
seglength = 100;

% Drainage area (m^2) used as the initial minimum value. Minimum is later
% changed with the ginput function. Critical drainage areas are very
% unlikely to be less than 1E+05 m^2, so it is safe to leave it there.
Threshold_area_m2_initial = 1E+06;

% Elevation interval (m) used to resample profile and potentially avoid the
% effects of "step" artifacts in the DEM.
contour_interval_m = 2.5;

% Enables (== 1) or diables (~= 1) the resampling of the long profile. As
% the code stands currently, I would not use the resampled profile.
perform_resampling = 0;

% This variable enables (== 1) or disables (~= 1) a figure showing the DEM.
% This can help when selecting a reach, but it can also take a long time on
% a large DEM file.
DEM_figure = 0;

% This enables (== 1) or disables (~= 1) the creation of a figure showing
% the DEM with steepness values superimposed on it. As stated above for
% DEM_figure, this can be helpful when selecting a reach but it can take a
% long time (or fail) for a large DEM file.
ksn_figure = 1;

% This option enables (== 1) or disables (~= 1) the option to treat the
% selected reaches as relict and transient sections of a transient profile.
% The relict reach is then projected forward to estimate a depth of
% incision. This incision could be caused by uplift and/or base-level fall.
% The exact reach number (Reach_num_ksn) used to select the relict and
% adjusted reaches is not enforced; you should be systematic in using the
% same number for each reach (e.g., 1 for adjusted, 2 for relict). The
% number of reaches for the steepness section (Reach_num_ksn) must be 2 for
% these analyses.
transience_projections = 0;

%% INITIAL DEM PROCESSING

cd(Input_location)

File_name = ['Set' num2str(Set, '%.0f') '_Scenario' num2str(Scenario, '%.0f') ...
    '_t_' num2str(time_Myr * 1000, '%.0f') '_kyr.mat'];

load(File_name)

DEM = H1;

if fillsinks_option == 1
    
    DEM = fillsinks(DEM);
    
end

FlwDir = FLOWobj(DEM);

FlwAcc = flowacc(FlwDir);

FlwAcc.Z(1,:) = 0;
FlwAcc.Z(end,:) = 0;
FlwAcc.Z(:,1) = 0;
FlwAcc.Z(:,end) = 0;

Grad = gradient8(DEM, 'per');

FlwDst = flowdistance(FlwDir);

Stream = STREAMobj(FlwDir, FlwAcc >= (Threshold_area_m2_initial / (DEM.cellsize ^ 2)));

ksn_grid = (Grad ./ 100) ./ ((FlwAcc .* (FlwAcc.cellsize ^ 2)) .^ (-theta_ref_ksn_pnt));

if Initialization_run == 1
    
    Kw_grid = GRIDobj(H1);
    
    Kw_grid.Z = ((Kw_grid.Z .* 0) + 1) * KS_ref;
    
end

%% CREATES A MAP THE DEM

if DEM_figure == 1
    
    figure(2)
    
    imageschs(DEM, DEM, 'colormap', 'gray')
    
    h = colorbar;
    ylabel(h, 'Elevation ASL (m)', 'fontsize', 12)
    title({['DEM for ' char(Basin_name_figures)]}, 'fontsize', 14)
    
    set(figure(2), 'Position', [25 50 550 400])
    
end

clear DEM_figure

%% CREATES A MAP OF KSN VALUES OVERLAYING A HILLSHADE

if ksn_figure == 1
    
    Kw_colors = gray(200);
    Kw_colors = Kw_colors(100:end,:);
    
    % This is located here instead of the top b/c you probably don't need
    % to change it.
    Threshold_area_ksn_map = 1E+06;
    Stream_map = STREAMobj(FlwDir, FlwAcc >= (Threshold_area_ksn_map / (DEM.cellsize ^ 2)));
    
    Map_Struct = STREAMobj2mapstruct(Stream_map, 'seglength', 1000, 'attributes',...
        {'ksn' ksn_grid @mean});
    symbolspec = makesymbolspec('line', {'ksn' [min([Map_Struct.ksn]) max([Map_Struct.ksn])] ...
        'color' jet(100) 'linewidth' 1.5});
    
    ksn_grid_map = ksn_grid;
    ksn_grid_map.Z(FlwAcc.Z <= (Threshold_area_ksn_map / (DEM.cellsize ^ 2))) = NaN;
    
    figure(10)
    
    imageschs(DEM, Kw_grid, 'colormap', Kw_colors)
    
    mapshow(Map_Struct, 'SymbolSpec', symbolspec);
    
    xlabel('X (km)', 'fontsize', 12)
    ylabel('Y (km)', 'fontsize', 12)
    
    h = colorbar;
    colormap(jet)
    ylabel(h, 'k_{sn} (m)', 'fontsize', 12)
    caxis([min([Map_Struct.ksn]) max([Map_Struct.ksn])])
    title({['Steepness for ' char(Basin_name_figures)]}, 'fontsize', 14)
    
    set(gcf, 'renderer', 'Painters')
    
    set(figure(10), 'Position', [25 25 550 400])
    
    clear Threshold_area_ksn_map Stream_map ksn_grid_map Map_Struct symbolspec
    
end

clear ksn_figure

%% ALLOWS THE USER TO SELECT A SPECIFIC REACH

Stream_mod = modify(Stream, 'interactive', 'reachselect');

close(figure(gcf))

close(figure(2))

% close(figure(3))

clear Stream

%% ACQUIRES THE NODE ATTRIBUTE LIST (NAL) FOR SELECTED STREAM

nal = getnal(Stream_mod, FlwDst, DEM, FlwAcc, Grad, Kw_grid, 'struct');

[L_sort, I] = sort(nal.FlwDst, 'descend');

Z_sort = zeros(max(I), 1);

DA_sort = zeros(max(I), 1);

Kw_sort = zeros(max(I), 1);

% These slopes are calculated manually with the smoothed profile, and are
% used in the SA plot. The first and last points do not have values though,
% so these are not used to calculate ksn.

S_sort = zeros(max(I), 1);

for i = 1:1:length(I)
    
    Z_sort(i,1) = nal.DEM(I(i,1),1);
    
    DA_sort(i,1) = nal.FlwAcc(I(i,1),1) .* (DEM.cellsize ^ 2);
    
    Kw_sort(i,1) = nal.Kw_grid(I(i,1),1);
    
end

ksn_sort = zeros(max(I),1);
ksn_sort(1,1) = NaN;

for i = 2:1:length(Z_sort)
    
    S_sort(i,1) = (Z_sort(i-1,1) - Z_sort(i,1)) / (L_sort(i-1,1) - L_sort(i,1));
    
    ksn_sort(i,1) = S_sort(i,1) / (DA_sort(i,1) ^ -theta_ref);
    
end

clear nal

%% DEALS WITH REACHES SELECTED TO END UPSTREAM OF OUTLET

% If you have the selected reach end at a specific point rather than the
% watershed outlet, the Stream_mod.distance values will represent the
% dstance upstream of that point while the L_sort values will still
% represent distance upstream of the outlet (from FlwDst). I check if the
% minimum Stream_mod.distance is less than the minimum L_sort value and
% adjust L_sort if it is.

L_adjustment = min(L_sort);

if min(Stream_mod.distance) < min(L_sort)
    
    L_sort = L_sort - L_adjustment;
    
end

%% RESAMPLES TRUNK OF THE BASIN AT A SET ELEVATION INTERVAL

if perform_resampling == 1
    
    % Dummy placeholder variables
    L_resampled = ones( floor((max(Z_sort) - min(Z_sort)) / contour_interval_m), 1) .* -1;
    Z_resampled = ones( floor((max(Z_sort) - min(Z_sort)) / contour_interval_m), 1) .* -1;
    DA_resampled = ones( floor((max(Z_sort) - min(Z_sort)) / contour_interval_m), 1) .* -1;
    S_resampled = ones( floor((max(Z_sort) - min(Z_sort)) / contour_interval_m), 1) .* -1;
    ksn_resampled = ones( floor((max(Z_sort) - min(Z_sort)) / contour_interval_m), 1) .* -1;
    
    L_resampled(1,1) = L_sort(max(I),1);
    Z_resampled(1,1) = Z_sort(max(I),1);
    
    if isnan(DA_sort(max(I),1)) == 0
        
        DA_resampled(1,1) = DA_sort(max(I),1);
        
    elseif isnan(DA_sort(max(I),1)) == 1
        
        i = max(I) - 1;
        
        while isnan(DA_sort(i,1)) == 1
            
            i = i - 1;
            
        end
        
        DA_resampled(1,1) = DA_sort(i,1);
        
    end
    
    resampling = 1;
    
    i = max(I);
    
    j = 1;
    
    while resampling == 1
        
        i = i - 1;
        
        if i > 0 && Z_sort(i,1) >= Z_resampled(j,1) + contour_interval_m && isnan(DA_sort(i,1)) == 0
            
            L_resampled(j + 1,1) = L_sort(i,1);
            
            Z_resampled(j + 1,1) = Z_sort(i,1);
            
            DA_resampled(j + 1,1) = DA_sort(i,1);
            
            j = j + 1;
            
        elseif i <= 0
            
            resampling = 0;
            
        end
        
    end
    
    % Gets rid of the excess placeholder values
    L_resampled = L_resampled(L_resampled ~= -1);
    Z_resampled = Z_resampled(Z_resampled ~= -1);
    DA_resampled = DA_resampled(DA_resampled ~= -1);
    
    L_resampled = flipud(L_resampled);
    Z_resampled = flipud(Z_resampled);
    DA_resampled = flipud(DA_resampled);
    
    for i = 2:1:length(L_resampled)
        
        S_resampled(i-1,1) = ( Z_resampled(i-1,1) - Z_resampled(i,1) ) / ( L_resampled(i-1,1) - L_resampled(i,1) );
        
        ksn_resampled(i-1,1) = S_resampled(i-1,1) / (DA_resampled(i-1,1) ^ -theta_ref);
        
    end
    
    % Gets rid of the excess placeholder values
    S_resampled = S_resampled(S_resampled ~= -1);
    ksn_resampled = ksn_resampled(S_resampled ~= -1);
    
    % Adding one extra NaN row makes these the same length as the other
    % resampled variables
    S_resampled = [NaN; S_resampled];
    ksn_resampled = [NaN; ksn_resampled];
    
end

clear I i j DEM FlwDir

%% RAW AND RESAMPLED LONG PROFILES

figure(3)

set(figure(3), 'defaultAxesColorOrder', [[0 0 0]; [0 0 0]]);

yyaxis left

h1 = plot(L_sort, Z_sort, 'linewidth', 3, 'color', 'k', 'marker', 'none', ...
    'linestyle', '-');
hold on

yyaxis right

h2 = plot(L_sort, ksn_sort, 'linewidth', 1, 'color', 'r', 'marker', 'none', ...
    'linestyle', ':');
hold on

if perform_resampling == 1
    
    yyaxis left
    
    h3 = plot(L_resampled, Z_resampled, 'linewidth', 2, 'color', 'b', 'marker', 'none', ...
        'linestyle', '--');
    hold on
    
    yyaxis right
    
    h4 = plot(L_resampled, ksn_resampled, 'linewidth', 2, 'color', 'y', 'marker', 'none', ...
        'linestyle', '--');
    hold off
    
    lgnd = legend([h1; h2; h3; h4], 'Z_{raw}', 'k_{sn raw}', 'Z_{resamp}', 'k_{sn resamp}', 'location', 'southeast');
    
else
    
    lgnd = legend([h1; h2], 'Z_{raw}', 'k_{sn raw}', 'location', 'southeast');
    
end

set(gca, 'Fontsize', 12)
lgnd.FontSize = 8;
xlabel('Distance Upstream (m)', 'FontSize', 14)

yyaxis left

ylabel('Elevation ASL (m)', 'FontSize', 14)

yyaxis right

if (2 * theta_ref) == 1
    
    ylabel('k_{sn} (m)', 'fontsize', 14)
    
else
    
    ylabel(['k_{sn} (m^{' num2str(2 * theta_ref, '%.2f') '})'], 'fontsize', 14)
    
end

title({['Long Profile for ' char(Basin_name_figures)]}, 'fontsize', 14)

set(figure(3), 'Position', [25 215 550 430])

clear lgnd h1 h2 h3 h4

%% CALCULATE CHI USING REFERENCE THETA

% Uses the approach of Perron and Royden (2013) to transform distance
% upstream using drainage area values. A steady-state reach with uniform
% properties should have a linear relationship between elevation and Chi.

% First, determine the average dx along the profile.
dx_all = zeros(length(L_sort), 1);

for i = 1:1:length(L_sort)
    
    if i == length(L_sort)
        
        dx_all(i,1) = L_sort(i,1);
        
    else
        
        dx_all(i,1) = (L_sort(i,1) - L_sort(i+1,1));
        
    end
    
end

dx_avg = mean(dx_all);

clear dx_all

% A0 is a reference drainage area. See Perron and Royden (2013).
A0 = 1E+06;

Chi = zeros(length(L_sort), 1);

Chi_partial = (A0 ./ DA_sort) .^ (theta_ref);

for i = length(Chi):-1:1
    
    Chi(i,1) = sum(Chi_partial(i:end,1)) * (dx_avg);
    
end

clear Chi_partial

%% USE DrEICH METHOD TO FIND CRITICAL DRAINAGE AREA

if use_DrEICH == 1
    
    % This variable is used to loop through DrEICH if you want to redo your
    % manual selection of a lower limit.
    continue_test_DrEICH = 0;
    
    while continue_test_DrEICH == 0
        
        figure(2)
        
        plot(Chi, Z_sort, 'color', 'k', 'linewidth', 2)
        
        set(gca, 'fontsize', 12)
        xlabel('\chi (m)', 'fontsize', 14)
        ylabel('Elevation ASL (m)', 'fontsize', 14)
        title(['\chi-Plot for ' char(Basin_name_figures)], 'fontsize', 14)
        
        set(figure(2), 'Position', [600 215 550 430])
        
        figure(4)
        
        loglog(DA_sort(S_sort ~= 0), S_sort(S_sort ~= 0), 'ko', 'linewidth', 1, 'linestyle', 'none')
        
        set(gca, 'fontsize', 12)
        xlabel('Drainage Area (m^2)', 'fontsize', 14)
        ylabel('Slope', 'fontsize', 14)
        title(['Slope-Area Plot for ' char(Basin_name_figures)], 'fontsize', 14)
        
        set(figure(4), 'Position', [1175 215 550 430])
        
        continue_test = 0;
        
        % Before using DrEICH, gives you the opportunity to clip out higher
        % drainage areas. For example, if you're dealing with a transient
        % profile, DrEICH would start at the bottom of the profile looking for
        % the start of nonlinearity in the chi plot. It would be confused,
        % however, by the nonlinearity you can get in a knickzone, and DrEICH
        % isn't meant to be used in that scenario.
        while continue_test == 0
            
            prompt = ['Set a lower limit for the fluvial section before identifying the critical drainage ' ...
                'area with DrEICH (''y'' or ''n'')? Necessary for a transient profile.'];
            
            prompt_title = 'Limit Fluvial Section?';
            
            Acr_response = inputdlg(prompt, prompt_title);
            
            Positive_response = 'y';
            
            Negative_response = 'n';
            
            if strcmp(char(Acr_response), Positive_response)
                
                DrEICH_clip = 1;
                
                continue_test = 1;
                
            end
            
            if strcmp(char(Acr_response), Negative_response)
                
                DrEICH_clip = 0;
                
                A_clip_DrEICH = 0;
                
                S_clip_DrEICH = 0;
                
                continue_test = 1;
                
            end
            
        end
        
        clear Positive_response Negative_response Acr_response continue_test
        
        %% CLIP FOR DrEICH, IF USER CHOOSES TO
        
        if DrEICH_clip == 1
            
            continue_test = 0;
            
            while continue_test == 0
                
                figure(2)
                
                plot(Chi, Z_sort, 'color', 'k', 'linewidth', 2)
                hold on
                
                set(gca, 'fontsize', 12)
                xlabel('\chi (m)', 'fontsize', 14)
                ylabel('Elevation ASL (m)', 'fontsize', 14)
                title(['\chi-Plot for ' char(Basin_name_figures)], 'fontsize', 14)
                
                set(figure(2), 'Position', [25 430 550 430])
                
                figure(4)
                
                loglog(DA_sort(S_sort ~= 0), S_sort(S_sort ~= 0), 'ko', 'linewidth', 1, 'linestyle', 'none')
                hold on
                
                set(gca, 'fontsize', 12)
                xlabel('Drainage Area (m^2)', 'fontsize', 14)
                ylabel('Slope', 'fontsize', 14)
                title('Select Lower Limit for Fluvial Section', 'fontsize', 14)
                
                set(figure(4), 'Position', [650 230 1150 630])
                
                [DA_clip_ginput, S_clip_ginput] = ginput(1);
                
                % Find the minimum distance between the clicked location and the data
                % points
                DA_sort_clip = DA_sort(S_sort ~= 0 & isnan(S_sort) == 0);
                S_sort_clip = S_sort(S_sort ~= 0 & isnan(S_sort) == 0);
                
                % Here, I transform the slope-area data so that distance from all
                % points can be calculated (i.e., distance isn't the same in a log-log
                % plot).
                DA_clip_diff = abs(log(DA_sort_clip) - log(DA_clip_ginput));
                S_clip_diff = abs(log(S_sort_clip) - log(S_clip_ginput));
                Dist_clip = ((DA_clip_diff .^ 2) + (S_clip_diff .^ 2)) .^ 0.5;
                
                A_clip_DrEICH = DA_sort_clip(Dist_clip == min(Dist_clip));
                S_clip_DrEICH = S_sort_clip(Dist_clip == min(Dist_clip));
                
                % Plot selection on SA Plot
                figure(4)
                
                loglog(A_clip_DrEICH, S_clip_DrEICH, 'marker', 'o', 'color', 'k', 'markerfacecolor', 'r', ...
                    'linewidth', 0.5, 'linestyle', 'none')
                hold on
                
                loglog([A_clip_DrEICH, A_clip_DrEICH], [min(S_sort_clip), max(S_sort_clip)], 'color', 'r', ...
                    'linewidth', 2, 'linestyle', '--')
                hold off
                
                % Plot selection on Chi Plot
                figure(2)
                
                plot(Chi(DA_sort == A_clip_DrEICH), Z_sort(DA_sort == A_clip_DrEICH), 'marker', 'o', ...
                    'color', 'r', 'linewidth', 1.5, 'linestyle', 'none')
                hold off
                
                % Evaluate
                prompt = 'Reselect the upper drainage area limit for the fluvial section (''y'' or ''n'')? ';
                
                prompt_title = 'Reselect?';
                
                Acr_response = inputdlg(prompt, prompt_title);
                
                Correct_response = 'n';
                
                if strcmp(char(Acr_response), Correct_response)
                    
                    continue_test = 1;
                    
                end
                
            end
            
        end
        
        close(figure(2), figure(4))
        
        clear Acr_response Correct_response DA_clip_ginput S_clip_ginput DA_clip_diff ...
            S_clip_diff
        
        % Uses these "DrEICH" verisons of main variables in the following
        % DrEICH section
        if DrEICH_clip == 1
            
            Z_sort_DrEICH = Z_sort(DA_sort <= A_clip_DrEICH);
            Chi_DrEICH = Chi(DA_sort <= A_clip_DrEICH) - min(Chi(DA_sort <= A_clip_DrEICH));
            DA_sort_DrEICH = DA_sort(DA_sort <= A_clip_DrEICH);
            
        elseif DrEICH_clip == 0
            
            Z_sort_DrEICH = Z_sort;
            Chi_DrEICH = Chi;
            DA_sort_DrEICH = DA_sort;
            
        end
        
        if smooth_for_DrEICH == 1
            
            Z_sort_DrEICH = smooth(Z_sort_DrEICH, smooth_window);
            
        end
        
        %%
        
        % Uses the Drainage Extraction by Identifying Channel Head (DrEICH) method.
        % Loop through different selections of the Chi plot, with channel and
        % hillslope segments varying in size. The R^2 of the channel segment and
        % the Durbin-Watson statistic (d) of the hllslope are used to calculate a
        % test value (t). The channel head / critical drainage area occurs at the
        % maximum t values. See Clubb et al. (2014).
        
        Channel_segments_R2 = zeros(length(Z_sort_DrEICH), 1);
        Channel_segments_ksn = zeros(length(Z_sort_DrEICH), 1);
        
        Hillslope_segments_d = zeros(length(Z_sort_DrEICH), 1);
        Hillslope_segments_ksn = zeros(length(Z_sort_DrEICH), 1);
        
        t_value_DrEICH = zeros(length(Z_sort_DrEICH), 1);
        
        % Start at 2 from the bottom, so the first channel segment has at least
        % 3 points. The last hillslope segment also needs to have 3 points, but
        % the last i value needs to be 4 b/c the hillslope segment is 1:1:i-1.
        % The channel segment is end:-1:i.
        for i = (length(Z_sort_DrEICH) - 2):-1:4
            
            % First, fit linear model to channel segment
            Channel_Chi = Chi_DrEICH(end:-1:i,1);
            
            Channel_Z = Z_sort_DrEICH(end:-1:i, 1);
            
            Chi_fitlm = fitlm(Channel_Chi - min(Channel_Chi), Channel_Z - min(Channel_Z), ...
                'linear', 'intercept', false);
            
            Coeffs_temp = table2array(Chi_fitlm.Coefficients);
            
            Channel_segments_ksn(i,1) = Coeffs_temp(1,1) * (A0 ^ theta_ref);
            
            Channel_segments_R2(i,1) = Chi_fitlm.Rsquared.Adjusted;
            
            % In case the data level off and R^2 is NaN, set it to a very small
            % number
            if isnan(Channel_segments_R2(i,1)) == 1
                
                Channel_segments_R2(i,1) = 1E-6;
                
            end
            
            
            % Now fit linear model to hillslope segment
            Hillslope_Chi = Chi_DrEICH(1:i-1, 1);
            
            Hillslope_Z = Z_sort_DrEICH(1:i-1, 1);
            
            Hillslope_fitlm = fitlm(Hillslope_Chi - min(Hillslope_Chi), Hillslope_Z - min(Hillslope_Z), ...
                'linear', 'intercept', false);
            
            Coeffs_temp = table2array(Hillslope_fitlm.Coefficients);
            
            Hillslope_segments_ksn(i,1) = Coeffs_temp(1,1) * (A0 ^ theta_ref);
            
            Hillslope_Z_pred = (Hillslope_Chi .* Hillslope_segments_ksn(i,1) .* A0 ^ (-theta_ref)) - ...
                min((Hillslope_Chi .* Hillslope_segments_ksn(i,1) .* A0 ^ (-theta_ref))) + min(Hillslope_Z);
            
            Hillslope_Z_residuals = Hillslope_Z_pred - Hillslope_Z;
            
            [~, Hillslope_segments_d(i,1)] = dwtest(Hillslope_Z_residuals, Hillslope_Chi - min(Hillslope_Chi));
            
            % Calculate test value t
            t_value_DrEICH(i,1) = Channel_segments_R2(i,1) - ((Hillslope_segments_d(i,1) - 2) / 2);
            
        end
        
        [row_max_t, ~] = find(t_value_DrEICH == max(t_value_DrEICH));
        
        Channel_Z = Z_sort_DrEICH(end:-1:row_max_t, 1);
        Channel_Chi = Chi_DrEICH(end:-1:row_max_t, 1);
        Channel_ksn = Channel_segments_ksn(row_max_t, 1);
        
        Hillslope_Z = Z_sort_DrEICH(1:(row_max_t-1), 1);
        Hillslope_Chi = Chi_DrEICH(1:(row_max_t-1), 1);
        Hillslope_ksn = Hillslope_segments_ksn(row_max_t, 1);
        
        A_cr = DA_sort_DrEICH(row_max_t, 1);
        
        clear Channel_segments_R2 Channel_segments_ksn Hillslope_segments_d ...
            Hillslope_segments_ksn Hillslope_Z_pred t_value_DrEICH Chi_fitlm ...
            Hillslope_fitlm row_max_t
        
        % Plot results from DrEICH
        figure(2)
        
        plot(Chi_DrEICH, Z_sort_DrEICH, 'color', 'k', 'linewidth', 2)
        hold on
        
        plot([min(Channel_Chi), max(Channel_Chi)], ((Channel_ksn * (A0 ^ -theta_ref)) .* ...
            [min(Channel_Chi), max(Channel_Chi)]) - ((Channel_ksn * (A0 ^ -theta_ref)) * ...
            min(Channel_Chi)) + min(Channel_Z), 'color', 'b', 'linewidth', 2, 'linestyle', '--')
        hold on
        
        plot([min(Hillslope_Chi), max(Hillslope_Chi)], ((Hillslope_ksn * (A0 ^ -theta_ref)) .* ...
            [min(Hillslope_Chi), max(Hillslope_Chi)]) - ((Hillslope_ksn * (A0 ^ -theta_ref)) * ...
            min(Hillslope_Chi)) + min(Hillslope_Z), 'color', 'r', 'linewidth', 2, 'linestyle', '--')
        hold on
        
        plot(max(Channel_Chi), max(Channel_Z), 'color', 'b', 'linewidth', 2, ...
            'marker', 'o', 'linestyle', 'none')
        hold on
        
        hold off
        
        ylim([min(Z_sort_DrEICH) (ceil((max(Z_sort_DrEICH) / 50)) * 50)])
        
        xaxis_limits = xlim;
        yaxis_limits = ylim;
        
        h = text(max(Channel_Chi) + ((xaxis_limits(1,2) - xaxis_limits(1,1)) / 25), max(Channel_Z) - ...
            ((yaxis_limits(1,2) - yaxis_limits(1,1)) / 15), ['k_{sn}: ' num2str(Channel_ksn, ...
            '%.2E') ' m^{' num2str(2 * theta_ref, '%.2f') '}']);
        set(h, 'fontsize', 10)
        set(h, 'color', 'b')
        set(h, 'rotation', 20)
        
        set(gca, 'fontsize', 12)
        xlabel('\chi (m)', 'fontsize', 14)
        ylabel('Elevation ASL (m)', 'fontsize', 14)
        title({['DrEICH Results for ' char(Basin_name_figures) ', A_{cr} = ' ...
            num2str(A_cr, '%.2E') ' m^2']}, 'fontsize', 12)
        lgnd = legend('All Data', 'Best-Fit Channel Segment', 'Best-Fit Hillslope Segment', ...
            ['Channel Head, ' num2str(max(Channel_Z), '%.2E') ' m ASL'], 'Location', 'Southeast');
        lgnd.FontSize = 10;
        
        set(figure(2), 'Position', [600 215 550 430])
        
        continue_test = 0;
        
        while continue_test == 0
            
            % Allows the user to judge the DrEICH estimate for critical drainage
            % area and manually select a value if desired.
            prompt = 'Use the critical drainage area derived from the DrEICH method (''y'' or ''n'')? Or redo selection of a lower limit (''r'')?';
            
            prompt_title = 'Use DrEICH?';
            
            Acr_response = inputdlg(prompt, prompt_title);
            
            Positive_response = 'y';
            
            Negative_response = 'n';
            
            Redo_response = 'r';
            
            if strcmp(char(Acr_response), Positive_response)
                
                continue_test = 1;
                
                continue_test_DrEICH = 1;
                
            elseif strcmp(char(Acr_response), Redo_response)
                
                continue_test = 1;
                
            elseif strcmp(char(Acr_response), Negative_response)
                
                use_DrEICH = 0;
                
                continue_test = 1;
                
                continue_test_DrEICH = 1;
                
            end
            
        end
        
        clear Hillslope_Chi Hillslope_ksn Hillslope_Z Hillslope_Z_residuals ...
            Channel_Chi Channel_ksn Channel_Z Channel_Z_residuals Acr_response ...
            Negative_response Z_sort_DrEICH Chi_DrEICH DA_sort_DrEICH Redo_response ...
            Positive_response
        
    end
    
end

clear continue_test continue_test_DrEICH

%% IF NOT USING DrEICH, MANUALLY PICK CRITICAL DRAINAGE AREA

% Note that you might want to pick a higher drainage area to avoid, for
% example, a different lithology at lower drainage areas. In that case, you
% wouldn't have a critical drainage area, but only a minimum.
if manual_Acr_option == 1
    
    if use_DrEICH == 0
        
        % If you choose to avoid DrEICH from the start, these variables won't
        % be set and the code will freak out when they are written in the
        % summary table.
        A_clip_DrEICH = 0;
        DrEICH_clip = 0;
        
        continue_test = 0;
        
        while continue_test == 0
            
            figure(1)
            
            loglog(DA_sort(S_sort ~= 0), S_sort(S_sort ~= 0), 'ko', 'linewidth', 1, ...
                'linestyle', 'none')
            hold on
            
            if perform_resampling == 1
                
                loglog(DA_resampled, S_resampled, 'bo', 'linewidth', 1, 'linestyle', 'none')
                hold on
                
            end
            
            set(gca, 'fontsize', 12)
            xlabel('Drainage Area (m^2)', 'fontsize', 16)
            ylabel('Slope', 'fontsize', 16)
            title({'Click to Select the Critical (or Minimum) Drainage Area'}, 'FontSize', 16)
            
            if perform_resampling == 1
                
                legend('Raw Data', 'Resampled Data', 'Location', 'SouthWest')
                
            end
            
            set(figure(1), 'Position', [650 230 1150 630])
            
            [DA_crit_ginput, S_crit_ginput] = ginput(1);
            
            % Find the minimum distance between the clicked location and the data
            % points
            DA_sort_clip = DA_sort(S_sort ~= 0 & isnan(S_sort) == 0);
            S_sort_clip = S_sort(S_sort ~= 0 & isnan(S_sort) == 0);
            
            % Here, I transform the slope-area data so that distance from all
            % points can be calculated (i.e., distance isn't the same in a log-log
            % plot).
            DA_crit_diff = abs(log(DA_sort_clip) - log(DA_crit_ginput));
            S_crit_diff = abs(log(S_sort_clip) - log(S_crit_ginput));
            Dist_crit = ((DA_crit_diff .^ 2) + (S_crit_diff .^ 2)) .^ 0.5;
            
            A_cr = DA_sort_clip(Dist_crit == min(Dist_crit));
            S_A_cr = S_sort_clip(Dist_crit == min(Dist_crit));
            
            figure(1)
            
            loglog(A_cr, S_A_cr, 'marker', 'o', 'color', 'k', 'markerfacecolor', 'r', ...
                'linewidth', 0.5, 'linestyle', 'none')
            hold on
            
            loglog([A_cr, A_cr], [min(S_sort_clip), max(S_sort_clip)], 'color', 'r', ...
                'linewidth', 2, 'linestyle', '--')
            hold off
            
            prompt = 'Reselect the critical (or minimum) drainage area (''y'' or ''n'')? ';
            
            prompt_title = 'Reselect?';
            
            Acr_response = inputdlg(prompt, prompt_title);
            
            Correct_response = 'n';
            
            if strcmp(char(Acr_response), Correct_response)
                
                figure(1)
                
                set(gca, 'fontsize', 12)
                xlabel('Drainage Area (m^2)', 'fontsize', 14)
                ylabel('Slope', 'fontsize', 14)
                title({[ 'Manual selection for ' char(Basin_name_figures) ': A_{cr} = ' num2str(A_cr, '%.2E')]}, 'FontSize', 12)
                
                set(figure(1), 'Position', [600 215 550 430])
                
                continue_test = 1;
                
            end
            
        end
        
        clear DA_crit_ginput S_crit_ginput DA_crit_diff S_crit_diff DA_sort_clip S_sort_clip ...
            Acr_response Correct_response Dist_crit S_A_cr
        
    end
    
elseif manual_Acr_option ~= 1
    
    A_cr = A_cr_ref;
    
end

%% CLIP PROFILE USING CRITICAL DRAINAGE AREA

% Storing these in case they are needed
L_sort_all = L_sort;
Z_sort_all = Z_sort;
DA_sort_all = DA_sort;
S_sort_all = S_sort;
Chi_all = Chi;
ksn_sort_all = ksn_sort;
Kw_sort_all = Kw_sort;

L_sort = L_sort(DA_sort >= A_cr);
Z_sort = Z_sort(DA_sort >= A_cr);
DA_sort = DA_sort(DA_sort >= A_cr);
S_sort = S_sort(DA_sort >= A_cr);
ksn_sort = ksn_sort(DA_sort >= A_cr);
Kw_sort = Kw_sort(DA_sort >= A_cr);

if perform_resampling == 1
    
    % Storing these in case they are needed
    L_resampled_all = L_resampled;
    Z_resampled_all = Z_resampled;
    DA_resampled_all = DA_resampled;
    S_resampled_all = S_resampled;
    ksn_resampled_all = ksn_resampled;
    
    L_resampled = L_resampled(DA_resampled >= A_cr);
    Z_resampled = Z_resampled(DA_resampled >= A_cr);
    DA_resampled = DA_resampled(DA_resampled >= A_cr);
    S_resampled = S_resampled(DA_resampled >= A_cr);
    ksn_resampled = ksn_resampled(DA_resampled >= A_cr);
    
end

% Number_pixels_remove = length(Stream_mod.distance(Stream_mod.distance > max(L_sort)));
%
% % Shortens the stream to account for the critical (or minimum) drainage
% % area.
% Stream_mod = modify(Stream_mod, 'shrinkfromtop', Number_pixels_remove);
%
% clear Number_pixels_remove

%% DECIDE BETWEEN RAW AND RESAMPLED PROFILES

if perform_resampling == 1
    
    continue_test = 0;
    
    while continue_test == 0
        
        prompt = 'Use the resampled profile instead of the raw profile (''y'' or ''n'')? ';
        
        prompt_title = 'Select Profile';
        
        Profile_response = inputdlg(prompt, prompt_title);
        
        Correct_response1 = 'n';
        
        Correct_response2 = 'y';
        
        if strcmp(char(Profile_response), Correct_response1)
            
            % Now, smooth elevation data. Didn't smooth before so the raw slope-area
            % data could be interpreted.
            
            Z_sort = smooth(Z_sort, smooth_window);
            
            for i = 2:1:length(Z_sort)
                
                S_sort(i,1) = (Z_sort(i-1,1) - Z_sort(i,1)) / (L_sort(i-1,1) - L_sort(i,1));
                
            end
            
            use_resampled = 0;
            
            continue_test = 1;
            
        elseif strcmp(char(Profile_response), Correct_response2)
            
            L_sort = L_resampled;
            Z_sort = Z_resampled;
            DA_sort = DA_resampled;
            S_sort = S_resampled;
            ksn_sort = ksn_resampled;
            
            use_resampled = 1;
            
            clear L_resampled Z_resampled DA_resampled S_resampled ksn_resampled
            
            continue_test = 1;
            
        end
        
    end
    
else
    
    % Now, smooth elevation data. Didn't smooth before so the raw slope-area
    % data could be interpreted.
    Z_sort = smooth(Z_sort, smooth_window);
    
    for i = 2:1:length(Z_sort)
        
        S_sort(i,1) = (Z_sort(i-1,1) - Z_sort(i,1)) / (L_sort(i-1,1) - L_sort(i,1));
        
    end
    
    use_resampled = 0;
    
end

%% UPDATE LONG PROFILE

close(figure(3))

figure(3)

set(figure(3), 'defaultAxesColorOrder', [[0 0 0]; [0 0 0]]);

yyaxis left

h1 = plot(L_sort(Kw_sort == min(Kw_sort)) ./ 1000, Z_sort(Kw_sort == min(Kw_sort)), 'linewidth', 2, 'color', 'k', 'marker', 'none', ...
    'linestyle', '-');
hold on

plot(L_sort(Kw_sort == max(Kw_sort)) ./ 1000, Z_sort(Kw_sort == max(Kw_sort)), 'linewidth', 2, 'color', [0.5 0.5 0.5], 'marker', 'none', ...
    'linestyle', '-');
hold on

yyaxis right

if smooth_ksn_data == 1
    
    h2 = plot(L_sort ./ 1000, smooth(ksn_sort, ksn_smooth_window), 'linewidth', 1, 'color', 'r', 'marker', 'none', ...
        'linestyle', ':');
    hold on
    
elseif smooth_ksn_data ~= 1
    
    h2 = plot(L_sort ./ 1000, ksn_sort, 'linewidth', 1, 'color', 'r', 'marker', 'none', ...
        'linestyle', ':');
    hold on
    
end

if use_resampled ~= 1
    
    if smooth_window == 1
        
        lgnd = legend([h1; h2], 'Z_{raw}', 'k_{sn raw}', 'location', 'southeast');
        
    else
        
        lgnd = legend([h1; h2], 'Z_{smooth}', 'k_{sn smooth}', 'location', 'southeast');
        
    end
    
elseif use_resampled == 1
    
    lgnd = legend([h1; h2], 'Z_{resamp}', 'k_{sn resamp}', 'location', 'southeast');
    
end

set(gca, 'Fontsize', 12)
lgnd.FontSize = 8;
xlabel('Distance Upstream (km)', 'FontSize', 14)

yyaxis left

hold off

ylabel('Elevation ASL (m)', 'FontSize', 14)

yyaxis right

hold off

if (2 * theta_ref) == 1
    
    ylabel('k_{sn} (m)', 'fontsize', 14)
    
else
    
    ylabel(['k_{sn} (m^{' num2str(2 * theta_ref, '%.2f') '})'], 'fontsize', 14)
    
end

title({['Long Profile for ' char(Basin_name_figures)]}, 'fontsize', 14)

set(figure(3), 'Position', [25 215 550 430])

clear lgnd h1 h2

if enforce_longprofile_lims == 1
    
    xlim([longprofile_xmin_km, longprofile_xmax_km])
    
    yyaxis left
    ylim([longprofile_zmin_m, longprofile_zmax_m])
    
    yyaxis right
    ylim([longprofile_ksnmin_val, longprofile_ksnmax_val])
    
end

%% SELECT REACH(ES) FROM SLOPE-AREA, USE TO CALCULATE CONCAVITY

if use_ref_theta ~= 1
    
    continue_test = 0;
    
    while continue_test == 0
        
        figure(4)
        
        loglog(DA_sort(S_sort ~= 0), S_sort(S_sort ~= 0), 'ko', 'linewidth', 1, ...
            'linestyle', 'none')
        hold on
        
        set(gca, 'fontsize', 12)
        xlabel('Drainage Area (m^2)', 'fontsize', 16)
        ylabel('Slope', 'fontsize', 16)
        
        set(figure(4), 'Position', [650 230 1150 630])
        
        % Now, I allow the steepness and concavity of a number of
        % reaches to be assessed (e.g., variable rock uplift rates).
        Reach_num_continue = 0;
        
        while Reach_num_continue == 0
            
            prompt = 'Specify the number of reaches you would like to assess. ';
            
            prompt_title = 'Number of Reaches';
            
            Reach_num_input = inputdlg(prompt, prompt_title);
            
            Reach_num_char = cell2mat(Reach_num_input);
            
            Reach_num_concav = str2double(Reach_num_char);
            
            % The number of reaches has to be a positive integer,
            % otherwise the code won't continue here.
            
            if mod(Reach_num_concav, 1) == 0 && Reach_num_concav > 0
                
                Reach_num_continue = 1;
                
            end
            
        end
        
        clear Reach_num_continue Reach_num_input Reach_num_char
        
        % This variable holds the drainage area bounds selected by the
        % user. First column is lower DA bound, second is higher.
        
        DA_m2_bounds_fromSA = zeros(Reach_num_concav, 2);
        S_bounds_fromSA = zeros(Reach_num_concav, 2);
        L_m_bounds_fromSA = zeros(Reach_num_concav, 2);
        Z_m_bounds_fromSA = zeros(Reach_num_concav, 2);
        Concav_hsv_ref = hsv(Reach_num_concav);
        
        for D = 1:1:Reach_num_concav
            
            % First, obtain the lower drainage area
            
            title({'Click Once to Select the Lowest',['Drainage Area For Reach #' num2str(D)]}, ...
                'FontSize', 16)
            
            [DA_bound_ginput, S_bound_ginput] = ginput(1);
            
            % Find the minimum distance between the clicked location and the data
            % points
            
            DA_sort_clip = DA_sort(S_sort ~= 0 & isnan(S_sort) == 0);
            S_sort_clip = S_sort(S_sort ~= 0 & isnan(S_sort) == 0);
            L_sort_clip = L_sort(S_sort ~= 0 & isnan(S_sort) == 0);
            Z_sort_clip = Z_sort(S_sort ~= 0 & isnan(S_sort) == 0);
            
            % Here, I transform the slope-area data so that distance from all
            % points can be calculated (i.e., distance isn't the same in a log-log
            % plot).
            DA_diff = abs(log(DA_sort_clip) - log(DA_bound_ginput(1,1)));
            S_diff = abs(log(S_sort_clip) - log(S_bound_ginput(1,1)));
            Dist = ((DA_diff .^ 2) + (S_diff .^ 2)) .^ 0.5;
            
            DA_m2_bounds_fromSA(D,1) = DA_sort_clip(Dist == min(Dist));
            S_bounds_fromSA(D,1) = S_sort_clip(Dist == min(Dist));
            L_m_bounds_fromSA(D,1) = L_sort_clip(Dist == min(Dist));
            Z_m_bounds_fromSA(D,1) = Z_sort_clip(Dist == min(Dist));
            
            % Now plot the results
            
            figure(4)
            
            loglog(DA_m2_bounds_fromSA(D,1), S_bounds_fromSA(D,1), 'marker', 'o', 'color', 'k', ...
                'markerfacecolor', Concav_hsv_ref(D,:), 'linewidth', 0.5, 'linestyle', 'none')
            hold on
            
            loglog([DA_m2_bounds_fromSA(D,1), DA_m2_bounds_fromSA(D,1)], [min(S_sort_clip), max(S_sort_clip)], 'color', ...
                Concav_hsv_ref(D,:), 'linewidth', 1, 'linestyle', '-')
            hold on
            
            % Now, obtain the higher drainage area
            
            title({'Click Once to Select the Highest',['Drainage Area For Reach #' num2str(D)]}, ...
                'FontSize', 16)
            
            [DA_bound_ginput, S_bound_ginput] = ginput(1);
            
            % Find the minimum distance between the clicked location and the data
            % points
            DA_sort_clip = DA_sort(S_sort ~= 0);
            S_sort_clip = S_sort(S_sort ~= 0);
            L_sort_clip = L_sort(S_sort ~= 0);
            Z_sort_clip = Z_sort(S_sort ~= 0);
            
            % Here, I transform the slope-area data so that distance from all
            % points can be calculated (i.e., distance isn't the same in a log-log
            % plot).
            DA_diff = abs(log(DA_sort_clip) - log(DA_bound_ginput(1,1)));
            S_diff = abs(log(S_sort_clip) - log(S_bound_ginput(1,1)));
            Dist = ((DA_diff .^ 2) + (S_diff .^ 2)) .^ 0.5;
            
            DA_m2_bounds_fromSA(D,2) = DA_sort_clip(Dist == min(Dist));
            S_bounds_fromSA(D,2) = S_sort_clip(Dist == min(Dist));
            L_m_bounds_fromSA(D,2) = L_sort_clip(Dist == min(Dist));
            Z_m_bounds_fromSA(D,2) = Z_sort_clip(Dist == min(Dist));
            
            % Now plot the results
            
            figure(4)
            
            loglog(DA_m2_bounds_fromSA(D,2), S_bounds_fromSA(D,2), 'marker', 'o', 'color', 'k', ...
                'markerfacecolor', Concav_hsv_ref(D,:), 'linewidth', 0.5, 'linestyle', 'none')
            hold on
            
            loglog([DA_m2_bounds_fromSA(D,2), DA_m2_bounds_fromSA(D,2)], [min(S_sort_clip), max(S_sort_clip)], 'color', ...
                Concav_hsv_ref(D,:), 'linewidth', 1, 'linestyle', '-')
            hold on
            
            yaxis_limits = ylim;
            
            % I'm writing the reach labels at these slope values to try and
            % avoid the text being on top of the data. Might not work
            % perfectly...
            text_y_val = logspace(log10(yaxis_limits(1,1)), log10(yaxis_limits(1,2)), 10);
            text_y_val = text_y_val(1,end-1);
            
            h = text(DA_m2_bounds_fromSA(D,1), text_y_val, ['#' num2str(D)]);
            set(h, 'fontsize', 12)
            set(h, 'color', Concav_hsv_ref(D,:))
            set(h, 'rotation', -20)
            
        end
        
        % I pause here so you can see the last point go down
        pause(0.1)
        
        clear DA_bound_ginput S_bound_ginput DA_diff S_diff DA_sort_clip ...
            S_sort_clip DA_bound_response Correct_response yaxis_limits text_y_val ...
            L_sort_clip Z_sort_clip
        
        figure(4)
        
        hold off
        
        loglog(DA_sort(S_sort ~= 0), S_sort(S_sort ~= 0), 'ko', 'linewidth', 1, ...
            'linestyle', 'none')
        hold on
        
        % These variables store the results from slope-area methods, with
        % reaches selected from a slope-area plot
        Reach_theta_fromSA = zeros(Reach_num_concav, 1);
        Reach_ks_fromSA = zeros(Reach_num_concav, 1);
        
        % Here, I calculate the best-fit concavity value using the chi method
        % of Perron and Royden (2013).
        min_concav_range = 0;
        delta_concav = 0.01;
        max_concav_range = 2;
        concav_range = min_concav_range:delta_concav:max_concav_range;
        
        % A0 is a reference drainage area. See Perron and Royden (2013).
        A0 = 1E+06;
        
        % These variables store the results from Chi methods, but with
        % reaches selected from a slope-area plot
        Reach_theta_fromChi = zeros(Reach_num_concav, 1);
        Reach_ks_fromChi = zeros(Reach_num_concav, 1);
        
        for D = 1:1:Reach_num_concav
            
            DA_m2_temp = DA_sort(DA_sort >= DA_m2_bounds_fromSA(D,1) & DA_sort <= DA_m2_bounds_fromSA(D,2));
            
            S_temp = S_sort(DA_sort >= DA_m2_bounds_fromSA(D,1) & DA_sort <= DA_m2_bounds_fromSA(D,2));
            
            Z_sort_temp = Z_sort(DA_sort >= DA_m2_bounds_fromSA(D,1) & DA_sort <= DA_m2_bounds_fromSA(D,2)) - ...
                min(Z_sort(DA_sort >= DA_m2_bounds_fromSA(D,1) & DA_sort <= DA_m2_bounds_fromSA(D,2)));
            
            %% CALCULATE CONCAVITY WITH SLOPE-AREA METHODS
            
            % First, slope-area methods
            SA_fit = fitlm(log(DA_m2_temp(S_temp ~= 0 & isnan(S_temp) == 0)), log(S_temp(S_temp ~= 0 & isnan(S_temp) == 0)));
            Coeffs_temp = table2array(SA_fit.Coefficients);
            Reach_theta_fromSA(D,1) = -Coeffs_temp(2,1);
            Reach_ks_fromSA(D,1) = exp(Coeffs_temp(1,1));
            
            figure(4)
            
            loglog([min(DA_m2_temp), max(DA_m2_temp)], Reach_ks_fromSA(D,1) .* ([min(DA_m2_temp), max(DA_m2_temp)] .^ ...
                -Reach_theta_fromSA(D,1)), 'color', Concav_hsv_ref(D,:), 'linestyle', '--', 'linewidth', 2)
            hold on
            
            yaxis_limits = ylim;
            
            % I'm writing the reach labels at these slope values to try and
            % avoid the text being on top of the data. Might not work
            % perfectly...
            text_y_val = logspace(log10(yaxis_limits(1,1)), log10(yaxis_limits(1,2)), 10);
            text_y_val = text_y_val(1,end-1);
            
            h = text(min(DA_m2_temp), text_y_val, ['k_{s ' num2str(D) '}: ' num2str(Reach_ks_fromSA(D,1), ...
                '%.2E') ' m^{' num2str(2 * Reach_theta_fromSA(D,1), '%.2f') '}']);
            set(h, 'fontsize', 10)
            set(h, 'color', Concav_hsv_ref(D,:))
            set(h, 'rotation', -20)
            
            %% CALCULATE CONCAVITY WITH CHI METHODS
            
            % Before calculating the best-fit concavity, I calculate the average dx
            % value (i.e., with a 10 meter DEM used, the flow length can be 10 m or
            % 14.1421 m across a cell).
            dx_all = zeros(length(L_sort), 1);
            
            for i = 1:1:length(L_sort)
                
                if i == length(L_sort)
                    
                    dx_all(i,1) = L_sort(i,1);
                    
                else
                    
                    dx_all(i,1) = (L_sort(i,1) - L_sort(i+1,1));
                    
                end
                
            end
            
            dx_avg = mean(dx_all);
            
            Chi_temp = zeros(length(DA_m2_temp), 1);
            R2_Chi = zeros(1, length(concav_range));
            Reach_ks_Chi_temp = zeros(1, length(concav_range));
            
            for r = 1:1:length(concav_range)
                
                Chi_partial = (A0 ./ DA_m2_temp) .^ (concav_range(1,r));
                
                for i = length(Chi_temp):-1:1
                    
                    Chi_temp(i,1) = sum(Chi_partial(i:end,1)) * (dx_avg);
                    
                end
                
                Chi_fitlm = fitlm(Chi_temp, Z_sort_temp, 'intercept', false);
                
                Coeffs_temp = table2array(Chi_fitlm.Coefficients);
                
                Reach_ks_Chi_temp(1,r) = Coeffs_temp(1,1) * (A0 ^ concav_range(1,r));
                
                R2_Chi(1,r) = Chi_fitlm.Rsquared.Adjusted;
                
            end
            
            if max(R2_Chi) > 0 && length(max(R2_Chi)) == 1
                
                Reach_theta_fromChi(D,1) = concav_range(R2_Chi == max(R2_Chi));
                
                Reach_ks_fromChi(D,1) = Reach_ks_Chi_temp(R2_Chi == max(R2_Chi));
                
                Chi_partial = (A0 ./ DA_m2_temp) .^ (concav_range(R2_Chi == max(R2_Chi)));
                
                for i = length(Chi_temp):-1:1
                    
                    Chi_temp(i,1) = sum(Chi_partial(i:end,1)) * dx_avg;
                    
                end
                
            elseif max(R2_Chi) > 0 && length(max(R2_Chi)) > 1
                
                % This is a problematic situation, in which the R^2 values
                % level off so there is no single maximum value. If this
                % happens, you shouldn't use the value selected! Look at
                % figure 6 to assess R^2 values.
                Reach_theta_fromChi(D,1) = concav_range(concav_range == max(concav_range(R2_Chi == max(R2_Chi))));
                
                Reach_ks_fromChi(D,1) = Reach_ks_Chi_temp(concav_range == max(concav_range(R2_Chi == max(R2_Chi))));
                
                Chi_partial = (A0 ./ DA_m2_temp) .^ (concav_range(concav_range == max(concav_range(R2_Chi == max(R2_Chi)))));
                
                for i = length(Chi_temp):-1:1
                    
                    Chi_temp(i,1) = sum(Chi_partial(i:end,1)) * dx_avg;
                    
                end
                
            else
                
                beep
                
                disp('Error in obtaining a concavity through the integral method, check input.')
                
            end
            
            % This is a figure showing the R^2 values for predicting the observed Z
            % values with a given theta value. The theta value selected is the one
            % with the highest R^2 value.
            
            figure(6)
            
            plot(concav_range, R2_Chi, 'color', Concav_hsv_ref(D,:), 'linewidth', 1)
            hold on
            
            plot(Reach_theta_fromChi(D,1), max(R2_Chi), 'marker', 'o', 'color', 'k', 'markerfacecolor', ...
                Concav_hsv_ref(D,:), 'linewidth', 0.5, 'linestyle', 'none')
            hold on
            
        end
        
        clear R2_Chi r min_concav_range max_concav_range delta_concav concav_range ...
            Chi_temp dx_avg Reach_ks_Chi_temp
        
        %% ADD LABELS TO GRAPHS
        
        for D = 1:1:Reach_num_concav
            
            if D == 1
                
                title_string_SA = ['\theta_{1}: ' num2str(Reach_theta_fromSA(D,1), '%.2f') ', '];
                
                title_string_Chi = ['\theta_{1}: ' num2str(Reach_theta_fromChi(D,1), '%.2f') ', '];
                
            elseif D ~= 1 && D ~= Reach_num_concav
                
                title_string_SA = [title_string_SA '\theta_' num2str(D) ':' num2str(Reach_theta_fromSA(D,1), '%.2f') ', '];
                
                title_string_Chi = [title_string_Chi ' \theta_' num2str(D) ':' num2str(Reach_theta_fromChi(D,1), '%.2f') ', '];
                
            elseif D == Reach_num_concav
                
                title_string_SA = [title_string_SA ' \theta_' num2str(D) ':' num2str(Reach_theta_fromSA(D,1), '%.2f') ];
                
                title_string_Chi = [title_string_Chi ' \theta_' num2str(D) ':' num2str(Reach_theta_fromChi(D,1), '%.2f') ];
                
            end
            
            % I add ks values to figure 6 here, b/c the limits of the graph
            % can be changing in the for loop above.
            figure(6)
            
            xaxis_limits = xlim;
            yaxis_limits = ylim;
            
            h = text(xaxis_limits(1,1) + ((xaxis_limits(1,2) - xaxis_limits(1,1)) / 20), yaxis_limits(1,1) + ...
                (D * (yaxis_limits(1,2) - yaxis_limits(1,1)) / 15), ['k_{s ' num2str(D) '}: ' num2str(Reach_ks_fromChi(D,1), ...
                '%.2E') ' m^{' num2str(2 * Reach_theta_fromChi(D,1), '%.2f') '}']);
            set(h, 'fontsize', 10)
            set(h, 'color', Concav_hsv_ref(D,:))
            
        end
        
        clear xaxis_limits yaxis_limits
        
        figure(4)
        
        hold off
        
        set(gca, 'fontsize', 12)
        xlabel('Drainage Area (m^2)', 'fontsize', 14)
        ylabel('Slope', 'fontsize', 14)
        
        title({['Slope-Area Data for ' char(Basin_name_figures) ','], ...
            char(title_string_SA)}, 'fontsize', 16)
        
        figure(6)
        
        hold off
        
        set(gca, 'fontsize', 12)
        xlabel('\theta', 'fontsize', 16)
        ylabel('R^2', 'fontsize', 14)
        title({['\chi Plot for ' char(Basin_name_figures) ','], ...
            char(title_string_Chi)}, 'fontsize', 12)
        
        set(figure(6), 'Position', [25 105 550 430])
        
        % Ask if the user wants to redo selection
        prompt = 'Reselect the reaches (''y'' or ''n'')? ';
        
        prompt_title = 'Reselect?';
        
        Reselect_response = inputdlg(prompt, prompt_title);
        
        Negative_response = 'n';
        
        if strcmp(char(Reselect_response), Negative_response)
            
            continue_test = 1;
            
        end
        
    end
    
end

if use_ref_theta ~= 1
    
    figure(4)
    
    title({['Slope-Area Data for ' char(Basin_name_figures) ','], ...
        char(title_string_SA)}, 'fontsize', 12)
    
    set(figure(4), 'Position', [25 160 550 430])
    
    clear Reselect_response continue_test continue_test_steady_state Negative_response ...
        DA_m2_temp S_temp SA_fit
    
end

%% CHOOSE CONCAVITY

if use_ref_theta == 0
    
    continue_test = 0;
    
    while continue_test == 0
        
        prompt = ['Enter the reach number corresponding with the concavity you want to use. ' ...
            'Otherwise, you can enter (''e'') a specific concavity value. '];
        
        prompt_title = 'Choose Concavity';
        
        Concavity_response = inputdlg(prompt, prompt_title);
        
        Correct_response_1 = 'e';
        
        if strcmp(char(Concavity_response), Correct_response_1)
            
            prompt = 'Enter the desired concavity value. ';
            
            prompt_title = 'Enter Concavity';
            
            theta_response = inputdlg(prompt, prompt_title);
            
            theta_response_char = cell2mat(theta_response);
            
            theta_ref = str2double(theta_response_char);
            
            continue_test = 1;
            
            clear theta_response theta_response_char
            
        else
            
            Reach_num_char = cell2mat(Concavity_response);
            
            Reach_num_selected = str2double(Reach_num_char);
            
            if mod(Reach_num_selected, 1) == 0 && Reach_num_selected > 0 && Reach_num_selected <= Reach_num_concav
                
                prompt = 'Use the value from slope-area (''s'') or slope-integral (''x'') methods? ';
                
                prompt_title = 'Choose Concavity';
                
                Concavity_response = inputdlg(prompt, prompt_title);
                
                Correct_response_1 = 's';
                
                Correct_response_2 = 'x';
                
                if strcmp(char(Concavity_response), Correct_response_1)
                    
                    theta_ref = Reach_theta_fromSA(Reach_num_selected, 1);
                    
                    continue_test = 1;
                    
                elseif strcmp(char(Concavity_response), Correct_response_2)
                    
                    theta_ref = Reach_theta_fromChi(Reach_num_selected, 1);
                    
                    continue_test = 1;
                    
                end
                
            end
            
            clear Reach_num_char Reach_num_selected
            
        end
        
    end
    
    clear Concavity_response Correct_response_1
    
end

%% CALCULATE CHI VALUES GIVEN CALCULATED OR ENFORCED CONCAVITY

% Uses the approach of Perron and Royden (2013) to transform distance
% upstream using drainage area values. A steady-state reach with uniform
% properties should have a linear relationship between elevation and Chi.

% First, determine the average dx along the profile. Although this is
% inluded above, it is here in case the use_ref_concav option was set to 1.
dx_all = zeros(length(L_sort), 1);

for i = 1:1:length(L_sort)
    
    if i == length(L_sort)
        
        dx_all(i,1) = L_sort(i,1);
        
    else
        
        dx_all(i,1) = (L_sort(i,1) - L_sort(i+1,1));
        
    end
    
end

dx_avg = mean(dx_all);

clear dx_all

% Now, calculate chi values

A0 = 1E+6;

Chi_partial = (A0 ./ DA_sort) .^ (theta_ref);

Chi = zeros(length(DA_sort),1);

for i = length(Chi):-1:1
    
    Chi(i,1) = sum(Chi_partial(i:end,1)) * dx_avg;
    
end

clear Chi_partial dx_avg

%% REQUESTS USER INPUT FOR ASSESSING THE CHI PLOT AND ADJUSTING THE REACHES

input_end_test = 0;

first_time_check = 0;

second_time_check = 0;

while input_end_test == 0
    
    if first_time_check ~= 0
        
        Affirmative_response = 'y';
        
        Negative_response = 'n';
        
        if second_time_check == 1
            
            prompt = ['Would you like to change the bounds of the selected ' ...
                'reaches (''y'' or ''n'')? '];
            
            prompt_title = 'Assess Selection(s)';
            
            User_input = inputdlg(prompt, prompt_title);
            
        else
            
            User_input = Affirmative_response;
            
        end
        
        % I allow the user to select multiple reaches on the Chi plot.
        % Here, the user specifies how many they want to select. Redundant,
        % b/c this is done above, but the slope-area and chi plot are
        % different (chi plot is usually more clear), and if you use a
        % reference concavity then the slope-areaabove is skipped.
        
        if strcmp(char(User_input), Affirmative_response) == 1
            
            Reach_num_continue = 0;
            
            while Reach_num_continue == 0
                
                prompt = 'Specify the number of reaches you would like to assess. ';
                
                prompt_title = 'Number of Reaches';
                
                Reach_num_input = inputdlg(prompt, prompt_title);
                
                Reach_num_char = cell2mat(Reach_num_input);
                
                Reach_num_ksn = str2double(Reach_num_char);
                
                % The number of reaches has to be a positive integer,
                % otherwise the code won't continue here.
                
                if mod(Reach_num_ksn, 1) == 0 && Reach_num_ksn > 0
                    
                    if Reach_num_ksn == 2
                        
                        ksn_jet_ref(1,:) = [0 0 1];
                        
                        ksn_jet_ref(2,:) = [1 0 0];
                        
                    else
                        
                        ksn_jet_ref = jet(Reach_num_ksn);
                        
                    end
                    
                    % There store the boundaries of the selected reach(es).
                    % First column is lower boundary, second is higher.
                    Chi_m_bounds_fromChiplot = zeros(Reach_num_ksn, 2);
                    Z_m_bounds_fromChiplot = zeros(Reach_num_ksn, 2);
                    L_m_bounds_fromChiplot = zeros(Reach_num_ksn, 2);
                    DA_m2_bounds_fromChiplot = zeros(Reach_num_ksn, 2);
                    S_bounds_fromChiplot = zeros(Reach_num_ksn, 2);
                    Reach_ks_Chi_fromChiplot = zeros(Reach_num_ksn, 1);
                    Reach_ks_SA_fromChiplot = zeros(Reach_num_ksn, 1);
                    
                    Reach_num_continue = 1;
                    
                end
                
            end
            
            clear Reach_num_input Reach_num_continue Reach_num_char
            
            % Now actually select the reach(es)
            
            figure(5)
            
            plot(Chi(Kw_sort == min(Kw_sort)), Z_sort(Kw_sort == min(Kw_sort)), ...
                'linewidth', 2, 'color', 'k')
            hold on
            
            plot(Chi(Kw_sort == max(Kw_sort)), Z_sort(Kw_sort == max(Kw_sort)), ...
                'linewidth', 3, 'color', [0.5, 0.5, 0.5])
            hold on
            
            set(gca, 'fontsize', 12)
            xlabel('\chi (m)', 'fontsize', 16)
            ylabel('Elevation (m)', 'fontsize', 16)
            
            for D = 1:1:Reach_num_ksn
                
                figure(5)
                
                title({'Click Once to Select the Lower',['Boundary of Reach #' num2str(D)]}, ...
                    'fontsize', 16)
                
                % Here, I use ginput to have the user select two points for the
                % lower and upper boundaris of each reach
                
                % First, lower boundary of the reach
                
                [Chi_bound_lower_ginput, Z_bound_lower_ginput] = ginput(1);
                
                % Find the closest point to the click for the lower boundary
                
                Chi_diff = abs(Chi - Chi_bound_lower_ginput(1,1));
                
                Z_diff = abs(Z_sort - Z_bound_lower_ginput(1,1));
                
                Chi_dist = ((Chi_diff .^ 2) + (Z_diff .^ 2)) .^ 0.5;
                
                Chi_m_bounds_fromChiplot(D,1) = Chi(Chi_dist == min(Chi_dist));
                Z_m_bounds_fromChiplot(D,1) = Z_sort(Chi_dist == min(Chi_dist));
                L_m_bounds_fromChiplot(D,1) = L_sort(Chi_dist == min(Chi_dist));
                DA_m2_bounds_fromChiplot(D,1) = DA_sort(Chi_dist == min(Chi_dist));
                S_bounds_fromChiplot(D,1) = S_sort(Chi_dist == min(Chi_dist));
                
                % Plot it
                
                plot(Chi_m_bounds_fromChiplot(D,1), Z_m_bounds_fromChiplot(D,1), 'color', ksn_jet_ref(D,:), 'marker', 'o', ...
                    'linewidth', 2, 'markersize', 8)
                hold on
                
                % Find the upper boundary
                
                title({'Click Once to Select the Upper',['Boundary of Reach #' num2str(D)]}, ...
                    'fontsize', 16)
                
                [Chi_bound_upper_ginput, Z_bound_upper_ginput] = ginput(1);
                
                % Find the closest point to the click for the upper boundary
                
                Chi_diff = abs(Chi - Chi_bound_upper_ginput(1,1));
                
                Z_diff = abs(Z_sort - Z_bound_upper_ginput(1,1));
                
                Chi_dist = ((Chi_diff .^ 2) + (Z_diff .^ 2)) .^ 0.5;
                
                Chi_m_bounds_fromChiplot(D,2) = Chi(Chi_dist == min(Chi_dist));
                Z_m_bounds_fromChiplot(D,2) = Z_sort(Chi_dist == min(Chi_dist));
                L_m_bounds_fromChiplot(D,2) = L_sort(Chi_dist == min(Chi_dist));
                DA_m2_bounds_fromChiplot(D,2) = DA_sort(Chi_dist == min(Chi_dist));
                S_bounds_fromChiplot(D,2) = S_sort(Chi_dist == min(Chi_dist));
                
                plot(Chi_m_bounds_fromChiplot(D,2), Z_m_bounds_fromChiplot(D,2), 'color', ksn_jet_ref(D,:), 'marker', 'o', ...
                    'linewidth', 2, 'markersize', 8)
                hold on
                
            end
            
            % I pause here so you can see the last point go down
            pause(0.1)
            
            hold off
            
            second_time_check = 1;
            
        elseif strcmp(char(User_input), Negative_response) == 1
            
            input_end_test = 1;
            
            break
            
        end
        
    end
    
    clear continue_test Chi_bound_lower_ginput Z_bound_lower_ginput Chi_bound_upper_ginput Z_bound_upper_ginput
    
    %% MAKES A CHI PLOT, SLOPE-AREA PLOT, AND GETS THE STEEPENSS FOR EACH REACH
    
    if first_time_check == 1 && input_end_test ~= 1
        
        figure(7)
        
        loglog(DA_sort(S_sort ~= 0), S_sort(S_sort ~= 0), 'ko', 'linewidth', 1, 'linestyle', 'none')
        hold on
        
        figure(5)
        
        hold off
        
        plot(Chi(Kw_sort == min(Kw_sort)), Z_sort(Kw_sort == min(Kw_sort)), ...
            'linewidth', 2, 'color', 'k')
        hold on
        
        plot(Chi(Kw_sort == max(Kw_sort)), Z_sort(Kw_sort == max(Kw_sort)), ...
            'linewidth', 3, 'color', [0.5, 0.5, 0.5])
        hold on
        
        for D = 1:1:Reach_num_ksn
            
            Z_sort_temp = Z_sort .* NaN;
            
            Z_sort_temp(Chi >= Chi_m_bounds_fromChiplot(D,1) & Chi <= Chi_m_bounds_fromChiplot(D,2)) = ...
                Z_sort(Chi >= Chi_m_bounds_fromChiplot(D,1) & Chi <= Chi_m_bounds_fromChiplot(D,2));
            
            % Here, I create a linear regression of the reach's elevations
            % with Chi. The resulting slope reflects steepness.
            Chi_regression = fitlm(Chi - min(Chi(isnan(Z_sort_temp) == 0)), ...
                Z_sort_temp - min(Z_sort_temp), 'linear', 'intercept', false);
            
            Coeffs_temp = table2array(Chi_regression.Coefficients);
            
            Reach_ks_Chi_fromChiplot(D,1) = Coeffs_temp(1,1) * (A0 ^ theta_ref);
            
            if isreal(Reach_ks_Chi_fromChiplot(D,1)) ~= 1
                
                Reach_ks_Chi_fromChiplot(D,1) = 0;
                
                beep
                
                disp('Selection resulted in an imaginary steepness value, check selection. ')
                
            end
            
            % Now, I calculate ks with slope-area methods applied to the entire
            % reach (as opposed to at-a-point).
            DA_m2_temp = DA_sort(Chi >= Chi_m_bounds_fromChiplot(D,1) & Chi <= Chi_m_bounds_fromChiplot(D,2));
            S_temp = S_sort(Chi >= Chi_m_bounds_fromChiplot(D,1) & Chi <= Chi_m_bounds_fromChiplot(D,2));
            
            Reach_ks_SA_fromChiplot(D,1) = exp(mean(log(S_temp(isnan(S_temp) == 0 & S_temp ~= 0)) - ...
                (-theta_ref .* log(DA_m2_temp(isnan(S_temp) == 0 & S_temp ~= 0)))));
            
            if isreal(Reach_ks_SA_fromChiplot(D,1)) ~= 1
                
                Reach_ks_SA_fromChiplot(D,1) = 0;
                
                beep
                
                disp('Selection resulted in an imaginary steepness value, check selection. ')
                
            end
            
            % Add reach to the Slope-Area plot
            
            figure(7)
            
            loglog([min(DA_m2_temp(S_temp ~= 0)), max(DA_m2_temp(S_temp ~= 0))], Reach_ks_SA_fromChiplot(D,1) .* ...
                ([min(DA_m2_temp(S_temp ~= 0)), max(DA_m2_temp(S_temp ~= 0))] .^ -theta_ref), 'color', ...
                ksn_jet_ref(D,:), 'linewidth', 2, 'linestyle', '--')
            hold on
            
            yaxis_limits = ylim;
            
            % I'm writing the reach labels at these slope values to try and
            % avoid the text being on top of the data. Might not work
            % perfectly...
            text_y_val = logspace(log10(yaxis_limits(1,1)), log10(yaxis_limits(1,2)), 10);
            text_y_val = text_y_val(1,end-1);
            
            h = text(min(DA_m2_temp(S_temp ~= 0)), text_y_val, ['k_{sn ' num2str(D) '}: ' num2str(Reach_ks_SA_fromChiplot(D,1), ...
                '%.2E') ' m^{' num2str(2 * theta_ref, '%.2f') '}']);
            set(h, 'fontsize', 10)
            set(h, 'color', ksn_jet_ref(D,:))
            set(h, 'rotation', -20)
            
            % Add reach to the Chi plot
            
            figure(5)
            
            plot(Chi(isnan(Z_sort_temp) == 0), ((Reach_ks_Chi_fromChiplot(D,1) / (A0 ^ theta_ref)) .* Chi(isnan(Z_sort_temp) == 0)) ...
                + min(Z_sort_temp) - min((Reach_ks_Chi_fromChiplot(D,1) / (A0 ^ theta_ref)) .* Chi(isnan(Z_sort_temp) == 0)), ...
                'color', ksn_jet_ref(D,:), 'linewidth', 2, 'linestyle', '--')
            hold on
            
            xaxis_limits = xlim;
            yaxis_limits = ylim;
            
            if min(Z_sort_temp(isnan(Z_sort_temp) == 0)) >= ( 4 * (yaxis_limits(1,2) - yaxis_limits(1,1)) / 5)
                
                h = text(min(Chi(isnan(Z_sort_temp) == 0)), min(Z_sort_temp(isnan(Z_sort_temp) == 0)) - ...
                    ((yaxis_limits(1,2) - yaxis_limits(1,1)) / 10), ['k_{sn ' num2str(D) '}: ' num2str(Reach_ks_Chi_fromChiplot(D,1), ...
                    '%.2E') ' m^{' num2str(2 * theta_ref, '%.2f') '}']);
                set(h, 'fontsize', 10)
                set(h, 'color', ksn_jet_ref(D,:))
                set(h, 'rotation', 20)
                
            else
                
                h = text(min(Chi(isnan(Z_sort_temp) == 0)) + ((xaxis_limits(1,2) - xaxis_limits(1,1)) / 10), min(Z_sort_temp(isnan(Z_sort_temp) == 0)) + ...
                    ((yaxis_limits(1,2) - yaxis_limits(1,1)) / 15), ['k_{sn ' num2str(D) '}: ' num2str(Reach_ks_Chi_fromChiplot(D,1), ...
                    '%.2E') ' m^{' num2str(2 * theta_ref, '%.2f') '}']);
                set(h, 'fontsize', 10)
                set(h, 'color', ksn_jet_ref(D,:))
                set(h, 'rotation', 20)
                
            end
            
            clear h
            
            if D == 1
                
                title_string_SA = ['k_{sn 1}: ' num2str(Reach_ks_SA_fromChiplot(D,1), '%.2E') ...
                    ' m^{' num2str(2 * theta_ref, '%.2f') '}, '];
                
                title_string_Chi = ['k_{sn 1}: ' num2str(Reach_ks_Chi_fromChiplot(D,1), '%.2E') ...
                    ' m^{' num2str(2 * theta_ref, '%.2f') '}, '];
                
            elseif D ~= 1 && D ~= Reach_num_ksn
                
                title_string_SA = [title_string_SA 'k_{sn ' num2str(D) '}: ' num2str(Reach_ks_SA_fromChiplot(D,1), '%.2E') ...
                    ' m^{' num2str(2 * theta_ref, '%.2f') '}, '];
                
                title_string_Chi = [title_string_Chi 'k_{sn ' num2str(D) '}: ' num2str(Reach_ks_Chi_fromChiplot(D,1), '%.2E') ...
                    ' m^{' num2str(2 * theta_ref, '%.2f') '}, '];
                
            elseif D == Reach_num_ksn
                
                title_string_SA = [title_string_SA 'k_{sn ' num2str(D) '}: ' num2str(Reach_ks_SA_fromChiplot(D,1), '%.2E') ...
                    ' m^{' num2str(2 * theta_ref, '%.2f') '} '];
                
                title_string_Chi = [title_string_Chi 'k_{sn ' num2str(D) '}: ' num2str(Reach_ks_Chi_fromChiplot(D,1), '%.2E') ...
                    ' m^{' num2str(2 * theta_ref, '%.2f') '} '];
                
            end
            
        end
        
        clear xaxis_limits yaxis_limits text_y_val
        
        figure(7)
        
        hold off
        
        set(gca, 'fontsize', 12)
        % Hopefully title_string isn't too long
        title({['Slope-Area Data for ' char({Basin_name_figures}) ' with \theta_{ref} = ' num2str(theta_ref, '%.2f') ...
            ','], char(title_string_SA)}, 'fontsize', 10)
        
        xlabel('Drainage Area (m^2)', 'fontsize', 14)
        ylabel('Slope', 'fontsize', 14)
        
        set(figure(7), 'Position', [600 55 550 430])
        
        drawnow
        
        figure(5)
        
        hold off
        
        set(gca, 'fontsize', 12)
        
        % Hopefully title_string isn't too long
        title({['\chi Plot for ' char({Basin_name_figures}) ' with \theta_{ref} = ' num2str(theta_ref, '%.2f') ...
            ','], char(title_string_Chi)}, 'fontsize', 10)
        
        xlabel('\chi (m)', 'fontsize', 14)
        ylabel('Elevation ASL (m)', 'fontsize', 14)
        
        set(figure(5), 'Position', [25 55 550 430])
        
        drawnow
        
        clear Chi_regression Coeffs_temp Z_sort_temp DA_m2_temp S_temp title_string_SA title_string_Chi Chi_fitlm
        
    end
    
    %% MAKES CHI PLOT FOR THE FIRST ITERATION THROUGH THE WHILE LOOP
    
    if first_time_check ~= 1 && input_end_test ~= 1
        
        figure(5)
        
        plot(Chi(Kw_sort == min(Kw_sort)), Z_sort(Kw_sort == min(Kw_sort)), ...
            'linewidth', 2, 'color', 'k')
        hold on
        
        plot(Chi(Kw_sort == max(Kw_sort)), Z_sort(Kw_sort == max(Kw_sort)), ...
            'linewidth', 3, 'color', [0.5, 0.5, 0.5])
        hold on
        
        set(gca, 'fontsize', 12)
        xlabel('\chi (m)', 'fontsize', 16)
        ylabel('Elevation ASL (m)', 'fontsize', 16)
        title({['\chi Plot for ' char({Basin_name_figures})]}, 'fontsize', 18)
        
        set(figure(5), 'Position', [25 55 550 430])
        
    end
    
    % The first time through this loop, first_time_check is still set to
    % zero. At the end of the first loop, it is set to one instead,
    % enabling many of the if statements that were ignored the first time.
    first_time_check = 1;
    
end

clear User_input Chi_diff Z_diff Chi_dist Affirmative_response Negative_response ...
    first_time_check second_time_check input_end_test Correct_response1 Correct_response2

%% PROJECTS RELICT AND ADJUSTED ELEVATIONS ON CHI PLOT

% This section is used to project the relict and adjusted elevations for a
% transient profile, in scenarios where transience is caused by an increase
% in rock-uplift rates. If you are not dealing with suh a scenario, then
% ignore this graph and the .txt output file relating to it!
if Reach_num_ksn == 2 && transience_projections == 1
    
    % The order of relict and adjusted, in terms of the reach numbers, is
    % not enforced here. You should, for example, be systematic in having
    % reach 1 be adjusted and reach 2 relict.
    
    if Reach_ks_Chi_fromChiplot(1,1) > Reach_ks_Chi_fromChiplot(2,1)
        
        Z_adj = ((Reach_ks_Chi_fromChiplot(1,1) / (A0 ^ theta_ref)) .* Chi) - ...
            min(((Reach_ks_Chi_fromChiplot(1,1) / (A0 ^ theta_ref)) .* Chi(Z_sort == Z_m_bounds_fromChiplot(1,1)))) + ...
            Z_m_bounds_fromChiplot(1,1);
        
        Z_adj(DA_sort > DA_m2_bounds_fromChiplot(1,1) | DA_sort < DA_m2_bounds_fromChiplot(1,2)) = 0;
        
        Z_rel = ((Reach_ks_Chi_fromChiplot(2,1) / (A0 ^ theta_ref)) .* Chi) - min(((Reach_ks_Chi_fromChiplot(2,1) / ...
            (A0 ^ theta_ref)) .* Chi(Z_sort >= Z_m_bounds_fromChiplot(2,1)))) + Z_m_bounds_fromChiplot(2,1);
        
        Z_rel(DA_sort > DA_m2_bounds_fromChiplot(2,1) | DA_sort < DA_m2_bounds_fromChiplot(2,2)) = 0;
        
        Z_rel_projected = ((Reach_ks_Chi_fromChiplot(2,1) / (A0 ^ theta_ref)) .* Chi) + Z_m_bounds_fromChiplot(2,1) - ...
            min(((Reach_ks_Chi_fromChiplot(2,1) / (A0 ^ theta_ref)) .* Chi(Z_sort >= Z_m_bounds_fromChiplot(2,1))));
        
        Z_rel_projected(DA_sort < DA_m2_bounds_fromChiplot(2,1)) = 0;
        
        Incision_projected = min(Z_rel_projected(Z_rel_projected ~= 0)) - Z_sort(Z_rel_projected == min(Z_rel_projected(Z_rel_projected ~= 0)));
        
    elseif Reach_ks_Chi_fromChiplot(1,1) < Reach_ks_Chi_fromChiplot(2,1)
        
        Z_adj = ((Reach_ks_Chi_fromChiplot(2,1) / (A0 ^ theta_ref)) .* Chi) - ...
            min(((Reach_ks_Chi_fromChiplot(1,1) / (A0 ^ theta_ref)) .* Chi(Z_sort == Z_m_bounds_fromChiplot(2,1)))) + ...
            Z_m_bounds_fromChiplot(2,1);
        
        Z_adj(DA_sort > DA_m2_bounds_fromChiplot(2,1) | DA_sort < DA_m2_bounds_fromChiplot(2,2)) = 0;
        
        Z_rel = ((Reach_ks_Chi_fromChiplot(1,1) / (A0 ^ theta_ref)) .* Chi) - min(((Reach_ks_Chi_fromChiplot(1,1) / ...
            (A0 ^ theta_ref)) .* Chi(Z_sort >= Z_m_bounds_fromChiplot(1,1)))) + Z_m_bounds_fromChiplot(1,1);
        
        Z_rel(DA_sort > DA_m2_bounds_fromChiplot(1,1) | DA_sort < DA_m2_bounds_fromChiplot(1,2)) = 0;
        
        Z_rel_projected = ((Reach_ks_Chi_fromChiplot(1,1) / (A0 ^ theta_ref)) .* Chi) + Z_m_bounds_fromChiplot(1,1) - ...
            min(((Reach_ks_Chi_fromChiplot(1,1) / (A0 ^ theta_ref)) .* Chi(Z_sort >= Z_m_bounds_fromChiplot(1,1))));
        
        Z_rel_projected(DA_sort < DA_m2_bounds_fromChiplot(1,1)) = 0;
        
        Incision_projected = min(Z_rel_projected(Z_rel_projected ~= 0)) - Z_sort(Z_rel_projected == min(Z_rel_projected(Z_rel_projected ~= 0)));
        
    end
    
    % Now Chi-Plot showing incision
    figure(8)
    
    plot(L_sort, Z_sort, 'linewidth', 2, 'color', 'k')
    hold on
    
    plot(L_sort(Z_rel ~= 0), Z_rel(Z_rel ~= 0), 'linewidth', 2, 'color', 'b', 'linestyle', '--')
    hold on
    
    plot(L_sort(Z_rel_projected ~= 0), Z_rel_projected(Z_rel_projected ~= 0), 'linewidth', 2, ...
        'color', 'b', 'linestyle', ':')
    hold on
    
    plot(L_sort(Z_adj ~= 0), Z_adj(Z_adj ~= 0), 'linewidth', 2, 'color', 'r', 'linestyle', '--')
    hold off
    
    set(gca, 'fontsize', 12)
    
    lgnd = legend('Observed', 'Relict', 'Projected Relict', 'Adjusted', 'location', 'southeast');
    lgnd.FontSize = 8;
    
    xlabel('Distance Upstream (m)', 'fontsize', 14)
    ylabel('Elevation ASL (m)', 'fontsize', 14)
    title({'Relict and Adjusted Reaches for', [char({Basin_name_figures}) ', Projected Incision: ' ...
        num2str(Incision_projected, '%.2E') ' m']}, 'fontsize', 12)
    
    set(figure(8), 'Position', [1175 55 550 430])
    
    % Now Chi-Plot showing incision
    figure(9)
    
    plot(Chi, Z_sort, 'linewidth', 2, 'color', 'k')
    hold on
    
    plot(Chi(Z_rel ~= 0), Z_rel(Z_rel ~= 0), 'linewidth', 2, 'color', 'b', 'linestyle', '--')
    hold on
    
    plot(Chi(Z_rel_projected ~= 0), Z_rel_projected(Z_rel_projected ~= 0), 'linewidth', 2, ...
        'color', 'b', 'linestyle', ':')
    hold on
    
    plot(Chi(Z_adj ~= 0), Z_adj(Z_adj ~= 0), 'linewidth', 2, 'color', 'r', 'linestyle', '--')
    hold off
    
    set(gca, 'fontsize', 12)
    
    lgnd = legend('Observed', 'Relict', 'Projected Relict', 'Adjusted', 'location', 'southeast');
    lgnd.FontSize = 8;
    
    xlabel('\chi (m)', 'fontsize', 14)
    ylabel('Elevation ASL (m)', 'fontsize', 14)
    title({'Relict and Adjusted Reaches for', [char({Basin_name_figures}) ', Projected Incision: ' ...
        num2str(Incision_projected, '%.2E') ' m']}, 'fontsize', 12)
    
    set(figure(9), 'Position', [1175 215 550 430])
    
elseif Reach_num_ksn ~= 2 || transience_projections ~= 1
    
    Incision_projected = 0;
    
end

clear lgnd D

%% CONSOLIDATE AND OUTPUT RESULTS

% Ask if the user wants to save the results

continue_test = 0;

while continue_test == 0
    
    prompt = 'Do you wish to save these results (''y'' or ''n'')?';
    
    prompt_title = 'Save Results?';
    
    save_response = inputdlg(prompt, prompt_title);
    
    Positive_response = 'y';
    
    Negative_response = 'n';
    
    if strcmp(char(save_response), Positive_response)
        
        write_output = 1;
        
        continue_test = 1;
        
    elseif strcmp(char(save_response), Negative_response)
        
        write_output = 0;
        
        continue_test = 1;
        
    end
    
end

clear Positive_response Negative_response prompt prompt_title save_response continue_test

% If the user does want to save results, a text file and shapefile are
% created.

if write_output == 1
    
    cd(Output_location)
    
    prompt = 'Enter a name for the output (do not include a file extension like .shp).';
    
    prompt_title = 'Enter Name';
    
    Output_name = inputdlg(prompt, prompt_title);
    
    % Create a new folder for just these output files. Otherwise, a
    % huge mess of output files becomes confusing.
    
    mkdir([char(Output_name) '-' date]);
    
    New_location = [Output_location '\' char(Output_name) '-' date];
    
    cd(New_location)
    
    % Shapefile output
    
    %% POLYLINE SHAPEFILE
    
    if make_shapefiles_option == 1
        
        % Grid of the reference concavity used.
        theta_ref_grid = ksn_grid;
        theta_ref_grid.Z(isnan(theta_ref_grid.Z) == 0) = theta_ref;
        
        % Grid of steepness values from slope-area methods applied to
        % entire reaches (as opposed to 'at-a-point' calculations).
        SA_ks_grid = ksn_grid;
        
        % First, set all non NaN values to zero
        SA_ks_grid.Z(isnan(SA_ks_grid.Z) == 0) = 0;
        
        % Grid of steepness values from Chi plot
        Chi_ks_grid = ksn_grid;
        
        % First, set all non NaN values to zero
        Chi_ks_grid.Z(isnan(Chi_ks_grid.Z) == 0) = 0;
        
        for D = 1:1:Reach_num_ksn
            
            % Now give the reaches included within the selected reach (on the Chi
            % plot) the steepness values from Slope-Area methods.
            SA_ks_grid.Z(isnan(SA_ks_grid.Z) == 0 & (FlwAcc.Z .* (FlwAcc.cellsize ^ 2)) >= ...
                DA_m2_bounds_fromChiplot(D,2) & (FlwAcc.Z .* (FlwAcc.cellsize ^ 2)) <= ...
                DA_m2_bounds_fromChiplot(D,1)) = Reach_ks_SA_fromChiplot(D,1);
            
            % Now give the reaches included within the selected reach (on the Chi
            % plot) the steepness values from Chi methods.
            Chi_ks_grid.Z(isnan(Chi_ks_grid.Z) == 0 & (FlwAcc.Z .* (FlwAcc.cellsize ^ 2)) >= ...
                DA_m2_bounds_fromChiplot(D,2) & (FlwAcc.Z .* (FlwAcc.cellsize ^ 2)) <= ...
                DA_m2_bounds_fromChiplot(D,1)) = Reach_ks_Chi_fromChiplot(D,1);
            
        end
        
        % Gridded matrix of steepness values from slope-integral methods applied to
        % entire reaches (as opposed to 'at-a-point' calculations).
        Chi_theta_grid = ksn_grid;
        
        % First, set all non NaN values to zero
        Chi_theta_grid.Z(isnan(Chi_theta_grid.Z) == 0) = 0;
        
        % Gridded matrix of concavity values from slope-area methods. Note
        % that the maximum value is a segment (see seglength) is used.
        SA_theta_grid = ksn_grid;
        
        % First, set all non NaN values to zero
        SA_theta_grid.Z(isnan(SA_theta_grid.Z) == 0) = 0;
        
        if use_ref_theta ~= 1
            
            % Now, measured concavities
            for D = 1:1:Reach_num_concav
                
                % Now give the reaches included within the selected reach (on the
                % SA plot) the convavity values from Chi methods.
                Chi_theta_grid.Z(isnan(Chi_theta_grid.Z) == 0 & (FlwAcc.Z .* (FlwAcc.cellsize ^ 2)) >= ...
                    DA_m2_bounds_fromSA(D,1) & (FlwAcc.Z .* (FlwAcc.cellsize ^ 2)) <= ...
                    DA_m2_bounds_fromSA(D,2)) = Reach_theta_fromChi(D,1);
                
                % Now give the reaches included within the selected reach (on the
                % SA plot) the convavity values from slope-area methods.
                SA_theta_grid.Z(isnan(SA_theta_grid.Z) == 0 & (FlwAcc.Z .* (FlwAcc.cellsize ^ 2)) >= ...
                    DA_m2_bounds_fromSA(D,1) & (FlwAcc.Z .* (FlwAcc.cellsize ^ 2)) <= ...
                    DA_m2_bounds_fromSA(D,2)) = Reach_theta_fromSA(D,1);
                
            end
            
        end
        
        % Drainage area grid
        DA_m2_grid = FlwAcc;
        DA_m2_grid.Z = DA_m2_grid.Z .* (FlwAcc.cellsize ^ 2);
        
        % Slope Grid
        Slope_grid = Grad;
        Slope_grid.Z = Slope_grid.Z ./ 100;
        
        % Distance Grid
        L_grid = FlwDst;
        L_grid.Z = L_grid.Z - L_adjustment;
        
        clear FlwAcc Grad FlwDst
        
        % Now create map structure, use to create an ArcGIS shapefile
        Map_Struct_Output = STREAMobj2mapstruct(Stream_mod, 'seglength', seglength, ...
            'attributes', {'L_m' L_grid @mean 'DA_m2' DA_m2_grid @mean 'Slope' Slope_grid @mean ...
            'ksn_pnt' ksn_grid @mean 'theta_ref' theta_ref_grid @mean ...
            'ksn_SA' SA_ks_grid @max 'ksn_Chi' Chi_ks_grid @max 'theta_SA' SA_theta_grid @max ...
            'theta_Chi' Chi_theta_grid @max});
        
        shapewrite(Map_Struct_Output, [char(Output_name) '.shp']);
        
    end
    
    %% CHANNEL HEAD POINT SHAPEFILE
    if make_shapefiles_option == 1
        
        % Now, point shapefile of the channel head. It's up to the user to know
        % if this channel head represents the critical drainage area or not
        % (e.g., if you limit the upstream extent of the profile due to a
        % change in lithology).
        Channel_head_X = Stream_mod.x(Stream_mod.distance == max(Stream_mod.distance));
        
        Channel_head_Y = Stream_mod.y(Stream_mod.distance == max(Stream_mod.distance));
        
        Channel_head_Z = max(Z_sort);
        
        Incision_shapefile = mappoint(Channel_head_X, Channel_head_Y, 'Elev_m', Channel_head_Z, 'DA_m2', A_cr);
        
        
        shapewrite(Incision_shapefile, [char(Output_name) '_channel_head.shp']);
        
    end
    
    %% CHANNEL HEAD TEXT SUMMARY
    
    %     % First, write the summary table with details about the options used
    %     % and overall results.
    %     Channel_head_textoutput = [Channel_head_X, Channel_head_Y, Channel_head_Z, A_cr];
    %
    %     file = fopen([char(Output_name) '_channel_head.txt'], 'w');
    %
    %     fprintf(file, ['Channel_head_X_in_projection_used \t Channel_head_Y_in_projection_used \t ' ...
    %         'Channel_head_Z_m \t Channel_head_DA_m2 \r']);
    %
    %     for i = 1:1:length(Channel_head_textoutput)
    %
    %         if i ~= length(Channel_head_textoutput)
    %
    %             fprintf(file, '%E \t', Channel_head_textoutput(1,i));
    %
    %         elseif i == length(Channel_head_textoutput)
    %
    %             fprintf(file, '%E \r', Channel_head_textoutput(1,i));
    %
    %         end
    %
    %     end
    %
    %     fclose('all');
    
    %% OPTIONS TEXT SUMMARY
    
    %     % Write a summary table with details about the options used.
    %     Results_summary_textoutput = [use_DrEICH, A_clip_DrEICH, DrEICH_clip, theta_ref, use_ref_theta, ...
    %         smooth_window, smooth_for_DrEICH, Threshold_area_m2_initial, contour_interval_m, use_resampled, ...
    %         transience_projections, fillsinks_option, seglength];
    %
    %     file = fopen([char(Output_name) '_options_summary.txt'], 'w');
    %
    %     fprintf(file, ['use_DrEICH_option \t A_clip_DrEICH_m2 \t DrEICH_clip_option \t theta_ref \t use_ref_theta_option \t '...
    %         'smooth_window \t smooth_for_DrEICH_option \t Threshold_A_m2_initial \t contour_interval_m \t use_resampled \t ' ...
    %         'transience_projections_option \t fillsinks_option \t seglength_m_for_polyline \r']);
    %
    %     for i = 1:1:length(Results_summary_textoutput)
    %
    %         if i ~= length(Results_summary_textoutput)
    %
    %             fprintf(file, '%E \t', Results_summary_textoutput(1,i));
    %
    %         elseif i == length(Results_summary_textoutput)
    %
    %             fprintf(file, '%E \r', Results_summary_textoutput(1,i));
    %
    %         end
    %
    %     end
    %
    %     fclose('all');
    
    %% TEXT FILE FOR KSN RESULTS
    
    % Now, write a file with the ksn results obtained by selecting
    % reaches on the Chi plot (figure 5).
    Results_ksn_textoutput = [(1:1:Reach_num_ksn)', ones(Reach_num_ksn, 1) .* theta_ref, Reach_ks_Chi_fromChiplot, ...
        Reach_ks_SA_fromChiplot, DA_m2_bounds_fromChiplot, L_m_bounds_fromChiplot, Z_m_bounds_fromChiplot];
    
    file = fopen([char(Output_name) '_ksn.txt'], 'w');
    
    fprintf(file, ['Reach_number \t theta_ref \t ksn_Chi \t ksn_SA \t DA_m2_bound_bottom \t ' ...
        'DA_m2_bound_top \t L_m_bound_bottom \t L_m_bound_top \t Z_m_bound_bottom \t Z_m_bound_top \r']);
    
    columns = size(Results_ksn_textoutput);
    columns = columns(1,2);
    
    for i = 1:1:Reach_num_ksn
        
        for j = 1:1:columns
            
            if j == 1
                
                fprintf(file, '%.0f \t', Results_ksn_textoutput(i,j));
                
            elseif j > 1 && j < columns
                
                fprintf(file, '%E \t', Results_ksn_textoutput(i,j));
                
            elseif j == columns
                
                fprintf(file, '%E \r', Results_ksn_textoutput(i,j));
                
            end
            
        end
        
    end
    
    fclose('all');
    
    %% TEXT FILE FOR CONCAVITY RESULTS
    
    if use_ref_theta ~= 1
        
        % Now, write a file with the concavity results obtained by selecting
        % reaches on the slope-area plot (figure 4).
        Results_theta_textoutput = [(1:1:Reach_num_concav)', Reach_theta_fromSA, Reach_theta_fromChi, ...
            Reach_ks_fromSA, Reach_ks_fromChi, DA_m2_bounds_fromSA, L_m_bounds_fromSA, Z_m_bounds_fromSA];
        
        file = fopen([char(Output_name) '_theta.txt'], 'w');
        
        fprintf(file, ['Reach_number \t theta_SA \t theta_Chi \t ks_SA \t ks_Chi \t DA_m2_bound_bottom \t DA_m2_bound_top \t ' ...
            'L_m_bound_bottom \t L_m_bound_top \t Z_m_bound_bottom \t Z_m_bound_top \r']);
        
        columns = size(Results_theta_textoutput);
        columns = columns(1,2);
        
        for i = 1:1:Reach_num_concav
            
            for j = 1:1:columns
                
                if j == 1
                    
                    fprintf(file, '%.0f \t', Results_theta_textoutput(i,j));
                    
                elseif j > 1 && j < columns
                    
                    fprintf(file, '%E \t', Results_theta_textoutput(i,j));
                    
                elseif j == columns
                    
                    fprintf(file, '%E \r', Results_theta_textoutput(i,j));
                    
                end
                
            end
            
        end
        
        fclose('all');
        
    end
    
    %% TEXT FILE FOR SMOOTHED PROFILE
    
    % Here, write the smoothed profile data.
    profile_data = [L_sort, Z_sort, DA_sort, S_sort, Chi];
    
    file = fopen([char(Output_name) '_profile_smoothed_clipped.txt'], 'w');
    
    fprintf(file, 'L_m \t Z_m \t DA_m2 \t S \t Chi_m \r');
    
    size_output = size(profile_data);
    
    rows = size_output(1,1);
    columns = size_output(1,2);
    
    for i = 1:1:rows
        
        for j = 1:1:columns
            
            if j ~= columns
                
                fprintf(file, '%E \t', profile_data(i,j));
                
            elseif j == columns
                
                fprintf(file, '%E \r', profile_data(i,j));
                
            end
            
        end
        
    end
    
    fclose('all');
    
    clear columns size_output i j
    
    %% TEXT FILE FOR RAW PROFILE
    
    % Here, write the unsmoothed, unclipped profile data
    profile_data_raw = [L_sort_all, Z_sort_all, DA_sort_all, S_sort_all, Chi_all];
    
    file = fopen([char(Output_name) '_profile_raw.txt'], 'w');
    
    fprintf(file, 'L_m \t Z_m \t DA_m2 \t S \t Chi_m \r');
    
    size_output = size(profile_data_raw);
    
    rows = size_output(1,1);
    columns = size_output(1,2);
    
    for i = 1:1:rows
        
        for j = 1:1:columns
            
            if j ~= columns
                
                fprintf(file, '%E \t', profile_data_raw(i,j));
                
            elseif j == columns
                
                fprintf(file, '%E \r', profile_data_raw(i,j));
                
            end
            
        end
        
    end
    
    fclose('all');
    
    %% SAVE FIGURES
    
    % Now, saving some of the main figures created. The .epsc figures
    % will have screwed up text, but you can edit that in Adobe
    % Illustrator.
    if manual_Acr_option == 1
        
        if use_DrEICH == 0
            
            saveas(figure(1), [char(Output_name) '_Fig1.fig'])
            saveas(figure(1), [char(Output_name) '_Fig1'], 'epsc')
            
        elseif use_DrEICH == 1
            
            saveas(figure(2), [char(Output_name) '_Fig2.fig'])
            saveas(figure(2), [char(Output_name) '_Fig2'], 'epsc')
            
        end
        
    end
    
    figure(3)
    set(gcf, 'renderer', 'Painters')
    
    saveas(figure(3), [char(Output_name) '_Fig3.fig'])
    saveas(figure(3), [char(Output_name) '_Fig3'], 'epsc')
    
    if use_ref_theta ~= 1
        
        figure(4)
        set(gcf, 'renderer', 'Painters')
        
        saveas(figure(4), [char(Output_name) '_Fig4.fig'])
        saveas(figure(4), [char(Output_name) '_Fig4'], 'epsc')
        
        figure(6)
        set(gcf, 'renderer', 'Painters')
        
        saveas(figure(6), [char(Output_name) '_Fig6.fig'])
        saveas(figure(6), [char(Output_name) '_Fig6'], 'epsc')
        
    end
    
    figure(5)
    set(gcf, 'renderer', 'Painters')
    
    saveas(figure(5), [char(Output_name) '_Fig5.fig'])
    saveas(figure(5), [char(Output_name) '_Fig5'], 'epsc')
    
    figure(7)
    set(gcf, 'renderer', 'Painters')
    
    saveas(figure(7), [char(Output_name) '_Fig7.fig'])
    saveas(figure(7), [char(Output_name) '_Fig7'], 'epsc')
    
    if Reach_num_ksn == 2 && transience_projections == 1
        
        figure(8)
        set(gcf, 'renderer', 'Painters')
        
        saveas(figure(8), [char(Output_name) '_Fig8.fig'])
        saveas(figure(8), [char(Output_name) '_Fig8'], 'epsc')
        
        figure(9)
        set(gcf, 'renderer', 'Painters')
        
        saveas(figure(9), [char(Output_name) '_Fig9.fig'])
        saveas(figure(9), [char(Output_name) '_Fig9'], 'epsc')
        
    end
    
    figure(10)
    set(gcf, 'renderer', 'Painters')
    
    saveas(figure(10), [char(Output_name) '_Fig10.fig'])
    saveas(figure(10), [char(Output_name) '_Fig10'], 'epsc')
    
end
