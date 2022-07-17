
% Navigate to the correct folder.
resultsdir = [pwd '\'];
fileprefix = 'Set4_ini';

AreaThresh = 1E+06;

min_time = 5E+07;
max_time = 2E+08;
time_interval = 1E+06;

Check_Streams_for_SS = 1;

make_video_option = 1;

enforce_clim_in_video_option = 0;

min_c = 0;
max_c = 2700;

make_elev_plot = 1;

% Threshold for % change in mean z for the consideration of steady state
percent_z_change_thresh = 1E-03;

% Threshold for time since changes to streams (Myr) for the consideration 
% of steady state
Compiled_time_since_S_change_threshold = 1E+07;

%% Movie and overview figure
% Create simulation movie of the results and generate an overview figure
% illustration of the results

if Check_Streams_for_SS == 1 || make_elev_plot == 1
    
    time_since_S_change = 0;
    
    time_range = min_time:time_interval:max_time;

    Compiled_time_since_S_change = zeros(1,length(time_range));
    
    Compiled_min_elev = zeros(1,length(time_range));
    Compiled_max_elev = zeros(1,length(time_range));
    Compiled_mean_elev = zeros(1,length(time_range));
    
end

if make_video_option == 1
    
    scrsz = get(0,'ScreenSize');
    movieTopo = figure('OuterPosition',[0.2*scrsz(4) 0.2*scrsz(4) scrsz(4) .5*scrsz(4)],'Name','Topo','color','white');
    
    ui = 0;
    
    % Create video file
    Vid = VideoWriter([fileprefix '_Simulation_z_over_t.avi']);
    
    open(Vid)
    
end

t_ref = 0;

for t = min_time:time_interval:max_time
    
    t_ref = t_ref + 1;
    
    load([resultsdir fileprefix '_t_' num2str(round(t / 1e3)) '_kyr.mat'],'H1');
    
    if t == min_time
        
        BORDER = GRIDobj(H1,'logical');
        
        BORDER.Z(:,1) = 1;
        BORDER.Z(:,end) = 1;
        BORDER.Z = BORDER.Z * 10000;
        
    end
    
    FD  = FLOWobj(H1+BORDER,'mex',true,'preprocess','c');
    S = STREAMobj(FD,flowacc(FD)>(AreaThresh / H1.cellsize^2));
    
    if make_video_option == 1
        
        ui = ui + 1;
        figure(movieTopo)
        
        if enforce_clim_in_video_option == 1
            
            imageschs(H1,H1,'caxis', [min_c max_c], 'colormap', 'gray', ...
                'ticksToKm', true, 'colorBarLabel', '\bfm');
            hold on
            
        elseif enforce_clim_in_video_option ~= 1
            
            imageschs(H1,H1, 'colormap', 'gray', ...
                'ticksToKm', true, 'colorBarLabel', '\bfm');
            hold on
            
        end
        
        S.x=S.x*1e-3;
        S.y=S.y*1e-3;
        plot(S,'b','linewidth',1);
        hold off
        xlabel('\bfX Coordinate, km')
        ylabel('\bfY Coordinate, km')
        title([num2str((t*1e-6), '%.2f') ' Myr']);
        
        
        
        set(gcf, 'renderer', 'Painters')
        
        movie1(ui)=getframe(gcf); %#ok<SAGROW>
        writeVideo(Vid,movie1(ui))
        
    end
    
    if Check_Streams_for_SS == 1
        
        if t ~= min_time
            
            % Just a thorough way of making sure it's the same stream
            % object
            if length(S.x) == length(last_S.x) && length(S.y) == length(last_S.y)
                
                if min(S.x == last_S.x) == 1 && min(S.y == last_S.y) == 1
                    
                    time_since_S_change = time_since_S_change + time_interval;
                
                end
                
            else
                
                time_since_S_change = 0;
                
            end
            
            Compiled_time_since_S_change(1,t_ref) = time_since_S_change;
            
        end
        
        % Store the current stream object so it can be compared with the next
        last_S = S;
        
    end
    
    if make_elev_plot == 1
        
        Compiled_min_elev(1,t_ref) = min(min(H1.Z));
        Compiled_max_elev(1,t_ref) = max(max(H1.Z));
        Compiled_mean_elev(1,t_ref) = mean(H1.Z, 'all');
        
    end
    
end

if make_elev_plot == 1
    
    figure(2)
    
    h_min = plot(time_range ./ 1E+06, Compiled_min_elev, 'color', 'b', 'marker', 's', ...
        'linestyle', 'none', 'linewidth', 1);
    hold on
    
    h_max = plot(time_range ./ 1E+06, Compiled_max_elev, 'color', 'r', 'marker', '^', ...
        'linestyle', 'none', 'linewidth', 1);
    hold on
    
    h_mean = plot(time_range ./ 1E+06, Compiled_mean_elev, 'color', 'k', 'marker', 'o', ...
        'linestyle', 'none', 'linewidth', 1);
    hold on
    
    set(gca, 'fontsize', 10)
    
    xlabel('t (Myr)', 'fontsize', 12)
    ylabel('z (m)', 'fontsize', 12)
    
    h_lgnd = legend([h_min, h_mean, h_max], 'Minimum', 'Mean', 'Maximum', ...
        'location', 'best');
    
    set(h_lgnd, 'fontsize', 10)
    
    set(gcf, 'renderer', 'Painters')
    
    saveas(figure(2),[fileprefix '_min_max_mean_z_over_t.fig'])
    saveas(figure(2),[fileprefix '_min_max_mean_over_t.png'])
    
    
    
    max_elev_change_percent = zeros(1,length(time_range));
    mean_elev_change_percent = zeros(1,length(time_range));
    
    for i = 2:1:length(time_range)
        
        max_elev_change_percent(1,i) = ((Compiled_max_elev(1,i) - Compiled_max_elev(1,i-1)) / Compiled_max_elev(1,i-1)) * 100;
        mean_elev_change_percent(1,i) = ((Compiled_mean_elev(1,i) - Compiled_mean_elev(1,i-1)) / Compiled_mean_elev(1,i-1)) * 100;
        
    end
    
    figure(3)
    
    h_max = semilogy(time_range ./ 1E+06, abs(max_elev_change_percent) .* 1E+03, 'color', 'r', 'marker', '^', ...
        'linestyle', 'none', 'linewidth', 1);
    hold on
    
    h_mean = semilogy(time_range ./ 1E+06, abs(mean_elev_change_percent) .* 1E+03, 'color', 'k', 'marker', 'o', ...
        'linestyle', 'none', 'linewidth', 1);
    hold on
    
    h_thresh = semilogy([min_time, max_time] ./ 1E+06, ...
        [percent_z_change_thresh, percent_z_change_thresh] .* 1E+03, 'color', 'k', 'marker', 'none', ...
        'linestyle', '--', 'linewidth', 1);
    hold on
    
    xlabel('t (Myr)', 'fontsize', 12)
    ylabel('Percent change in z', 'fontsize', 12)
    
    h_lgnd = legend([h_mean, h_max, h_thresh], 'Mean z', 'Maximum z', 'Threshold for S.S.', ...
        'location', 'southwest');
    
    set(h_lgnd, 'fontsize', 10)
    
    set(gcf, 'renderer', 'Painters')
    
    saveas(figure(3),[fileprefix '_Rate_of_max_mean_z_Changes_over_t.fig'])
    saveas(figure(3),[fileprefix '_Rate_of_max_mean_z_Changes_over_t.png'])
    
    save([fileprefix '_Elev_plot_Workspace.mat'])
    
end

if make_video_option == 1
    
    close(Vid)
    
end

if Check_Streams_for_SS == 1
    
    figure(4)
    
    plot(time_range ./ 1E+06, Compiled_time_since_S_change ./ 1E+06, ...
        'marker', 'o', 'color', 'b', ...
        'linewidth', 1, 'linestyle', 'none')
    hold on
    
    h_thresh = plot([min_time, max_time] ./ 1E+06, ...
        [Compiled_time_since_S_change_threshold, Compiled_time_since_S_change_threshold] ./ 1E+06, ...
        'marker', 'none', 'color', 'k', ...
        'linewidth', 1, 'linestyle', '--');
    hold on
    
    legend(h_thresh, 'Threshold for S.S.', 'location', 'northwest')
    
    xlabel('t (Myr)', 'fontsize', 12)
    ylabel('Time since changes to streams (Myr)', 'fontsize', 12)
    
    set(gcf, 'renderer', 'Painters')
    
    saveas(figure(4),[fileprefix '_Check_Streams_for_SS.fig'])
    saveas(figure(4),[fileprefix '_Check_Streams_for_SS.png'])
    
    save([fileprefix '_Check_Streams_Workspace.mat'])
    
end
