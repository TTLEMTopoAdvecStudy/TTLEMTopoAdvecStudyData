
Original_location = pwd;

Mapview_Fig_location = 'D:\Users\nmitch10\D_TTLEM_Sims_no_implCLF\Mapview_Figs_May13';

% First, load the desired time step

Set_title = 4;

Scenario_title = 0;

% Time in years
t_ref = 2E+08;

% use_ksn_lims = 0;
% 
% min_ksn = 50;
% max_ksn = 430;

use_BORDER_option = 0;

save_option = 1;

Relief_window_1_m = 2500;

use_Relief_lims = 1;

min_Relief_m = 575;
max_Relief_m = 1675;

plot_only_South_Basins = 0;

plot_only_strong_Relief = 1;

% m^2
A_cr = 1E+06;

% m^2
A0 = 1E+06;

% m yr^-1
U1 = 0.0001;
U2 = U1 * 2;

% m ^(1 - (2 * n * theta_ref)) yr^-1
Ks = 2.5E-09;
Kw_over_Ks = 10;
Kw = Ks * Kw_over_Ks;

Vy = U1 * 5;

n = 2;
m = 1;

theta_ref = 0.5;

ksn1 = (U1 / Ks) ^ (1 / n);
ksn2 = (U2 / Ks) ^ (1 / n);

thresh_DivOrder_factor = 0.33;

KW_colors = gray(200);
KW_colors = KW_colors(100:end,:);

Relief_colors = cbrewer('seq', 'YlOrRd', 100, 'pchip');

%%

BORDER = GRIDobj(H1);
BORDER.Z(:,:) = 0;

if use_BORDER_option == 1
    
    BORDER.Z(:,1) = 1000;
    BORDER.Z(:,end) = 1000;
    
end

FlwDir = FLOWobj(H1 + BORDER,'mex',true,'preprocess','carve');

FlwAcc = flowacc(FlwDir);

FlwAcc.Z(1,:) = 0;
FlwAcc.Z(end,:) = 0;
FlwAcc.Z(:,1) = 0;
FlwAcc.Z(:,end) = 0;

Stream = STREAMobj(FlwDir, FlwAcc > (A_cr / H1.cellsize^2));

Grad = gradient8(H1, 'per');

ksn_grid = (Grad ./ 100) ./ ((FlwAcc .* (FlwAcc.cellsize ^ 2)) .^ (-theta_ref));

% DRAINAGE DIVIDE NETWORK
DrnDiv = DIVIDEobj(FlwDir, Stream);
DrnDiv = cleanedges(DrnDiv,FlwDir);

DrnDiv = sort(DrnDiv);

DrnDiv = divorder(DrnDiv,'topo');

% Only the divide network over the threshold
thresh_DivOrder = (max(DrnDiv.order) * thresh_DivOrder_factor);

DrnDiv_plot = DrnDiv;
DrnDiv_plot.order(DrnDiv_plot.order < thresh_DivOrder) = NaN;

Relief_m = localtopography(H1, Relief_window_1_m);

[x,y] = getcoordinates(H1,'GRIDobj');

N_Basins = drainagebasins(FlwDir, x.Z(1,:), y.Z(1,:));
N_Basins.Z(N_Basins.Z ~= 0) = 1;

% figure(4)
%
% imageschs(H1, N_Basins)

S_Basins = drainagebasins(FlwDir, x.Z(1,:), y.Z(end,:));
S_Basins.Z(S_Basins.Z ~= 0) = 1;

%%

check_Kw_grid = exist('Kw_grid','var');

if exist('Kw_grid','var') == 0
    
    Kw_grid = GRIDobj(H1);
    Kw_grid.Z = Ks;
    
end

figure(2)

% imageschs(H1, Kw_grid, 'colormap', KW_colors, ...
%     'ticksToKm', false, 'colorBarLabel', '\bfm');
% hold on

if plot_only_South_Basins == 1
    
    Relief_m.Z(S_Basins.Z ~= 1) = NaN;
    
end

if plot_only_strong_Relief == 1
    
    Relief_m.Z(Kw_grid.Z == max(max(Kw_grid.Z))) = NaN;
    
end

if use_Relief_lims == 1
    
    imageschs(H1, Relief_m, 'colormap', Relief_colors, ...
        'caxis', [min_Relief_m, max_Relief_m], 'ticksToKm', false, 'colorBarLabel', '\bfm');
    hold on
    
elseif use_Relief_lims ~= 1
    
    imageschs(H1, Relief_m, 'colormap', Relief_colors, ...
        'ticksToKm', false, 'colorBarLabel', '\bfm');
    hold on
    
end

% plot(DrnDiv_plot, 'color', [1 1 1])
% hold on

% plot(DrnDiv_plot, 'color', [1 1 1])
% hold on
% 
% Map_Struct = STREAMobj2mapstruct(Stream, 'seglength', 1000, 'attributes',...
%     {'ksn' ksn_grid @mean});
% 
% if use_ksn_lims == 1
%     
% symbolspec = makesymbolspec('line', {'ksn' [min_ksn max_ksn] ...
%     'color' flipud(cbrewer('div', 'RdYlBu', 100, 'pchip')) 'linewidth' 1.5});
% 
% elseif use_ksn_lims ~= 1
%     
%     symbolspec = makesymbolspec('line', {'ksn' [min([Map_Struct.ksn]) max([Map_Struct.ksn])] ...
%     'color' flipud(cbrewer('div', 'RdYlBu', 100, 'pchip')) 'linewidth' 1.5});
% 
% end
% 
% mapshow(Map_Struct, 'SymbolSpec', symbolspec);

h = colorbar;
colormap(Relief_colors)
ylabel(h, 'Relief (m)', 'fontsize', 10)

if use_Relief_lims == 1
    
    caxis([min_Relief_m max_Relief_m])
    
elseif use_Relief_lims ~= 1
    
    caxis([min(min(Relief_m.Z)) max(max(Relief_m.Z))])
    
end

min(min(Relief_m.Z)) 
max(max(Relief_m.Z))

% Stream.x=Stream.x*1e-3;
% Stream.y=Stream.y*1e-3;
plot(Stream,'b','linewidth',1.5);
hold off

xlabel('\bfX Coordinate, km')
ylabel('\bfY Coordinate, km')
title(['Set ' num2str(Set_title, '%.0f')  ' Scenario ' num2str(Scenario_title, '%.0f') ...
    ', t = ' num2str((t_ref * 1e-6), '%.2f') ' Myr']);

set(gcf, 'renderer', 'Painters')



if save_option == 1
    
    cd(Mapview_Fig_location)
    
    saveas(figure(2),['Set' num2str(Set_title, '%.0f') ...
        '_Scenario' num2str(Scenario_title, '%.0f') '_Relief_t0.fig'])
    saveas(figure(2),['Set' num2str(Set_title, '%.0f') ...
        '_Scenario' num2str(Scenario_title, '%.0f') '_Relief_t0'], 'epsc')
    
    cd(Original_location)
    
end
