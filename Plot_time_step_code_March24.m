
Original_location = pwd;

Mapview_Fig_location = 'D:\Users\USERNAME\D_TTLEM_Sims_no_implCLF\Mapview_figures_April7';

% First, load the desired time step

Set_title = 5;

Scenario_title = 11;

% Time in years
t_ref = 2.5E+07;

use_ksn_lims = 1;

min_ksn = 50;
max_ksn = 300;

use_BORDER_option = 0;

save_option = 1;

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

ksn_colors = cbrewer('div', 'RdYlBu', 100, 'pchip');
ksn_colors = flipud(ksn_colors);

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

%%

figure(2)

check_Kw_grid = exist('Kw_grid','var');

if exist('Kw_grid','var') == 0
    
    Kw_grid = GRIDobj(H1);
    Kw_grid.Z = Ks;
    
end

imageschs(H1, Kw_grid, 'colormap', KW_colors, ...
    'ticksToKm', false, 'colorBarLabel', '\bfm');
hold on

plot(DrnDiv_plot, 'color', [1 1 1])
hold on

Map_Struct = STREAMobj2mapstruct(Stream, 'seglength', 1000, 'attributes',...
    {'ksn' ksn_grid @mean});

if use_ksn_lims == 1
    
symbolspec = makesymbolspec('line', {'ksn' [min_ksn max_ksn] ...
    'color' flipud(cbrewer('div', 'RdYlBu', 100, 'pchip')) 'linewidth' 1.5});

elseif use_ksn_lims ~= 1
    
    symbolspec = makesymbolspec('line', {'ksn' [min([Map_Struct.ksn]) max([Map_Struct.ksn])] ...
    'color' flipud(cbrewer('div', 'RdYlBu', 100, 'pchip')) 'linewidth' 1.5});

end

mapshow(Map_Struct, 'SymbolSpec', symbolspec);

h = colorbar;
colormap(ksn_colors)
ylabel(h, 'k_{sn} (m)', 'fontsize', 10)

if use_ksn_lims == 1
    
    caxis([min_ksn max_ksn])
    
elseif use_ksn_lims ~= 1
    
    caxis([min([Map_Struct.ksn]) max([Map_Struct.ksn])])
    
end

% Stream.x=Stream.x*1e-3;
% Stream.y=Stream.y*1e-3;
% plot(Stream,'b','linewidth',1);
% hold off

xlabel('\bfX Coordinate, km')
ylabel('\bfY Coordinate, km')
title(['Set ' num2str(Set_title, '%.0f')  ' Scenario ' num2str(Scenario_title, '%.0f') ...
    ', t = ' num2str((t_ref * 1e-6), '%.2f') ' Myr']);

set(gcf, 'renderer', 'Painters')



if save_option == 1
    
    cd(Mapview_Fig_location)
    
    saveas(figure(2),['Set' num2str(Set_title, '%.0f') ...
        'Scenario' num2str(Scenario_title, '%.0f') '.fig'])
    saveas(figure(2),['Set' num2str(Set_title, '%.0f') ...
        'Scenario' num2str(Scenario_title, '%.0f')], 'epsc')
    
    cd(Original_location)
    
end
