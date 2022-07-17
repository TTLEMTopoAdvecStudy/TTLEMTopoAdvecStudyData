
% 1 for ksn, 0 for Relief
ksn_or_Relief_option = 0;

y_min_ksn = 0.2;
y_max_ksn = 2.2;

y_min_Relief = 0.75;
y_max_Relief = 1.25;

x_min = 0;
x_max = 100;

if ksn_or_Relief_option == 1
    
    ylim([y_min_ksn y_max_ksn])
    
elseif ksn_or_Relief_option == 0
    
    ylim([y_min_Relief y_max_Relief])
    
end

xlim([x_min x_max])

set(gcf, 'renderer', 'Painters')

saveas(figure(1), 'FORMATTED_Relief_in_S_over_time_Fig_LowDiff', 'epsc')
