
KW_over_KS = 10;

%mean_ksn_S_i = mean_ksn_in_strong_S_over_time(1,1);

% m per yr
U = 1E-04;

KS_ref = 2.5E-09;

KW = KS_ref * KW_over_KS;

n = 2;

K_star = (KW_over_KS) ^ (1 / (1 - n));

ksn_W = (U / KW) ^ (1 / n);
ksn_S = (U / KS_ref) ^ (1 / n);

ksn_S_amplified_pred = (1 / K_star) * ksn_W;

ksn_S_star_amplified_pred = ksn_S_amplified_pred / ksn_S

ksn_S_star_amplified_pred = (KS_ref / KW) ^ (1 / (n * (1 - n)))

%ksn_S_amplified_star_pred = ksn_S_amplified_pred / mean_ksn_S_i
