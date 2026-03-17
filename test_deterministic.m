clear
clc

% Old 7D initial condition:
% [v; n; h; alpha; vollung; PO2lung; PO2blood]
old = [-61.2345;
        0.000316;
        0.71447;
        4.9339e-05;
        2.0128;
        116.7726;
        115.4154];

% Lift scalar n into 5-state K occupancy
n_scalar = old(2);

Kinit = [0.000316;
         0.000316;
         0.000316;
         0.000316;
         0.000316];

% New 11D initial condition:
% [v; n0; n1; n2; n3; n4; h; alpha; vollung; PO2lung; PO2blood]
inits = [old(1);
         Kinit;
         old(3);
         old(4);
         old(5);
         old(6);
         old(7)];

tf   = 200000;    % ms
dt   = 0.05;    % ms
Area = 100000;     % try something moderate first

[t,U] = simulate_closedloop_Kfull(tf,dt,inits,Area);

v        = U(:,1);
n0       = U(:,2);
n1       = U(:,3);
n2       = U(:,4);
n3       = U(:,5);
n4       = U(:,6);
h        = U(:,7);
alpha    = U(:,8);
vollung  = U(:,9);
PO2lung  = U(:,10);
PO2blood = U(:,11);

time_s = t/1000;
gtonic = 0.3*(1 - tanh((PO2blood - 85)./30));
colors = lines(5);

% Figure 1: 2x3 closed-loop grid
figure('Position',[100 100 900 600],'Color','w');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile(1); hold on
plot(time_s, v, 'k', 'LineWidth', 1.3);
ylabel('V (mV)');
title('CPG (V)');
xlim([0 tf/1000]); box on; grid on

nexttile(2)
plot(time_s, alpha, 'k', 'LineWidth', 1.5);
ylabel('\alpha');
title('Motor Pool');
xlim([0 tf/1000]); box on; grid on

nexttile(3)
plot(time_s, vollung, 'k', 'LineWidth', 1.5);
ylabel('Vol_{lung}');
title('Lung Volume');
xlim([0 tf/1000]); box on; grid on

nexttile(4); hold on
plot(time_s, gtonic, 'k', 'LineWidth', 1.3);
ylabel('g_{tonic}');
title('Chemosensation');
xlim([0 tf/1000]); box on; grid on

nexttile(5)
plot(time_s, PO2blood, 'k', 'LineWidth', 1.5);
ylabel('P_{aO_2} (mmHg)');
title('Blood Oxygen');
xlim([0 tf/1000]); box on; grid on

nexttile(6)
plot(time_s, PO2lung, 'k', 'LineWidth', 1.5);
ylabel('P_{O_2}^{lung} (mmHg)');
title('Lung Oxygen');
xlabel('Time (s)');
xlim([0 tf/1000]); box on; grid on

% Figure 2: potassium occupancy states
figure('Position',[150 150 900 450],'Color','w'); hold on
plot(time_s, n0, 'Color', colors(1,:), 'LineWidth', 1.4);
plot(time_s, n1, 'Color', colors(2,:), 'LineWidth', 1.4);
plot(time_s, n2, 'Color', colors(3,:), 'LineWidth', 1.4);
plot(time_s, n3, 'Color', colors(4,:), 'LineWidth', 1.4);
plot(time_s, n4, 'Color', colors(5,:), 'LineWidth', 1.4);

xlabel('Time (s)');
ylabel('Occupancy Fraction');
title('Potassium Channel State Occupancies');
legend({'n0','n1','n2','n3','n4'}, 'Location', 'best');
xlim([0 tf/1000]); box on; grid on