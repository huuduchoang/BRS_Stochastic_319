clear; clc; close all;

%% 1. Initial condition for 12D system
% Old reduced state:
% [v; n; hp; alpha; vollung; PO2lung; PO2blood]
old7 = [-50.3697203300953;
        0.00462275707173923;
        0.750556122453219;
        3.93159280387281e-05;
        2.00894512526427;
        94.4645088086363;
        93.3277058316367];

v0       = old7(1);
n_scalar = old7(2);
hp0      = old7(3);
alpha0   = old7(4);
vollung0 = old7(5);
PO2lung0 = old7(6);
PO2blood0= old7(7);

% Lift scalar n into 5-state K occupancy
Kinit = [
    (1 - n_scalar)^4;
    4 * n_scalar * (1 - n_scalar)^3;
    6 * n_scalar^2 * (1 - n_scalar)^2;
    4 * n_scalar^3 * (1 - n_scalar);
    n_scalar^4
];

% Fast sodium factor consistent with your 12D model
hf0 = 1 - n_scalar;

% 12D initial condition:
% [v; n0; n1; n2; n3; n4; hp; hf; alpha; vollung; PO2lung; PO2blood]
inits12 = [
    v0;
    Kinit;
    hp0;
    hf0;
    alpha0;
    vollung0;
    PO2lung0;
    PO2blood0
];

disp(['Sum of Kinit = ' num2str(sum(Kinit), '%.16f')]);

%% 2. Simulation settings
tf = 100000;              % ms
dt = 0.05;               % ms
NK_target = 1e4;         % target K-channel count
Area = NK_target / 18;   % because solver uses NK = round(18*Area)

params_noisy = struct();
params_noisy.M = 8e-6;
params_noisy.noise_on = true;

params_clean = struct();
params_clean.M = 8e-6;
params_clean.noise_on = false;

%% 3. Run non-noisy 12D simulation
disp('Running non-noisy 12D simulation...');
[t_clean, U_clean] = simulate_closedloop_Kfull(tf, dt, inits12, Area, params_clean);

%% 4. Run noisy 12D simulation
disp('Running noisy 12D simulation...');
rng(1);   % reproducible noisy run
[t_noisy, U_noisy] = simulate_closedloop_Kfull(tf, dt, inits12, Area, params_noisy);

%% 5. Extract variables
% Non-noisy
time_clean    = t_clean / 1000;
v_clean       = U_clean(:,1);
n0_clean      = U_clean(:,2);
n1_clean      = U_clean(:,3);
n2_clean      = U_clean(:,4);
n3_clean      = U_clean(:,5);
n4_clean      = U_clean(:,6);
hp_clean      = U_clean(:,7);
hf_clean      = U_clean(:,8);
alpha_clean   = U_clean(:,9);
vollung_clean = U_clean(:,10);
PO2lung_clean = U_clean(:,11);
PO2blood_clean= U_clean(:,12);
gtonic_clean  = 0.3 * (1 - tanh((PO2blood_clean - 85) ./ 30));

% Noisy
time_noisy    = t_noisy / 1000;
v_noisy       = U_noisy(:,1);
n0_noisy      = U_noisy(:,2);
n1_noisy      = U_noisy(:,3);
n2_noisy      = U_noisy(:,4);
n3_noisy      = U_noisy(:,5);
n4_noisy      = U_noisy(:,6);
hp_noisy      = U_noisy(:,7);
hf_noisy      = U_noisy(:,8);
alpha_noisy   = U_noisy(:,9);
vollung_noisy = U_noisy(:,10);
PO2lung_noisy = U_noisy(:,11);
PO2blood_noisy= U_noisy(:,12);
gtonic_noisy  = 0.3 * (1 - tanh((PO2blood_noisy - 85) ./ 30));

%% 6. Plot superimposed traces
figure('Position',[80 80 1300 800],'Color','w');
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

c_noisy = [0 0.4470 0.7410];
c_clean = 'k';

% V
nexttile; hold on
plot(time_noisy, v_noisy, 'Color', c_noisy, 'LineWidth', 1.3, 'DisplayName', 'Noisy');
plot(time_clean, v_clean, 'Color', c_clean, 'LineWidth', 1.1, 'DisplayName', 'Non-noisy');
title('Voltage');
ylabel('V (mV)');
grid on; box on; legend('Location','best');

% n4
nexttile; hold on
plot(time_noisy, n4_noisy, 'Color', c_noisy, 'LineWidth', 1.3);
plot(time_clean, n4_clean, 'Color', c_clean, 'LineWidth', 1.1);
title('Open K fraction');
ylabel('n_4');
grid on; box on;

% hp
nexttile; hold on
plot(time_noisy, hp_noisy, 'Color', c_noisy, 'LineWidth', 1.3);
plot(time_clean, hp_clean, 'Color', c_clean, 'LineWidth', 1.1);
title('Persistent Na inactivation');
ylabel('h_p');
grid on; box on;

% hf
nexttile; hold on
plot(time_noisy, hf_noisy, 'Color', c_noisy, 'LineWidth', 1.3);
plot(time_clean, hf_clean, 'Color', c_clean, 'LineWidth', 1.1);
title('Fast Na factor');
ylabel('h_f');
grid on; box on;

% alpha
nexttile; hold on
plot(time_noisy, alpha_noisy, 'Color', c_noisy, 'LineWidth', 1.3);
plot(time_clean, alpha_clean, 'Color', c_clean, 'LineWidth', 1.1);
title('Motor pool');
ylabel('\alpha');
grid on; box on;

% lung volume
nexttile; hold on
plot(time_noisy, vollung_noisy, 'Color', c_noisy, 'LineWidth', 1.3);
plot(time_clean, vollung_clean, 'Color', c_clean, 'LineWidth', 1.1);
title('Lung volume');
ylabel('Vol_{lung}');
grid on; box on;

% lung oxygen
nexttile; hold on
plot(time_noisy, PO2lung_noisy, 'Color', c_noisy, 'LineWidth', 1.3);
plot(time_clean, PO2lung_clean, 'Color', c_clean, 'LineWidth', 1.1);
title('Lung oxygen');
ylabel('PO2_{lung}');
xlabel('Time (s)');
grid on; box on;

% blood oxygen
nexttile; hold on
plot(time_noisy, PO2blood_noisy, 'Color', c_noisy, 'LineWidth', 1.3);
plot(time_clean, PO2blood_clean, 'Color', c_clean, 'LineWidth', 1.1);
title('Blood oxygen');
ylabel('PO2_{blood}');
xlabel('Time (s)');
grid on; box on;

% chemosensory drive
nexttile; hold on
plot(time_noisy, gtonic_noisy, 'Color', c_noisy, 'LineWidth', 1.3);
plot(time_clean, gtonic_clean, 'Color', c_clean, 'LineWidth', 1.1);
title('Chemosensory drive');
ylabel('g_{tonic}');
xlabel('Time (s)');
grid on; box on;

sgtitle('12D Closed-Loop: Noisy vs Non-noisy');

% %% 7. Optional: save output
% save('compare_noisy_vs_nonoisy_12D.mat', ...
%     't_clean', 'U_clean', 't_noisy', 'U_noisy', ...
%     'time_clean', 'time_noisy', ...
%     'v_clean', 'v_noisy', ...
%     'n0_clean', 'n1_clean', 'n2_clean', 'n3_clean', 'n4_clean', ...
%     'n0_noisy', 'n1_noisy', 'n2_noisy', 'n3_noisy', 'n4_noisy', ...
%     'hp_clean', 'hp_noisy', 'hf_clean', 'hf_noisy', ...
%     'alpha_clean', 'alpha_noisy', ...
%     'vollung_clean', 'vollung_noisy', ...
%     'PO2lung_clean', 'PO2lung_noisy', ...
%     'PO2blood_clean', 'PO2blood_noisy', ...
%     'gtonic_clean', 'gtonic_noisy', ...
%     'inits12', 'params_clean', 'params_noisy', 'tf', 'dt', 'Area', ...
%     '-v7.3');