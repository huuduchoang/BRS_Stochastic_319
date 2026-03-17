clear; clc;

%% 1. Initial Conditions Setup
% Deterministic 7D initial condition:
% [v; n; hp; alpha; vollung; PO2lung; PO2blood]
old = [-61.2345;
        0.000316;
        0.71447;
        4.9339e-05;
        2.0128;
        116.7726;
        115.4154];

% Lift scalar n into 5-state K occupancy
n_scalar = old(2);
Kinit = [
    (1 - n_scalar)^4;
    4 * n_scalar * (1 - n_scalar)^3;
    6 * n_scalar^2 * (1 - n_scalar)^2;
    4 * n_scalar^3 * (1 - n_scalar);
    n_scalar^4
];

% Initialize hp and hf for stochastic model
hp0 = old(3);
hf0 = 1 - n_scalar;

% New 12D initial condition:
% [v; n0; n1; n2; n3; n4; hp; hf; alpha; vollung; PO2lung; PO2blood]
inits = [old(1);
         Kinit;
         hp0;
         hf0;
         old(4);
         old(5);
         old(6);
         old(7)];

disp(['Sum of Kinit = ' num2str(sum(Kinit), '%.16f')])

tf = 100000;   % ms
dt = 0.05;    % ms

% If your solver uses NK = round(18*Area), then this gives about 1e4 K channels
NK_target = 1e4;
Area = NK_target / 18;

%% 2. Run Stochastic Model
disp('Running Stochastic Simulation...');
[t_sim, U_sim] = simulate_closedloop_Kfull(tf, dt, inits, Area);

% Extract Stochastic Data
v        = U_sim(:,1);
n0       = U_sim(:,2);
n1       = U_sim(:,3);
n2       = U_sim(:,4);
n3       = U_sim(:,5);
n4       = U_sim(:,6);
hp       = U_sim(:,7);
hf       = U_sim(:,8);
alpha    = U_sim(:,9);
vollung  = U_sim(:,10);
PO2lung  = U_sim(:,11);
PO2blood = U_sim(:,12);

time_s = t_sim / 1000;
gtonic = 0.3 * (1 - tanh((PO2blood - 85) ./ 30));

%% 3. Run Deterministic Model (7D)
disp('Running Deterministic ODE...');
[t_det, U_det] = ode15s(@closedloop_modified, [0 tf], old);

% Extract Deterministic Data
time_s_det   = t_det / 1000;
v_det        = U_det(:,1);
n_det        = U_det(:,2);
hp_det       = U_det(:,3);
alpha_det    = U_det(:,4);
vollung_det  = U_det(:,5);
PO2lung_det  = U_det(:,6);
PO2blood_det = U_det(:,7);

% Deterministic fast sodium factor for comparison
hf_det = 1 - n_det;

% Calculate Chemosensation for Deterministic
gtonic_det = 0.3 * (1 - tanh((PO2blood_det - 85) ./ 30));

%% 4. Figure 1: 2x3 Closed-loop Grid (Superimposed)
figure('Position',[100 100 1000 600],'Color','w');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% Colors
stoch_color = [0 0.4470 0.7410];
det_color   = 'k';

% 1. CPG (V)
nexttile(1); hold on
plot(time_s, v, 'Color', stoch_color, 'LineWidth', 1.5, 'DisplayName', 'Stochastic');
plot(time_s_det, v_det, 'Color', det_color, 'LineWidth', 1.2, 'DisplayName', 'Deterministic');
ylabel('V (mV)');
title('CPG (V)');
xlim([0 tf/1000]);
box on; grid on
legend('Location','best');

% 2. Motor Pool
nexttile(2); hold on
plot(time_s, alpha, 'Color', stoch_color, 'LineWidth', 1.5);
plot(time_s_det, alpha_det, 'Color', det_color, 'LineWidth', 1.2);
ylabel('\alpha');
title('Motor Pool');
xlim([0 tf/1000]);
box on; grid on

% 3. Lung Volume
nexttile(3); hold on
plot(time_s, vollung, 'Color', stoch_color, 'LineWidth', 1.5);
plot(time_s_det, vollung_det, 'Color', det_color, 'LineWidth', 1.2);
ylabel('Vol_{lung}');
title('Lung Volume');
xlim([0 tf/1000]);
box on; grid on

% 4. Chemosensation
nexttile(4); hold on
plot(time_s, gtonic, 'Color', stoch_color, 'LineWidth', 1.5);
plot(time_s_det, gtonic_det, 'Color', det_color, 'LineWidth', 1.2);
ylabel('g_{tonic}');
title('Chemosensation');
xlim([0 tf/1000]);
box on; grid on

% 5. Blood Oxygen
nexttile(5); hold on
plot(time_s, PO2blood, 'Color', stoch_color, 'LineWidth', 1.5);
plot(time_s_det, PO2blood_det, 'Color', det_color, 'LineWidth', 1.2);
ylabel('P_{aO_2} (mmHg)');
title('Blood Oxygen');
xlabel('Time (s)');
xlim([0 tf/1000]);
box on; grid on

% 6. Lung Oxygen
nexttile(6); hold on
plot(time_s, PO2lung, 'Color', stoch_color, 'LineWidth', 1.5);
plot(time_s_det, PO2lung_det, 'Color', det_color, 'LineWidth', 1.2);
ylabel('P_{O_2}^{lung} (mmHg)');
title('Lung Oxygen');
xlabel('Time (s)');
xlim([0 tf/1000]);
box on; grid on

%% Optional: comparison plots for K / sodium inactivation variables
figure('Position',[120 120 1000 700],'Color','w');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile; hold on
plot(time_s, n4, 'Color', stoch_color, 'LineWidth', 1.5, 'DisplayName', 'Stochastic n_4');
plot(time_s_det, n_det.^4, 'Color', det_color, 'LineWidth', 1.2, 'DisplayName', 'Deterministic n^4');
ylabel('K activation');
title('n_4 vs n^4');
xlim([0 tf/1000]);
box on; grid on
legend('Location','best');

nexttile; hold on
plot(time_s, hp, 'Color', stoch_color, 'LineWidth', 1.5, 'DisplayName', 'Stochastic h_p');
plot(time_s_det, hp_det, 'Color', det_color, 'LineWidth', 1.2, 'DisplayName', 'Deterministic h_p');
ylabel('h_p');
title('Persistent sodium inactivation');
xlim([0 tf/1000]);
box on; grid on
legend('Location','best');

nexttile; hold on
plot(time_s, hf, 'Color', stoch_color, 'LineWidth', 1.5, 'DisplayName', 'Stochastic h_f');
plot(time_s_det, hf_det, 'Color', det_color, 'LineWidth', 1.2, 'DisplayName', 'Deterministic 1-n');
ylabel('h_f');
xlabel('Time (s)');
title('Fast sodium factor');
xlim([0 tf/1000]);
box on; grid on
legend('Location','best');

%% Save all variables
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
baseName  = ['closedloop_compare_' timestamp];

save([baseName '.mat'], '-v7.3');

disp(['Saved all variables to: ' baseName '.mat']);

%% --- LOCAL FUNCTION FOR DETERMINISTIC MODEL ---
function z = closedloop_modified(~,u)
    % State variables
    % u = [v; n; hp; alpha; vollung; PO2lung; PO2blood]
    v        = u(1);
    n        = u(2);
    hp       = u(3);
    alpha    = u(4);
    vollung  = u(5);
    PO2lung  = u(6);
    PO2blood = u(7);

    z = zeros(7,1);

    %% CPG
    C = 21;
    gnap = 2.8;
    gna  = 28;
    gk   = 11.2;
    gl   = 2.8;

    Ena  = 50;
    Ek   = -85;
    El   = -65;
    Esyn = 0;

    % persistent sodium
    theta_mp = -40;
    sigma_mp = -6;
    theta_hp = -48;
    sigma_hp = 6;
    taumax_hp = 10000;

    mp_inf = 1/(1 + exp((v - theta_mp)/sigma_mp));
    hp_inf = 1/(1 + exp((v - theta_hp)/sigma_hp));
    tau_hp = taumax_hp / cosh((v - theta_hp)/(2*sigma_hp));

    Inap = gnap * mp_inf * hp * (v - Ena);

    % transient sodium
    theta_m = -34;
    sigma_m = -5;
    m_inf = 1/(1 + exp((v - theta_m)/sigma_m));

    Ina = gna * (m_inf^3) * (1 - n) * (v - Ena);

    % potassium
    theta_n = -29;
    sigma_n = -4;
    taumax_n = 10;

    Ik = gk * (n^4) * (v - Ek);
    n_inf = 1/(1 + exp((v - theta_n)/sigma_n));
    tau_n = taumax_n / cosh((v - theta_n)/(2*sigma_n));

    % leak
    Il = gl * (v - El);

    %% Motor pool
    r = 0.001;
    Tmax = 1;
    VT = 2;
    Kp = 5;

    NT = Tmax / (1 + exp(-(v - VT)/Kp));

    %% Lung volume
    E1 = 0.0025;
    E2 = 0.4;
    Vol0 = 2;
    dvolrhs = max(0, -E1*(vollung - Vol0) + E2*alpha);

    %% Lung oxygen
    PO2ext = (760 - 47)*0.21;
    R = 62.364;
    Temp = 310;
    taulb = 500;

    %% Blood oxygen
    M = 8e-6;
    Hb = 150;
    volblood = 5;
    eta = Hb*1.36;
    gamma = volblood/22400;
    betaO2 = 0.03;

    c = 2.5;
    K = 26;

    SaO2 = (PO2blood^c)/(PO2blood^c + K^c);
    CaO2 = eta*SaO2 + betaO2*PO2blood;
    partial = (c*PO2blood^(c-1)) * ...
              (1/(PO2blood^c + K^c) - (PO2blood^c)/((PO2blood^c + K^c)^2));

    Jlb = (1/taulb) * (PO2lung - PO2blood) * (vollung/(R*Temp));
    Jbt = M * CaO2 * gamma;

    %% Chemosensory feedback
    gtonic = 0.3 * (1 - tanh((PO2blood - 85)/30));
    Itonic = gtonic * (v - Esyn);

    %% Differential equations
    z(1) = (-Inap - Ina - Ik - Il - Itonic) / C;
    z(2) = (n_inf - n) / tau_n;
    z(3) = (hp_inf - hp) / tau_hp;
    z(4) = r*NT*(1 - alpha) - r*alpha;
    z(5) = -E1*(vollung - Vol0) + E2*alpha;
    z(6) = (1/vollung)*(PO2ext - PO2lung)*dvolrhs - Jlb*(R*Temp/vollung);
    z(7) = (Jlb - Jbt) / (gamma*(betaO2 + eta*partial));
end