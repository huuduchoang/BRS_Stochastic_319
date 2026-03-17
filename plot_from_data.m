%% 4. Plot only from 140000 ms to the end
t_start_ms = 0;
t_start_s  = t_start_ms / 1000;
% Calculate eigenvalues of K for different voltages
%% Figure 1: 2x3 Closed-loop Grid (Superimposed)
fig1 = figure('Position',[100 100 1000 600],'Color','w');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% Colors
stoch_color = [0 0.4470 0.7410]; % Blue
det_color   = 'k';               % Solid Black

% 1. CPG (V)
nexttile(1); hold on
plot(time_s, v, 'Color', stoch_color, 'LineWidth', 1.5, 'DisplayName', 'Stochastic');
plot(time_s_det, v_det, 'Color', det_color, 'LineWidth', 1.2, 'DisplayName', 'Deterministic');
ylabel('V (mV)'); title('CPG (V)');
xlim([t_start_s tf/1000]); box on; grid on
legend('Location','best');

% 2. Motor Pool
nexttile(2); hold on
plot(time_s, alpha, 'Color', stoch_color, 'LineWidth', 1.5);
plot(time_s_det, alpha_det, 'Color', det_color, 'LineWidth', 1.2);
ylabel('\alpha'); title('Motor Pool');
xlim([t_start_s tf/1000]); box on; grid on

% 3. Lung Volume
nexttile(3); hold on
plot(time_s, vollung, 'Color', stoch_color, 'LineWidth', 1.5);
plot(time_s_det, vollung_det, 'Color', det_color, 'LineWidth', 1.2);
ylabel('Vol_{lung}'); title('Lung Volume');
xlim([t_start_s tf/1000]); box on; grid on

% 4. Chemosensation
nexttile(4); hold on
plot(time_s, gtonic, 'Color', stoch_color, 'LineWidth', 1.5);
plot(time_s_det, gtonic_det, 'Color', det_color, 'LineWidth', 1.2);
ylabel('g_{tonic}'); title('Chemosensation');
xlim([t_start_s tf/1000]); box on; grid on

% 5. Blood Oxygen
nexttile(5); hold on
plot(time_s, PO2blood, 'Color', stoch_color, 'LineWidth', 1.5);
plot(time_s_det, PO2blood_det, 'Color', det_color, 'LineWidth', 1.2);
ylabel('P_{aO_2} (mmHg)'); title('Blood Oxygen');
xlabel('Time (s)');
xlim([t_start_s tf/1000]); box on; grid on

% 6. Lung Oxygen
nexttile(6); hold on
plot(time_s, PO2lung, 'Color', stoch_color, 'LineWidth', 1.5);
plot(time_s_det, PO2lung_det, 'Color', det_color, 'LineWidth', 1.2);
ylabel('P_{O_2}^{lung} (mmHg)'); title('Lung Oxygen');
xlabel('Time (s)');
xlim([t_start_s tf/1000]); box on; grid on

%% 5. Figure 2: K activation comparison (n4 vs n^4)
fig2 = figure('Position',[150 150 1000 400],'Color','w');
hold on
plot(time_s, n4, 'Color', stoch_color, 'LineWidth', 1.5, 'DisplayName', 'Stochastic n_4');
plot(time_s_det, n_det.^4, 'Color', det_color, 'LineWidth', 1.2, 'DisplayName', 'Deterministic n^4');
ylabel('K activation');
xlabel('Time (s)');
title('Stochastic n_4 vs Deterministic n^4');
xlim([t_start_s tf/1000]);
box on
grid on
legend('Location','best')