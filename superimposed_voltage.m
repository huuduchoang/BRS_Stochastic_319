%% Results Generator: Stochastic Convergence
clear; clc;

% --- Settings ---
Area_test = 1e4; 
tf = 20000; % Time to run (ms)
dt_list = [0.1, 0.05, 0.01];
colors = {'#0072BD', '#D95319', '#7E2F8E'}; % Blue, Orange, Purple

% --- Initial Conditions ---
% u = [v; n; h; alpha; vollung; PO2lung; PO2blood]
old = [-61.2345; 0.000316; 0.71447; 4.9339e-05; 2.0128; 116.7726; 115.4154];

% Lifting 7D initial conditions to 11D for the Stochastic Model
n_scalar = old(2);
Kinit = [(1-n_scalar)^4; 4*n_scalar*(1-n_scalar)^3; 6*n_scalar^2*(1-n_scalar)^2; 4*n_scalar^3*(1-n_scalar); n_scalar^4];
inits_11D = [old(1); Kinit; old(3); old(4); old(5); old(6); old(7)];

% Create Figure (Adjusted width for a single plot)
figure('Position', [100, 100, 800, 500]);
hold on;

%% Plot: Voltage Convergence (Multiple dt)
for i = 1:length(dt_list)
    dt_val = dt_list(i);
    
    % Call your custom solver directly for each dt
    [t_v, u_v] = simulate_closedloop_Kfull(tf, dt_val, inits_11D, Area_test);
    
    % u_v(:, 1) extracts just the voltage column
    plot(t_v, u_v(:, 1), 'Color', colors{i}, 'LineWidth', 1, ...
        'DisplayName', ['dt = ' num2str(dt_val)]);
end

title('Voltage: Time-Step Convergence');
xlabel('Time (s)'); 
ylabel('Membrane Potential (mV)');
legend('Location', 'northeast'); 
grid on;
