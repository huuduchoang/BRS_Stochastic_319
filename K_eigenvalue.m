% Script to calculate eigenvalues of the K-state matrix A_K(V)
% for voltages between -85 mV and 50 mV,
% and find the largest eigenvalue at each voltage.

clear; clc;

%% Voltage range
Vvals = linspace(-85, 50, 500);
nV = numel(Vvals);

%% Preallocate storage
eigvals_K   = zeros(5, nV);
lambda_max  = zeros(1, nV);

%% Compute eigenvalues at each voltage
for i = 1:nV
    V = Vvals(i);
    A = AK_matrix(V);

    lam = eig(A);

    % Sort by real part, descending
    [~, idx] = sort(abs(real(lam)), 'descend');
    lam = lam(idx);

    eigvals_K(:, i) = lam;
    lambda_max(i)   = lam(1);   % largest eigenvalue by real part
end


%% Plot largest eigenvalue
figure;
plot(Vvals, abs(real(lambda_max)), 'LineWidth', 2);
xlabel('Voltage V (mV)');
ylabel('Largest eigenvalue');
title('Largest Eigenvalue of A_K(V) vs Voltage');
grid on;
box on;

%% Display a few sample values
sample_V = [-85, -65, -45, -25, 0, 25, 50];
disp('Largest eigenvalue of A_K(V) at sample voltages:');
for k = 1:numel(sample_V)
    V = sample_V(k);
    A = AK_matrix(V);
    lam = eig(A);
    [~, idx] = max(real(lam));
    fprintf('V = %6.1f mV   lambda_max = %.12g\n', V, lam(idx));
end

%% ----- Local functions -----

function A = AK_matrix(V)
    an = alphan_pt(V);
    bn = betan_pt(V);

    A = [ -4*an,          bn,              0,              0,      0;
            4*an, -3*an - bn,           2*bn,              0,      0;
               0,        3*an, -2*an - 2*bn,           3*bn,      0;
               0,           0,           2*an,    -an - 3*bn,   4*bn;
               0,           0,              0,             an, -4*bn ];
end

function a = alphan_pt(V)
    theta_n  = -29;   % mV
    sigma_n  = -4;    % mV
    taun_bar = 10;    % ms

    a = (1/(2*taun_bar)) * exp(-(V - theta_n)/(2*sigma_n));
end

function b = betan_pt(V)
    theta_n  = -29;   % mV
    sigma_n  = -4;    % mV
    taun_bar = 10;    % ms

    b = (1/(2*taun_bar)) * exp((V - theta_n)/(2*sigma_n));
end