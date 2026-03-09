% Script to reproduce Figure 1
% calls 'closedloop.m' 

clear all

%% Initial conditions [V n h alpha vollung PO2lung PO2blood]

inits = [-61.2345 0.000316 0.71447 4.9339e-05 2.0128 116.7726 115.4154];

tf = 360000;

options = odeset('RelTol',1e-9,'AbsTol',1e-9);

[t,u] = ode15s('closedloop_modified',[0 tf],inits,options);

v = u(:,1);
alpha = u(:,4);
vollung = u(:,5);
PO2lung = u(:,6);
PO2blood = u(:,7);
gtonic = 0.3*(1-tanh((PO2blood-85)./30));


%% make plots

set(0,'DefaultAxesFontSize',24)

figure(1)
clf

lw = 3;

subplot(2,3,1)
plot(t/1000, v, 'k', 'LineWidth', lw)
set(gca, 'Box', 'off', 'TickDir', 'out', ...
         'XLimMode', 'auto', 'YLimMode', 'auto', ...
         'XTickMode', 'auto', 'YTickMode', 'auto')
ylabel('$V$', 'Interpreter', 'latex')
xlabel('$t$ (s)', 'Interpreter', 'latex')

subplot(2,3,2)
plot(t/1000, alpha, 'k', 'LineWidth', lw)
set(gca, 'Box', 'off', 'TickDir', 'out', ...
         'XLimMode', 'auto', 'YLimMode', 'auto', ...
         'XTickMode', 'auto', 'YTickMode', 'auto')
ylabel('$\alpha$', 'Interpreter', 'latex')
xlabel('$t$ (s)', 'Interpreter', 'latex')

subplot(2,3,3)
plot(t/1000, vollung, 'k', 'LineWidth', lw)
set(gca, 'Box', 'off', 'TickDir', 'out', ...
         'XLimMode', 'auto', 'YLimMode', 'auto', ...
         'XTickMode', 'auto', 'YTickMode', 'auto')
ylabel('$\mathrm{vol}_\mathrm{L}$', 'Interpreter', 'latex')
xlabel('$t$ (s)', 'Interpreter', 'latex')

subplot(2,3,4)
plot(t/1000, gtonic, 'k', 'LineWidth', lw)
set(gca, 'Box', 'off', 'TickDir', 'out', ...
         'XLimMode', 'auto', 'YLimMode', 'auto', ...
         'XTickMode', 'auto', 'YTickMode', 'auto')
ylabel('$g_\mathrm{tonic}$', 'Interpreter', 'latex')
xlabel('$t$ (s)', 'Interpreter', 'latex')

subplot(2,3,5)
plot(t/1000, PO2blood, 'k', 'LineWidth', lw)
set(gca, 'Box', 'off', 'TickDir', 'out', ...
         'XLimMode', 'auto', 'YLimMode', 'auto', ...
         'XTickMode', 'auto', 'YTickMode', 'auto')
ylabel('$P_\mathrm{A}\mathrm{O}_2$', 'Interpreter', 'latex')
xlabel('$t$ (s)', 'Interpreter', 'latex')

subplot(2,3,6)
plot(t/1000, PO2lung, 'k', 'LineWidth', lw)
set(gca, 'Box', 'off', 'TickDir', 'out', ...
         'XLimMode', 'auto', 'YLimMode', 'auto', ...
         'XTickMode', 'auto', 'YTickMode', 'auto')
ylabel('$P_\mathrm{a}\mathrm{O}_2$', 'Interpreter', 'latex')
xlabel('$t$ (s)', 'Interpreter', 'latex')

set(gcf, 'Position', get(0, 'screensize'))