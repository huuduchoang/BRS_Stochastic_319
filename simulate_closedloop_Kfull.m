function [t,U] = simulate_closedloop_Kfull(tf,dt,inits,Area,params)
% Full K-channel noise:
%   - Euler-Maruyama for 5-state potassium occupancy vector
%   - ode15s for the remaining deterministic variables
%
% State:
%   u = [v; n0; n1; n2; n3; n4; hp; hf; alpha; vollung; PO2lung; PO2blood]

    t  = (0:dt:tf).';
    nt = length(t);

    % number of K channels
    NK = round(18*Area);

    U      = zeros(nt,12);
    U(1,:) = inits(:).';

    opts = odeset('RelTol',1e-9,'AbsTol',1e-9);

    for k = 1:nt-1
        uk = U(k,:).';
        tk = t(k);
        % hk = t(k+1) - t(k);

        % Current voltage and K-state vector
        v_old = uk(1);
        K_old = uk(2:6);

        %------------------------------------------
        % Step 1: Euler-Maruyama for full K block
        %------------------------------------------
        xi = randn(8,1);

        % K_new = K_old ...
        %       + hk * AK_matrix(v_old) * K_old ...
        %       + sqrt(hk) * DKfull_PT(v_old, K_old, NK, xi);
        if params.noise_on
         K_new = K_old ...
              + dt * AK_matrix(v_old) * K_old + sqrt(dt) * DKfull_PT(v_old, K_old, NK, xi);
        else
            K_new = K_old ...
              + dt * AK_matrix(v_old) * K_old;
        end

        % % Numerical projection back to simplex
        % K_new = max(K_new,0);
        % sK = sum(K_new);
        % 
        % if sK <= 0
        %     K_new = [1;0;0;0;0];
        % else
        %     K_new = K_new / sK;
        % end

        %------------------------------------------
        % Step 2: deterministic substep for rest
        %------------------------------------------
        % x = [v; hp; hf; alpha; vollung; PO2lung; PO2blood]
        x0 = uk([1 7 8 9 10 11 12]);

        rhs = @(tt,x) closedloop_rest_rhs_Kfull(tt, x, K_new, params);
        % rhs = @(tt,x) closedloop_rest_rhs_Kfull(tt, x, (K_new + K_old)/2);

        [~,Xsol] = ode15s(rhs,[tk t(k+1)],x0,opts);
        x_new = Xsol(end,:).';

        %------------------------------------------
        % Reassemble
        %------------------------------------------
        u_new = zeros(12,1);
        u_new(1)    = x_new(1);   % v
        u_new(2:6)  = K_new;      % n0..n4
        u_new(7)    = x_new(2);   % hp
        u_new(8)    = x_new(3);   % hf
        u_new(9)    = x_new(4);   % alpha
        u_new(10)   = x_new(5);   % vollung
        u_new(11)   = x_new(6);   % PO2lung
        u_new(12)   = x_new(7);   % PO2blood

        U(k+1,:) = u_new.';
    end
end

function z = closedloop_rest_rhs_Kfull(~,x,Kstates, params)
% x = [v; hp; hf; alpha; vollung; PO2lung; PO2blood]
% Kstates = [n0;n1;n2;n3;n4], frozen during this substep

    v        = x(1);
    hp       = x(2);
    hf       = x(3);
    alpha    = x(4);
    vollung  = x(5);
    PO2lung  = x(6);
    PO2blood = x(7);

    n4 = Kstates(5);   % conducting K state

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

    mp_inf = 1/(1+exp((v-theta_mp)/sigma_mp));
    hp_inf = 1/(1+exp((v-theta_hp)/sigma_hp));
    tau_hp = taumax_hp/cosh((v-theta_hp)/(2*sigma_hp));

    Inap = gnap * mp_inf * hp * (v-Ena);

    % transient sodium
    theta_m = -34;
    sigma_m = -5;
    m_inf   = 1/(1+exp((v-theta_m)/sigma_m));

    % hf should behave like (1 - n) from the original model
    theta_n  = -29;
    sigma_n  = -4;
    taumax_n = 10;

    n_inf  = 1/(1+exp((v-theta_n)/sigma_n));
    hf_inf = 1 - n_inf;
    tau_hf = taumax_n/cosh((v-theta_n)/(2*sigma_n));

    Ina = gna * (m_inf^3) * hf * (v-Ena);

    % potassium: use open-channel fraction n4
    Ik = gk * n4 * (v-Ek);

    % leak
    Il = gl * (v-El);

    %% Motor pool
    r    = 0.001;
    Tmax = 1;
    VT   = 2;
    Kp   = 5;

    NT = Tmax/(1+exp(-(v-VT)/Kp));

    %% Lung volume
    E1   = 0.0025;
    E2   = 0.4;
    Vol0 = 2;

    dvolrhs = -E1*(vollung-Vol0)+E2*alpha;

    %% Lung oxygen
    PO2ext = (760-47)*0.21;
    R      = 62.364;
    Temp   = 310;
    taulb  = 500;

    %% Blood oxygen
    Hb       = 150;
    volblood = 5;
    eta      = Hb*1.36;
    gamma    = volblood/22400;
    betaO2   = 0.03;

    c = 2.5;
    K = 26;

    SaO2 = (PO2blood^c)/(PO2blood^c + K^c);
    CaO2 = eta*SaO2 + betaO2*PO2blood;
    partial = (c*PO2blood^(c-1)) * ...
              (1/(PO2blood^c+K^c) - (PO2blood^c)/((PO2blood^c+K^c)^2));

    Jlb = (1/taulb)*(PO2lung-PO2blood)*(vollung/(R*Temp));
    Jbt = params.M*CaO2*gamma;

    %% Chemosensory feedback
    gtonic = 0.3*(1-tanh((PO2blood-85)/30));
    Itonic = gtonic*(v-Esyn);

    %% ODEs
    z(1) = (-Inap - Ina - Ik - Il - Itonic)/C;
    z(2) = (hp_inf - hp)/tau_hp;
    z(3) = (hf_inf - hf)/tau_hf;
    z(4) = r*NT*(1-alpha) - r*alpha;
    z(5) = -E1*(vollung-Vol0)+E2*alpha;
    z(6) = (1/vollung)*(PO2ext-PO2lung)*max(0,dvolrhs) - Jlb*(R*Temp/vollung);
    z(7) = (Jlb-Jbt)/(gamma*(betaO2+eta*partial));
end

function A = AK_matrix(V)
    an = alphan_pt(V);
    bn = betan_pt(V);

    A = [ -4*an,          bn,              0,              0,      0;
            4*an, -3*an - bn,           2*bn,              0,      0;
               0,        3*an, -2*an - 2*bn,           3*bn,      0;
               0,           0,           2*an,    -an - 3*bn,   4*bn;
               0,           0,              0,             an, -4*bn ];
end

function D = DKfull_PT(V,Y,Nch,M)
% Full K diffusion term from the 5-state/8-edge channel model
%
% Inputs:
%   V   = voltage
%   Y   = [n0;n1;n2;n3;n4]
%   Nch = total number of K channels
%   M   = randn(8,1)
%
% Output:
%   D   = 5x1 stochastic increment coefficient

    D = zeros(5,1);

    an = alphan_pt(V);
    bn = betan_pt(V);

    % numerical safeguard, following the HHSS14D style
    Y = abs(Y);

    N8 = zeros(8,1);

    % directed edges
    N8(1) = sqrt(4*an*Y(1))*M(1);   % 0 -> 1
    N8(2) = sqrt(bn*Y(2))*M(2);     % 1 -> 0
    N8(3) = sqrt(3*an*Y(2))*M(3);   % 1 -> 2
    N8(4) = sqrt(2*bn*Y(3))*M(4);   % 2 -> 1
    N8(5) = sqrt(2*an*Y(3))*M(5);   % 2 -> 3
    N8(6) = sqrt(3*bn*Y(4))*M(6);   % 3 -> 2
    N8(7) = sqrt(an*Y(4))*M(7);     % 3 -> 4
    N8(8) = sqrt(4*bn*Y(5))*M(8);   % 4 -> 3

    % net node balances
    D(1) = -N8(1) + N8(2);
    D(2) = -N8(2) - N8(3) + N8(1) + N8(4);
    D(3) = -N8(4) - N8(5) + N8(3) + N8(6);
    D(4) = -N8(6) - N8(7) + N8(5) + N8(8);
    D(5) = -N8(8) + N8(7);

    D = D / sqrt(Nch);
end

function a = alphan_pt(V)
    theta_n  = -29;   % mV
    sigma_n  = -4;    % mV
    taun_bar = 10;    % ms

    a = (1/(2*taun_bar)) * exp(-(V - theta_n)/(2*sigma_n));
    % equivalent to: 0.05 * exp((V + 29)/8)
end

function b = betan_pt(V)
    theta_n  = -29;   % mV
    sigma_n  = -4;    % mV
    taun_bar = 10;    % ms

    b = (1/(2*taun_bar)) * exp((V - theta_n)/(2*sigma_n));
    % equivalent to: 0.05 * exp(-(V + 29)/8
end
