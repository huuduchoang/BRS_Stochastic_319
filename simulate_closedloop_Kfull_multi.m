function [t, U] = simulate_closedloop_Kfull_multi(tf, dt, inits, Area, params)
% simulate_closedloop_Kfull_multi
% Multi-neuron Full K-channel noise:
%   - Euler-Maruyama for 5-state potassium occupancy vector for N neurons
%   - ode15s for the remaining deterministic variables
%
% State mapping per time step (Length: 9*N + 4):
%   v        (1:N)
%   K-states (N+1 : 6*N) -> 5 states per neuron
%   hp       (6*N+1 : 7*N)
%   hf       (7*N+1 : 8*N)
%   s        (8*N+1 : 9*N)  -> synaptic gating
%   alpha    (9*N+1)
%   vollung  (9*N+2)
%   PO2lung  (9*N+3)
%   PO2blood (9*N+4)

    t  = (0:dt:tf).';
    nt = length(t);
    N  = params.N;

    % number of K channels per neuron
    NK = round(18*Area);

    U      = zeros(nt, 9*N + 4);
    U(1,:) = inits(:).';

    % Relaxed tolerances slightly for speed in population models
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9); 

    for k = 1:nt-1
        uk = U(k,:).';
        tk = t(k);

        % Extract current voltage and K-state vectors
        v_old = uk(1:N);
        K_old = reshape(uk(N+1 : 6*N), 5, N); % 5xN matrix
        
        % Extract deterministic states for the ODE step
        hp_old = uk(6*N+1 : 7*N);
        hf_old = uk(7*N+1 : 8*N);
        s_old  = uk(8*N+1 : 9*N);
        globals_old = uk(9*N+1 : end); 

        %------------------------------------------
        % Step 1: Euler-Maruyama for full K block (N neurons)
        %------------------------------------------
        K_new = zeros(5, N);
        for i = 1:N
            xi = randn(8,1);
            if isfield(params, 'noise_on') && params.noise_on
                K_new(:,i) = K_old(:,i) ...
                   + dt * AK_matrix(v_old(i)) * K_old(:,i) ...
                   + sqrt(dt) * DKfull_PT(v_old(i), K_old(:,i), NK, xi);
            else
                K_new(:,i) = K_old(:,i) ...
                   + dt * AK_matrix(v_old(i)) * K_old(:,i);
            end
            
            % % Simple projection to keep probabilities valid and sum to 1
            % K_new(:,i) = max(K_new(:,i), 0);
            % sK = sum(K_new(:,i));
            % if sK > 0
            %     K_new(:,i) = K_new(:,i) / sK;
            % else
            %     K_new(:,i) = [1; 0; 0; 0; 0];
            % end
        end

        %------------------------------------------
        % Step 2: deterministic substep for rest
        %------------------------------------------
        % x0 = [v(1:N); hp(1:N); hf(1:N); s(1:N); alpha; vollung; PO2lung; PO2blood]
        % Total size of ODE state: 4*N + 4
        x0 = [v_old; hp_old; hf_old; s_old; globals_old];

        rhs = @(tt,x) closedloop_rest_rhs_Kfull_multi(tt, x, K_new, params);

        [~,Xsol] = ode15s(rhs, [tk t(k+1)], x0, opts);
        x_new = Xsol(end,:).';

        %------------------------------------------
        % Reassemble into full state vector
        %------------------------------------------
        u_new = zeros(9*N + 4, 1);
        
        u_new(1:N)         = x_new(1:N);           % v
        u_new(N+1 : 6*N)   = K_new(:);             % K states (flattened to vector)
        u_new(6*N+1 : 7*N) = x_new(N+1 : 2*N);     % hp
        u_new(7*N+1 : 8*N) = x_new(2*N+1 : 3*N);   % hf
        u_new(8*N+1 : 9*N) = x_new(3*N+1 : 4*N);   % s
        u_new(9*N+1 : end) = x_new(4*N+1 : end);   % globals (alpha, vollung, PO2lung, PO2blood)

        U(k+1,:) = u_new.';
    end
end

% =========================================================================
% ODE Right-Hand Side Function
% =========================================================================
function z = closedloop_rest_rhs_Kfull_multi(~, x, Kstates, params)
% x = [v(1:N); hp(1:N); hf(1:N); s(1:N); alpha; vollung; PO2lung; PO2blood]
% Kstates = 5xN matrix, frozen during this substep

    N = params.N;

    v        = x(1 : N);
    hp       = x(N+1 : 2*N);
    hf       = x(2*N+1 : 3*N);
    s        = x(3*N+1 : 4*N);
    alpha    = x(4*N+1);
    vollung  = x(4*N+2);
    PO2lung  = x(4*N+3);
    PO2blood = x(4*N+4);

    n4 = Kstates(5, :).';   % conducting K state for all N neurons (column vector)

    z = zeros(4*N + 4, 1);

    %% Unpack Heterogeneous Parameters
    El = params.Eleak;
    gnap = params.gnap;
    gsyn = params.gsyn;
    
    phi_vec     = params.phi;
    thetaO2_vec = params.thetaO2;
    sigmaO2_vec = params.sigmaO2;

    %% CPG Parameters
    C = 21;
    gna  = 28;
    gk   = 11.2;
    gl   = 2.8;

    Ena  = 50;
    Ek   = -85;
    Esyn = 0;

    % persistent sodium
    theta_mp = -40; sigma_mp = -6;
    theta_hp = -48; sigma_hp = 6;
    taumax_hp = 10000;

    mp_inf = 1 ./ (1+exp((v-theta_mp)/sigma_mp));
    hp_inf = 1 ./ (1+exp((v-theta_hp)/sigma_hp));
    tau_hp = taumax_hp ./ cosh((v-theta_hp)/(2*sigma_hp));

    Inap = gnap .* mp_inf .* hp .* (v-Ena);

    % transient sodium
    theta_m = -34; sigma_m = -5;
    m_inf   = 1 ./ (1+exp((v-theta_m)/sigma_m));

    theta_n  = -29; sigma_n  = -4; taumax_n = 10;
    n_inf  = 1 ./ (1+exp((v-theta_n)/sigma_n));
    hf_inf = 1 - n_inf;
    tau_hf = taumax_n ./ cosh((v-theta_n)/(2*sigma_n));

    Ina = gna .* (m_inf.^3) .* hf .* (v-Ena);

    % potassium: use open-channel fraction n4
    Ik = gk .* n4 .* (v-Ek);

    % leak
    Il = gl .* (v-El);

    % synaptic gating
    tau_s = 5; k_r = 1;
    thetas = -10; sigmas = -5;
    s_inf = 1 ./ (1 + exp((v - thetas)/sigmas));
    g_in = gsyn.' * s; 
    Isyn = g_in .* (v - Esyn);

    %% Chemosensory feedback (Affects all N neurons)
    gtonic_vec = phi_vec .* (1 - tanh((PO2blood - thetaO2_vec) ./ sigmaO2_vec));
    Itonic = gtonic_vec .* (v - Esyn);

    %% Motor pool (Pools from all N neurons)
    r = 0.001; Tmax = 1; VT = 2; Kp = 5;
    NT = 1 ./ (1+exp(-(v-VT)/Kp));
    Tpop = Tmax * mean(NT); 

    %% Lung volume
    E1 = 0.0025; E2 = 0.4; Vol0 = 2;
    dvolrhs = -E1*(vollung-Vol0)+E2*alpha;

    %% Lung oxygen
    PO2ext = (760-47)*0.21; R = 62.364; Temp = 310; taulb = 500;

    %% Blood oxygen
    Hb = 150; volblood = 5; eta = Hb*1.36; gamma = volblood/22400; betaO2 = 0.03;
    c = 2.5; K = 26;

    SaO2 = (PO2blood^c)/(PO2blood^c + K^c);
    CaO2 = eta*SaO2 + betaO2*PO2blood;
    partial = (c*PO2blood^(c-1)) * (1/(PO2blood^c+K^c) - (PO2blood^c)/((PO2blood^c+K^c)^2));

    Jlb = (1/taulb)*(PO2lung-PO2blood)*(vollung/(R*Temp));
    Jbt = params.M*CaO2*gamma;

    %% ODEs
    z(1:N)       = (-Inap - Ina - Ik - Il - Itonic - Isyn) ./ C;
    z(N+1:2*N)   = (hp_inf - hp) ./ tau_hp;
    z(2*N+1:3*N) = (hf_inf - hf) ./ tau_hf;
    z(3*N+1:4*N) = ((1 - s) .* s_inf  -  k_r .* s) ./ tau_s;
    z(4*N+1)     = r*Tpop*(1-alpha) - r*alpha;
    z(4*N+2)     = -E1*(vollung-Vol0) + E2*alpha;
    z(4*N+3)     = (1/vollung)*(PO2ext-PO2lung)*max(0,dvolrhs) - Jlb*(R*Temp/vollung);
    z(4*N+4)     = (Jlb-Jbt)/(gamma*(betaO2+eta*partial));
end

% =========================================================================
% Helper Functions for K-Channel Noise
% =========================================================================
function A = AK_matrix(V)
    an = alphan_pt(V);
    bn = betan_pt(V);

    A = [ -4*an,          bn,              0,              0,      0;
            4*an, -3*an - bn,           2*bn,              0,      0;
               0,       3*an, -2*an - 2*bn,           3*bn,      0;
               0,          0,           2*an,    -an - 3*bn,   4*bn;
               0,          0,              0,             an, -4*bn ];
end

function D = DKfull_PT(V,Y,Nch,M)
% Full K diffusion term from the 5-state/8-edge channel model
%
% Inputs:
%   V   = voltage (scalar)
%   Y   = [n0;n1;n2;n3;n4] for a single neuron
%   Nch = total number of K channels
%   M   = randn(8,1)
%
% Output:
%   D   = 5x1 stochastic increment coefficient

    D = zeros(5,1);

    an = alphan_pt(V);
    bn = betan_pt(V);

    % numerical safeguard
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
end

function b = betan_pt(V)
    theta_n  = -29;   % mV
    sigma_n  = -4;    % mV
    taun_bar = 10;    % ms
    b = (1/(2*taun_bar)) * exp((V - theta_n)/(2*sigma_n));
end