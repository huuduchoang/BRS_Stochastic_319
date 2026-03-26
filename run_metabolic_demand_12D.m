function results = run_metabolic_demand_12D(varargin)
% RUN_FIG8_12D_CLOSEDLOOP
% Closed-loop-only driver for the 12D model, analogous to the old Figure 8 experiment.
%
% Assumes the simulator has signature:
%   [t,U] = simulate_closedloop_Kfull(tf, dt, inits, Area, params)
%
% with params containing at least:
%   params.M
%   params.noise_on
%
% 12D state ordering:
%   u = [v; n0; n1; n2; n3; n4; hp; hf; alpha; vollung; PO2lung; PO2blood]
%
% Output:
%   results struct with fields:
%       .Mvals
%       .avgPO2blood          (mean across realizations)
%       .avgPO2blood_all      (nRuns x nM)
%       .u_final_all          (12 x nM x nRuns final states)
%       .settings
%
% Example:
%   results = run_fig8_12D_closedloop();
%
%   figure;
%   plot(results.Mvals, results.avgPO2blood, 'k', 'LineWidth', 2);
%   xlabel('M'); ylabel('Average P_{aO_2}');
%   grid on;

    %% ---------------------------
    %  User-adjustable defaults
    %  ---------------------------
    p = inputParser;

    addParameter(p, 'M0',         8e-6);
    addParameter(p, 'Mvals',      2e-6:0.1e-6:18e-6);
    addParameter(p, 'tf_seg',     6e4);         % ms, same spirit as old script
    addParameter(p, 'dt',         0.05);        % ms
    addParameter(p, 'Area',       1e4/18);      % gives about 1e4 K channels if NK = round(18*Area)
    addParameter(p, 'noise_on',   false);       % start deterministic first
    addParameter(p, 'nRuns',      1);           % >1 useful only when noise_on = true
    addParameter(p, 'rng_seed',   1);
    addParameter(p, 'save_results', false);
    addParameter(p, 'save_name',  'fig8_12D_closedloop_results.mat');

    parse(p, varargin{:});
    S = p.Results;

    %% ---------------------------
    %  Initial condition
    %  ---------------------------
    % Old 7D deterministic IC:
    % [v; n; hp; alpha; vollung; PO2lung; PO2blood]
   old7 = [-50.3697203300953;
        0.00462275707173923;
        0.750556122453219;
        3.93159280387281e-05;
        2.00894512526427;
        94.4645088086363;
        93.3277058316367];

    inits12 = build_initial_12D_from_old7(old7);

    %% ---------------------------
    %  Storage
    %  ---------------------------
    Mvals = S.Mvals(:).';
    nM    = numel(Mvals);
    nRuns = S.nRuns;

    avgPO2blood_all = nan(nRuns, nM);
    u_final_all     = nan(12, nM, nRuns);

    %% ---------------------------
    %  Run experiment
    %  ---------------------------
    for r = 1:nRuns
        if S.noise_on
            rng(S.rng_seed + r - 1, 'twister');
        end

        params = struct();
        params.noise_on = S.noise_on;
        params.M        = S.M0;

        % ------------------------
        % Baseline burn-in at M0
        % ------------------------
        [~, U0] = simulate_closedloop_Kfull(S.tf_seg, S.dt, inits12, S.Area, params);
        inits1 = U0(end,:).';

        [~, U1] = simulate_closedloop_Kfull(S.tf_seg, S.dt, inits1, S.Area, params);
        initsM = U1(end,:).';

        % ------------------------
        % Sweep over M values
        % ------------------------
        for ix = 1:nM
            params.M = Mvals(ix);

            fprintf('Run %d/%d, M = %.6g\n', r, nRuns, params.M);

            % Replicate the old "chain several long segments" approach
            [~, U2] = simulate_closedloop_Kfull(S.tf_seg, S.dt, initsM,         S.Area, params);
            [~, U3] = simulate_closedloop_Kfull(S.tf_seg, S.dt, U2(end,:).',    S.Area, params);
            [~, U4] = simulate_closedloop_Kfull(S.tf_seg, S.dt, U3(end,:).',    S.Area, params);
            [~, U5] = simulate_closedloop_Kfull(S.tf_seg, S.dt, U4(end,:).',    S.Area, params);
            [t6, U6] = simulate_closedloop_Kfull(S.tf_seg, S.dt, U5(end,:).',   S.Area, params);

            % Average PO2blood over the final segment only
            PO2blood6 = U6(:,12);
            avgPO2blood_all(r, ix) = trapz(t6, PO2blood6) / (t6(end) - t6(1));

            % Save final state for continuation and diagnostics
            u_final_all(:, ix, r) = U6(end,:).';
            initsM = U6(end,:).';
        end
    end

    %% ---------------------------
    %  Aggregate results
    %  ---------------------------
    results = struct();
    results.Mvals            = Mvals;
    results.avgPO2blood      = mean(avgPO2blood_all, 1, 'omitnan');
    results.avgPO2blood_all  = avgPO2blood_all;
    results.u_final_all      = u_final_all;
    results.settings         = S;

    %% ---------------------------
    %  Optional save
    %  ---------------------------
    if S.save_results
        save(S.save_name, 'results', '-v7.3');
        fprintf('Saved results to %s\n', S.save_name);
    end
end


function inits12 = build_initial_12D_from_old7(old7)
% old7 = [v; n; hp; alpha; vollung; PO2lung; PO2blood]

    v0       = old7(1);
    n_scalar = old7(2);
    hp0      = old7(3);
    alpha0   = old7(4);
    vollung0 = old7(5);
    PO2lung0 = old7(6);
    PO2blood0= old7(7);

    % Binomial lift of scalar n to 5-state K occupancy
    Kinit = [
        (1 - n_scalar)^4;
        4 * n_scalar * (1 - n_scalar)^3;
        6 * n_scalar^2 * (1 - n_scalar)^2;
        4 * n_scalar^3 * (1 - n_scalar);
        n_scalar^4
    ];

    % Fast sodium factor consistent with old reduced model
    hf0 = 1 - n_scalar;

    % 12D state:
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
end