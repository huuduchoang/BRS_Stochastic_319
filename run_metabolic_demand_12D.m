function results = run_metabolic_demand_12D(varargin)
% RUN_METABOLIC_DEMAND_12D
% Closed-loop-only driver for the 12D model, with progress diagnostics
% and checkpoint saving.
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

    %% ---------------------------
    %  User-adjustable defaults
    %  ---------------------------
    p = inputParser;

    addParameter(p, 'M0',              8e-6);
    addParameter(p, 'Mvals',           2e-6:0.1e-6:18e-6);
    addParameter(p, 'tf_seg',          6e4);              % ms
    addParameter(p, 'dt',              0.05);             % ms
    addParameter(p, 'Area',            1e4/18);           % gives about 1e4 K channels if NK = round(18*Area)
    addParameter(p, 'noise_on',        false);
    addParameter(p, 'nRuns',           1);
    addParameter(p, 'rng_seed',        1);
    addParameter(p, 'save_results',    false);
    addParameter(p, 'save_name',       'fig8_12D_closedloop_results.mat');
    addParameter(p, 'checkpoint_name', 'fig8_12D_closedloop_checkpoint.mat');
    addParameter(p, 'use_parallel',    true);
    addParameter(p, 'save_every_M',    true);
    addParameter(p, 'verbose',         true);

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
    run_completed   = false(nRuns, nM);

    total_jobs = nRuns * nM;
    t_global   = tic;

    if S.verbose
        fprintf('\n=========================================\n');
        fprintf('Starting 12D metabolic-demand experiment\n');
        fprintf('nRuns = %d, nM = %d, total jobs = %d\n', nRuns, nM, total_jobs);
        fprintf('noise_on = %d, use_parallel = %d\n', S.noise_on, S.use_parallel);
        fprintf('Checkpoint file: %s\n', S.checkpoint_name);
        fprintf('=========================================\n\n');
    end

    %% ---------------------------
    %  Run experiment
    %  ---------------------------
    if S.use_parallel && nRuns > 1
        dq = parallel.pool.DataQueue;
        afterEach(dq, @update_from_packet);

        parfor r = 1:nRuns
            run_one_realization(r, S, inits12, Mvals, nM, dq, true);
        end
    else
        for r = 1:nRuns
            run_one_realization(r, S, inits12, Mvals, nM, @update_from_packet, false);
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
    results.run_completed    = run_completed;
    results.settings         = S;
    results.elapsed_seconds  = toc(t_global);

    %% ---------------------------
    %  Optional final save
    %  ---------------------------
    if S.save_results
        save(S.save_name, 'results', '-v7.3');
        if S.verbose
            fprintf('\nSaved final results to %s\n', S.save_name);
        end
    end

    %% ---------------------------
    %  Nested callback for progress/checkpoint
    %  ---------------------------
    function update_from_packet(packet)
        switch packet.type
            case 'run_start'
                if S.verbose
                    fprintf('Starting run %d of %d\n', packet.r, nRuns);
                end

            case 'M_done'
                avgPO2blood_all(packet.r, packet.ix) = packet.avgPO2blood;
                u_final_all(:, packet.ix, packet.r)  = packet.ufinal;
                run_completed(packet.r, packet.ix)   = true;

                completed_jobs = nnz(run_completed);
                pct_complete   = 100 * completed_jobs / total_jobs;
                elapsed_sec    = toc(t_global);

                if S.verbose
                    fprintf(['Run %d/%d, M index %d/%d, M = %.6g, ' ...
                             'avgPO2blood = %.6f, progress = %d/%d (%.1f%%), elapsed = %.1f s\n'], ...
                             packet.r, nRuns, packet.ix, nM, packet.M, ...
                             packet.avgPO2blood, completed_jobs, total_jobs, pct_complete, elapsed_sec);
                end

                if S.save_every_M
                    checkpoint = struct();
                    checkpoint.Mvals            = Mvals;
                    checkpoint.avgPO2blood_all  = avgPO2blood_all;
                    checkpoint.u_final_all      = u_final_all;
                    checkpoint.run_completed    = run_completed;
                    checkpoint.settings         = S;
                    checkpoint.elapsed_seconds  = elapsed_sec;

                    save(S.checkpoint_name, 'checkpoint', '-v7.3');
                end

            case 'run_done'
                if S.verbose
                    fprintf('Completed run %d of %d\n', packet.r, nRuns);
                end
        end
    end
end


function run_one_realization(r, S, inits12, Mvals, nM, reporter, useDataQueue)
    if S.noise_on
        rng(S.rng_seed + r - 1, 'twister');
    end

    send_or_call(reporter, useDataQueue, struct('type', 'run_start', 'r', r));

    params = struct();
    params.noise_on = S.noise_on;
    params.M        = S.M0;

    % ------------------------
    % Baseline burn-in at M0
    % ------------------------
    [~, U0] = simulate_closedloop_Kfull(1e5, S.dt, inits12, S.Area, params);
    inits1 = U0(end,:).';

    [~, U1] = simulate_closedloop_Kfull(1e5, S.dt, inits1, S.Area, params);
    initsM = U1(end,:).';

    % ------------------------
    % Sweep over M values
    % ------------------------
    for ix = 1:nM
        params.M = Mvals(ix);

        [~, U2]  = simulate_closedloop_Kfull(S.tf_seg, S.dt, initsM,       S.Area, params);
        [~, U3]  = simulate_closedloop_Kfull(S.tf_seg, S.dt, U2(end,:).',  S.Area, params);
        [~, U4]  = simulate_closedloop_Kfull(S.tf_seg, S.dt, U3(end,:).',  S.Area, params);
        [~, U5]  = simulate_closedloop_Kfull(S.tf_seg, S.dt, U4(end,:).',  S.Area, params);
        [t6, U6] = simulate_closedloop_Kfull(S.tf_seg, S.dt, U5(end,:).',  S.Area, params);

        % Average PO2blood over final segment only
        PO2blood6    = U6(:,12);
        avgPO2blood  = trapz(t6, PO2blood6) / (t6(end) - t6(1));
        ufinal       = U6(end,:).';

        % Send partial result back to client
        packet = struct();
        packet.type         = 'M_done';
        packet.r            = r;
        packet.ix           = ix;
        packet.M            = params.M;
        packet.avgPO2blood  = avgPO2blood;
        packet.ufinal       = ufinal;

        send_or_call(reporter, useDataQueue, packet);

        % continuation in M
        initsM = ufinal;
    end

    send_or_call(reporter, useDataQueue, struct('type', 'run_done', 'r', r));
end


function send_or_call(reporter, useDataQueue, packet)
    if useDataQueue
        send(reporter, packet);
    else
        reporter(packet);
    end
end


function inits12 = build_initial_12D_from_old7(old7)
% old7 = [v; n; hp; alpha; vollung; PO2lung; PO2blood]

    v0        = old7(1);
    n_scalar  = old7(2);
    hp0       = old7(3);
    alpha0    = old7(4);
    vollung0  = old7(5);
    PO2lung0  = old7(6);
    PO2blood0 = old7(7);

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