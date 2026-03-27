function results = run_metabolic_demand_12D(varargin)
% RUN_METABOLIC_DEMAND_12D
% Closed-loop-only driver for the 12D model, with checkpoint saving and progress display.
%
% Assumes the simulator has signature:
%   [t,U] = simulate_closedloop_Kfull(tf, dt, inits, Area, params)
%
% with params containing at least:
%   params.M
%   params.noise_on

    %% ---------------------------
    %  User-adjustable defaults
    %  ---------------------------
    p = inputParser;

    addParameter(p, 'M0',              8e-6);
    addParameter(p, 'Mvals',           2e-6:0.1e-6:18e-6);
    addParameter(p, 'tf_seg',          6e4);         % ms
    addParameter(p, 'dt',              0.05);        % ms
    addParameter(p, 'Area',            1e4/18);      % gives about 1e4 K channels if NK = round(18*Area)
    addParameter(p, 'noise_on',        false);
    addParameter(p, 'nRuns',           1);
    addParameter(p, 'rng_seed',        1);
    addParameter(p, 'save_results',    false);
    addParameter(p, 'save_name',       'fig8_12D_closedloop_results.mat');
    addParameter(p, 'checkpoint_name', 'fig8_12D_closedloop_checkpoint.mat');
    addParameter(p, 'save_every_M',    true);

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

    run_completed = false(nRuns, nM);

    total_jobs = nRuns * nM;
    job_counter = 0;
    t_global = tic;

    fprintf('\nStarting 12D metabolic-demand experiment\n');
    fprintf('nRuns = %d, nM = %d, total jobs = %d\n', nRuns, nM, total_jobs);
    fprintf('Checkpoint file: %s\n\n', S.checkpoint_name);

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

        fprintf('=============================\n');
        fprintf('Starting run %d of %d\n', r, nRuns);
        fprintf('=============================\n');

        % ------------------------
        % Baseline burn-in at M0
        % ------------------------
        fprintf('  Burn-in 1 at M0 = %.6g ... ', S.M0);
        t_step = tic;
        [~, U0] = simulate_closedloop_Kfull(12e4, S.dt, inits12, S.Area, params);
        fprintf('done (%.2f s)\n', toc(t_step));
        inits1 = U0(end,:).';

        fprintf('  Burn-in 2 at M0 = %.6g ... ', S.M0);
        t_step = tic;
        [~, U1] = simulate_closedloop_Kfull(12e4, S.dt, inits1, S.Area, params);
        fprintf('done (%.2f s)\n', toc(t_step));
        initsM = U1(end,:).';

        % ------------------------
        % Sweep over M values
        % ------------------------
        for ix = 1:nM
            params.M = Mvals(ix);
            job_counter = job_counter + 1;

            fprintf('\nRun %d/%d, M index %d/%d, M = %.6g\n', r, nRuns, ix, nM, params.M);
            fprintf('Overall progress: %d / %d (%.1f%%)\n', ...
                job_counter, total_jobs, 100*job_counter/total_jobs);

            t_job = tic;

            % Replicate the old "chain several long segments" approach
            fprintf('  Segment 2 ... ');
            t_step = tic;
            [~, U2] = simulate_closedloop_Kfull(S.tf_seg, S.dt, initsM,       S.Area, params);
            fprintf('done (%.2f s)\n', toc(t_step));

            fprintf('  Segment 3 ... ');
            t_step = tic;
            [~, U3] = simulate_closedloop_Kfull(S.tf_seg, S.dt, U2(end,:).',  S.Area, params);
            fprintf('done (%.2f s)\n', toc(t_step));

            fprintf('  Segment 4 ... ');
            t_step = tic;
            [~, U4] = simulate_closedloop_Kfull(S.tf_seg, S.dt, U3(end,:).',  S.Area, params);
            fprintf('done (%.2f s)\n', toc(t_step));

            fprintf('  Segment 5 ... ');
            t_step = tic;
            [~, U5] = simulate_closedloop_Kfull(S.tf_seg, S.dt, U4(end,:).',  S.Area, params);
            fprintf('done (%.2f s)\n', toc(t_step));

            fprintf('  Segment 6 ... ');
            t_step = tic;
            [t6, U6] = simulate_closedloop_Kfull(S.tf_seg, S.dt, U5(end,:).', S.Area, params);
            fprintf('done (%.2f s)\n', toc(t_step));

            % Average PO2blood over final segment only
            PO2blood6 = U6(:,12);
            avgPO2blood_all(r, ix) = trapz(t6, PO2blood6) / (t6(end) - t6(1));

            % Save final state for continuation and diagnostics
            u_final_all(:, ix, r) = U6(end,:).';
            run_completed(r, ix)  = true;
            initsM = U6(end,:).';

            fprintf('  avgPO2blood = %.6f\n', avgPO2blood_all(r, ix));
            fprintf('  Finished this M in %.2f s\n', toc(t_job));

            % ------------------------
            % Save checkpoint
            % ------------------------
            if S.save_every_M
                checkpoint = struct();
                checkpoint.Mvals            = Mvals;
                checkpoint.avgPO2blood_all  = avgPO2blood_all;
                checkpoint.u_final_all      = u_final_all;
                checkpoint.run_completed    = run_completed;
                checkpoint.settings         = S;
                checkpoint.current_run      = r;
                checkpoint.current_M_index  = ix;
                checkpoint.elapsed_seconds  = toc(t_global);

                save(S.checkpoint_name, 'checkpoint', '-v7.3');
                fprintf('  Checkpoint saved to %s\n', S.checkpoint_name);
            end
        end

        fprintf('\nCompleted run %d of %d\n\n', r, nRuns);
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

    fprintf('Total elapsed time: %.2f s\n', results.elapsed_seconds);

    %% ---------------------------
    %  Optional final save
    %  ---------------------------
    if S.save_results
        save(S.save_name, 'results', '-v7.3');
        fprintf('Saved final results to %s\n', S.save_name);
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