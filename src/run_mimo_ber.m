function results = run_mimo_ber(Nt, Nr, detector, rho, varargin)
%RUN_MIMO_BER  BER vs. SNR of a spatially correlated Rayleigh MIMO link (BPSK, ZF/MMSE).
%
%   results = RUN_MIMO_BER(Nt, Nr, detector, rho)
%   results = RUN_MIMO_BER(Nt, Nr, detector, rho, Name, Value, ...)
%
%   Simulates BPSK spatial multiplexing over a flat Rayleigh channel with
%   Kronecker/Toeplitz antenna correlation, compares the ideal channel
%   (rho = 0) against the requested correlation factor, and evaluates the
%   average channel condition number for a range of correlation values.
%
%   Inputs
%     Nt, Nr    number of transmit / receive antennas (use Nr >= Nt for ZF)
%     detector  'zf' | 'mmse' | 'both'   (or 1 | 2 | 3)
%     rho       spatial correlation factor, 0 <= rho < 1  (report uses 0.7)
%
%   Name-Value options
%     'SNR_dB'          SNR grid in dB                      (default 0:2:20)
%     'NumBits'         bits simulated per SNR/rho point    (default 1e4)
%     'CondTrials'      channel draws per rho for COND(H)   (default 200)
%     'RhoRange'        rho grid for the condition number   (default 0:0.1:0.9)
%     'Spectrum'        stream Tx/Rx to a spectrumAnalyzer  (default false)
%                       (needs DSP System Toolbox)
%     'SpectrumSNR_dB'  SNR at which the spectrum is shown  (default 10)
%     'Seed'            RNG seed for reproducibility        (default [] = none)
%     'Plot'            draw the figures                    (default true)
%     'SaveDir'         if non-empty, save PNGs and a .mat here (default '')
%
%   Output (struct)
%     Nt, Nr, rho_values, SNR_dB, numBits,
%     BER_ZF, BER_MMSE   [numel(rho_values) x numel(SNR_dB)], NaN if unused
%     rho_range, avg_cond
%
%   Signal model:  y = H*x + n,  x in {-1,+1}^Nt,  n ~ CN(0, noiseVar*I),
%   noiseVar = 10^(-SNR_dB/10)   (i.e. SNR = Es/N0 per receive antenna branch
%   of a single transmit stream, unit-energy symbols, E|h_ij|^2 = 1).
%
%   Example
%     results = run_mimo_ber(2, 2, 'both', 0.7, 'Seed', 42);
%
%   See also MAIN, DETECT_ZF, DETECT_MMSE, CORRELATION_FACTORS, GENERATE_CHANNEL.

    %% ---- parse inputs -------------------------------------------------
    p = inputParser;
    addParameter(p, 'SNR_dB',         0:2:20);
    addParameter(p, 'NumBits',        1e4);
    addParameter(p, 'CondTrials',     200);
    addParameter(p, 'RhoRange',       0:0.1:0.9);
    addParameter(p, 'Spectrum',       false);
    addParameter(p, 'SpectrumSNR_dB', 10);
    addParameter(p, 'Seed',           []);
    addParameter(p, 'Plot',           true);
    addParameter(p, 'SaveDir',        '');
    parse(p, varargin{:});
    o = p.Results;

    if isnumeric(detector)
        names    = {'zf', 'mmse', 'both'};
        detector = names{detector};
    end
    detector = lower(detector);
    useZF    = any(strcmp(detector, {'zf',   'both'}));
    useMMSE  = any(strcmp(detector, {'mmse', 'both'}));
    if ~(useZF || useMMSE)
        error('run_mimo_ber:detector', "detector must be 'zf', 'mmse' or 'both'.");
    end
    if Nr < Nt
        warning('run_mimo_ber:antennas', ...
            'Nr < Nt: the channel is rank deficient, so ZF is not well defined.');
    end
    if ~isempty(o.Seed)
        rng(o.Seed);
    end

    %% ---- set-up -------------------------------------------------------
    rho_values = [0, rho];                    % ideal vs. correlated
    nRho       = numel(rho_values);
    nSNR       = numel(o.SNR_dB);
    numTrials  = floor(o.NumBits / Nt);       % each trial carries Nt bits
    totalBits  = numTrials * Nt;

    BER_ZF   = nan(nRho, nSNR);
    BER_MMSE = nan(nRho, nSNR);

    % Optional spectrum analyser (Tx vs Rx)
    useSpec = false;
    nCols   = 1;                              % samples per stream per trial
    if o.Spectrum
        try
            spec = spectrumAnalyzer('NumInputPorts', 2, 'ChannelNames', {'Tx', 'Rx'});
            [~, specIdx] = min(abs(o.SNR_dB - o.SpectrumSNR_dB));
            useSpec = true;
            nCols   = 20;                     % repeat symbols to get a visible spectrum
        catch
            warning('run_mimo_ber:spectrum', ...
                'spectrumAnalyzer unavailable (needs DSP System Toolbox). Skipping spectrum.');
        end
    end

    %% ---- Monte-Carlo simulation --------------------------------------
    for r = 1:nRho
        C = correlation_factors(Nt, Nr, rho_values(r));

        for s = 1:nSNR
            noiseVar = 10^(-o.SNR_dB(s)/10);
            errZF    = 0;
            errMMSE  = 0;
            showSpec = useSpec && (r == nRho) && (s == specIdx);

            for k = 1:numTrials
                bits = randi([0 1], Nt, 1);
                x    = (2*bits - 1) * ones(1, nCols);          % BPSK: 0->-1, 1->+1
                H    = generate_channel(C);
                n    = sqrt(noiseVar/2) * (randn(Nr, 1) + 1j*randn(Nr, 1));
                y    = H*x + n*ones(1, nCols);

                if showSpec && mod(k, 50) == 0
                    spec(x(:), y(:));
                end

                if useZF
                    xhat  = detect_zf(H, y(:, 1));
                    errZF = errZF + sum(bits ~= (real(xhat) > 0));
                end
                if useMMSE
                    xhat    = detect_mmse(H, y(:, 1), noiseVar);
                    errMMSE = errMMSE + sum(bits ~= (real(xhat) > 0));
                end
            end

            if useZF,   BER_ZF(r, s)   = errZF   / totalBits; end
            if useMMSE, BER_MMSE(r, s) = errMMSE / totalBits; end
        end
    end
    if useSpec
        release(spec);
    end

    %% ---- channel conditioning ----------------------------------------
    avg_cond = avg_condition_number(Nt, Nr, o.RhoRange, o.CondTrials);

    %% ---- pack results -------------------------------------------------
    results = struct( ...
        'Nt', Nt, 'Nr', Nr, 'rho_values', rho_values, ...
        'SNR_dB', o.SNR_dB, 'numBits', totalBits, ...
        'BER_ZF', BER_ZF, 'BER_MMSE', BER_MMSE, ...
        'rho_range', o.RhoRange, 'avg_cond', avg_cond);

    %% ---- plots / export ----------------------------------------------
    if o.Plot || ~isempty(o.SaveDir)
        vis = 'off';
        if o.Plot, vis = 'on'; end
        tag = sprintf('%dx%d', Nt, Nr);
        figs = {};  files = {};

        if useZF
            figs{end+1}  = plot_ber(o.SNR_dB, BER_ZF, rho_values, 'ZF Detector BER Performance', vis);
            files{end+1} = [tag '_zf_ber.png'];
        end
        if useMMSE
            figs{end+1}  = plot_ber(o.SNR_dB, BER_MMSE, rho_values, 'MMSE Detector BER Performance', vis);
            files{end+1} = [tag '_mmse_ber.png'];
        end
        figs{end+1}  = plot_cond(o.RhoRange, avg_cond, vis);
        files{end+1} = [tag '_condition_number.png'];

        if ~isempty(o.SaveDir)
            if ~exist(o.SaveDir, 'dir'), mkdir(o.SaveDir); end
            for i = 1:numel(figs)
                save_figure(figs{i}, fullfile(o.SaveDir, files{i}));
            end
            save(fullfile(o.SaveDir, ['results_' tag '.mat']), 'results');
            fprintf('Saved %d figures and results_%s.mat to %s\n', numel(figs), tag, o.SaveDir);
        end
    end
end

%% ======================= local helpers ==================================
function fig = plot_ber(SNR_dB, BER, rho_values, ttl, vis)
    fig = figure('Visible', vis);
    semilogy(SNR_dB, BER(1, :), 'b-o', 'LineWidth', 1.5); hold on;
    semilogy(SNR_dB, BER(2, :), 'r-s', 'LineWidth', 1.5);
    grid on; xlabel('SNR (dB)'); ylabel('BER'); title(ttl);
    legend(sprintf('\\rho = %g (ideal)', rho_values(1)), ...
           sprintf('\\rho = %g (correlated)', rho_values(2)), 'Location', 'southwest');
end

function fig = plot_cond(rho_range, avg_cond, vis)
    fig = figure('Visible', vis);
    plot(rho_range, avg_cond, 'm-o', 'LineWidth', 1.5);
    grid on; xlabel('Correlation Factor \rho'); ylabel('Average Condition Number');
    title('Condition Number vs Correlation');
end

function save_figure(fig, path)
    try
        exportgraphics(fig, path, 'Resolution', 200);   % R2020a+
    catch
        saveas(fig, path);
    end
end
