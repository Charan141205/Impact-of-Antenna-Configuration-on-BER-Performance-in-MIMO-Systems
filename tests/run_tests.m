function run_tests()
%RUN_TESTS  Lightweight self-checks (no test framework needed; MATLAB or Octave).
%
%   Run from the repository root:   >> run_tests   (or  >> addpath tests; run_tests)

    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
    seed_rng(1);
    n = 0;  failed = 0;

    % 1. rho = 0 -> identity correlation
    C = correlation_factors(3, 3, 0);
    [n, failed] = check(norm(C.sqrtRt - eye(3)) < 1e-12 && norm(C.sqrtRr - eye(3)) < 1e-12, ...
        'rho = 0 gives identity correlation', n, failed);

    % 2. square root is a true square root of the Toeplitz matrix
    C = correlation_factors(3, 4, 0.9);
    [n, failed] = check(norm(C.sqrtRt*C.sqrtRt - toeplitz(0.9.^(0:2))) < 1e-10 && ...
                        norm(C.sqrtRr*C.sqrtRr - toeplitz(0.9.^(0:3))) < 1e-10, ...
        'sqrtm factors reproduce Toeplitz matrices', n, failed);

    % 3. channel dimensions
    H = generate_channel(correlation_factors(2, 3, 0.5));
    [n, failed] = check(isequal(size(H), [3 2]), 'channel is Nr-by-Nt', n, failed);

    % 4. noiseless ZF recovers the transmitted symbols exactly
    C = correlation_factors(2, 2, 0.3);
    H = generate_channel(C);  x = [1; -1];
    [n, failed] = check(norm(detect_zf(H, H*x) - x) < 1e-6, 'ZF is exact without noise', n, failed);

    % 5. MMSE -> ZF as noise variance -> 0
    [n, failed] = check(norm(detect_mmse(H, H*x, 0) - detect_zf(H, H*x)) < 1e-6, ...
        'MMSE equals ZF when noiseVar = 0', n, failed);

    % 6. MMSE shrinks the estimate when noise is large (regularisation)
    y = H*x;
    [n, failed] = check(norm(detect_mmse(H, y, 10)) < norm(detect_zf(H, y)), ...
        'MMSE regularises (smaller norm) at high noise', n, failed);

    % 7. condition number grows with correlation
    ac = avg_condition_number(2, 2, [0 0.9], 300);
    [n, failed] = check(ac(2) > ac(1), 'cond(H) increases with rho', n, failed);

    % 8. end-to-end smoke test
    res = run_mimo_ber(2, 2, 'both', 0.7, 'SNR_dB', [0 20], 'NumBits', 4000, ...
                       'CondTrials', 20, 'Plot', false);
    [n, failed] = check(all(res.BER_ZF(:, 2)   < res.BER_ZF(:, 1)) && ...
                        all(res.BER_MMSE(:, 2) < res.BER_MMSE(:, 1)), ...
        'BER falls as SNR rises', n, failed);
    [n, failed] = check(res.BER_ZF(2, 2) > res.BER_ZF(1, 2), ...
        'correlation hurts ZF BER at 20 dB', n, failed);
    [n, failed] = check(res.BER_MMSE(1, 1) <= res.BER_ZF(1, 1) + 0.02, ...
        'MMSE no worse than ZF at 0 dB (rho = 0)', n, failed);

    fprintf('\n%d/%d checks passed.\n', n - failed, n);
    if failed > 0
        error('run_tests:failed', '%d check(s) failed.', failed);
    end
end

function [n, failed] = check(cond, name, n, failed)
    n = n + 1;
    if cond
        fprintf('  PASS  %s\n', name);
    else
        fprintf('  FAIL  %s\n', name);
        failed = failed + 1;
    end
end

function seed_rng(s)
    try
        rng(s);
    catch
        randn('state', s);  rand('state', s); %#ok<RAND>
    end
end
