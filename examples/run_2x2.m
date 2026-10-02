%% RUN_2X2  Reproduce the 2x2 results of the report (rho = 0.7, ZF + MMSE).
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
outDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'generated');

results = run_mimo_ber(2, 2, 'both', 0.7, ...
    'Seed', 42, 'Spectrum', true, 'SaveDir', outDir);

% Tip: for smoother curves use more bits, e.g. 'NumBits', 1e5
