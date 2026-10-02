%% RUN_3X3  Reproduce the 3x3 results of the report (rho = 0.7, ZF + MMSE).
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
outDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'generated');

results = run_mimo_ber(3, 3, 'both', 0.7, ...
    'Seed', 42, 'Spectrum', true, 'SaveDir', outDir);
