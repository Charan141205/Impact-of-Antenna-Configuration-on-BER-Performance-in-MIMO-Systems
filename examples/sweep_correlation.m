%% SWEEP_CORRELATION  BER at a fixed SNR for every correlation factor (2x2 and 3x3).
%  Not part of the original report - a quick extension showing how BER
%  grows with rho for both detectors.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

rhos   = 0:0.1:0.9;
snr    = 10;                                   % dB
configs = [2 2; 3 3];

figure; hold on; grid on;
for c = 1:size(configs, 1)
    Nt = configs(c, 1);  Nr = configs(c, 2);
    berZF = zeros(size(rhos));  berMMSE = zeros(size(rhos));
    for i = 1:numel(rhos)
        res = run_mimo_ber(Nt, Nr, 'both', rhos(i), ...
            'SNR_dB', snr, 'NumBits', 2e4, 'CondTrials', 1, 'Plot', false, 'Seed', i);
        berZF(i)   = res.BER_ZF(2, 1);
        berMMSE(i) = res.BER_MMSE(2, 1);
    end
    semilogy(rhos, berZF,   '-o', 'DisplayName', sprintf('ZF %dx%d',   Nt, Nr));
    semilogy(rhos, berMMSE, '-s', 'DisplayName', sprintf('MMSE %dx%d', Nt, Nr));
end
set(gca, 'YScale', 'log');
xlabel('Correlation factor \rho'); ylabel(sprintf('BER @ %d dB SNR', snr));
title('BER vs spatial correlation'); legend('Location', 'northwest');
