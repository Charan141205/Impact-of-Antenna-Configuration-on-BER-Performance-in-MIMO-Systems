%% MAIN  Interactive entry point (mirrors the workflow in the PBL report).
%
%  Prompts for the antenna configuration, detector and correlation factor,
%  then runs the BER-vs-SNR simulation, the condition-number study and the
%  Tx/Rx spectrum plot.
%
%  Requires MATLAB; the spectrum plot additionally needs DSP System Toolbox.

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));

Nt              = input('Nt (number of transmit antennas): ');
Nr              = input('Nr (number of receive antennas): ');
detector_choice = input('Detector  [1-ZF, 2-MMSE, 3-Both]: ');
rho_user        = input('Correlation factor rho (0-0.9): ');

results = run_mimo_ber(Nt, Nr, detector_choice, rho_user, ...
                       'Spectrum', true, ...
                       'SaveDir',  fullfile(fileparts(mfilename('fullpath')), 'results', 'generated'));
