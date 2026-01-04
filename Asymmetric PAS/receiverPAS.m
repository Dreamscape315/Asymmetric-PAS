%==========================================================================
% True PAS Receiver with Sign-Magnitude Representation
%==========================================================================
% Description:
%   Processes received 4-PAM signal with true PAS (sign-magnitude). 
%   Includes synchronization, Volterra equalization, separate amplitude
%   and sign LLR calculation, LDPC decoding, and CCDM deshaping.
%
% Processing Chain:
%   Synchronization -> Downsampling -> Volterra Equalization ->
%   Amplitude Correction -> LLR Calculation (amp+sign) ->
%   LDPC Decoding -> CCDM Decoding
%
% Representation:
%   Amplitude: {1, 3} (2 levels, shaped via CCDM)
%   Sign: {-1, +1} (uniform via LDPC parity)
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

clc;
clear;
close all;

%----------------------------------------------------------------------
% Step 1: Load Configuration and Data
%----------------------------------------------------------------------
knowledge_file = 'example/Knowledge_truePAS.mat';
rx_file = 'example/RigolDS2.csv';

fprintf('Loading transmitter configuration: %s\n', knowledge_file);
load(knowledge_file, 'finalSerialSignal', 'constellationSymbolsAll', ...
    'originalSrcBitsAll_I', 'originalExtraBitsAll_I', 'p_quant_I', ...
    'pVector', 'k', 'n', 'n_i', 'numExtraBits', 'cfgDec', ...
    'rateLDPC', 'ldpcBlockLength', 'totalBlocks', 'numIter_pcs', 'numZero');

if ~exist('A', 'var')
    A = [-3, -1, 1, 3];
end

originalSignal = constellationSymbolsAll(:);
bitsMatrix = reshape(originalSrcBitsAll_I, k, totalBlocks);
extraBitsMatrix = reshape(originalExtraBitsAll_I, numExtraBits, totalBlocks);
srcBitStream = originalSrcBitsAll_I(:);
pAmp = double(p_quant_I(:).');
pAmp = pAmp / sum(pAmp);

if exist(rx_file, 'file') == 2
    fprintf('Loading received waveform: %s\n', rx_file);
    rxOsc = readmatrix(rx_file);
else
    warning('Received waveform not found, using ideal TX signal');
    rxOsc = finalSerialSignal;
end

rxOsc = rxOsc(:);
txSerial = finalSerialSignal(:);

% Sampling parameters
upsampleRate = 5;
fprintf('Using %dx oversampling\n', upsampleRate);

L = upsampleRate * length(txSerial);
originalSignalUpsampled = repelem(txSerial, upsampleRate);

[cfgEnc_tmp, cfgDec_tmp] = getProtoMatrix(ldpcBlockLength, ldpcBlockLength * rateLDPC);
if ~exist('cfgDec', 'var') || isempty(cfgDec)
    cfgDec = cfgDec_tmp;
end

%----------------------------------------------------------------------
% Step 2: Signal Synchronization (Cross-Correlation)
%----------------------------------------------------------------------
rx = rxOsc(:);
ref = originalSignalUpsampled(:);
[c, lags] = xcorr(rx, ref);
[~, idx_max] = max(abs(c));
syncIndex = lags(idx_max) + 1;

s = max(1, syncIndex);
e = min(numel(rx), syncIndex + L - 1);
syncSignal = rx(s:e);

fprintf('Synchronization complete, sync index: %d\n', syncIndex);

%----------------------------------------------------------------------
% Step 3: Remove Preamble Sequence
%----------------------------------------------------------------------
syncSignal = syncSignal(numZero * upsampleRate + 1:end);
fprintf('Preamble removed: %d samples\n', numZero * upsampleRate);

%----------------------------------------------------------------------
% Step 4: Downsampling
%----------------------------------------------------------------------
U = upsampleRate;
K = floor(length(syncSignal) / U);
mid = floor(U / 2) + 1;
signal_downsample = syncSignal(mid : U : mid + (K - 1) * U);

power_rx_downsample = mean(signal_downsample.^2);
fprintf('\n=== Received Signal Power ===\n');
fprintf('Downsampled signal power: %.6f\n', power_rx_downsample);

%----------------------------------------------------------------------
% Step 5: Volterra Nonlinear Equalization
%----------------------------------------------------------------------
m = 2;
if ~exist('rateLDPC', 'var')
    rateLDPC = 2/3;
end
totalSymbols = ldpcBlockLength / 2 * totalBlocks;
trainLen = floor(totalSymbols * 0.1);
trainStart = floor(totalSymbols * 0.1);

fprintf('\n=== Training Configuration ===\n');
fprintf('Total symbols: %d\n', totalSymbols);
fprintf('Training length: %d (%.1f%%)\n', trainLen, trainLen / totalSymbols * 100);
fprintf('Training start: %d\n', trainStart + 1);
fprintf('Training range: [%d, %d]\n', trainStart + 1, trainStart + trainLen);

txTrain = originalSignal(trainStart + 1 : trainStart + trainLen);
rxTrain = signal_downsample(trainStart + 1 : trainStart + trainLen);

Nsys_sym = (ldpcBlockLength * rateLDPC) / m;
Npar_sym = (ldpcBlockLength / m) - Nsys_sym;
p_parity = [0.25, 0.25, 0.25, 0.25];
p_data_levels = [pAmp(2)/2, pAmp(1)/2, pAmp(1)/2, pAmp(2)/2];
p_effective = (Nsys_sym * p_data_levels + Npar_sym * p_parity) / (Nsys_sym + Npar_sym);

power_tx = mean(originalSignal.^2);
fprintf('Transmitted signal power: %.6f\n', power_tx);

Ntap = 15;
delta = floor((Ntap - 1) / 2);
Nr = numel(rxTrain) - Ntap + 1;
if Nr <= 0
    error('Insufficient training samples for %d-tap Volterra filter', Ntap);
end

Phi = zeros(Nr, 2 * Ntap + 1);
for i = 1:Nr
    seg = rxTrain(i + Ntap - 1 : -1 : i);
    Phi(i, 1:Ntap) = seg;
    Phi(i, Ntap + 1:2 * Ntap) = seg.^2;
    Phi(i, end) = 1;
end
d = txTrain((1:Nr) + delta);

lambda = 1e-3;
theta = (Phi' * Phi + lambda * eye(size(Phi, 2))) \ (Phi' * d);

fprintf('\n=== Volterra Filter Parameters ===\n');
fprintf('Number of taps: %d\n', Ntap);
fprintf('Delay: %d\n', delta);
fprintf('Training samples: %d\n', Nr);
fprintf('Regularization lambda: %.2e\n', lambda);

Nr_full = numel(signal_downsample) - Ntap + 1;
Phi_full = zeros(Nr_full, 2 * Ntap + 1);
for i = 1:Nr_full
    seg = signal_downsample(i + Ntap - 1 : -1 : i);
    Phi_full(i, 1:Ntap) = seg;
    Phi_full(i, Ntap + 1:2 * Ntap) = seg.^2;
    Phi_full(i, end) = 1;
end
y_hat = Phi_full * theta;

signal_before_correction = zeros(size(signal_downsample));
signal_before_correction(delta + 1 : delta + Nr_full) = y_hat;
signal_before_correction(1:delta) = y_hat(1);
signal_before_correction(delta + Nr_full + 1:end) = y_hat(end);

%----------------------------------------------------------------------
% Step 6: Amplitude and DC Offset Correction
%----------------------------------------------------------------------
rxTrain_eq = signal_before_correction(trainStart + 1 : trainStart + trainLen);
X = [rxTrain_eq(:), ones(trainLen, 1)];
theta_lin = X \ txTrain(:);
a_scale = theta_lin(1);
b_offset = theta_lin(2);

fprintf('\n=== Amplitude Correction Parameters ===\n');
fprintf('Scaling factor a: %.6f\n', a_scale);
fprintf('DC offset b: %.6f\n', b_offset);

signal_equalized = a_scale * signal_before_correction + b_offset;

power_equalized = mean(signal_equalized.^2);
fprintf('Equalized signal power: %.6f\n', power_equalized);
fprintf('Power ratio (equalized/TX): %.4f (%.2f dB)\n', power_equalized/power_tx, 10*log10(power_equalized/power_tx));

%----------------------------------------------------------------------
% Step 7: Noise Variance Estimation
%----------------------------------------------------------------------
rxTrain_eq_var = signal_equalized(trainStart + 1 : trainStart + trainLen);
txTrain_var = originalSignal(trainStart + 1 : trainStart + trainLen);

errorVec_global = rxTrain_eq_var(:) - double(txTrain_var(:));
variance_global = mean(abs(errorVec_global).^2);
fprintf('\n=== Noise Statistics ===\n');
fprintf('Global variance: %.6f, std: %.4f\n', variance_global, sqrt(variance_global));

snrLinear = mean(double(txTrain_var).^2) / max(variance_global, 1e-12);
snrDb = 10 * log10(snrLinear);
fprintf('Estimated SNR: %.2f dB\n', snrDb);

rxMatrix = reshape(signal_equalized, ldpcBlockLength / 2, totalBlocks);
recoveredBitsMat = zeros(k, totalBlocks);
recoveredExtraBitsMat = zeros(numExtraBits, totalBlocks);
LLR_CLIP = 40;

%----------------------------------------------------------------------
% Step 8: PAS LLR Calculation and LDPC Decoding
%----------------------------------------------------------------------
fprintf('\n=== PAS LLR + LDPC Decoding ===\n');
fprintf('Amplitude prior: [%.6f, %.6f]\n', pAmp);
fprintf('Using global noise variance\n');
fprintf('Decoding noise variance: %.6e\n', variance_global);

for blk = 1:totalBlocks
    rxBlock = rxMatrix(:, blk);
    
    llr = llr_pas_true(rxBlock, pAmp, variance_global, n, numExtraBits, ldpcBlockLength, LLR_CLIP);
    
    decodedInfo = ldpcDecode(llr, cfgDec, numIter_pcs, 'OutputFormat', 'info');
    decodedInfo = decodedInfo(:);

    ampBits = decodedInfo(1:n);
    extraBits = decodedInfo(n + (1:numExtraBits));

    recoveredBits = ccdm.decode(ampBits', n_i, k)';
    recoveredBitsMat(:, blk) = recoveredBits;
    recoveredExtraBitsMat(:, blk) = extraBits(:);
end

%----------------------------------------------------------------------
% Step 9: Performance Evaluation
%----------------------------------------------------------------------
blkErr_amp = sum(recoveredBitsMat ~= bitsMatrix, 1);
blkErr_sign = sum(recoveredExtraBitsMat ~= extraBitsMatrix, 1);
numErrorBlocks_amp = sum(blkErr_amp > 0);
numErrorBlocks_sign = sum(blkErr_sign > 0);

BER_amp = mean(originalSrcBitsAll_I ~= recoveredBitsMat(:));
BER_sign = mean(originalExtraBitsAll_I ~= recoveredExtraBitsMat(:));

demodSymbols = pamdemod(signal_equalized, 4, 0, 'bin');
remodSymbols = pammod(demodSymbols, 4, 0, 'bin');
SER_uncoded = mean(originalSignal ~= remodSymbols(:));

fprintf('\n=== Performance Statistics ===\n');
fprintf('Source bit BER: %g (%.4f%%)\n', BER_amp, BER_amp * 100);
fprintf('Sign bit BER: %g (%.4f%%)\n', BER_sign, BER_sign * 100);
fprintf('Pre FEC SER: %g (%.4f%%)\n', SER_uncoded, SER_uncoded * 100);
fprintf('Error blocks (amp): %d / %d\n', numErrorBlocks_amp, totalBlocks);
fprintf('Error blocks (sign): %d / %d\n', numErrorBlocks_sign, totalBlocks);
fprintf('\n=== Receiver Processing Complete ===\n');
