%==========================================================================
% APAS Receiver with Gray Decoding and CCDM Deshaping
%==========================================================================
% Description:
%   Processes received 4-PAM signal with APAS (Asymmetric Probabilistic
%   Amplitude Shaping). Includes synchronization, Volterra equalization,
%   LDPC decoding, and CCDM deshaping.
%
% Processing Chain:
%   Synchronization -> Downsampling -> Volterra Equalization ->
%   Amplitude Correction -> Deinterleaving -> LLR Calculation ->
%   LDPC Decoding -> Gray Remapping -> CCDM Decoding
%
% Gray Code Mapping:
%   00 -> -3,  01 -> -1,  11 -> +1,  10 -> +3
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
knowledge_file = 'example/Knowledge_gray.mat';
rx_file = 'example/RigolDS3.csv';

fprintf('Loading transmitter configuration: %s\n', knowledge_file);
load(knowledge_file);

fprintf('Loading received waveform: %s\n', rx_file);
rxOsc = readmatrix(rx_file);

% Sampling parameters
upsampleRate = 5;
fprintf('Using %dx oversampling\n', upsampleRate);
fprintf('Gray mapping: 00->-3, 01->-1, 11->+1, 10->+3\n');

L = upsampleRate * length(finalSerialSignal);

[cfgEnc, cfgDec] = getProtoMatrix(ldpcBlockLength, ldpcBlockLength * rateLDPC);

originalSignalUpsampled = repelem(finalSerialSignal, upsampleRate);

%----------------------------------------------------------------------
% Step 2: Signal Synchronization (Cross-Correlation)
%----------------------------------------------------------------------
rx = rxOsc(:);
ref = originalSignalUpsampled(:);

[c,lags] = xcorr(rx, ref);
[~, idx] = max(abs(c));
syncIndex = lags(idx) + 1;

s = max(1, syncIndex);
e = min(numel(rx), syncIndex + L - 1);
syncSignal = rx(s:e);

fprintf('Synchronization complete, sync index: %d\n', syncIndex);

%----------------------------------------------------------------------
% Step 3: Remove Preamble Sequence
%----------------------------------------------------------------------
syncSignal = syncSignal(numZero * upsampleRate + 1:end);

%----------------------------------------------------------------------
% Step 4: Downsampling
%----------------------------------------------------------------------
U  = upsampleRate;
K  = floor(length(syncSignal)/U);
mid = floor(U/2) + 1;
signal_downsample = syncSignal(mid : U : mid + (K-1)*U);

power_rx_downsample = mean(signal_downsample.^2);
fprintf('\n=== Received Signal Power ===\n');
fprintf('Downsampled signal power: %.6f\n', power_rx_downsample);

%----------------------------------------------------------------------
% Step 5: Volterra Nonlinear Equalization
%----------------------------------------------------------------------
if ~exist('rateLDPC', 'var')
    rateLDPC = 2/3;
end

totalSymbols = ldpcBlockLength / 2 * totalBlocks;
trainLen = floor(totalSymbols * 0.1);
trainStart = floor(totalSymbols * 0.1);

fprintf('\n=== Training Configuration ===\n');
fprintf('Training length: %d (%.1f%%), start position: %d\n', trainLen, trainLen/totalSymbols*100, trainStart+1);

% Calculate effective probability distribution
m = 2;
Nsys_sym = (ldpcBlockLength * rateLDPC) / m;
Npar_sym = (ldpcBlockLength / m) - Nsys_sym;
p_parity = [0.25, 0.25, 0.25, 0.25];

total_sym_per_block = Nsys_sym + Npar_sym;
p_effective = (Nsys_sym * pVector + Npar_sym * p_parity) / total_sym_per_block;

power_tx = mean(originalSignal.^2);
fprintf('Transmitted signal power: %.6f\n', power_tx);

txTrain = originalSignal(trainStart+1 : trainStart+trainLen);
rxTrain = signal_downsample(trainStart+1 : trainStart+trainLen);

Ntap = 15;
delta = floor((Ntap-1)/2);
Nr = numel(rxTrain) - Ntap + 1;

if Nr <= 0
    error('Insufficient training samples for %d-tap Volterra equalizer', Ntap);
end

Phi = zeros(Nr, 2*Ntap + 1);
for i = 1:Nr
    seg = rxTrain(i+Ntap-1:-1:i);
    Phi(i, 1:Ntap) = seg;
    Phi(i, Ntap+1:2*Ntap) = seg.^2;
    Phi(i, end) = 1;
end
d = txTrain( (1:Nr) + delta );

lambda = 1e-3;
theta = (Phi' * Phi + lambda * eye(size(Phi,2))) \ (Phi' * d);

fprintf('\n=== Volterra Equalizer ===\n');
fprintf('Number of taps: %d, training samples: %d\n', Ntap, Nr);

Nr_full = numel(signal_downsample) - Ntap + 1;
Phi_full = zeros(Nr_full, 2*Ntap + 1);
for i = 1:Nr_full
    seg = signal_downsample(i+Ntap-1:-1:i);
    Phi_full(i, 1:Ntap) = seg;
    Phi_full(i, Ntap+1:2*Ntap) = seg.^2;
    Phi_full(i, end) = 1;
end
y_hat = Phi_full * theta;

signal_before_correction = zeros(size(signal_downsample));
signal_before_correction(delta+1 : delta+Nr_full) = y_hat;
signal_before_correction(1:delta) = y_hat(1);
signal_before_correction(delta+Nr_full+1:end) = y_hat(end);

%----------------------------------------------------------------------
% Step 6: Amplitude and DC Offset Correction
%----------------------------------------------------------------------
txTrain_aligned = txTrain;
rxTrain_eq = signal_before_correction(trainStart+1 : trainStart+trainLen);

X = [rxTrain_eq(:), ones(trainLen, 1)];
theta = X \ txTrain_aligned(:);
a_scale = theta(1);
b_offset = theta(2);

fprintf('\n=== Amplitude Correction ===\n');
fprintf('Scaling factor: %.6f, DC offset: %.6f\n', a_scale, b_offset);

signal_equalized = a_scale * signal_before_correction + b_offset;

power_equalized = mean(signal_equalized.^2);
fprintf('Equalized power: %.6f (relative to TX: %.2f dB)\n', power_equalized, 10*log10(power_equalized/power_tx));

%----------------------------------------------------------------------
% Step 7: Noise Variance Estimation
%----------------------------------------------------------------------
rxTrain_eq_var = signal_equalized(trainStart+1 : trainStart+trainLen);
txTrain_var    = originalSignal(trainStart+1 : trainStart+trainLen);

errorVec_global = rxTrain_eq_var(:) - double(txTrain_var(:));
variance_global = mean(abs(errorVec_global).^2);
fprintf('\nGlobal noise variance: %.6f, std: %.4f\n', variance_global, sqrt(variance_global));

rxMatrix = reshape(signal_equalized, ldpcBlockLength/2, totalBlocks);

recoveredBitsMat = zeros(k, totalBlocks);

%----------------------------------------------------------------------
% Step 8: Demodulation and Decoding
%----------------------------------------------------------------------
fprintf('\n=== Starting Demodulation and Decoding ===\n');
fprintf('Using global variance: %.6f\n', variance_global);

for blk = 1:totalBlocks
    rxBlock = rxMatrix(:, blk);
    rxBlock = ldpc_deinterleave(rxBlock, rateLDPC);
    
    llr = llr_pam_gray(rxBlock, p_effective, sqrt(variance_global), A, ldpcBlockLength, rateLDPC);
    
    decodedBits = ldpcDecode(llr', cfgDec, numIter_pcs, 'OutputFormat', 'info');
    
    remapped = bits_remaping_gray(decodedBits);
    
    recoveredBitsMat(:, blk) = ccdm.decode(remapped, n_i, k)';
end

%----------------------------------------------------------------------
% Step 9: Bit Error Rate Statistics
%----------------------------------------------------------------------
blkErr = sum(recoveredBitsMat ~= bitsMatrix, 1);
numErrorBlocks = sum(blkErr > 0);

BER = mean(srcBitStream ~= recoveredBitsMat(:));

%----------------------------------------------------------------------
% Step 10: Uncoded SER Reference (Hard Decision)
%----------------------------------------------------------------------
hardDecision = zeros(size(signal_equalized));
for i = 1:length(signal_equalized)
    [~, idx] = min(abs(signal_equalized(i) - A));
    hardDecision(i) = A(idx);
end

SER_uncoded = mean(originalSignal ~= hardDecision(:));

fprintf('\n=== Performance Statistics ===\n');
fprintf('BER: %g (%.4f%%)\n', BER, BER*100);
fprintf('Pre FEC SER: %g (%.4f%%)\n', SER_uncoded, SER_uncoded*100);
fprintf('Error blocks: %d / %d\n', numErrorBlocks, totalBlocks);
fprintf('\n=== Receiver Processing Complete ===\n');
