%==========================================================================
% 4-PAM Receiver with Gray Decoding and LDPC FEC
%==========================================================================
% Description:
%   Processes received 4-PAM signal with synchronization, Volterra
%   nonlinear equalization, amplitude correction, and LDPC decoding.
%
% Processing Chain:
%   Synchronization -> Downsampling -> Volterra Equalization -> 
%   Amplitude Correction -> Deinterleaving -> LLR Calculation -> 
%   LDPC Decoding
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
knowledge_file = 'example/Knowledge_4pam_gray.mat';
rx_file = 'example/RigolDS1.csv';

fprintf('Loading transmitter configuration: %s\n', knowledge_file);
load(knowledge_file);

fprintf('Loading received waveform: %s\n', rx_file);
rxOsc = readmatrix(rx_file);

% Sampling parameters
upsampleRate = 5;  % Fixed 5x oversampling
fprintf('Using %dx oversampling\n', upsampleRate);
fprintf('Gray mapping: 00->-3, 01->-1, 11->+1, 10->+3\n');

% Calculate expected signal length
L = upsampleRate * length(finalSerialSignal);

[cfgEnc, cfgDec] = getProtoMatrix(ldpcBlockLength, ldpcBlockLength * rateLDPC);

% Upsample reference signal to match receiver sampling rate
originalSignalUpsampled = repelem(finalSerialSignal, upsampleRate);

%----------------------------------------------------------------------
% Step 2: Signal Synchronization (Cross-Correlation)
%----------------------------------------------------------------------
rx = rxOsc(:); 
ref = originalSignalUpsampled(:);

% Find optimal synchronization point via cross-correlation
[c,lags] = xcorr(rx, ref);
[~, idx] = max(abs(c));
syncIndex = lags(idx) + 1;

% Extract synchronized segment
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

% Calculate received signal power
power_rx_downsample = mean(signal_downsample.^2);
fprintf('\n=== Received Signal Power ===\n');
fprintf('Downsampled signal power: %.6f\n', power_rx_downsample);

%----------------------------------------------------------------------
% Step 5: Volterra Nonlinear Equalization
%----------------------------------------------------------------------
totalSymbols = ldpcBlockLength / 2 * totalBlocks;
trainLen = floor(totalSymbols * 0.1);
trainStart = floor(totalSymbols * 0.1);

fprintf('\n=== Training Configuration ===\n');
fprintf('Training length: %d (%.1f%%), start position: %d\n', trainLen, trainLen/totalSymbols*100, trainStart+1);

% Uniform distribution statistics
p_effective = [0.25, 0.25, 0.25, 0.25];
power_tx = mean(originalSignal.^2);
fprintf('Transmitted signal power: %.6f\n', power_tx);

% Extract training sequences
txTrain = originalSignal(trainStart+1 : trainStart+trainLen);
rxTrain = signal_downsample(trainStart+1 : trainStart+trainLen);

% Volterra equalizer parameters
Ntap = 15;
delta = floor((Ntap-1)/2);  % Group delay
Nr = numel(rxTrain) - Ntap + 1;

if Nr <= 0
    error('Insufficient training samples for %d-tap Volterra equalizer', Ntap);
end

% Build Volterra feature matrix (linear + quadratic + bias terms)
Phi = zeros(Nr, 2*Ntap + 1);
for i = 1:Nr
    seg = rxTrain(i+Ntap-1:-1:i);
    Phi(i, 1:Ntap) = seg;           % Linear terms
    Phi(i, Ntap+1:2*Ntap) = seg.^2; % Quadratic terms
    Phi(i, end) = 1;                % Bias term
end
d = txTrain( (1:Nr) + delta );

% Ridge regression to solve for coefficients (regularization to prevent ill-conditioning)
lambda = 1e-3;
theta = (Phi' * Phi + lambda * eye(size(Phi,2))) \ (Phi' * d);

fprintf('\n=== Volterra Equalizer ===\n');
fprintf('Number of taps: %d, training samples: %d\n', Ntap, Nr);

% Apply equalizer to full received sequence
Nr_full = numel(signal_downsample) - Ntap + 1;
Phi_full = zeros(Nr_full, 2*Ntap + 1);
for i = 1:Nr_full
    seg = signal_downsample(i+Ntap-1:-1:i);
    Phi_full(i, 1:Ntap) = seg;
    Phi_full(i, Ntap+1:2*Ntap) = seg.^2;
    Phi_full(i, end) = 1;
end
y_hat = Phi_full * theta;

% Compensate for group delay and handle edge effects
signal_before_correction = zeros(size(signal_downsample));
signal_before_correction(delta+1 : delta+Nr_full) = y_hat;
signal_before_correction(1:delta) = y_hat(1);
signal_before_correction(delta+Nr_full+1:end) = y_hat(end);

%----------------------------------------------------------------------
% Step 6: Amplitude and DC Offset Correction
%----------------------------------------------------------------------
txTrain_aligned = txTrain;
rxTrain_eq = signal_before_correction(trainStart+1 : trainStart+trainLen);

% Least squares estimation of linear transformation parameters
X = [rxTrain_eq(:), ones(trainLen, 1)];
theta = X \ txTrain_aligned(:);
a_scale = theta(1);    % Scaling factor
b_offset = theta(2);   % DC offset

fprintf('\n=== Amplitude Correction ===\n');
fprintf('Scaling factor: %.6f, DC offset: %.6f\n', a_scale, b_offset);

% Apply correction to entire sequence
signal_equalized = a_scale * signal_before_correction + b_offset;

% Calculate equalized signal power
power_equalized = mean(signal_equalized.^2);
fprintf('Equalized power: %.6f (relative to TX: %.2f dB)\n', power_equalized, 10*log10(power_equalized/power_tx));

%----------------------------------------------------------------------
% Step 7: Noise Variance Estimation
%----------------------------------------------------------------------
rxTrain_eq_var = signal_equalized(trainStart+1 : trainStart+trainLen);
txTrain_var    = originalSignal(trainStart+1 : trainStart+trainLen);
errorVec       = rxTrain_eq_var(:) - double(txTrain_var(:));
variance       = mean(abs(errorVec).^2);

% Reshape received signal into block matrix structure
rxMatrix = reshape(signal_equalized, ldpcBlockLength/2, totalBlocks);

% Initialize recovered bit matrix
recoveredBitsMat = zeros(k, totalBlocks);

%----------------------------------------------------------------------
% Step 8: Demodulation and Decoding
%----------------------------------------------------------------------
% Block-wise processing: Deinterleaving -> LLR calculation -> LDPC decoding
llr_all = [];
for blk = 1:totalBlocks
    rxBlock = rxMatrix(:, blk);
    rxBlock = ldpc_deinterleave(rxBlock, rateLDPC);
    
    % Calculate log-likelihood ratios (LLR) using Gray code mapping
    llr = llr_pam_gray(rxBlock, p_effective, sqrt(variance), A, ldpcBlockLength, rateLDPC);
    
    % Collect LLR from first block for analysis
    if blk == 1
        llr_all = llr;
    end
    
    % LDPC soft-decision decoding
    decodedBits = ldpcDecode(llr', cfgDec, numIter_pcs, 'OutputFormat', 'info');
    
    % Store decoded bits (no CCDM decoding for uniform scheme)
    recoveredBitsMat(:, blk) = decodedBits;
end

%----------------------------------------------------------------------
% Step 9: Bit Error Rate Statistics
%----------------------------------------------------------------------
blkErr = sum(recoveredBitsMat ~= bitsMatrix, 1);
numErrorBlocks = sum(blkErr > 0);

% Calculate overall coded BER
BER = mean(srcBitStream ~= recoveredBitsMat(:));

%----------------------------------------------------------------------
% Step 10: Uncoded SER Reference (Hard Decision)
%----------------------------------------------------------------------
% Hard decision with Gray code mapping
hardDecision = zeros(size(signal_equalized));
for i = 1:length(signal_equalized)
    [~, idx] = min(abs(signal_equalized(i) - A));
    hardDecision(i) = A(idx);
end

% Calculate uncoded hard-decision symbol error rate
SER_uncoded = mean(originalSignal ~= hardDecision(:));

fprintf('\n=== Performance Statistics ===\n');
fprintf('BER: %g (%.4f%%)\n', BER, BER*100);
fprintf('Pre FEC SER: %g (%.4f%%)\n', SER_uncoded, SER_uncoded*100);
fprintf('Error blocks: %d / %d\n', numErrorBlocks, totalBlocks);
fprintf('\n=== Receiver Processing Complete ===\n');
