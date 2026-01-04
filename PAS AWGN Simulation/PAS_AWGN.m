%==========================================================================
% 4-PAM PAS LDPC-Coded AWGN Channel Simulation
%==========================================================================
% Description:
%   Monte Carlo simulation of 4-PAM modulation with Probabilistic Amplitude
%   Shaping (PAS) and LDPC forward error correction over AWGN channel.
%   Uses CCDM (Constant Composition Distribution Matching) for amplitude
%   shaping with non-uniform amplitude distribution.
%
% Processing Chain:
%   Random Bits -> CCDM Encoding -> LDPC Encoding ->
%   4-PAM Modulation -> AWGN Channel ->
%   PAS-Aware LLR Calculation -> LDPC Decoding -> CCDM Decoding
%
% Performance Metric:
%   Bit Error Rate (BER) vs. Eb/N0
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

clear;
%close all;
clc;

rng(42);  % Fixed seed for reproducibility

%----------------------------------------------------------------------
% Simulation Parameters
%----------------------------------------------------------------------
numTrials = 10000;                              % Monte Carlo trials per SNR point
EbN0Range = 1:0.5:10;                          % Eb/N0 range in dB for BER simulation
[cfgEnc, cfgDec] = getProtoMatrix(1944, 1296); % LDPC code configuration (rate 2/3)
numTasks = 500;                                % Parallel tasks for parfor
numIter = 25;                                  % LDPC decoding iterations
M = 4;                                         % PAM order (4-PAM)
m = 2;                                         % Amplitude bits per symbol (|x| in {1,3})
v = 0.814;                                     % Shaping parameter (0=uniform, 1=max)

%----------------------------------------------------------------------
% PAS Configuration
%----------------------------------------------------------------------
pVector = generatePAS(v);  % Returns [P(|x|=1), P(|x|=3)]
fprintf('=== PAS Configuration ===\n');
fprintf('Shaping parameter v: %.2f\n', v);
fprintf('Amplitude probabilities: [%.4f, %.4f]\n', pVector);

% LDPC parameters
infoLength = cfgEnc.NumInformationBits;
blockLength = cfgDec.BlockLength;
rateLDPC = infoLength / blockLength;

% CCDM setup
n = blockLength / m;                    % Amplitude symbols per block
rateGamma = rateLDPC * m - (m - 1);     % Sign bit rate
numExtraBits = round(rateGamma * n);    % Systematic sign bits

% Initialize CCDM
[p_quant, k, n_i] = ccdm.initialize(pVector, n);
fprintf('CCDM rate k/n: %d/%d = %.6f bits/symbol\n', k, n, k/n);

%----------------------------------------------------------------------
% Rate Analysis
%----------------------------------------------------------------------
entropy = -sum(pVector .* log2(pVector + eps));
ccdmRate = k / n;
signRate = numExtraBits / n;
transmissionRate = ccdmRate + signRate;
channel_per_use = entropy + rateGamma;

fprintf('\n=== Rate Metrics ===\n');
fprintf('Amplitude entropy H(p):   %.6f bits/use\n', entropy);
fprintf('CCDM amplitude rate:      %.6f bits/use\n', ccdmRate);
fprintf('Sign bit rate:            %.6f bits/use\n', signRate);
fprintf('Total transmission rate:  %.6f bits/use\n', transmissionRate);
fprintf('Theoretical rate:         %.6f bits/use\n', channel_per_use);
fprintf('Rate gap:                 %.6f bits/use\n', channel_per_use - transmissionRate);

%----------------------------------------------------------------------
% Power Analysis
%----------------------------------------------------------------------
% 4-PAM constellation: {-3, -1, +1, +3}
A = [-3, -1, 1, 3];
symbolEnergy = mean(abs(A).^2);

% Average power under PAS distribution (uniform signs)
p_amplitudes = [pVector(2)/2, pVector(1)/2, pVector(1)/2, pVector(2)/2];
avgPower_PAS = sum(A.^2 .* p_amplitudes);
fprintf('\n=== Power Metrics ===\n');
fprintf('Uniform 4-PAM power:      %.4f\n', symbolEnergy);
fprintf('Shaped PAS power:         %.4f\n', avgPower_PAS);
fprintf('Shaping gain (power):     %.4f dB\n', 10*log10(symbolEnergy/avgPower_PAS));

trialsPerTask = numTrials / numTasks;
berKnown = zeros(1, length(EbN0Range));

%----------------------------------------------------------------------
% Parameter Sanity Check
%----------------------------------------------------------------------
fprintf('\n=== Debug: Parameter Check ===\n');
fprintf('n (amplitude symbols): %d\n', n);
fprintf('k (source bits): %d\n', k);
fprintf('numExtraBits: %d\n', numExtraBits);
fprintf('infoLength (LDPC info): %d\n', infoLength);
fprintf('blockLength (LDPC total): %d\n', blockLength);
fprintf('Parity length: %d\n', blockLength - infoLength);
fprintf('Expected: n + numExtraBits should equal infoLength\n');
fprintf('Check: %d + %d = %d (infoLength = %d)\n', n, numExtraBits, n+numExtraBits, infoLength);
if n + numExtraBits ~= infoLength
    warning('Configuration mismatch: n + numExtraBits ~= infoLength');
end

%----------------------------------------------------------------------
% Main Monte Carlo Simulation Loop
%----------------------------------------------------------------------
for EbN0Index = 1:length(EbN0Range)
    tic;
    EbN0 = EbN0Range(EbN0Index);
    
    % Calculate noise variance from Eb/N0
    EbN0_linear = 10^(EbN0 / 10);
    Eb = avgPower_PAS / transmissionRate;
    N0 = Eb / EbN0_linear;
    variance = N0 / 2;
    
    taskErrorsKnown = zeros(1, numTasks);
    
    % Parallel Monte Carlo trials
    parfor taskId = 1:numTasks
        localErrorsKnown = 0;
        
        for trial = 1:trialsPerTask
            % Step 1: Generate source bits for CCDM
            srcBits = randi([0, 1], k, 1);
            
            % Step 2: Generate systematic sign bits
            extraBits = randi([0, 1], numExtraBits, 1);
            
            % Step 3: CCDM encoding (source bits -> amplitude bits)
            ccdmBits = ccdm.encode(srcBits', n_i);
            ccdmBits = ccdmBits(:);
            
            % Step 4: LDPC encoding (amplitude bits + sign bits)
            ldpcInput = [ccdmBits; extraBits];
            ldpcParity = ldpcEncode(ldpcInput, cfgEnc, "OutputFormat", "parity");
            
            % Step 5: Map amplitude bits to {1, 3}
            amplitudes = ccdmBits * 2 + 1;
            
            % Step 6: Map sign bits to {-1, +1}
            signBits = [extraBits; ldpcParity];
            signs = signBits * 2 - 1;
            
            % Step 7: Generate 4-PAM symbols
            txSymbols = amplitudes .* signs;
            
            % Step 8: AWGN channel
            noise = sqrt(variance) * randn(size(txSymbols));
            rxSymbols = txSymbols + noise;
            
            % Step 9: PAS-aware LLR calculation
            llr = llr_pas_true(rxSymbols, pVector, variance, n, numExtraBits, blockLength, 40);
            
            % Step 10: LDPC decoding (information bits only)
            decodedInfo = ldpcDecode(llr, cfgDec, numIter, "OutputFormat", "info");
            decodedInfo = decodedInfo(:);
            
            % Step 11: Split amplitude and sign bits
            decodedAmpBits = decodedInfo(1:n);
            decodedExtraBits = decodedInfo(n+1:end);
            
            % Step 12: CCDM decoding (amplitude bits -> source bits)
            recoveredSrcBits = ccdm.decode(double(decodedAmpBits)', n_i, k)';
            
            % Step 13: Count errors
            srcErrors = nnz(double(srcBits) ~= double(recoveredSrcBits));
            extraErrors = nnz(double(extraBits) ~= double(decodedExtraBits));
            localErrorsKnown = localErrorsKnown + srcErrors + extraErrors;
        end
        taskErrorsKnown(taskId) = localErrorsKnown;
    end
    
    totalErrors = sum(taskErrorsKnown);
    totalInfoBits = numTrials * (k + numExtraBits);
    berKnown(1, EbN0Index) = totalErrors / totalInfoBits;
    fprintf('Eb/N0 = %.1f dB, BER = %.2e, elapsed: ', EbN0, berKnown(EbN0Index));
    toc;
end

fprintf('\nPAS simulation complete (v = %.2f)\n', v);

%----------------------------------------------------------------------
% Plot BER Curve
%----------------------------------------------------------------------
figure;
semilogy(EbN0Range, berKnown, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate');
title(sprintf('PAS %d-PAM over AWGN (v=%.2f, Rate %.3f)', M, v, transmissionRate));

%----------------------------------------------------------------------
% Display Simulation Summary
%----------------------------------------------------------------------
fprintf('\n=== Simulation Summary ===\n');
fprintf('CCDM source bits: %d\n', numTrials * k);
fprintf('Extra sign bits: %d\n', numTrials * numExtraBits);
fprintf('Total information bits: %d\n', numTrials * (k + numExtraBits));
fprintf('Modulation: %d-PAM with PAS\n', M);
fprintf('Shaping parameter: v = %.2f\n', v);
fprintf('LDPC rate: %.3f\n', rateLDPC);
fprintf('Transmission rate: %.6f bits/symbol\n', transmissionRate);

%----------------------------------------------------------------------
% Save Results
%----------------------------------------------------------------------
EbN0 = EbN0Range;
BER = berKnown;
save('results_PAS_AWGN.mat', 'EbN0', 'BER');
fprintf('Results saved to: results_PAS_AWGN.mat\n');
