%==========================================================================
% 4-PAM LDPC-Coded AWGN Channel Simulation (Uniform Signaling)
%==========================================================================
% Description:
%   Monte Carlo simulation of 4-PAM modulation with LDPC forward error
%   correction over AWGN channel. Uses uniform symbol distribution and
%   Gray coding for bit mapping.
%
% Processing Chain:
%   Random Bits -> LDPC Encoding -> Gray PAM Modulation ->
%   AWGN Channel -> Soft Demodulation -> LDPC Decoding
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
[cfgEnc, cfgDec] = getProtoMatrix(1944, 972);  % LDPC code configuration (rate 1/2)
numTasks = 100;                                % Parallel tasks for parfor
numIter = 25;                                  % LDPC decoding iterations
M = 4;                                         % PAM order (4-PAM)
bitsPerSymbol = log2(M);                       % Bits per symbol (2 for 4-PAM)

%----------------------------------------------------------------------
% Derived Parameters
%----------------------------------------------------------------------
infoLength = cfgEnc.NumInformationBits;
blockLength = cfgDec.BlockLength;
numBits = numTrials * infoLength;
rate = infoLength / blockLength;
trialsPerTask = numTrials / numTasks;
berKnown = zeros(1, length(EbN0Range));

% 4-PAM constellation and average energy
constellation = pammod(0:M-1, M, 0, 'gray');
symbolEnergy = mean(abs(constellation).^2);
constellation = constellation(:);

% Natural bit weights (MSB -> LSB)
bitWeights = 2.^(bitsPerSymbol - 1:-1:0).';

%----------------------------------------------------------------------
% Main Monte Carlo Simulation Loop
%----------------------------------------------------------------------
for EbN0Index = 1:length(EbN0Range)
    tic;
    EbN0 = EbN0Range(EbN0Index);
    
    % Calculate noise variance from Eb/N0
    EbN0_linear = 10^(EbN0 / 10);
    EsN0_linear = bitsPerSymbol * rate * EbN0_linear;
    variance = symbolEnergy / (2 * EsN0_linear);
    sqrtVariance = sqrt(variance);
    
    taskErrorsKnown = zeros(1, numTasks);
    
    % Parallel Monte Carlo trials
    parfor taskId = 1:numTasks
        localErrorsKnown = 0;
        
        for trial = 1:trialsPerTask
            % Step 1: Generate random information bits
            infoBits = randi([0, 1], infoLength, 1);
            
            % Step 2: LDPC encoding
            codeword = ldpcEncode(infoBits, cfgEnc, "OutputFormat", "whole");
            
            % Step 3: Group bits into PAM symbols (natural labeling)
            bitMatrix = reshape(codeword, bitsPerSymbol, []).';
            symbolIndices = bitMatrix * bitWeights;
            
            % Step 4: 4-PAM modulation (Gray mapping)
            codewordMod = pammod(symbolIndices, M, 0, 'gray');
            
            % Step 5: AWGN channel (real-valued)
            noise = sqrtVariance * randn(size(codewordMod));
            receive = codewordMod + noise;
            
            % Step 6: Soft-demodulation (LLR calculation)
            llr = pamLLR(receive, M, variance, constellation);
            
            % Step 7: LDPC decoding (information bits only)
            decodedBits = ldpcDecode(llr, cfgDec, numIter, "OutputFormat", "info");
            
            % Step 8: Count bit errors
            localErrorsKnown = localErrorsKnown + nnz(infoBits ~= decodedBits);
        end
        taskErrorsKnown(taskId) = localErrorsKnown;
    end
    
    totalErrors = sum(taskErrorsKnown);
    berKnown(1, EbN0Index) = totalErrors / numBits;
    fprintf('Eb/N0 = %.1f dB, BER = %.2e, elapsed: ', EbN0, berKnown(EbN0Index));
    toc;
end

fprintf('LDPC simulation complete. Rate: %.3f\n', rate);

%----------------------------------------------------------------------
% Plot BER Curve
%----------------------------------------------------------------------
figure;
semilogy(EbN0Range, berKnown, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate');
title(sprintf('LDPC-coded %d-PAM over AWGN (Rate %.3f, Gray Mapping)', M, rate));

%----------------------------------------------------------------------
% Display Simulation Summary
%----------------------------------------------------------------------
fprintf('\n=== Simulation Summary ===\n');
fprintf('Total information bits: %d\n', numBits);
fprintf('Modulation: %d-PAM (Gray)\n', M);
fprintf('LDPC rate: %.3f\n', rate);

%----------------------------------------------------------------------
% Save Results
%----------------------------------------------------------------------
EbN0 = EbN0Range;
BER = berKnown;
save('results_PAM_LDPC_AWGN.mat', 'EbN0', 'BER');
fprintf('Results saved to: results_PAM_LDPC_AWGN.mat\n');

