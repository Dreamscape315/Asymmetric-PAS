%==========================================================================
% 4-PAM Transmitter with Gray Coding and LDPC FEC
%==========================================================================
% Description:
%   Generates 4-PAM modulated signal with uniform amplitude distribution,
%   Gray bit mapping, and LDPC forward error correction (rate 1/2).
%
% Gray Code Mapping:
%   00 -> -3,  01 -> -1,  11 -> +1,  10 -> +3
%   (Adjacent constellation points differ by only 1 bit)
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

clc;
clear;
rng(42);

%----------------------------------------------------------------------
% Configuration Parameters
%----------------------------------------------------------------------
numZero = 3000;                 % Preamble length for synchronization
numIter_pcs = 25;               % LDPC decoding iterations
rateLDPC = 1/2;                 % LDPC code rate
ldpcBlockLength = 1944;         % LDPC codeword length
modulationOrder = 4;            % 4-PAM modulation
A = [-3 -1 1 3];                % Constellation points
totalBlocks = 150;              % Number of blocks to transmit

% Initialize LDPC encoder/decoder configuration
[cfgEnc, cfgDec] = getProtoMatrix(ldpcBlockLength, ldpcBlockLength * rateLDPC);

% Calculate number of information bits per block
k = ldpcBlockLength * rateLDPC;

%----------------------------------------------------------------------
% Generate Random Bit Stream
%----------------------------------------------------------------------
srcBitStream = randi([0 1], totalBlocks * k, 1);
bitsMatrix = reshape(srcBitStream, k, totalBlocks);
txAllSig = zeros(ldpcBlockLength/2, totalBlocks);

%----------------------------------------------------------------------
% Block-wise Processing: LDPC Encoding -> Modulation
%----------------------------------------------------------------------
for blk = 1:totalBlocks
    blockBits = bitsMatrix(:, blk);
    
    % Step 1: LDPC encoding (add redundancy for error protection)
    codedBits = ldpcEncode(blockBits, cfgEnc, 'OutputFormat', 'whole');
    
    % Step 2: Gray code remapping and modulation to PAM symbols
    txIdx = double(bits_remaping_gray(codedBits));
    txSig = A(txIdx + 1);
    txSig = ldpc_interleave(txSig, rateLDPC);
    txAllSig(:,blk) = txSig.';
end

% Flatten signal to column vector
originalSignal = txAllSig(:);

%----------------------------------------------------------------------
% Add Preamble Sequence (for synchronization)
%----------------------------------------------------------------------
mu_preamble = mean(originalSignal);
preambleSignal = mu_preamble * ones(numZero, 1);
finalSerialSignal = [preambleSignal; txAllSig(:)];

fprintf('\n=== Preamble Information ===\n');
fprintf('Preamble length: %d samples\n', numZero);
fprintf('Preamble level: %.4f (signal mean)\n', mu_preamble);
fprintf('Preamble energy: %.4f\n', mu_preamble^2);

%----------------------------------------------------------------------
% Save Waveform and Configuration
%----------------------------------------------------------------------
% Create output directory if not exists
outputDir = 'example';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

waveformFile = fullfile(outputDir, '4pam_tx_waveform_gray.csv');
knowledgeFile = fullfile(outputDir, 'Knowledge_4pam_gray.mat');

writematrix(finalSerialSignal, waveformFile);

% Save all parameters for receiver
pVector = [0.25, 0.25, 0.25, 0.25];  % Uniform distribution
save(knowledgeFile, 'finalSerialSignal', 'numZero', 'totalBlocks', ...
    'modulationOrder', 'srcBitStream', 'txAllSig', 'originalSignal', ...
    'k', 'ldpcBlockLength', 'rateLDPC', 'pVector', 'A', 'cfgDec','numIter_pcs',...
    'bitsMatrix', 'mu_preamble');

fprintf('\n=== File Save Status ===\n');
fprintf('Waveform saved to: %s\n', waveformFile);
fprintf('Knowledge saved to: %s\n', knowledgeFile);
fprintf('Using Gray mapping: 00->-3, 01->-1, 11->+1, 10->+3\n');

%----------------------------------------------------------------------
% Statistical Analysis
%----------------------------------------------------------------------
% Count symbol occurrences
counts = zeros(1, 4);
for i = 1:4
    counts(i) = sum(originalSignal == A(i));
end
p_actual = counts / length(originalSignal);

% Calculate entropy
H_actual = -sum(p_actual(p_actual>0) .* log2(p_actual(p_actual>0)));

% Calculate spectral efficiency (excluding preamble)
num_data_symbols = length(originalSignal);
total_info_bits = totalBlocks * k;
symbols_per_block = ldpcBlockLength / 2;
spectral_efficiency = k / symbols_per_block;  % bits/symbol

%----------------------------------------------------------------------
% Display Statistics
%----------------------------------------------------------------------
fprintf('\n=== Transmitter Statistics ===\n');
fprintf('Actual probability: [%.4f, %.4f, %.4f, %.4f]\n', p_actual);
fprintf('Entropy: %.6f bits/symbol\n', H_actual);
fprintf('Spectral efficiency: %.6f bits/symbol (4-PAM + LDPC %.0f%%)\n', spectral_efficiency, rateLDPC*100);
fprintf('Symbol probabilities: [%s]\n', sprintf('%.4f ', p_actual));
fprintf('\n=== Transmitter Configuration Complete ===\n');
