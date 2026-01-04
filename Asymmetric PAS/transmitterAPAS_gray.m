%==========================================================================
% APAS Transmitter with Gray Coding and CCDM Shaping
%==========================================================================
% Description:
%   Generates 4-PAM modulated signal with CCDM (Constant Composition
%   Distribution Matching) shaping, Gray bit mapping, and LDPC FEC.
%   Achieves non-uniform amplitude distribution for power efficiency.
%
% Processing Chain:
%   CCDM Shaping -> Gray Bit Mapping -> LDPC Encoding -> 
%   Gray Remapping -> PAM Modulation -> Interleaving
%
% Gray Code Mapping:
%   00 -> -3,  01 -> -1,  11 -> +1,  10 -> +3
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
rateLDPC = 3/4;                 % LDPC code rate
ldpcBlockLength = 1944;         % LDPC codeword length
modulationOrder = 4;            % 4-PAM modulation
v = 0.0001;                     % Shaping parameter (0=uniform, 1=max shaping)

pVector = generatePAS(modulationOrder, v);
A = [-3 -1 1 3];                % Constellation points
totalBlocks = 150;              % Number of blocks to transmit

% Initialize LDPC encoder/decoder configuration
[cfgEnc, cfgDec] = getProtoMatrix(ldpcBlockLength, ldpcBlockLength * rateLDPC);

% Initialize CCDM (Constant Composition Distribution Matching) parameters
[p_quant_I, k, n_i] = ccdm.initialize(pVector, ldpcBlockLength * rateLDPC/2);

%----------------------------------------------------------------------
% Generate Random Bit Stream
%----------------------------------------------------------------------
srcBitStream = randi([0 1], totalBlocks * k, 1);
bitsMatrix = reshape(srcBitStream, k, totalBlocks);
txAllSig = zeros(ldpcBlockLength/2, totalBlocks);

%----------------------------------------------------------------------
% Block-wise Processing: CCDM Shaping -> Bit Mapping -> LDPC -> Modulation
%----------------------------------------------------------------------
for blk = 1:totalBlocks
    blockBits = bitsMatrix(:, blk);
    
    % Step 1: CCDM encoding - shape information bits to non-uniform distribution
    shapedSymbols = ccdm.encode(blockBits, n_i).';

    % Step 2: Bit mapping - convert symbol indices to bit sequence (Gray code)
    mappedBits = bits_maping_gray(shapedSymbols);
    
    % Step 3: LDPC encoding - add redundancy for error protection
    codedBits = ldpcEncode(mappedBits, cfgEnc, 'OutputFormat', 'whole');
    
    % Step 4: Gray code remapping and modulation to PAM symbols
    txIdx = double(bits_remaping_gray(codedBits));
    txSig = A(txIdx + 1);
    txSig = ldpc_interleave(txSig, rateLDPC);
    txAllSig(:,blk) = txSig.';
end

% Flatten signal to column vector
originalSignal = txAllSig(:);

%----------------------------------------------------------------------
% Add Preamble Sequence
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

waveformFile = fullfile(outputDir, 'pas_tx_waveform_extreme_gray.csv');
knowledgeFile = fullfile(outputDir, 'Knowledge_gray.mat');

writematrix(finalSerialSignal, waveformFile);

save(knowledgeFile, 'finalSerialSignal', 'numZero', 'totalBlocks', ...
    'modulationOrder', 'srcBitStream', 'txAllSig', 'originalSignal', ...
    'k', 'ldpcBlockLength', 'rateLDPC', 'pVector', 'p_quant_I', 'A', 'cfgDec','numIter_pcs',...
    'n_i', 'bitsMatrix', 'mu_preamble');

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

% Calculate spectral efficiency
num_data_symbols = length(originalSignal);
total_info_bits = totalBlocks * k;
symbols_per_block = ldpcBlockLength / 2;
spectral_efficiency_pas = k / symbols_per_block;  % bits/symbol
spectral_efficiency_4pam = 2 * rateLDPC;          % bits/symbol (without PAS)

%----------------------------------------------------------------------
% Display Statistics
%----------------------------------------------------------------------
fprintf('\n=== Transmitter Statistics ===\n');
fprintf('Actual probability: [%.4f, %.4f, %.4f, %.4f]\n', p_actual);
fprintf('Entropy: %.6f bits/symbol\n', H_actual);
fprintf('Spectral efficiency: %.6f bits/symbol\n', spectral_efficiency_pas);
fprintf('  (Regular 4-PAM+LDPC: %.6f bits/symbol)\n', spectral_efficiency_4pam);
fprintf('\n=== Transmitter Configuration Complete ===\n');
