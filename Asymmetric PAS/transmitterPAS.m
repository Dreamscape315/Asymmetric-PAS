%==========================================================================
% True PAS Transmitter with Sign-Magnitude Representation
%==========================================================================
% Description:
%   Generates 4-PAM modulated signal using sign-magnitude representation
%   with probabilistic amplitude shaping (PAS). Amplitude bits are shaped
%   via CCDM while sign bits remain uniform for optimal rate.
%
% Processing Chain:
%   CCDM Shaping (amplitude) -> LDPC Encoding (amplitude+sign) ->
%   Sign Assignment -> Constellation Mapping
%
% Representation:
%   Amplitude: {1, 3} (2 levels, non-uniform via CCDM)
%   Sign: {-1, +1} (uniform via LDPC parity)
%   Constellation: {-3, -1, +1, +3}
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
numZero = 3000;                 % Preamble length
numIter_pcs = 25;               % LDPC iterations
rateLDPC = 2/3;                 % LDPC code rate
ldpcBlockLength = 1944;         % LDPC codeword length
m = 2;                          % Number of amplitude levels
v = 0.8;                        % Shaping parameter
totalBlocks = 150;              % Number of blocks

%----------------------------------------------------------------------
% Output Configuration
%----------------------------------------------------------------------
outputDir = 'example';
waveformFile = fullfile(outputDir, 'pas_true_tx_waveform.csv');
knowledgeFile = fullfile(outputDir, 'Knowledge_truePAS.mat');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%----------------------------------------------------------------------
% Distribution Matching Setup
%----------------------------------------------------------------------
pVector = generatePAS(m, v);
[cfgEnc, cfgDec] = getProtoMatrix(ldpcBlockLength, ldpcBlockLength * rateLDPC);
n = cfgEnc.BlockLength / m;
rateGamma = rateLDPC * m - (m - 1);
numExtraBits = round(rateGamma * n);
[p_quant_I, k, n_i] = ccdm.initialize(pVector, n);

%----------------------------------------------------------------------
% Rate Metrics
%----------------------------------------------------------------------
entropy = -sum(pVector .* log2(pVector + eps));
channel_per_use = entropy + rateGamma;
ccdmRate = k / n;
signRate = numExtraBits / n;
transmissionRate = ccdmRate + signRate;
rateGap = channel_per_use - transmissionRate;

fprintf('\n=== Rate Metrics ===\n');
fprintf('Target amplitude entropy (H(p)): %.6f bits/use\n', entropy);
fprintf('CCDM amplitude rate:            %.6f bits/use\n', ccdmRate);
fprintf('Sign bit rate:                  %.6f bits/use\n', signRate);
fprintf('Total transmission rate:        %.6f bits/use\n', transmissionRate);
fprintf('Theoretical rate (H(p)+gamma):  %.6f bits/use\n', channel_per_use);
fprintf('Rate gap (theory - actual):     %.6f bits/use\n', rateGap);

%----------------------------------------------------------------------
% Pre-allocation
%----------------------------------------------------------------------
constellationSymbolsAll = zeros(totalBlocks * n, 1);
originalSrcBitsAll_I = zeros(totalBlocks * k, 1);
originalExtraBitsAll_I = zeros(totalBlocks * numExtraBits, 1);

%----------------------------------------------------------------------
% Main Processing Loop
%----------------------------------------------------------------------
for blk = 1:totalBlocks
    % Generate random bit streams for this block
    srcBitStream_I = randi([0 1], 1, k + numExtraBits);

    % Separate source bits (for CCDM) and extra bits (for signs)
    srcBits_I = srcBitStream_I(1:k);
    extraBits_I = srcBitStream_I(k + 1:k + numExtraBits);

    % Store original source and extra bits for BER evaluation
    idxSrc = (blk - 1) * k + (1:k);
    originalSrcBitsAll_I(idxSrc) = srcBits_I;

    idxExtra = (blk - 1) * numExtraBits + (1:numExtraBits);
    originalExtraBitsAll_I(idxExtra) = extraBits_I;

    % CCDM encoding to determine amplitude levels
    ccdmSrcBits_I = ccdm.encode(srcBits_I, n_i);

    % Prepare input for LDPC encoder
    ldpcInput_I = [ccdmSrcBits_I, extraBits_I];

    % LDPC encoding to get parity bits (for signs)
    ldpcParity_I = ldpcEncode(ldpcInput_I', cfgEnc, 'OutputFormat', 'parity');

    % Constellation mapping using CCDM and LDPC outputs
    amplitude_I = ccdmSrcBits_I * 2 + 1;              % Map bits to {1, 3}
    signParity_I = ldpcParity_I * 2 - 1;              % Map parity to {-1, +1}
    signExtrabits_I = extraBits_I * 2 - 1;            % Map extra bits to {-1, +1}
    signVector = [signExtrabits_I'; signParity_I]';   % Combine sign sources

    constellationBlock = amplitude_I .* signVector;

    % Store this block's symbols
    idxSym = (blk - 1) * n + (1:n);
    constellationSymbolsAll(idxSym) = constellationBlock(:);
end

%----------------------------------------------------------------------
% Distribution Inspection
%----------------------------------------------------------------------
[symbolValues, ~, symbolIdx] = unique(constellationSymbolsAll);
symbolCounts = accumarray(symbolIdx, 1);
symbolProb = symbolCounts / numel(constellationSymbolsAll);
symbolDistribution = table(symbolValues, symbolCounts, symbolProb);

amplitudeValues = unique(abs(constellationSymbolsAll));
amplitudeCounts = arrayfun(@(val) sum(abs(constellationSymbolsAll) == val), amplitudeValues);
amplitudeProb = amplitudeCounts / numel(constellationSymbolsAll);
amplitudeDistribution = table(amplitudeValues, amplitudeCounts, amplitudeProb);

fprintf('\n=== Constellation Distribution ===\n');
disp(symbolDistribution);
fprintf('\n=== Amplitude Distribution (magnitude only) ===\n');
disp(amplitudeDistribution);

%----------------------------------------------------------------------
% Add Preamble Sequence
%----------------------------------------------------------------------
originalSignal = constellationSymbolsAll;
mu_preamble = mean(originalSignal);
preambleSignal = mu_preamble * ones(numZero, 1);
finalSerialSignal = [preambleSignal; constellationSymbolsAll];

%----------------------------------------------------------------------
% Save Results
%----------------------------------------------------------------------
writematrix(finalSerialSignal, waveformFile);
save(knowledgeFile, 'finalSerialSignal', 'constellationSymbolsAll', ...
    'originalSrcBitsAll_I', 'originalExtraBitsAll_I', 'pVector', ...
    'p_quant_I', 'k', 'n', 'n_i', 'numExtraBits', 'ccdmRate', ...
    'signRate', 'transmissionRate', 'channel_per_use', 'rateGap', ...
    'cfgEnc', 'cfgDec', 'rateLDPC', 'ldpcBlockLength', 'totalBlocks', ...
    'm', 'v', 'numIter_pcs', 'numZero', 'mu_preamble', 'originalSignal');

fprintf('\n=== Output Files ===\n');
fprintf('Waveform CSV saved to: %s\n', waveformFile);
fprintf('Knowledge MAT saved to: %s\n', knowledgeFile);
fprintf('\n=== Transmitter Configuration Complete ===\n');
