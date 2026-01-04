%==========================================================================
% Log-Likelihood Ratio (LLR) Calculation for Gray-coded 4-PAM
%==========================================================================
% Description:
%   Computes bit-wise LLRs for Gray-coded 4-PAM with probabilistic
%   amplitude shaping (PAS). Supports non-uniform amplitude distributions
%   for systematic bits and uniform distribution for parity bits.
%
% Inputs:
%   y - Received amplitude sequence (1xN or Nx1)
%   p_data - Prior probability for data symbols (1x4, from CCDM/shaping)
%   noiseStd - Noise standard deviation (scalar or vector matching y)
%   A - Decision levels (1x4), e.g., [-3 -1 1 3]
%   ldpcBlockLength - LDPC codeword length (e.g., 648, 1296, 1944)
%   rateLDPC - LDPC code rate (e.g., 1/2, 2/3, 3/4, 5/6)
%
% Outputs:
%   llr_bits - Bit LLRs (1x(2*N)), format: [b1, b0, b1, b0, ...]
%              Systematic bits first, then parity bits
%
% Gray Code Mapping:
%   A(1) = -3 <-> [0,0]
%   A(2) = -1 <-> [0,1]
%   A(3) = +1 <-> [1,1]
%   A(4) = +3 <-> [1,0]
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

function [llr_bits] = llr_pam_gray(y, p_data, noiseStd, A, ldpcBlockLength, rateLDPC)
    y = y(:).';
    Ntot = numel(y);
    
    % Calculate number of systematic and parity symbols from LDPC parameters
    m = 2;  % PAM4: 2 bits per symbol
    Nsys_sym = (ldpcBlockLength * rateLDPC) / m;  % Systematic symbols
    Npar_sym = (ldpcBlockLength / m) - Nsys_sym;  % Parity symbols
    assert(Nsys_sym + Npar_sym == Ntot, 'Length of y does not match code rate/block length');
    
    % Prior probabilities (uniform for parity region)
    p_parity = [0.25 0.25 0.25 0.25];
    
    % Noise variance
    if isscalar(noiseStd)
        noiseVar = (noiseStd^2) * ones(1, Ntot);
        if noiseStd == 0, noiseVar(:) = 1e-10; end
    else
        noiseVar = (noiseStd(:).').^2;
        noiseVar(noiseVar == 0) = 1e-10;
        assert(numel(noiseVar) == Ntot, 'noiseStd length must match y');
    end
    
    % Gray code label definition
    % A(1)=-3: [0,0], A(2)=-1: [0,1], A(3)=+1: [1,1], A(4)=+3: [1,0]
    labels = [0 0; 0 1; 1 1; 1 0];  % Row index corresponds to A index 1..4
    
    idx_b1_0 = find(labels(:,1) == 0).';  % A(1), A(2) -> -3, -1
    idx_b1_1 = find(labels(:,1) == 1).';  % A(3), A(4) -> +1, +3
    idx_b0_0 = find(labels(:,2) == 0).';  % A(1), A(4) -> -3, +3
    idx_b0_1 = find(labels(:,2) == 1).';  % A(2), A(3) -> -1, +1
    
    % Partition indices (systematic vs parity)
    data_idx   = 1:Nsys_sym;                      % Systematic symbols (use p_data)
    parity_idx = (Nsys_sym+1):(Nsys_sym+Npar_sym);% Parity symbols (uniform)
    
    % Construct prior matrix 4xNtot (data region: p_data, parity region: p_parity)
    priors = repmat(p_parity(:), 1, Ntot);
    priors(:, data_idx) = repmat(p_data(:), 1, numel(data_idx));
    log_prior = log(priors + eps);  % 4xNtot
    
    % Log-likelihood (4xNtot)
    A = A(:).';
    dist = (y - A.').^2;              % 4xNtot
    log_like = -dist ./ (2*noiseVar); % 4xNtot
    
    log_post = log_like + log_prior;  % 4xNtot
    
    % Bit LLRs (log-sum-exp per bit position)
    lp_b1_0 = logsumexp_rows(log_post, idx_b1_0);  % 1xNtot
    lp_b1_1 = logsumexp_rows(log_post, idx_b1_1);
    lp_b0_0 = logsumexp_rows(log_post, idx_b0_0);
    lp_b0_1 = logsumexp_rows(log_post, idx_b0_1);
    
    llr_b1 = lp_b1_0 - lp_b1_1;  % 1xNtot (MSB)
    llr_b0 = lp_b0_0 - lp_b0_1;  % 1xNtot (LSB)
    
    % LLR clipping for numerical stability
    L = 40;
    llr_b1 = max(-L, min(L, llr_b1));
    llr_b0 = max(-L, min(L, llr_b0));
    
    % Concatenate [b1, b0] for each symbol to get 2*Ntot bits
    llr_bits = reshape([llr_b1; llr_b0], 1, []);  % 1x(2*Ntot)
end

function out = logsumexp_rows(M, row_idx)
    % Numerically stable log-sum-exp over selected rows
    S = M(row_idx, :);        % |row_idx| x N
    m = max(S, [], 1);        % 1 x N
    out = m + log(sum(exp(S - m), 1) + eps);
end
