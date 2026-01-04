%==========================================================================
% PAS-Aware LLR Calculation for 4-PAM with CCDM Shaping
%==========================================================================
% Description:
%   Computes PAS-aware log-likelihood ratios for 4-PAM constellation
%   {-3, -1, +1, +3} with non-uniform amplitude distribution. Uses
%   amplitude priors from CCDM shaping while assuming uniform sign bits.
%
% Inputs:
%   y - Received symbols (column vector, length n)
%   pAmp - Amplitude prior probabilities [P(|x|=1), P(|x|=3)]
%   noiseVar - AWGN noise variance (sigma^2)
%   n - Number of amplitude symbols per block
%   numExtraBits - Number of systematic sign bits
%   ldpcBlockLength - Total LDPC codeword length
%   llrClip - LLR clipping value for numerical stability (default: 40)
%
% Outputs:
%   llr_block - LLR values for entire LDPC block (column vector)
%               Structure: [amplitude LLRs; sign LLRs; parity sign LLRs]
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

function llr_block = llr_pas_true(y, pAmp, noiseVar, n, numExtraBits, ldpcBlockLength, llrClip)
    %----------------------------------------------------------------------
    % Input Validation and Initialization
    %----------------------------------------------------------------------
    if nargin < 7 || isempty(llrClip)
        llrClip = 40;
    end
    
    y = y(:)';
    sigma2 = noiseVar;
    if sigma2 <= 0
        sigma2 = 1e-8;
    end
    
    % Normalize amplitude priors
    pAmp = pAmp(:)';
    pAmp = pAmp / sum(pAmp);
    
    % Precompute log-priors
    logHalf = log(0.5);
    logAmp1 = log(pAmp(1));  % log P(|x|=1)
    logAmp3 = log(pAmp(2));  % log P(|x|=3)
    
    %----------------------------------------------------------------------
    % Amplitude Bit LLR: |x|=1 (bit 0) vs |x|=3 (bit 1)
    %----------------------------------------------------------------------
    % P(|x|=1|y) = P(x=-1|y) + P(x=+1|y)
    logNum1 = logAmp1 + logHalf - ((y - 1).^2) ./ (2 * sigma2);
    logNum2 = logAmp1 + logHalf - ((y + 1).^2) ./ (2 * sigma2);
    
    % P(|x|=3|y) = P(x=-3|y) + P(x=+3|y)
    logDen1 = logAmp3 + logHalf - ((y - 3).^2) ./ (2 * sigma2);
    logDen2 = logAmp3 + logHalf - ((y + 3).^2) ./ (2 * sigma2);
    
    ampLLR = logsumexp_pair(logNum1, logNum2) - logsumexp_pair(logDen1, logDen2);
    
    %----------------------------------------------------------------------
    % Sign Bit LLR: sign=-1 (bit 0) vs sign=+1 (bit 1)
    %----------------------------------------------------------------------
    % P(sign=-1|y) = P(x=-1|y) + P(x=-3|y)
    logNeg1 = logAmp1 - ((y + 1).^2) ./ (2 * sigma2);
    logNeg2 = logAmp3 - ((y + 3).^2) ./ (2 * sigma2);
    
    % P(sign=+1|y) = P(x=+1|y) + P(x=+3|y)
    logPos1 = logAmp1 - ((y - 1).^2) ./ (2 * sigma2);
    logPos2 = logAmp3 - ((y - 3).^2) ./ (2 * sigma2);
    
    signLLR = logsumexp_pair(logNeg1, logNeg2) - logsumexp_pair(logPos1, logPos2);
    
    %----------------------------------------------------------------------
    % LLR Clipping for Numerical Stability
    %----------------------------------------------------------------------
    ampLLR = max(-llrClip, min(llrClip, ampLLR));
    signLLR = max(-llrClip, min(llrClip, signLLR));
    
    %----------------------------------------------------------------------
    % Assemble Full LDPC Block LLRs
    %----------------------------------------------------------------------
    parityLen = ldpcBlockLength - (n + numExtraBits);
    if parityLen < 0
        error('Invalid LDPC block configuration: parity length < 0');
    end
    
    llr_block = zeros(1, ldpcBlockLength);
    llr_block(1:n) = ampLLR;                                        % Amplitude bits
    llr_block(n + (1:numExtraBits)) = signLLR(1:numExtraBits);     % Systematic sign bits
    llr_block(n + numExtraBits + (1:parityLen)) = signLLR(numExtraBits + 1 : numExtraBits + parityLen);  % Parity sign bits
    llr_block = llr_block(:);
end

%----------------------------------------------------------------------
% Helper Function: Numerically Stable Log-Sum-Exp for Pairs
%----------------------------------------------------------------------
function s = logsumexp_pair(a, b)
    % Computes log(exp(a) + exp(b)) avoiding overflow/underflow
    mx = max(a, b);
    mn = min(a, b);
    s = mx + log1p(exp(mn - mx));
end
