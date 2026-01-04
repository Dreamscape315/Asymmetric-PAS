%==========================================================================
% LLR Calculation for True PAS (Sign-Magnitude Representation)
%==========================================================================
% Description:
%   Computes PAS-aware LLRs for 4-PAM with sign-magnitude representation.
%   Amplitude bits use non-uniform priors from shaping, while sign bits
%   are uniformly distributed.
%
% Inputs:
%   y - Received amplitudes (column vector, length n)
%   pAmp - Amplitude prior probabilities [P(|x|=1), P(|x|=3)]
%   noiseVar - AWGN variance
%   n - Number of amplitude symbols per block
%   numExtraBits - Number of systematic sign bits
%   ldpcBlockLength - Total LDPC codeword length (amplitude + sign bits)
%   llrClip - Clipping value for numerical stability (default: 40)
%
% Outputs:
%   llr_block - LLR vector for entire LDPC block (column vector)
%               [amplitude LLRs; systematic sign LLRs; parity sign LLRs]
%
% Note:
%   For 4-PAM with |x| in {1,3} and sign in {-1,+1}:
%   Constellation: {-3, -1, +1, +3}
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

function llr_block = llr_pas_true(y, pAmp, noiseVar, n, numExtraBits, ldpcBlockLength, llrClip)
    if nargin < 7 || isempty(llrClip)
        llrClip = 40;
    end
    
    y = y(:)';
    sigma2 = noiseVar;
    if sigma2 <= 0
        sigma2 = 1e-8;
    end
    
    pAmp = pAmp(:)';
    pAmp = pAmp / sum(pAmp);
    
    logHalf = log(0.5);
    logAmp1 = log(pAmp(1));
    logAmp3 = log(pAmp(2));
    
    % Precompute terms for amplitude LLR (mixture of +/-1 vs +/-3)
    logNum1 = logAmp1 + logHalf - ((y - 1).^2) ./ (2 * sigma2);
    logNum2 = logAmp1 + logHalf - ((y + 1).^2) ./ (2 * sigma2);
    logDen1 = logAmp3 + logHalf - ((y - 3).^2) ./ (2 * sigma2);
    logDen2 = logAmp3 + logHalf - ((y + 3).^2) ./ (2 * sigma2);
    
    ampLLR = logsumexp_pair(logNum1, logNum2) - logsumexp_pair(logDen1, logDen2);
    
    % Sign-bit LLR: sign = -1 (bit 0) vs +1 (bit 1)
    logNeg1 = logAmp1 - ((y + 1).^2) ./ (2 * sigma2);
    logNeg2 = logAmp3 - ((y + 3).^2) ./ (2 * sigma2);
    logPos1 = logAmp1 - ((y - 1).^2) ./ (2 * sigma2);
    logPos2 = logAmp3 - ((y - 3).^2) ./ (2 * sigma2);
    
    signLLR = logsumexp_pair(logNeg1, logNeg2) - logsumexp_pair(logPos1, logPos2);
    
    % Apply LLR clipping
    ampLLR = max(-llrClip, min(llrClip, ampLLR));
    signLLR = max(-llrClip, min(llrClip, signLLR));
    
    % Compute parity length
    parityLen = ldpcBlockLength - (n + numExtraBits);
    if parityLen < 0
        error('Invalid LDPC block configuration: parity length < 0');
    end
    
    % Assemble LLR block: [amplitude; systematic signs; parity signs]
    llr_block = zeros(1, ldpcBlockLength);
    llr_block(1:n) = ampLLR;
    llr_block(n + (1:numExtraBits)) = signLLR(1:numExtraBits);
    llr_block(n + numExtraBits + (1:parityLen)) = signLLR(numExtraBits + 1 : numExtraBits + parityLen);
    llr_block = llr_block(:);
end

function s = logsumexp_pair(a, b)
    % Numerically stable computation of log(exp(a) + exp(b))
    mx = max(a, b);
    mn = min(a, b);
    s = mx + log1p(exp(mn - mx));
end
