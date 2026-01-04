%==========================================================================
% PAM LLR Calculation for AWGN Channel
%==========================================================================
% Description:
%   Computes bit-wise log-likelihood ratios (LLRs) for M-PAM modulation
%   over real-valued AWGN channel. Uses max-log approximation for
%   numerical stability.
%
% Inputs:
%   rxSymbols - Received symbols (real-valued vector)
%   M - PAM modulation order (2, 4, 8, ...)
%   noiseVariance - AWGN noise variance (sigma^2)
%   constellation - PAM constellation points (M×1 vector, ordered by index)
%
% Outputs:
%   llr - Log-likelihood ratios (column vector, length = bitsPerSymbol × numSymbols)
%         LLR definition: log(P(bit=0|rx) / P(bit=1|rx))
%
% Examples:
%   constellation = [-3, -1, 1, 3];
%   llr = pamLLR(rxSymbols, 4, 0.1, constellation);
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

function llr = pamLLR(rxSymbols, M, noiseVariance, constellation)
    bitsPerSymbol = log2(M);
    symbolOrder = (0:M-1).';
    
    llr = zeros(numel(rxSymbols) * bitsPerSymbol, 1);
    invTwoVar = 1 / (2 * noiseVariance);
    
    for symIdx = 1:numel(rxSymbols)
        rx = rxSymbols(symIdx);
        
        % Log-likelihood metric per hypothesis (up to an additive constant)
        metricLog = -((rx - constellation).^2) * invTwoVar;
        
        % Partition by transmitter (natural) bit labels
        for bitPos = 1:bitsPerSymbol
            currentBits = bitget(symbolOrder, bitsPerSymbol - bitPos + 1);
            lse0 = localLogSumExp(metricLog(currentBits == 0));
            lse1 = localLogSumExp(metricLog(currentBits == 1));
            llr((symIdx-1)*bitsPerSymbol + bitPos) = lse0 - lse1;
        end
    end
end

%----------------------------------------------------------------------
% Helper Function: Numerically Stable Log-Sum-Exp
%----------------------------------------------------------------------
function s = localLogSumExp(x)
    % Computes log(sum(exp(x))) avoiding overflow/underflow
    m = max(x);
    s = m + log(sum(exp(x - m)));
end