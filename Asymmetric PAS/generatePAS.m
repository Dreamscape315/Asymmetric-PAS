%==========================================================================
% Generate Probabilistic Amplitude Shaping (PAS) Distribution
%==========================================================================
% Description:
%   Generates probability distribution for amplitude shaping based on
%   shaping parameter v. Controls trade-off between power efficiency
%   and spectral efficiency.
%
% Inputs:
%   pam_order - PAM modulation order (2 or 4)
%   v - Shaping parameter (0-1 range)
%       v=0: uniform distribution (no shaping)
%       v=1: maximum shaping (exponential distribution)
%
% Outputs:
%   pVector - Probability distribution array
%
% Examples:
%   pVector = generatePAS(4, 0.5);  % Moderate shaping
%   pVector = generatePAS(4, 0.8);  % Strong shaping
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

function [pVector] = generatePAS(pam_order, v)
    % Limit v to [0,1] range
    v = max(0, min(1, v));
    
    switch pam_order
        case 2
            % 2-PAM: Two probability values
            % v=0: [0.5, 0.5] (uniform)
            % v=1: [0.9, 0.1] (maximum shaping)
            p_high = 0.5 + 0.4 * v;  % Range from 0.5 to 0.9
            p_low = 1 - p_high;      % Range from 0.5 to 0.1
            pVector = [p_high, p_low];
            
        case 4
            % 4-PAM: Four probability values
            if v == 0
                % Uniform distribution
                pVector = [0.25, 0.25, 0.25, 0.25];
            else
                % Exponential decay distribution
                % Larger v results in steeper decay
                alpha = 1 + 3*v;  % alpha ranges from 1 to 4
                weights = exp(-(0:3) * alpha);
                pVector = weights / sum(weights);
            end
            
        otherwise
            error('Only 2-PAM and 4-PAM are currently supported');
    end
end
