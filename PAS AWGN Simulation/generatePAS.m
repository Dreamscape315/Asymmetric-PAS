%==========================================================================
% Generate Probabilistic Amplitude Shaping (PAS) Distribution
%==========================================================================
% Description:
%   Generates probability distribution for amplitude shaping based on
%   shaping parameter v. Returns [P(|x|=1), P(|x|=3)] using a simple
%   linear shaping model for 2-amplitude PAM.
%
% Inputs:
%   v - Shaping parameter (0-1 range)
%       v=0: uniform distribution (no shaping)
%       v=1: maximum shaping (strongest non-uniform distribution)
%
% Outputs:
%   pVector - Probability distribution array [P(|x|=1), P(|x|=3)]
%
% Examples:
%   pVector = generatePAS(0);    % Uniform: [0.5, 0.5]
%   pVector = generatePAS(0.5);  % Moderate shaping
%   pVector = generatePAS(1);    % Maximum shaping: [0.9, 0.1]
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

function [pVector] = generatePAS(v)
    % Limit v to [0,1] range
    v = max(0, min(1, v));
    
    % Linear shaping model
    % v=0: [0.5, 0.5] (uniform)
    % v=1: [0.9, 0.1] (maximum shaping)
    p_high = 0.5 + 0.4 * v;  % P(|x|=1): Range from 0.5 to 0.9
    p_low = 1 - p_high;      % P(|x|=3): Range from 0.5 to 0.1
    pVector = [p_high, p_low];
end
