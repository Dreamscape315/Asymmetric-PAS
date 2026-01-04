%==========================================================================
% LDPC Symbol Deinterleaver (Inverse of ldpc_interleave)
%==========================================================================
% Description:
%   Deinterleaves LDPC-encoded symbols by separating alternated systematic
%   and parity symbols back to their original order. This is the inverse
%   operation of ldpc_interleave.
%
% Inputs:
%   y - Interleaved symbol sequence (row or column vector)
%       Format: Systematic and parity symbols alternated
%   rateLDPC - LDPC code rate (1/2, 2/3, 3/4, 5/6)
%
% Outputs:
%   x - Deinterleaved symbol sequence (same shape as input)
%       Format: [systematic1, systematic2, ..., parity1, parity2, ...]
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

function x = ldpc_deinterleave(y, rateLDPC)
    r = isrow(y); 
    y = y(:);
    N = numel(y);
    
    % Calculate number of systematic and parity symbols
    Nsys = floor(N * rateLDPC);  % Number of systematic symbols
    Npar = N - Nsys;             % Number of parity symbols
    
    % Determine interleaving pattern based on code rate
    sys_ratio = Nsys / Npar;
    if abs(sys_ratio - round(sys_ratio)) > 1e-6
        error('Code rate %.4f not supported for this interleaving scheme', rateLDPC);
    end
    sys_per_par = round(sys_ratio);  % Systematic symbols per group
    
    % Deinterleave: separate systematic and parity symbols from alternated sequence
    num_groups = Npar;
    sys_bits = zeros(Nsys, 1);
    par_bits = zeros(Npar, 1);
    
    for g = 1:num_groups
        % Current group starting position
        group_start = (g-1) * (sys_per_par + 1) + 1;
        
        % Extract systematic symbols
        sys_start = (g-1) * sys_per_par + 1;
        sys_end = min(g * sys_per_par, Nsys);
        num_sys = sys_end - sys_start + 1;
        sys_bits(sys_start : sys_end) = y(group_start : group_start + num_sys - 1);
        
        % Extract parity symbol
        par_bits(g) = y(group_start + num_sys);
    end
    
    % Restore sequence: [systematic; parity]
    x = [sys_bits; par_bits];
    
    % Restore original shape (row/column)
    if r
        x = x.'; 
    end
    
    % Safety check: ensure length unchanged
    assert(numel(x) == N, 'LDPC deinterleave changed sequence length');
end
