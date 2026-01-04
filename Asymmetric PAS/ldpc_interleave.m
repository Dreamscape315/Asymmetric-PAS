%==========================================================================
% LDPC Symbol Interleaver (Multiple Code Rates Supported)
%==========================================================================
% Description:
%   Interleaves LDPC-encoded symbols by alternating systematic and parity
%   symbols. Improves burst error resistance and decoding performance.
%
% Inputs:
%   x - Input symbol sequence (row or column vector)
%       Format: [systematic1, systematic2, ..., parity1, parity2, ...]
%   rateLDPC - LDPC code rate (1/2, 2/3, 3/4, 5/6)
%
% Outputs:
%   y - Interleaved symbol sequence (same shape as input)
%       Format: Systematic and parity symbols alternated
%
% Interleaving Pattern (based on code rate):
%   - 1/2: [S, P, S, P, S, P, ...]           1 systematic : 1 parity
%   - 2/3: [S, S, P, S, S, P, ...]           2 systematic : 1 parity
%   - 3/4: [S, S, S, P, S, S, S, P, ...]     3 systematic : 1 parity
%   - 5/6: [S, S, S, S, S, P, ...]           5 systematic : 1 parity
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

function y = ldpc_interleave(x, rateLDPC)
    r = isrow(x); 
    x = x(:); 
    N = numel(x);
    
    % Calculate number of systematic and parity symbols
    Nsys = floor(N * rateLDPC);  % Number of systematic symbols
    Npar = N - Nsys;             % Number of parity symbols
    
    % Extract systematic and parity bits
    sys_bits = x(1:Nsys);
    par_bits = x(Nsys+1:end);
    
    % Determine interleaving pattern based on code rate
    % Interleaving ratio = systematic : parity
    sys_ratio = Nsys / Npar;  % Systematic symbols per parity symbol
    
    % Verify ratio is integer (ensures perfect interleaving)
    if abs(sys_ratio - round(sys_ratio)) > 1e-6
        error('Code rate %.4f not supported (systematic/parity = %.2f is not integer)', rateLDPC, sys_ratio);
    end
    sys_per_par = round(sys_ratio);  % Systematic symbols per group
    
    % Interleave: reorganize systematic and parity bits
    % Each group contains sys_per_par systematic + 1 parity symbol
    num_groups = Npar;  % Number of groups = number of parity symbols
    y = zeros(N, 1);
    
    for g = 1:num_groups
        % Current group starting position
        group_start = (g-1) * (sys_per_par + 1) + 1;
        
        % Fill systematic symbols
        sys_start = (g-1) * sys_per_par + 1;
        sys_end = min(g * sys_per_par, Nsys);
        num_sys = sys_end - sys_start + 1;
        y(group_start : group_start + num_sys - 1) = sys_bits(sys_start : sys_end);
        
        % Fill parity symbol
        y(group_start + num_sys) = par_bits(g);
    end
    
    % Restore original shape (row/column)
    if r
        y = y.'; 
    end
    
    % Safety check: ensure length unchanged
    assert(numel(y) == N, 'LDPC interleave changed sequence length');
end
