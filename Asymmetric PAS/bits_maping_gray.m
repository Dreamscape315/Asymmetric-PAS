%==========================================================================
% Gray Code Bit Mapping (Symbol Index to Bits)
%==========================================================================
% Description:
%   Converts symbol indices (0-3) to bit sequences using Gray code mapping.
%   Gray code ensures adjacent constellation points differ by only 1 bit,
%   minimizing bit error rate in noisy channels.
%
% Inputs:
%   a - Symbol indices (row/column vector, elements in {0,1,2,3})
%
% Outputs:
%   bits - 2N×1 double column vector (MSB-first, [b1;b0] concatenated)
%
% Gray Code Mapping (1-bit difference between adjacent points):
%   Index 0 -> [0,0] -> A(1) = -3
%   Index 1 -> [0,1] -> A(2) = -1
%   Index 2 -> [1,1] -> A(3) = +1
%   Index 3 -> [1,0] -> A(4) = +3
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

function bits = bits_maping_gray(a)
    a = uint8(a(:));
    N = numel(a);

    % Gray code conversion: index XOR (index >> 1)
    gray = bitxor(a, bitshift(a, -1));
    
    bits_tmp = zeros(2*N, 1, 'uint8');
    bits_tmp(1:2:end) = bitshift(gray, -1);  % b1 (MSB)
    bits_tmp(2:2:end) = bitand(gray, 1);     % b0 (LSB)

    bits = double(bits_tmp);
end
