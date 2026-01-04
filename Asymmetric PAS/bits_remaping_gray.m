%==========================================================================
% Gray Code Bit Remapping (Bits to Symbol Index)
%==========================================================================
% Description:
%   Converts bit sequences back to symbol indices (0-3) using Gray code
%   inverse mapping. This is the inverse operation of bits_maping_gray.
%
% Inputs:
%   bits - Binary bit sequence (row/column vector, length = 2N, MSB-first)
%
% Outputs:
%   dec - N×1 double, Gray code corresponding decimal integers (0..3)
%
% Gray Code Mapping (1-bit difference between adjacent points):
%   [0,0] -> Index 0 -> A(1) = -3
%   [0,1] -> Index 1 -> A(2) = -1
%   [1,1] -> Index 2 -> A(3) = +1
%   [1,0] -> Index 3 -> A(4) = +3
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

function dec = bits_remaping_gray(bits)
    bits = uint8(bits(:));
    assert(mod(numel(bits), 2) == 0, 'Bit sequence length must be even');

    b1 = bits(1:2:end);   % MSB
    b0 = bits(2:2:end);   % LSB
    
    % Combine into Gray code
    gray = bitshift(b1, 1) + b0;
    
    % Gray to binary inverse transformation
    % For 2-bit: binary = gray XOR (gray >> 1)
    dec = double(bitxor(gray, bitshift(gray, -1)));
end
