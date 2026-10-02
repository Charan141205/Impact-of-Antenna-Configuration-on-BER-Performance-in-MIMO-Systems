function H = generate_channel(C)
%GENERATE_CHANNEL  One realisation of a spatially correlated Rayleigh MIMO channel.
%
%   H = GENERATE_CHANNEL(C) returns an Nr-by-Nt complex channel matrix
%
%       H = Rr^(1/2) * Hw * Rt^(1/2)
%
%   where Hw has i.i.d. CN(0,1) entries (Rayleigh fading) and Rr^(1/2),
%   Rt^(1/2) come from CORRELATION_FACTORS (Kronecker correlation model).
%
%   See also CORRELATION_FACTORS.

    Hw = (randn(C.Nr, C.Nt) + 1j*randn(C.Nr, C.Nt)) / sqrt(2);
    H  = C.sqrtRr * Hw * C.sqrtRt;
end
