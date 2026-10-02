function C = correlation_factors(Nt, Nr, rho)
%CORRELATION_FACTORS  Matrix square roots of the Toeplitz correlation matrices.
%
%   C = CORRELATION_FACTORS(Nt, Nr, rho) builds the exponential (Toeplitz)
%   correlation matrices used by the Kronecker model
%
%       Rt = toeplitz(rho.^(0:Nt-1))      transmit side
%       Rr = toeplitz(rho.^(0:Nr-1))      receive side
%
%   and returns a struct with their principal square roots, so that the
%   (expensive) SQRTM call is done once per correlation value instead of
%   once per channel realisation.
%
%   Fields of C:  Nt, Nr, rho, sqrtRt, sqrtRr
%
%   rho = 0 gives identity matrices (ideal, uncorrelated channel);
%   rho -> 1 gives increasingly correlated antennas.
%
%   See also GENERATE_CHANNEL.

    if rho < 0 || rho >= 1
        error('correlation_factors:rho', 'rho must satisfy 0 <= rho < 1.');
    end

    C.Nt     = Nt;
    C.Nr     = Nr;
    C.rho    = rho;
    C.sqrtRt = sqrtm(toeplitz(rho.^(0:Nt-1)));
    C.sqrtRr = sqrtm(toeplitz(rho.^(0:Nr-1)));
end
