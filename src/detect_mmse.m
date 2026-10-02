function xhat = detect_mmse(H, y, noiseVar)
%DETECT_MMSE  Linear MMSE detector (unit-energy symbols).
%
%   xhat = DETECT_MMSE(H, y, noiseVar) applies
%
%       W_MMSE = (H'H + noiseVar*I)^-1 H'
%
%   to the received vector(s) y. The regularisation term trades a small
%   amount of residual interference for much lower noise enhancement, so
%   MMSE degrades far more gracefully than ZF on correlated channels.
%
%   As noiseVar -> 0 the MMSE detector converges to ZF.
%
%   See also DETECT_ZF.

    Nt   = size(H, 2);
    xhat = (H'*H + noiseVar*eye(Nt)) \ (H' * y);
end
