function xhat = detect_zf(H, y)
%DETECT_ZF  Zero-Forcing detector.
%
%   xhat = DETECT_ZF(H, y) applies  W_ZF = (H'H)^-1 H'  to the received
%   vector(s) y. The pseudo-inverse is used for numerical robustness; it is
%   identical to (H'H)^-1 H' whenever H has full column rank (Nr >= Nt).
%
%   Removes inter-stream interference completely but amplifies noise when
%   H is ill-conditioned (e.g. at high spatial correlation).
%
%   See also DETECT_MMSE.

    xhat = pinv(H) * y;
end
