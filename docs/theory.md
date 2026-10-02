# Theory and implementation notes

This page summarises the model behind the simulation and where each piece lives in the code.

## 1. System model

A flat-fading MIMO link with $N_t$ transmit and $N_r$ receive antennas:

$$\mathbf{y} = \mathbf{H}\mathbf{x} + \mathbf{n}$$

| Symbol | Meaning | Code |
| --- | --- | --- |
| $\mathbf{x}\in\{-1,+1\}^{N_t}$ | BPSK symbols (bit 0 → −1, bit 1 → +1) | `run_mimo_ber.m` |
| $\mathbf{H}\in\mathbb{C}^{N_r\times N_t}$ | Channel matrix | `generate_channel.m` |
| $\mathbf{n}\sim\mathcal{CN}(\mathbf{0},\sigma^2\mathbf{I})$ | AWGN | `run_mimo_ber.m` |

## 2. Correlated Rayleigh channel (Kronecker model)

$$\mathbf{H} = \mathbf{R}_r^{1/2}\,\mathbf{H}_w\,\mathbf{R}_t^{1/2}$$

- $\mathbf{H}_w$ has i.i.d. $\mathcal{CN}(0,1)$ entries (Rayleigh fading).
- $\mathbf{R}_t$ and $\mathbf{R}_r$ are exponential (Toeplitz) correlation matrices:

$$[\mathbf{R}]_{ij} = \rho^{|i-j|}, \qquad \mathbf{R}=\mathrm{toeplitz}(1,\rho,\rho^2,\dots)$$

$\rho = 0$ → identity (ideal, uncorrelated); $\rho\to 1$ → highly correlated antennas.
Implemented in `correlation_factors.m` (square roots are computed once per ρ for speed).

## 3. Detectors

**Zero Forcing** — inverts the channel, removing inter-stream interference but amplifying noise:

$$\mathbf{W}_{ZF} = (\mathbf{H}^H\mathbf{H})^{-1}\mathbf{H}^H \quad (\texttt{detect\_zf.m})$$

**MMSE** — regularises the inverse with the noise variance:

$$\mathbf{W}_{MMSE} = (\mathbf{H}^H\mathbf{H} + \sigma^2\mathbf{I})^{-1}\mathbf{H}^H \quad (\texttt{detect\_mmse.m})$$

Hard decision: $\hat b_i = \mathbb{1}\{\Re(\hat x_i) > 0\}$.

## 4. Metrics

**BER** = (number of wrong bits) / (total bits), averaged over random channel and noise realisations.

**Condition number** $\kappa(\mathbf{H}) = \sigma_{\max}/\sigma_{\min}$ (`avg_condition_number.m`).
For ZF the noise on each stream is scaled by the diagonal of $(\mathbf{H}^H\mathbf{H})^{-1}$, which blows up as
$\sigma_{\min}\to 0$, i.e. as $\kappa$ grows. Correlation pushes $\sigma_{\min}$ down, which is why $\kappa$ rises
with ρ and ZF degrades first.

## 5. SNR definition used here

With unit-energy symbols and $\mathbb{E}|h_{ij}|^2 = 1$:

$$\sigma^2 = 10^{-\mathrm{SNR_{dB}}/10}$$

So "SNR" is the per-stream, per-receive-antenna ratio $E_s/N_0$. It is **not** normalised by $N_t$ or by
$E_b$, so curves are comparable between configurations only under this convention.

## 6. A note on 2×2 vs 3×3

For spatial multiplexing with a linear detector, the diversity order is $N_r - N_t + 1$. For square systems
($N_t = N_r$) that is 1 regardless of size, so going 2×2 → 3×3 does not steepen the BER slope; it transmits
one more stream per channel use, and a larger random matrix tends to have a worse condition number. To buy
real diversity you need $N_r > N_t$ (e.g. 2×3), which this code supports out of the box:

```matlab
run_mimo_ber(2, 3, 'both', 0.7);   % Nt = 2, Nr = 3
```

## 7. Known simplifications

- Flat fading (no frequency selectivity / OFDM), BPSK only, perfect CSI at the receiver
- Linear detectors only (no SIC, ML or sphere decoding)
- Same noise sample repeated across the 20 spectrum columns (only used for the spectrum plot; detection
  uses the first column)
