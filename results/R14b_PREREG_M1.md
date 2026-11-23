# PRE-REGISTRATION — collision primes of J_0(23)  (frozen 2026-11-19T17:23:03.502750Z, before any scan beyond p<=101)
Collision: p != 2,23 with a_p(f) = a_p(f^sigma), i.e. v_p = 0 where a_p(f) = u_p + v_p*phi, phi^2=phi+1.
Known before freezing: collisions <=101 are 13,19,29,43 (R14).

## T (theorem predictions, must hold 100%; one failure falsifies the proof or its imported inputs)
T1  every odd p != 23: v_p even.
T2  every p != 23: a_p(f) - 1 - p in prime P11 above 11 (one fixed root r in {4,8} of phi mod 11 for all p).
T3  every collision: a_p == p+1 (mod 11).
T4  every collision: a_p odd  <=>  p represented by 2x^2+xy+3y^2.

## M1 model (heuristic; Haar on model l-adic images, independence across l, SU(2)xSU(2) Sato-Tate)
F_2 = 22/9 (image mod 2 in GL2(F2), full preimage), F_5 ~ 1.03472 (ramified, full), F_11 ~ 1.09924 (Eisenstein, diag model),
generic Van Hirtum factors elsewhere (incl. l=3).  F_hat = 3.17290;  Prob(collision at p) = 1.916948 * p^(-1/2) ... i.e. K = 1.916948.
N_pred(x) = K * sum_{3<=p<=x, p!=23} p^(-1/2):
  x=   100000: N_pred = 132.5
  x=  1000000: N_pred = 336.3
  x= 10000000: N_pred = 883.2
  x= 20000000: N_pred = 1187.9
Success band for M1: |N_obs - N_pred| <= 3*sqrt(N_pred) at every x  (Poisson-type noise only).
Disagreement beyond the band => M1 wrong; diagnose with E-tests below (which l is off).

## E (local frequency tests over all good odd primes p<=X, binomial 3-sigma)
E2: P(4|v)=0.583333, P(8|v)=0.302083, P(16|v)=0.152344
E3: P(3|v)=0.362500, P(9|v)=0.121759
E5: P(5|v)=0.206667, P(25|v)=0.041387
E11: P(11|v)=0.099167, P(121|v)=0.009084

## H_gap (Van Hirtum's 10% excess): ratio sum p^-1/2 / (2 sqrt x/log x) = 1.174 at 1e7 (computed), decreasing; predicts any
fit of c*sqrt(x)/log x to data <=1e8 overshoots the asymptotic c by ~10-17%.
