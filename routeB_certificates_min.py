#!/usr/bin/env python3
# routeB_certificates_min.py
# Minimal reproducible script to compute S(t), N(t), M2(t) bounds and solve calibration T(a)
# Requires: mpmath
import mpmath as mp
import csv
mp.mp.dps = 50

# group-specific data for SU(2) as example: use integer l>=1, dim d_l = 2l+1, Casimir c2(l)=l(l+1)
def S_N_trunc(t, Lmax=200):
    S = mp.mpf('0')
    N = mp.mpf('0')
    for l in range(1, Lmax+1):
        dl = 2*l+1
        c2 = l*(l+1)
        term = dl * mp.e**(-t * c2)
        S += term
        N += dl * c2 * mp.e**(-t * c2)
    return S, N

def M2_bound(t, Lmax=200):
    S, N = S_N_trunc(t, Lmax)
    if S >= 1:
        return None, S, N
    return N / (1 - S), S, N

# Solve for T(a) given target alpha = 1 - c*a via root find on a(a) = 2*CD*M2(T) = 1 - c*a
def find_T_for_a(a, c=0.8, CD=24, t_guess=3.5):
    target = 1 - c*a
    f = lambda T: 2*CD*(M2_bound(T)[0]) - target
    # bracket and find root
    lo = 0.5; hi = 20.0
    # ensure S<1 in bracket
    if M2_bound(lo)[0] is None:
        lo = 1.0
    try:
        T = mp.findroot(f, t_guess)
    except Exception:
        # fallback bisection
        for _ in range(100):
            mid = (lo+hi)/2
            mval = M2_bound(mid)[0]
            if mval is None:
                lo = mid
                continue
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        T = (lo+hi)/2
    m2, S, N = M2_bound(T)
    alpha = 1 - c*a
    muB = -mp.log( alpha )
    KA = 1 - 12 * m2
    return float(T), float(m2), float(S), float(alpha), float(muB), float(muB / a), float(KA)

# Example table for a in [0.5,0.25,0.125]
rows = []
for a in [0.5, 0.25, 0.125]:
    T, m2, S, alpha, muB, muB_a, KA = find_T_for_a(a, c=0.8, CD=24, t_guess=3.6)
    rows.append( (a, T, m2, S, alpha, muB, muB_a, KA) )

# write CSV
with open('routeB_certificates_min.csv','w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['a','T(a)','M2(T)','S(T)','alpha','muB','muB/a','KA'])
    for r in rows:
        w.writerow(r)

print("Wrote routeB_certificates_min.csv with rows:")
for r in rows:
    print(r)
