
import mpmath as mp
mp.mp.dps = 80

def su2_tail_bounds(t, L):
    f_S  = lambda x: (2*x+1)**2 * mp.e**(-t*x*(x+1))
    f_N  = lambda x: (2*x+1)**2 * (x*(x+1)) * mp.e**(-t*x*(x+1))
    S_tail = mp.quad(f_S, [L, mp.inf])
    N_tail = mp.quad(f_N, [L, mp.inf])
    return S_tail, N_tail

def su2_series_parts(t, Lmax=200):
    S = mp.mpf('0'); N = mp.mpf('0')
    for ell in range(1, Lmax+1):
        d = 2*ell + 1
        c2 = ell*(ell+1)
        e = mp.e**(-t*c2)
        S += d**2 * e
        N += d**2 * c2 * e
    return S, N

def su2_M2_upper(t, Lmax=200):
    S_part, N_part = su2_series_parts(t, Lmax=Lmax)
    S_tail, N_tail = su2_tail_bounds(t, Lmax)
    S = S_part + S_tail
    N = N_part + N_tail
    denom = 1 - S
    if denom <= 0:
        return mp.inf, S, N
    return N/denom, S, N

def find_tau_for_alpha(alpha_target, CD=24, Lmax=250):
    target_M2 = mp.mpf(alpha_target) / (2*CD)
    tL, tR = mp.mpf('0.5'), mp.mpf('12.0')
    for _ in range(80):
        M2R, _, _ = su2_M2_upper(tR, Lmax=Lmax)
        if M2R < target_M2:
            break
        tR *= 1.25
    for _ in range(90):
        tm = (tL + tR)/2
        M2m, _, _ = su2_M2_upper(tm, Lmax=Lmax)
        if M2m == mp.inf:
            tL = tm; continue
        if M2m > target_M2:
            tL = tm
        else:
            tR = tm
    t = tR
    M2t, S, _ = su2_M2_upper(t, Lmax=Lmax)
    alpha = 2*CD*M2t
    return t, M2t, S, alpha

if __name__ == "__main__":
    a_list = [0.5, 0.25, 0.125]
    c_calib = 0.8
    CD = 24
    print("# a, tau(a), M2, S, alpha, muB, muB/a, kappa_A")
    for a in a_list:
        alpha_target = max(1 - c_calib*a, 0.2)
        t, M2t, S, alpha = find_tau_for_alpha(alpha_target, CD=CD, Lmax=250)
        muB = -mp.log(alpha)
        kappa_A = 1 - 12*M2t
        print(float(a), float(t), float(M2t), float(S), float(alpha), float(muB), float(muB/a), float(kappa_A))
