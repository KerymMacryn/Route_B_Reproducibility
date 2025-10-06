
import mpmath as mp
mp.mp.dps = 70

def su3_dim(p, q):
    return (p+1)*(q+1)*(p+q+2)//2

def su3_casimir(p, q):
    return (p*p + q*q + p*q + 3*p + 3*q)/3.0

def su3_series_parts(t, Pmax=40):
    S = mp.mpf('0'); N = mp.mpf('0')
    for p in range(Pmax+1):
        for q in range(Pmax+1):
            if p==0 and q==0: continue
            d = su3_dim(p,q)
            c2 = su3_casimir(p,q)
            e = mp.e**(-t*c2)
            S += (d*d)*e
            N += (d*d)*c2*e
    return S, N

def su3_tail_bound(t, Pmax):
    C1 = 8.0; C2 = 1.0/9.0
    R = float(Pmax)
    u = t*C2*R*R
    tail_S = 0.5 * mp.gammainc(3.5, u) / ((t*C2)**3.5)
    tail_N = 0.5 * mp.gammainc(4.5, u) / ((t*C2)**4.5)
    return C1*tail_S, C1*tail_N

def su3_M2_upper(t, Pmax=40):
    S_part, N_part = su3_series_parts(t, Pmax=Pmax)
    S_tail, N_tail = su3_tail_bound(t, Pmax)
    S = S_part + S_tail
    N = N_part + N_tail
    denom = 1 - S
    if denom <= 0:
        return mp.inf, S, N
    return N/denom, S, N

def find_tau_for_alpha_su3(alpha_target, CD=24, Pmax=40):
    target_M2 = mp.mpf(alpha_target) / (2*CD)
    tL, tR = mp.mpf('0.25'), mp.mpf('15.0')
    for _ in range(50):
        M2R, _, _ = su3_M2_upper(tR, Pmax=Pmax)
        if M2R < target_M2: break
        tR *= 1.3
    for _ in range(90):
        tm = (tL + tR)/2
        M2m, _, _ = su3_M2_upper(tm, Pmax=Pmax)
        if M2m == mp.inf:
            tL = tm; continue
        if M2m > target_M2:
            tL = tm
        else:
            tR = tm
    t = tR
    M2t, S, _ = su3_M2_upper(t, Pmax=Pmax)
    alpha = 2*CD*M2t
    return t, M2t, S, alpha

if __name__ == "__main__":
    a_list = [0.5, 0.25, 0.125]
    c_calib = 0.8
    CD = 24
    Pmax = 42
    print("# a, tau(a), M2, S, alpha, muB, muB/a, kappa_A")
    for a in a_list:
        alpha_target = max(1 - c_calib*a, 0.25)
        t, M2t, S, alpha = find_tau_for_alpha_su3(alpha_target, CD=CD, Pmax=Pmax)
        muB = -mp.log(alpha)
        kappa_A = 1 - 12*M2t
        print(float(a), float(t), float(M2t), float(S), float(alpha), float(muB), float(muB/a), float(kappa_A))
