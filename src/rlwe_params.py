#!/usr/bin/env python3
"""Ring-LWE parameter estimates for the BGV instantiation over R_0[Y]/(Y^h).

Noise model (Section "Parameter selection for the Ring-LWE instantiation"):
  B_e   = ceil(6 sigma)                      tail cutoff of the error
  delta = 2 sqrt(N_0 h)                      heuristic coefficient expansion
  F     = B_e (1 + 2 delta)                  noise of a fresh encryption of zero
  B_0   = F                                  noise of an issued tag
  B_{j+1} = d (B_j + F)                      a transition merges at most d-1 degraded
                                             tags into one tag and rerandomizes once;
                                             multiplication by Y adds no noise
  B_L   = F (d^L (2d-1) - d) / (d-1)         for d >= 2, and (L+1) F for d = 1
  log2 Q = ceil(1 + log2 t + log2 B_L)       decryption is correct while B_L < Q/(2t) - 1/2
N_0 is the smallest dimension whose 128-bit limit in the Homomorphic Encryption
Security Standard (ternary secret, sigma = 3.2) is at least log2 Q.
Tag size |w| = 2 N_0 h ceil(log2 Q) bits.
margin_bits = max_log2Q - log2Q; max_d and max_L are the largest fan-in (at
horizon L) and the largest horizon (at fan-in d) whose log2 Q still fits the
limit of the same N_0.

The plaintext modulus t is 65537 while 65537 > 2^(h-1) and 2 N_0 divides
65536; otherwise it is the smallest probable prime t > 2^(h-1) with
t = 1 mod 2 N_0 (Miller-Rabin, 40 random bases).

Usage: python3 rlwe_params.py [h ...] > file.csv   (default h = 16)
       python3 rlwe_params.py 16 > ../results_rlwe_params.csv
       python3 rlwe_params.py 10 > ../results_rlwe_params_h10.csv
       python3 rlwe_params.py 20 50 100 > ../results_rlwe_params_large_h.csv
"""
import math
import random
import sys

HESTD_128 = [(1024, 27), (2048, 54), (4096, 109), (8192, 218), (16384, 438), (32768, 881)]


def log2_bound(d, L, F):
    if d == 1:
        return math.log2((L + 1) * F)
    return math.log2(F) + math.log2((d ** L * (2 * d - 1) - d) / (d - 1))


def is_probable_prime(n, rounds=40):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d, r = d // 2, r + 1
    rng = random.Random(n)
    for _ in range(rounds):
        x = pow(rng.randrange(2, n - 1), d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def plaintext_modulus(h, n0):
    if 65537 > 2 ** (h - 1) and 65536 % (2 * n0) == 0:
        return 65537
    m = 2 * n0
    t = (2 ** (h - 1) // m + 1) * m + 1
    while not is_probable_prime(t):
        t += m
    return t


def estimate(h, L, d, sigma=3.2):
    be = math.ceil(6 * sigma)
    for n0, qmax in HESTD_128:
        t = plaintext_modulus(h, n0)
        delta = 2 * math.sqrt(n0 * h)
        F = be * (1 + 2 * delta)
        logq = math.ceil(1 + math.log2(t) + log2_bound(d, L, F))
        if logq <= qmax:
            return n0, qmax, logq, F, t
    return None


def log2q(n0, h, t, L, d, sigma=3.2):
    F = math.ceil(6 * sigma) * (1 + 4 * math.sqrt(n0 * h))
    return math.ceil(1 + math.log2(t) + log2_bound(d, L, F))


def main():
    rows = [("linear_chain", 1), ("two_state", 2), ("merge_heavy", 4)]
    hs = [int(a) for a in sys.argv[1:]] or [16]
    L = 32
    w = sys.stdout.write
    w("profile,h,t,L,d,N0,max_log2Q,log2Q,log2F,N,tag_bits,tag_bytes,tag_MiB,margin_bits,max_d,max_L\n")
    for h, (name, d) in [(h, r) for h in hs for r in rows]:
        n0, qmax, logq, F, t = estimate(h, L, d)
        bits = 2 * n0 * h * logq
        md = d
        while log2q(n0, h, t, L, md + 1) <= qmax:
            md += 1
        lo, hi = L, L
        while log2q(n0, h, t, hi, d) <= qmax:
            lo, hi = hi, 2 * hi
        while hi - lo > 1:
            mid = (lo + hi) // 2
            lo, hi = (mid, hi) if log2q(n0, h, t, mid, d) <= qmax else (lo, mid)
        ml = lo
        w(f"{name},{h},{t},{L},{d},{n0},{qmax},{logq},{math.log2(F):.3f},{n0 * h},{bits},{bits // 8},{bits / 8 / 2**20:.3f},{qmax - logq},{md},{ml}\n")


if __name__ == "__main__":
    main()
