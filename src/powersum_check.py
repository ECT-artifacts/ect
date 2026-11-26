#!/usr/bin/env python3
"""Check the power-sum encoding of Section "More live identifiers".

Part 1 checks the decoding of power sums over F_P, with P the first prime above 2^20: Berlekamp-Massey, the roots of the recurrence among the points of the issued identifiers, and the multiplicities by Gaussian elimination.  It decodes random multisets of at most K points, with K = 59 and K = 32, and once with all 2^20 identifiers issued; then it decodes depths that hold more than K points, which the tracer cannot recover.

Part 2 runs the Damgard-Jurik variant end to end with a 3072-bit modulus: real encryption, degradation, merging, decryption and tracing, on random histories of boards, ingestions and payments.  Every tag of every history is traced and compared with its live set, computed from the contributions.  Two configurations of Table "powersum" are run: both tags degraded at h = 6 with K = 32, and the capped policy at h = 10 with K = 59.

Part 3 fills one depth of a Damgard-Jurik tag to exactly K distinct points, through real Tag, Merge and Degrade, and traces it; then it traces a depth of K+1 points, which gives bottom, and the tag degraded once more, in which the K points have expired.

Modular exponentiation calls GMP through ctypes when libgmp is installed, and Python's pow otherwise (much slower).  The randomness is seeded, the encryption randomness included, so a run is reproducible; the script is a correctness check, not an implementation for use.

usage: python3 powersum_check.py [seed] > ../results_powersum_check.txt
"""
import ctypes
import ctypes.util
import math
import random
import sys
import time
from collections import Counter

# ---------------------------------------------------------------------------
# Modular exponentiation through GMP.

class _Mpz(ctypes.Structure):
    _fields_ = [("alloc", ctypes.c_int), ("size", ctypes.c_int), ("d", ctypes.c_void_p)]

_gmp = None
_lib = ctypes.util.find_library("gmp")
if _lib:
    try:
        _gmp = ctypes.CDLL(_lib)
    except OSError:
        _gmp = None

if _gmp:
    _P = ctypes.POINTER(_Mpz)
    _gmp.__gmpz_init.argtypes = [_P]
    _gmp.__gmpz_clear.argtypes = [_P]
    _gmp.__gmpz_import.argtypes = [_P, ctypes.c_size_t, ctypes.c_int, ctypes.c_size_t, ctypes.c_int, ctypes.c_size_t, ctypes.c_void_p]
    _gmp.__gmpz_export.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t), ctypes.c_int, ctypes.c_size_t, ctypes.c_int, ctypes.c_size_t, _P]
    _gmp.__gmpz_export.restype = ctypes.c_void_p
    _gmp.__gmpz_powm.argtypes = [_P, _P, _P, _P]
    _gmp.__gmpz_nextprime.argtypes = [_P, _P]
    _gmp.__gmpz_sizeinbase.argtypes = [_P, ctypes.c_int]
    _gmp.__gmpz_sizeinbase.restype = ctypes.c_size_t

    def _to_mpz(x, z):
        b = x.to_bytes((x.bit_length() + 7) // 8 or 1, "big")
        _gmp.__gmpz_import(ctypes.byref(z), len(b), 1, 1, 1, 0, b)

    def _from_mpz(z):
        n = (_gmp.__gmpz_sizeinbase(ctypes.byref(z), 2) + 7) // 8
        buf = ctypes.create_string_buffer(n + 1)
        cnt = ctypes.c_size_t(0)
        _gmp.__gmpz_export(buf, ctypes.byref(cnt), 1, 1, 1, 0, ctypes.byref(z))
        return int.from_bytes(buf.raw[:cnt.value], "big")

    def _call(f, *args):
        zs = [_Mpz() for _ in range(len(args) + 1)]
        for z in zs:
            _gmp.__gmpz_init(ctypes.byref(z))
        for z, a in zip(zs[1:], args):
            _to_mpz(a, z)
        f(*(ctypes.byref(z) for z in zs))
        out = _from_mpz(zs[0])
        for z in zs:
            _gmp.__gmpz_clear(ctypes.byref(z))
        return out

    def powm(b, e, m):
        return _call(_gmp.__gmpz_powm, b % m, e, m)

    def nextprime(x):
        return _call(_gmp.__gmpz_nextprime, x)
else:
    def powm(b, e, m):
        return pow(b, e, m)

    def _is_probable_prime(n, rounds=40):
        if n < 2:
            return False
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
            if n % p == 0:
                return n == p
        d, s = n - 1, 0
        while d % 2 == 0:
            d, s = d // 2, s + 1
        rng = random.Random(n)
        for _ in range(rounds):
            x = pow(rng.randrange(2, n - 1), d, n)
            if x in (1, n - 1):
                continue
            for _ in range(s - 1):
                x = x * x % n
                if x == n - 1:
                    break
            else:
                return False
        return True

    def nextprime(x):
        x += 1
        while not _is_probable_prime(x):
            x += 1
        return x

# ---------------------------------------------------------------------------
# Power sums over F_P (Definition "Power-sum encoding", Lemma "Decoding").

def first_prime_above(x):
    return nextprime(x)

def power_sums(mult, K, P):
    """s_j = sum_x c(x) x^j mod P for j = 0..2K-1."""
    s = [0] * (2 * K)
    for x, c in mult.items():
        y = 1
        for j in range(2 * K):
            s[j] = (s[j] + c * y) % P
            y = y * x % P
    return s

def berlekamp_massey(s, P):
    """Shortest linear recurrence of s over F_P: C = 1 + c_1 z + ... + c_L z^L."""
    C, B = [1], [1]
    L, m, b = 0, 1, 1
    for n in range(len(s)):
        d = s[n]
        for i in range(1, L + 1):
            d = (d + C[i] * s[n - i]) % P
        if d == 0:
            m += 1
            continue
        coef = d * pow(b, P - 2, P) % P
        T = C[:]
        if len(C) < len(B) + m:
            C = C + [0] * (len(B) + m - len(C))
        for i in range(len(B)):
            C[i + m] = (C[i + m] - coef * B[i]) % P
        if 2 * L <= n:
            L, B, b, m = n + 1 - L, T, d, 1
        else:
            m += 1
    return C[:L + 1] + [0] * (L + 1 - len(C)), L

def solve_mod(A, y, P):
    """Solve A c = y over F_P by Gaussian elimination; A is square and invertible."""
    n = len(A)
    M = [row[:] + [y[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = next((r for r in range(col, n) if M[r][col]), None)
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        inv = pow(M[col][col], P - 2, P)
        M[col] = [v * inv % P for v in M[col]]
        for r in range(n):
            if r != col and M[r][col]:
                f = M[r][col]
                M[r] = [(a - f * c) % P for a, c in zip(M[r], M[col])]
    return [M[r][n] for r in range(n)]

def decode(s, K, P, points):
    """Recover the multiset of points from its 2K power sums, in the three steps of Section "Preliminaries".  Returns {point: multiplicity}, or None for bottom."""
    C, L = berlekamp_massey(s, P)
    if L > K:
        return None
    lam = C[::-1]                     # Lambda(z) = z^L C(1/z), monic of degree L
    roots = []
    for x in points:
        v = 0
        for a in reversed(lam):       # Horner, highest coefficient first
            v = (v * x + a) % P
        if v == 0:
            roots.append(x)
    if len(roots) != L:
        return None
    if L == 0:
        return {}
    A = [[pow(x, j, P) for x in roots] for j in range(L)]
    c = solve_mod(A, s[:L], P)
    if c is None:
        return None
    return dict(zip(roots, c))

def check_multiset(mult, s_int, K, P):
    """The tests of Trace: at most K points, multiplicities in [1, K], and the slots equal to the power sums as integers."""
    if mult is None or len(mult) > K or any(not 1 <= c <= K for c in mult.values()):
        return False
    for j in range(2 * K):
        if s_int[j] != sum(c * pow(x, j, P) for x, c in mult.items()):
            return False
    return True

# ---------------------------------------------------------------------------
# Part 1: decoding.

def part1(rng, P, out):
    out("Part 1: decoding of power sums over F_P, P = %d (first prime above 2^20)" % P)
    universe = 1 << 20                # identifiers 1..2^20, point of identifier i is i
    for K in (59, 32):
        issued = rng.sample(range(1, universe + 1), 4000)
        ok = trials = 0
        t0 = time.time()
        for _ in range(300):
            k = rng.randint(0, K)
            pts = rng.sample(issued, k)
            mult = {x: rng.randint(1, K) for x in pts}
            s_int = [sum(c * pow(x, j, P) for x, c in mult.items()) for j in range(2 * K)]
            got = decode([v % P for v in s_int], K, P, issued)
            trials += 1
            ok += (got == mult and check_multiset(got, s_int, K, P)) or (k == 0 and got == {})
        out("  K = %d, 4000 issued identifiers, at most K points with multiplicities in [1, K]: %d of %d decoded exactly (%.1f s)" % (K, ok, trials, time.time() - t0))
        assert ok == trials
    K = 59
    t0 = time.time()
    pts = rng.sample(range(1, universe + 1), K)
    mult = {x: rng.randint(1, K) for x in pts}
    s_int = [sum(c * pow(x, j, P) for x, c in mult.items()) for j in range(2 * K)]
    got = decode([v % P for v in s_int], K, P, range(1, universe + 1))
    good = got == mult and check_multiset(got, s_int, K, P)
    out("  K = 59, all 2^20 identifiers issued, 59 points: decoded exactly: %s (%.1f s)" % (good, time.time() - t0))
    assert good
    for K in (59, 32):
        issued = rng.sample(range(1, universe + 1), 4000)
        outcome = Counter()
        for _ in range(100):
            k = rng.randint(K + 1, 2 * K)
            pts = rng.sample(issued, k)
            mult = {x: rng.randint(1, K) for x in pts}
            s_int = [sum(c * pow(x, j, P) for x, c in mult.items()) for j in range(2 * K)]
            got = decode([v % P for v in s_int], K, P, issued)
            if got is None or not check_multiset(got, s_int, K, P):
                outcome["bottom"] += 1
            elif got == mult:
                outcome["exact"] += 1
            else:
                outcome["wrong"] += 1
        out("  K = %d, more than K points (K+1 to 2K): %d bottom, %d wrong multiset, %d exact, of 100" % (K, outcome["bottom"], outcome["wrong"], outcome["exact"]))
        assert outcome["exact"] == 0

# ---------------------------------------------------------------------------
# Damgard-Jurik with s = h (Section "Damgard-Jurik instantiation").

class DJ:
    def __init__(self, rng, nbits, h):
        self.rng, self.h = rng, h
        while True:
            p = nextprime(rng.getrandbits(nbits // 2) | (3 << (nbits // 2 - 2)))
            q = nextprime(rng.getrandbits(nbits // 2) | (3 << (nbits // 2 - 2)))
            N = p * q
            if p != q and N.bit_length() == nbits and math.gcd(N, (p - 1) * (q - 1)) == 1:
                break
        self.N, self.p, self.q = N, p, q
        self.Np = [N ** k for k in range(h + 2)]
        self.lam = (p - 1) * (q - 1) // math.gcd(p - 1, q - 1)
        self.lam_inv = pow(self.lam, -1, self.Np[h])
        self.inv_fact = [pow(math.factorial(k), -1, self.Np[h + 1]) for k in range(h + 1)]

    def one_plus_N_pow(self, m, s):
        """(1+N)^m mod N^{s+1} by the binomial expansion."""
        mod = self.Np[s + 1]
        acc, falling = 1, 1
        for k in range(1, s + 1):
            falling = falling * (m - k + 1) % mod
            acc = (acc + falling * self.inv_fact[k] % mod * self.Np[k]) % mod
        return acc

    def rerandomizer(self):
        """rho^{N^h} mod N^{h+1} with rho uniform in Z*_N, as h successive N-th powers."""
        while True:
            x = self.rng.randrange(1, self.N)
            if math.gcd(x, self.N) == 1:
                break
        for s in range(1, self.h + 1):
            x = powm(x, self.N, self.Np[s + 1])
        return x

    def enc(self, m):
        return self.one_plus_N_pow(m % self.Np[self.h], self.h) * self.rerandomizer() % self.Np[self.h + 1]

    def degrade(self, w):
        return powm(w, self.N, self.Np[self.h + 1]) * self.rerandomizer() % self.Np[self.h + 1]

    def merge(self, ws):
        acc = self.rerandomizer()
        for w in ws:
            acc = acc * w % self.Np[self.h + 1]
        return acc

    def dec(self, c):
        s, N, Np = self.h, self.N, self.Np
        a = powm(c, self.lam, Np[s + 1])
        i = 0
        for j in range(1, s + 1):
            t1 = ((a % Np[j + 1]) - 1) // N
            t2 = i
            for k in range(2, j + 1):
                i -= 1
                t2 = t2 * i % Np[j]
                t1 = (t1 - t2 * Np[k - 1] * self.inv_fact[k]) % Np[j]
            i = t1
        return i * self.lam_inv % Np[s]

# ---------------------------------------------------------------------------
# The power-sum variant: encoding, Tag and Trace (paragraph "Damgard-Jurik" of Section "More live identifiers").

class Variant:
    def __init__(self, dj, K, P):
        self.dj, self.K, self.P = dj, K, P
        self.B = K * P
        assert self.B ** (2 * K) <= dj.N, "(KP)^{2K} must not exceed N"
        self.issued = {}              # point -> identifier, recorded by the authority

    def m(self, x):
        return sum(pow(x, j, self.P) * self.B ** j for j in range(2 * self.K))

    def tag(self, x, b, hprime=None):
        dj = self.dj
        hprime = dj.h if hprime is None else hprime
        if b:
            self.issued[x] = x
            return dj.enc(dj.Np[dj.h - hprime] * self.m(x))
        return dj.enc(0)

    def trace(self, w):
        dj, K, P, B = self.dj, self.K, self.P, self.B
        v = dj.dec(w)
        if v == 0:
            return None, 0
        S, most = set(), 0
        points = list(self.issued)
        for t in range(dj.h):
            v, d = divmod(v, dj.N)
            if d == 0:
                continue
            sig = []
            for _ in range(2 * K):
                d, r = divmod(d, B)
                sig.append(r)
            if d:
                return None, most
            mult = decode([x % P for x in sig], K, P, points)
            if not check_multiset(mult, sig, K, P):
                return None, most
            most = max(most, len(mult))
            S |= {self.issued[x] for x in mult}
        return S, most

# ---------------------------------------------------------------------------
# Histories: tags with their contributions, and the two payment policies.

class Tag:
    """A ciphertext with its contributions: Counter of (identifier or None for a dummy, depth); expired ones are dropped."""
    def __init__(self, w, contrib):
        self.w, self.contrib = w, contrib

def shift(contrib, k, h):
    out = Counter()
    for (i, d), c in contrib.items():
        if d + k < h:
            out[(i, d + k)] += c
    return out

def counts(contrib, h):
    cnt = [0] * h
    for (_, d), c in contrib.items():
        cnt[d] += c
    return cnt

def live_set(contrib):
    S = {i for (i, _), c in contrib.items() if i is not None and c}
    return S if S else None

def part2(rng, out, policy, h, K, P, accounts, epochs, payments, warrant):
    t0 = time.time()
    dj = DJ(rng, 3072, h)
    V = Variant(dj, K, P)
    out("Part 2: Damgard-Jurik with power sums, 3072-bit N, h = %d, K = %d, %s policy: %d accounts, %d epochs, %d payments per epoch, warrant probability %.2f (setup %.0f s)" % (h, K, "both tags degraded" if policy == "both" else "capped", accounts, epochs, payments, warrant, time.time() - t0))
    universe = 1 << 20
    fresh = rng.sample(range(1, universe + 1), accounts * epochs)
    acct = [None] * accounts          # current tag of each account
    last = [None] * accounts          # epoch of the last ingestion
    traced = mismatch = capped = 0
    most_points = 0
    max_cnt = [0] * h

    def check(tag):
        nonlocal traced, mismatch, most_points
        S, most = V.trace(tag.w)
        most_points = max(most_points, most)
        expect = live_set(tag.contrib)
        traced += 1
        if S != expect:
            mismatch += 1
            out("    MISMATCH: traced %s, live set %s" % (sorted(S) if S else S, sorted(expect) if expect else expect))
        cnt = counts(tag.contrib, h)
        for t in range(h):
            max_cnt[t] = max(max_cnt[t], cnt[t])
            if policy == "both":
                assert cnt[t] <= 2 ** t, "Lemma Contributions at a depth"
            else:
                assert cnt[t] <= K

    def degrade(tag, k=1):
        w = tag.w
        for _ in range(k):
            w = dj.degrade(w)
        return Tag(w, shift(tag.contrib, k, h))

    def merge(a, b):
        return Tag(dj.merge([a.w, b.w]), a.contrib + b.contrib)

    def ingest(a, e):
        k = 1 if last[a] is None else min(e_now - last[a], h)
        base = acct[a] if acct[a] is not None else Tag(V.tag(None, 0), Counter({(None, 0): 1}))
        new = merge(degrade(base, k), e)
        last[a] = e_now
        acct[a] = new
        check(new)

    for e_now in range(epochs):
        board = []
        for a in range(accounts):
            if rng.random() < warrant:
                x = fresh.pop()
                board.append(Tag(V.tag(x, 1), Counter({(x, 0): 1})))
            else:
                board.append(Tag(V.tag(None, 0), Counter({(None, 0): 1})))
        dormant = set(rng.sample(range(accounts), accounts // 4)) if e_now else set()
        active = [a for a in range(accounts) if a not in dormant]
        for a in active:
            ingest(a, board[a])
        for _ in range(payments):
            s, r = rng.sample(active, 2)
            ws = degrade(acct[s])
            if policy == "both":
                wr = degrade(acct[r])
            else:
                d = 0
                while True:
                    shifted = shift(acct[r].contrib, d, h)
                    if all(c <= K for c in counts(shifted + ws.contrib, h)):
                        break
                    d += 1
                wr = degrade(acct[r], d) if d else acct[r]
                capped += d > 0
            acct[s], acct[r] = ws, merge(ws, wr)
            check(acct[s])
            check(acct[r])
    out("  tags traced: %d, mismatches with the live set: %d" % (traced, mismatch))
    out("  largest number of contributions at each depth (dummies and copies included): %s" % max_cnt)
    if policy == "cap":
        out("  payments in which the cap degraded the recipient's tag: %d of %d" % (capped, epochs * payments))
    out("  largest number of distinct points decoded at one depth: %d" % most_points)
    out("  time: %.0f s" % (time.time() - t0))
    assert mismatch == 0

def part3(rng, out, h, K, P):
    """A depth filled to K distinct points, and to K+1, through real Tag, Merge and Degrade."""
    t0 = time.time()
    dj = DJ(rng, 3072, h)
    V = Variant(dj, K, P)
    pts = rng.sample(range(1, (1 << 20) + 1), K + 4)
    full, extra, other = pts[:K], pts[K], pts[K + 1:]
    tags = {x: V.tag(x, 1) for x in full + [extra] + other}
    # K distinct points, some with multiplicity above 1, degraded to depth h-1; three more points at depth 0.
    ws = [tags[x] for x in full] + [tags[x] for x in rng.sample(full, 5)]
    w = dj.merge(ws)
    for _ in range(h - 1):
        w = dj.degrade(w)
    w = dj.merge([w] + [tags[x] for x in other])
    S, most = V.trace(w)
    ok1 = S == set(full) | set(other) and most == K
    out("Part 3: Damgard-Jurik with power sums, 3072-bit N, h = %d, K = %d" % (h, K))
    out("  depth %d holds %d distinct points, depth 0 holds 3: traced exactly: %s" % (h - 1, K, ok1))
    w2 = dj.merge([tags[x] for x in full + [extra]])
    S2, _ = V.trace(w2)
    out("  a depth holds %d distinct points: Trace returns bottom: %s" % (K + 1, S2 is None))
    w3 = dj.degrade(w)
    S3, _ = V.trace(w3)
    out("  the tag degraded once more, so the K points expire: traced exactly the 3 points: %s" % (S3 == set(other)))
    out("  time: %.0f s" % (time.time() - t0))
    assert ok1 and S2 is None and S3 == set(other)

def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    rng = random.Random(seed)
    def out(line):
        print(line, flush=True)
    out("powersum_check.py, seed %d, modular exponentiation by %s" % (seed, "GMP (%s)" % _lib if _gmp else "Python"))
    P = first_prime_above(1 << 20)
    part1(rng, P, out)
    part2(rng, out, "both", h=6, K=32, P=P, accounts=8, epochs=6, payments=12, warrant=0.6)
    part2(rng, out, "cap", h=10, K=59, P=P, accounts=8, epochs=4, payments=20, warrant=0.6)
    part3(rng, out, h=6, K=32, P=P)
    part3(rng, out, h=10, K=59, P=P)
    out("all checks passed")

if __name__ == "__main__":
    main()
