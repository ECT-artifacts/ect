#!/usr/bin/env python3
"""Compare three payment policies by the reach of a trace.

sender: the sender's tag is degraded, the recipient's is not (the policy of the paper).
cap:    as sender, and the recipient's tag is degraded the fewest times for which every depth has at most K contributions.
both:   the sender's tag and the recipient's tag are both degraded on every payment.

Measured per traced board entry: the number of accounts whose tag holds it before it expires, and the longest chain of payments along which it travels.  Measured per tag: the contributions at each depth (dummies included, copies counted), and the number of rerandomisation noise terms in each block of a Ring-LWE tag.
"""
import random
import statistics
import sys


def degrade(tag):
    return [dict()] + tag[:-1]


def degrade_noise(nz):
    return [1] + [x + 1 for x in nz[:-1]]


def counts(tag):
    return [sum(d.values()) for d in tag]


def merge(a, b):
    out = []
    for da, db in zip(a, b):
        d = dict(da)
        for k, v in db.items():
            d[k] = d.get(k, 0) + v
        out.append(d)
    return out


def simulate(policy, K=None, h=10, n=2000, epochs=25, rate=3, f=0.005, zipf=1.0, seed=1):
    rng = random.Random(seed)
    tags = [[dict() for _ in range(h)] for _ in range(n)]
    noise = [[0] * h for _ in range(n)]
    path = [dict() for _ in range(n)]
    weights = [1 / (i + 1) ** zipf for i in range(n)]
    next_id = 0
    reach, born, longest = {}, {}, {}
    forced = []
    depth_max = [0] * h
    noise_max = 0
    distinct_max = 0
    for epoch in range(epochs):
        for a in range(n):
            entry = (next_id, rng.random() < f)
            next_id += 1
            t = degrade(tags[a])
            t[0][entry] = 1
            tags[a] = t
            noise[a] = [x + 2 for x in degrade_noise(noise[a])]  # the entry is a fresh encryption, and the merge adds another
            if entry[1]:
                reach[entry], born[entry], longest[entry] = {a}, epoch, 0
                path[a][entry] = 0
        for payee in rng.choices(range(n), weights, k=rate * n):
            payer = rng.randrange(n)
            if payer == payee:
                continue
            s = degrade(tags[payer])
            tags[payer] = s
            ns = degrade_noise(noise[payer])
            noise[payer] = ns
            r, nr = tags[payee], noise[payee]
            if policy == "both":
                r, nr = degrade(r), degrade_noise(nr)
            elif policy == "cap":
                cs, j = counts(s), 0
                while any(x + y > K for x, y in zip(cs, counts(r))):
                    r, nr = degrade(r), degrade_noise(nr)
                    j += 1
                forced.append(j)
            tags[payee] = merge(s, r)
            noise[payee] = [x + y + 1 for x, y in zip(ns, nr)]  # the merge adds a fresh encryption of zero
            for d in s:
                for entry in d:
                    if entry[1]:
                        reach[entry].add(payee)
                        hops = path[payer].get(entry, 0) + 1
                        if hops > path[payee].get(entry, -1):
                            path[payee][entry] = hops
                        if hops > longest[entry]:
                            longest[entry] = hops
        for a in range(n):
            for depth, d in enumerate(tags[a]):
                c = sum(d.values())
                if c > depth_max[depth]:
                    depth_max[depth] = c
                distinct_max = max(distinct_max, sum(1 for e in d if e[1]))
            noise_max = max(noise_max, max(noise[a]))
    old = [e for e in reach if born[e] < epochs - h]
    return {
        "reach": statistics.mean(len(reach[e]) for e in old),
        "chain": statistics.mean(longest[e] for e in old),
        "depth_max": depth_max,
        "noise_max": noise_max,
        "distinct_max": distinct_max,
        "forced": statistics.mean(forced) if forced else 0,
        "forced_any": sum(1 for x in forced if x) / len(forced) if forced else 0,
    }


def line(label, r, h):
    full = max(r["depth_max"])
    print(f"  {label:26s} reach {r['reach']:6.1f} accounts, chain {r['chain']:4.1f} payments, fullest depth {full:6d} contributions, {r['distinct_max']:3d} distinct traced, noise terms {r['noise_max']:7d}")
    sys.stdout.flush()


def policies(rate):
    print(f"2000 accounts, {rate} payments per account per epoch, 0.5% of entries traced, 25 epochs")
    line("sender, h=10", simulate("sender", rate=rate), 10)
    for K in (59, 1024):
        r = simulate("cap", K, rate=rate)
        line(f"cap K={K}, h=10", r, 10)
        print(f"      recipient degraded in {r['forced_any']:.1%} of payments")
    for h in (6, 7, 8, 10, 12, 14, 16):
        r = simulate("both", h=h, rate=rate)
        line(f"both, h={h}", r, h)
        ok = all(c <= 2 ** d for d, c in enumerate(r["depth_max"]))
        print(f"      depth t has at most 2^t contributions: {ok}; noise terms at most 3(2^h - 1): {r['noise_max'] <= 3 * (2 ** h - 1)}")
    for h in (6, 7, 8):
        line(f"sender, h={h}", simulate("sender", h=h, rate=rate), h)


if __name__ == "__main__":
    policies(3)
    policies(1)
