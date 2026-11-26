#!/usr/bin/env python3
"""paper_tables.py -> prints, from the results files in .., the LaTeX rows of
the evaluation tables of the paper and the derived numbers its text uses.

Usage (from src/src, after run_benchmarks.sh): python3 paper_tables.py

Every number printed here is read from a results file or computed from one by
the formula printed next to it; the paper's tables are these rows.
"""
import csv
import math
import os

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def rows(name, skip_to=None):
    with open(os.path.join(R, name)) as f:
        lines = f.read().splitlines()
    if skip_to:
        lines = lines[next(i for i, l in enumerate(lines) if l.startswith(skip_to)):]
    return list(csv.DictReader(lines))


def ms(x):
    """milliseconds, three significant digits, thin space in thousands"""
    x = float(x)
    if x < 1:
        s = f"{x:.2f}"
    elif x < 100:
        s = f"{x:.1f}"
    else:
        s = f"{x:.0f}"
    if "." not in s and len(s) > 4:
        s = f"{int(s):,}".replace(",", "\\,")
    return s


def kib(b):
    v = int(b) / 1024
    return f"{v:.3f}".rstrip("0").rstrip(".") if v < 100 else f"{v:.0f}"


def gib(b):
    return f"{b / 2**30:.2f}"


def pick(rs, label_prefix, op):
    for r in rs:
        if r.get("label", "").startswith(label_prefix) and r["operation"] == op:
            return r
    return None


def main():
    eg = {r["operation"]: r for r in rows("results_micro_bench_k28.csv")}
    egt = rows("results_bench_trace_k28.csv", skip_to="mult_max")
    egb = rows("results_eg_board.csv")
    dj = rows("results_dj_h10.csv")
    djl = rows("results_dj_large_h.csv") if os.path.exists(os.path.join(R, "results_dj_large_h.csv")) else []
    rl = rows("results_rlwe_bench.csv")
    rp = {(r["profile"], r["h"]): r for f in ("results_rlwe_params_h10.csv", "results_rlwe_params.csv",
                                              "results_rlwe_params_large_h.csv")
          for r in rows(f)}

    djl10 = "dj_h10_"
    prof = [("linear_chain", "2048", "36"), ("two_state", "4096", "65"), ("merge_heavy", "4096", "97")]

    def rlr(p, h, op):
        n0, lq = rp[(p, str(h))]["N0"], rp[(p, str(h))]["log2Q"]
        return pick(rl, f"rlwe_{p}_N0{n0}_h{h}_logQ{lq}", op)

    def djr(op):
        return pick(dj, djl10, op)

    print("% tab:costs (h = 10), mean ms")
    lines = [
        ("\\genwatermark{} (tag and signature)", eg["tag"], djr("tag"), [rlr(p, 10, "tag") for p, _, _ in prof]),
        ("$\\mathsf{Verify}$ (signature)", eg["verify_signature"], djr("verify_signature"),
         [rlr(p, 10, "verify_signature") for p, _, _ in prof]),
        ("\\degrade", eg["degrade"], djr("degrade"), [rlr(p, 10, "degrade") for p, _, _ in prof]),
        ("\\mergewatermarks, $k=2$", eg["merge_k2"], djr("merge_k2"), [rlr(p, 10, "merge_k2") for p, _, _ in prof]),
    ]
    for name, a, b, cs in lines:
        print(f"{name} & {ms(a['mean_ms'])} & {ms(b['mean_ms'])} & " + " & ".join(ms(c["mean_ms"]) for c in cs) + " \\\\")
    for name, o in (("\\degrade, online", "degrade_online"),
                    ("\\mergewatermarks, $k=2$, online", "merge_k2_online"),
                    ("Rerandomizer, precomputed offline", "rerandomizer")):
        print(f"{name} & --- & {ms(djr(o)['mean_ms'])} & --- & --- & --- \\\\")
    tn = [r for r in dj if r["operation"].startswith("trace_K") and r["operation"] != "trace_K1"][0]
    rt = [[r for r in rl if r["label"].startswith(f"rlwe_{p}_N0{rp[(p, '10')]['N0']}_h10_")
           and r["operation"].startswith("trace_K") and r["operation"] != "trace_K1"][0] for p, _, _ in prof]
    print(f"\\trace, every slot live & {ms(eg['trace_k28']['mean_ms'])} & {ms(tn['mean_ms'])} & "
          + " & ".join(ms(r["mean_ms"]) for r in rt) + " \\\\")
    print(f"Live identifiers $K$ & 28 & {tn['operation'][7:]} & " + " & ".join(r["operation"][7:] for r in rt) + " \\\\")
    print(f"Tag (KiB) & {kib(eg['degrade']['bytes'])} & {kib(djr('degrade')['bytes'])} & "
          + " & ".join(kib(rlr(p, 10, "degrade")["bytes"]) for p, _, _ in prof) + " \\\\")
    sig = [int(eg["tag"]["bytes"]) - int(eg["degrade"]["bytes"]),
           int(djr("tag")["bytes"]) - int(djr("degrade")["bytes"])] + \
          [int(rlr(p, 10, "tag")["bytes"]) - int(rlr(p, 10, "degrade")["bytes"]) for p, _, _ in prof]
    print("Signature (bytes) & " + " & ".join(str(s) for s in sig) + " \\\\")
    print(f"% runs: ElGamal {eg['tag']['runs']} (trace {eg['trace_k28']['runs']}), "
          f"DJ {djr('tag')['runs']} (trace {tn['runs']}), Ring-LWE {rlr('two_state', 10, 'tag')['runs']}")
    print(f"% ElGamal degradation proof: prove {ms(eg['degrade_and_prove']['mean_ms'])} ms, "
          f"verify {ms(eg['verify_degradation_proof']['mean_ms'])} ms, "
          f"{int(eg['degrade_and_prove']['bytes']) - int(eg['degrade']['bytes'])} bytes; merge k=10 {ms(eg['merge_k10']['mean_ms'])} ms")
    print("% ElGamal trace sweep (min of runs): " + "; ".join(
        f"M<={r['mult_max']} K={r['K']} {r['trace_ms']} ms" for r in egt))
    print(f"% trace one identifier: DJ {ms(djr('trace_K1')['mean_ms'])} ms, Ring-LWE "
          + ", ".join(ms(rlr(p, 10, 'trace_K1')['mean_ms']) for p, _, _ in prof) + " ms")

    print("\n% tab:ledger (h = 10, 10^6 people, Merkle path of 20 hashes of 32 bytes)")
    cols = [("eg", None), ("dj", None)] + [("rl", p) for p, _, _ in prof]

    def op(c, o):
        k, p = c
        return eg[o] if k == "eg" else djr(o) if k == "dj" else rlr(p, 10, o)
    print("\\multicolumn{6}{l}{\\emph{Per payment}}\\\\")
    print("Ciphertext $\\hat w_S$, pruned (KiB) & " + " & ".join(kib(op(c, "degrade")["bytes"]) for c in cols) + " \\\\")
    print("\\degrade{} and \\mergewatermarks{} (ms) & " + " & ".join(
        ms(float(op(c, "degrade")["mean_ms"]) + float(op(c, "merge_k2")["mean_ms"])) for c in cols) + " \\\\")
    print("\\quad online, after precomputation (ms) & --- & " + ms(float(djr("degrade_online")["mean_ms"]) + float(djr("merge_k2_online")["mean_ms"])) + " & --- & --- & --- \\\\")
    print("% DJ offline per payment: two rerandomizers, " + ms(2 * float(djr("rerandomizer")["mean_ms"])) + " ms")
    print("\\multicolumn{6}{l}{\\emph{Per person and epoch}}\\\\")
    print("Entry, signature, path (KiB) & " + " & ".join(
        f"{int(op(c, 'tag')['bytes']) + 640:d}" if False else kib(int(op(c, "tag")["bytes"]) + 20 * 32) for c in cols) + " \\\\")
    egb = {r["operation"]: r for r in rows("results_eg_board_1e6.csv")}
    bd = [egb["board_entry"], djr("board_entry")]
    bv = [egb["verify_board_entry"], djr("verify_board_entry")]
    eb = [int(eg["degrade"]["bytes"]), int(djr("degrade")["bytes"])]  # the entry itself
    print("Flag ciphertext and entry proof (KiB) & " + " & ".join(
        f"{(int(r['bytes']) - e) / 1024:.1f}" for r, e in zip(bd, eb)) + " & --- & --- & --- \\\\")
    print("Verifying the entry proof (s) & " + " & ".join(
        f"{float(r['mean_ms']) / 1000:.2f}" if float(r["mean_ms"]) < 1e4 else f"{float(r['mean_ms']) / 1000:.0f}" for r in bv) + " & --- & --- & --- \\\\")
    print("\\multicolumn{6}{l}{\\emph{Per epoch}}\\\\")
    print("Entries in the public store (GiB) & " + " & ".join(
        f"{v:.2f}" if v < 10 else f"{v:.0f}" for v in (int(op(c, "tag")["bytes"]) * 10**6 / 2**30 for c in cols)) + " \\\\")
    print("Issuing the entries (core-hours) & " + " & ".join(
        f"{v:.1f}" if v < 10 else f"{v:.0f}" for v in (float(op(c, "tag")["mean_ms"]) * 1e6 / 3.6e6 for c in cols)) + " \\\\")
    print("Entry proofs in the public store (GiB) & " + " & ".join(
        f"{(int(r['bytes']) - e) * 10**6 / 2**30:.1f}" for r, e in zip(bd, eb)) + " & --- & --- & --- \\\\")
    print("Proving the board (core-hours) & " + " & ".join(
        ms(float(r["mean_ms"]) * 1e6 / 3.6e6) for r in bd) + " & --- & --- & --- \\\\")
    print("% store: 10^6 * tag bytes (with the signature) / 2^30; hours: 10^6 * tag ms / 3.6e6")
    print("% board rows: bytes of board_entry minus the entry; one core, one entry, times 10^6")
    for r in rows("results_eg_board_1e6.csv"):
        if r["operation"].endswith("_t16"):
            print(f"% ElGamal {r['operation']}: {float(r['mean_ms']) / 3.6e6:.2f} h = {float(r['mean_ms']) / 6e4:.1f} min")

    print("\n% tab:deploy: h, slots, tag (KiB), degrade (ms), trace every slot (ms)")
    print(f"Exp.\\ ElGamal, $\\mu=1$ & 10 & 28 & {kib(eg['degrade']['bytes'])} & {ms(eg['degrade']['mean_ms'])} & {ms(eg['trace_k28']['mean_ms'])} \\\\")
    for r in [djr("degrade")] + [r for r in djl if r["operation"] == "degrade"]:
        lab = r["label"]
        h = int(lab.split("_")[1][1:])
        n = lab.split("_n")[-1]
        src = dj if lab.startswith(djl10) else djl
        tr = next(x for x in src if x["label"] == lab and x["operation"] == f"trace_K{n}")
        print(f"Damg{{\\aa}}rd--Jurik & {h} & {n} & {kib(r['bytes'])} & {ms(r['mean_ms'])} & {ms(tr['mean_ms'])} \\\\")
    for r in [djr("degrade_online")] + [r for r in djl if r["operation"] == "degrade_online"]:
        lab = r["label"]; src = dj if lab.startswith(djl10) else djl
        rr = next(x for x in src if x["label"] == lab and x["operation"] == "rerandomizer")
        mo = next(x for x in src if x["label"] == lab and x["operation"] == "merge_k2_online")
        print(f"% DJ {lab}: degrade online {ms(r['mean_ms'])} ms, merge k=2 online {ms(mo['mean_ms'])} ms, rerandomizer {ms(rr['mean_ms'])} ms")
    for h in (10, 16):
        for p, _, _ in prof:
            q = rp[(p, str(h))]
            d = rlr(p, h, "degrade")
            tr = [r for r in rl if r["label"] == d["label"] and r["operation"].startswith("trace_K") and r["operation"] != "trace_K1"][0]
            print(f"Ring-LWE, $d={q['d']}$ & {h} & {q['N0']} & {kib(d['bytes'])} & {ms(d['mean_ms'])} & {ms(tr['mean_ms'])} \\\\")
    for h in (20, 50, 100):
        q = rp.get(("two_state", str(h)))
        if q:
            print(f"% Ring-LWE d=2 h={h}: N0 {q['N0']}, tag {kib(q['tag_bytes'])} KiB (parameters only, not measured, not in the table)")
    for h in (20, 50, 100):
        for p, _, _ in prof:
            q = rp.get((p, str(h)))
            if q:
                print(f"% Ring-LWE estimate {p} h={h}: N0 {q['N0']}, log2Q {q['log2Q']}, log2 t {int(q['t']).bit_length()}, "
                      f"tag {kib(q['tag_bytes'])} KiB ({q['tag_MiB']} MiB), not measured")


def freq():
    with open(os.path.join(R, "results_cpu_frequency.txt")) as f:
        k = sorted(int(l.split()[1]) for l in f if l.strip())
    m = k[len(k) // 2] if len(k) % 2 else (k[len(k) // 2 - 1] + k[len(k) // 2]) / 2
    print(f"\n% paragraph Measurements: core frequency median {m / 1e6:.1f} GHz, "
          f"min {k[0] / 1e6:.1f} GHz, max {k[-1] / 1e6:.1f} GHz over {len(k)} samples")


if __name__ == "__main__":
    main()
    freq()
