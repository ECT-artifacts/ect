# Track me if you can: Ephemeral coin tracing

Prototype implementation of **Ephemeral Coin Tracing (ECT)**

The prototype implements the three instantiations of the paper: exponential ElGamal (`ect.c`, on RELIC), Damgård–Jurik (`dj.c`, on GMP) and Ring-LWE (`rlwe.c`, BGV over R_0[Y]/(Y^h) with its own negacyclic NTT), together with the board proofs of the public count for exponential ElGamal (`eg_board.c`) and Damgård–Jurik (`dj_board.c`).

The scheme operations and ledger checks are deliberately separated. `ect_merge`, `dj_merge` and `rlwe_merge` are the public ECT `Merge` algorithm: they combine the ciphertexts homomorphically and rerandomize the result with a fresh encryption of zero. Issuance-signature checks and the prototype degradation-proof checks model the host ledger's validity predicate; they are not performed by `Merge`.

Issuance signatures: Schnorr in the ElGamal group for exponential ElGamal, Ed25519 for Damgård–Jurik and ML-DSA-44 for Ring-LWE (both through OpenSSL), so that Ring-LWE issuance is post-quantum as well.



## Dependencies

-  [RELIC](https://github.com/relic-toolkit/relic) | v. 0.7.0 | Big integer arithmetic, modular exponentiation, SHA-256 
- [GMP](https://gmplib.org/) | v. 6.3.0 | Backend for RELIC multi-precision arithmetic 
- [SageMath](https://www.sagemath.org/) | v. 10.8  | Group parameter generation
- [OpenSSL](https://www.openssl.org/) | v. 3.5 | Randomness (RAND_bytes), SHA-256 / SHA-512 for the board proofs, Ed25519 and ML-DSA-44 signatures (ML-DSA needs OpenSSL 3.5 or later)

### System requirements

- GCC with C11 support (`-std=c11`)
- Linux (tested on Ubuntu 24.04 and Fedora 44)

---

## Building with Docker

The `Dockerfile` builds the prototype on Ubuntu 24.04, with OpenSSL 3.5.9 and RELIC 0.7.0 built from pinned commits (about ten minutes):

```bash
docker build -t ect .
docker run --rm ect make tests                            # the correctness tests
docker run --rm -v "$PWD:/ect" ect make benchmarks        # writes the results files into this repository
docker run --rm -v "$PWD:/ect" ect make everything
```

In a container, `run_benchmarks.sh` leaves the tuned profile alone; the core pinning and the frequency sampling work as on the host.

## Building RELIC

RELIC must be compiled with the following flags to support the 3072-bit group arithmetic used in this scheme:

```bash
git clone https://github.com/relic-toolkit/relic.git
cd relic
mkdir build && cd build

cmake \
  -DALLOC=DYNAMIC \
  -DBN_PRECI=8192 \
  ..

make -j$(nproc)
sudo make install
```

## Building the prototype

```bash
cd src
make
```

**NOTE**

If RELIC is not installed system-wide:

```bash
make RELIC_INC=/path/to/relic/include RELIC_LIB=/path/to/relic/lib
```
---

## Running

```bash
./ect_demo
```

For the focused measurements reported in the paper, build and run the
microbenchmark instead of the full simulation:

```bash
cd src
make micro-bench
```

The optional arguments are `./micro_bench RUNS TRACE_RUNS K` (defaults:
`10 3 N`, currently `10 3 28`). Output is CSV and each operation is timed in isolation after an
untimed warm-up.

## Damgård–Jurik, Ring-LWE and the board proofs

```bash
cd src
make tests        # test_eg_board, test_rlwe on the six parameter sets, test_dj at h = 10 (several minutes)
make benchmarks   # builds the binaries and runs run_benchmarks.sh (about 75 minutes; the Damgard-Jurik board proof and Damgard-Jurik at h = 100 dominate)
python3 paper_tables.py   # prints the rows of the paper's evaluation tables from the results files
make everything   # the tests, the benchmarks, the policy simulation, the power-sum check, the table rows (results_tables.txt) and the 10^6-entry board, one after the other (about eight hours)
```

From the paper repository, `make bench` and `make bench-tests` call these two targets.

### Rerunning the measurements of the paper

Run `make benchmarks` on an idle machine on AC power, with the performance governor. `run_benchmarks.sh` pins every benchmark to one core (`CPU=2` by default; `CPU=n make benchmarks` selects another, preferably a performance core), runs them one after the other, single-threaded, and samples the frequency of that core every five seconds. It regenerates every results file below. Then `python3 paper_tables.py` prints the table rows of Section "Evaluation" and the numbers its text quotes, each with the results file and the formula it comes from; the paper's numbers are replaced by these.

| File | Content | Paper |
|---|---|---|
| `results_environment.txt` | date, CPU, pinned core and its maximum frequency, governor, power, kernel, compiler, flags, library versions | Section "Implementation and measurements", paragraph "Measurements" |
| `results_cpu_frequency.txt` | frequency of the pinned core every five seconds (Unix time, kHz) | same paragraph (median, minimum, maximum) |
| `results_micro_bench_k28.csv` | exponential ElGamal operations, h = 10 (`./micro_bench 100 3 28`) | Tables "costs", "ledger" and "deploy" (ElGamal column and row); the degradation proof in paragraph "Exponential ElGamal" |
| `results_bench_trace_k28.csv` | ElGamal trace time against K for multiplicity bounds 1, 2 and 4, fastest of two runs (`./bench_trace 2`) | Section "Cost of the operations", last paragraph |
| `results_eg_board.csv` | ElGamal board proof per entry, count proof and a board of 100 entries (`./bench_eg_board 10 100 1`) | not used: the paper takes the per-entry costs from the board of 10^6 entries |
| `results_dj_h10.csv` | Damgård–Jurik at h = 10, with the rerandomizer and the online parts of degrade and merge timed separately, and the board proof for one entry (`./bench_dj 10 5 3 1 1`) | Tables "costs", "ledger" and "deploy", and the online costs in Sections "Cost of the operations", "Cost on the ledger" and "Parameters for a deployment" |
| `results_dj_large_h.csv` | Damgård–Jurik at h = 20, 50 and 100, one run each without warm-up, with the same split (`./bench_dj h 1 1 0 1`) | Table "deploy" and the trade-off paragraph of Section "Parameters for a deployment" |
| `results_rlwe_bench.csv` | Ring-LWE at the parameter sets of `results_rlwe_params.csv` (h = 16) and `results_rlwe_params_h10.csv` (h = 10) (`./bench_rlwe PROFILE N0 h log2Q PRIMES RUNS`) | Tables "costs", "ledger" and "deploy" |
| `results_tests_rlwe.txt` | `test_rlwe` at the six parameter sets: 32 transitions of fan-in d and the smallest noise margin | paragraph "Ring-LWE" of Section "Implementation and measurements" |
| `results_rlwe_params.csv` | Ring-LWE parameter sets at h = 16 (`python3 rlwe_params.py 16`), with the modulus margin, largest fan-in and largest horizon of each | Table "rlwe-params" and its caption |
| `results_rlwe_params_h10.csv` | the same at h = 10 (`python3 rlwe_params.py 10`) | Tables "instantiations" (tag size at h = 10) and "deploy" |
| `results_rlwe_params_large_h.csv` | estimates at h = 20, 50 and 100 with the smallest admissible t (`python3 rlwe_params.py 20 50 100`), not measured | Table "deploy" and the trade-off paragraph |
| `results_tests.txt` | output of `make tests` (`make everything` only) | not used |
| `results_tables.txt` | output of `python3 paper_tables.py` (`make everything` only) | the rows of the evaluation tables and the numbers the text quotes |
| `results_eg_board_1e6.csv` | ElGamal board of 10^6 entries on 16 threads (`make board-1e6`, about six hours; progress in `results_eg_board_1e6_progress.txt`) | Table "ledger" (the ElGamal board rows) and the times of the board in Section "Cost of a payment and of an epoch" |
| `results_policy_sim.txt` | simulation of the three payment policies (sender-side, capped, both tags degraded) on 2000 accounts over 25 epochs: reach of a traced entry, contributions per depth, Ring-LWE noise terms (`python3 policy_sim.py > ../results_policy_sim.txt`, about three minutes, seeded; part of `make everything`, not of `make benchmarks`) | Section "More live identifiers", paragraph "Parameters and cost" |
| `results_powersum_check.txt` | check of the power-sum encoding: decoding of power sums over F_P with P the first prime above 2^20, at K = 59 and 32 and once with all 2^20 identifiers issued; the Damgård–Jurik variant end to end at 3072 bits on random histories, with both tags degraded at h = 6, K = 32 and with the capped policy at h = 10, K = 59, every tag traced against its live set; a depth filled to K points and to K + 1 (`python3 powersum_check.py > ../results_powersum_check.txt` or `make powersum-check`, about fifteen minutes, seeded; Python, with GMP through ctypes when libgmp is installed; part of `make everything`) | Section "More live identifiers", paragraph "Parameters and cost" |

The CSV columns of the benchmarks are `label,operation,runs,mean_ms,sd_ms,min_ms,bytes`, where `bytes` is the size of the output of the operation (for a tag, the ciphertext and, for `tag`, its signature; for the board, the entries, signatures, bit ciphertexts and proofs). Ring-LWE sizes are 2 N_0 h ceil(log2 Q) bits per tag.

The run of 24 November 2026 at 16:47 behind the current numbers was on a 12th Gen Intel Core i7-1260P on AC power, Fedora 44 (Linux 7.3) with the performance governor and the throughput-performance tuned profile, gcc 16.2.1 with `-std=c11 -O2`, GMP 6.3.0, OpenSSL 3.5.9 and RELIC 0.7.0, pinned to core 2; the frequency of that core had a median of 3.5 GHz (3.2 to 4.0 GHz, `results_cpu_frequency.txt`) against a maximum of 4.7 GHz.

Implementation notes:
- Damgård–Jurik samples the rerandomizer rho in Z*_N, which gives the same distribution of rho^{N^h} mod N^{h+1} as sampling it in Z*_{N^{h+1}}, and computes rho^{N^h} as h successive N-th powers on the moduli N^2, ..., N^{h+1}.
- The Damgård–Jurik rerandomizer does not depend on the tag, so it can be precomputed: `dj_rerandomizer` computes it, and `dj_degrade_pre` and `dj_merge_pre` are the online parts of degrade and merge: w^N t and t prod_j w_j for a precomputed rerandomizer t. `dj_degrade` and `dj_merge` are `dj_rerandomizer` followed by these.
- Ring-LWE uses a ternary secret and a discrete Gaussian error of standard deviation 3.2 truncated at 20, sampled from a table (not constant time), and t = 65537.
- The board proofs are sent as their challenges and responses; the verifier recomputes the first messages and checks the hash (SHA-512 reduced mod q for ElGamal, SHA-256 with 256-bit challenges for Damgård–Jurik). Damgård–Jurik responses are sent mod N and must be units.


## Known build issues

### "cannot find -lgmp"

GMP development headers are required at link time. Install them before building the prototype:

#### Ubuntu/Debian
```bash
sudo apt install libgmp-dev
```

#### Fedora/RHEL
```bash
sudo dnf install gmp-devel
```
### "error while loading shared libraries: librelic.so"

After installing RELIC, the dynamic linker cache may need to be updated, especially if RELIC installed to /usr/local/lib:

```bash
echo "/usr/local/lib" | sudo tee /etc/ld.so.conf.d/relic.conf
sudo ldconfig
```
Alternatively, set LD_LIBRARY_PATH for the current session:

```bash
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```
