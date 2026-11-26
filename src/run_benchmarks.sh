#!/bin/sh
# run_benchmarks.sh -> runs the benchmarks of the paper's evaluation and
# writes the results next to src/.
# Usage (from src/src): ./run_benchmarks.sh [suite|board|all|everything],
# through make benchmarks (suite), make board-1e6 (board), make all-benchmarks
# (all) or make everything (everything), which build the binaries first.
# everything: the correctness tests (results_tests.txt), the suite, the policy
# simulation (results_policy_sim.txt), the power-sum check
# (results_powersum_check.txt), the table rows of paper_tables.py
# (results_tables.txt) and the board; about eight hours.
# suite: every operation, one after the other, single-threaded and pinned to
# one performance core (the i7-1260P has performance cores 0-7 and efficiency
# cores 8-15).  src/src/paper_tables.py then prints the paper's table rows.
# board: the ElGamal board of BOARD_M entries (default 10^6) on THREADS
# threads (default 16), unpinned, with its progress in
# results_eg_board_1e6_progress.txt; it takes about six hours.
# CPU=n selects another core.  During the suite, the frequency of that core is
# sampled every five seconds into results_cpu_frequency.txt.  The tuned profile is set to
# throughput-performance for the run and restored afterwards; PROFILE=name
# selects another.
set -e
MODE=${1:-suite}
case $MODE in suite|board|all|everything) ;; *) echo "usage: $0 [suite|board|all|everything]" >&2; exit 1;; esac
OUT=..
BOARD_M=${BOARD_M:-1000000}
THREADS=${THREADS:-16}
CPU=${CPU:-2}
PIN="taskset -c $CPU"
HDR="label,operation,runs,mean_ms,sd_ms,min_ms,bytes"
PROFILE=${PROFILE:-throughput-performance}
# Without tuned-adm (in the container, for example) the profile is left alone.
TUNED=$(command -v tuned-adm || true)
OLD_PROFILE=
if [ -n "$TUNED" ]; then
  OLD_PROFILE=$(tuned-adm active | sed 's/^Current active profile: //')
  tuned-adm profile "$PROFILE"
fi
restore() { [ -n "$TUNED" ] && tuned-adm profile "$OLD_PROFILE"; true; }
trap restore EXIT

{
  echo "date: $(date -Iseconds)"
  echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
  echo "pinned to: cpu $CPU, max $(cat /sys/devices/system/cpu/cpu$CPU/cpufreq/cpuinfo_max_freq) kHz"
  echo "governor: $(cat /sys/devices/system/cpu/cpu$CPU/cpufreq/scaling_governor)"
  echo "energy performance preference: $(cat /sys/devices/system/cpu/cpu$CPU/cpufreq/energy_performance_preference 2>/dev/null || echo unknown)"
  echo "tuned profile: $(tuned-adm active 2>/dev/null | sed 's/^Current active profile: //' || true)"
  echo "power: AC online $(cat /sys/class/power_supply/AC/online 2>/dev/null || echo unknown)"
  echo "kernel: $(uname -sr)"
  echo "compiler: $(gcc --version | head -1)"
  echo "flags: -std=c11 -Wall -Wextra -Wpedantic -O2"
  echo "gmp: $(rpm -q gmp 2>/dev/null || echo unknown)"
  echo "openssl: $(openssl version)"
  echo "relic: $(grep -m1 'RLC_VERSION' /usr/local/include/relic/relic_conf.h | awk '{print $3}')"
  echo "runs: one process at a time; suite single-threaded and pinned, board on $THREADS threads unpinned"
  echo "mode: $MODE"
} > $OUT/results_environment.txt

F=/sys/devices/system/cpu/cpu$CPU/cpufreq/scaling_cur_freq
SAMPLER=
trap '[ -n "$SAMPLER" ] && kill $SAMPLER; restore' EXIT

suite() {
( while :; do echo "$(date +%s) $(cat $F)"; sleep 5; done ) > $OUT/results_cpu_frequency.txt &
SAMPLER=$!

# Exponential ElGamal: operations (Merge rerandomizes) and board proof.
$PIN ./micro_bench 100 3 28 > $OUT/results_micro_bench_k28.csv
$PIN ./bench_trace 2 > $OUT/results_bench_trace_k28.csv 2>&1
{ echo "$HDR"; $PIN ./bench_eg_board 10 100 1; } > $OUT/results_eg_board.csv

# Ring-LWE parameter sets (no timing): h = 16, h = 10, and estimates at larger h.
python3 rlwe_params.py 16 > $OUT/results_rlwe_params.csv
python3 rlwe_params.py 10 > $OUT/results_rlwe_params_h10.csv
python3 rlwe_params.py 20 50 100 > $OUT/results_rlwe_params_large_h.csv

# Ring-LWE noise test: L = 32 transitions of fan-in d at each parameter set.
{
  $PIN ./test_rlwe 2048 16 36 1 1 32
  $PIN ./test_rlwe 4096 16 65 2 2 32
  $PIN ./test_rlwe 4096 16 97 2 4 32
  $PIN ./test_rlwe 2048 10 36 1 1 32
  $PIN ./test_rlwe 4096 10 65 2 2 32
  $PIN ./test_rlwe 4096 10 97 2 4 32
} > $OUT/results_tests_rlwe.txt 2>&1

# Ring-LWE at the parameter sets of results_rlwe_params.csv (h = 16) and
# results_rlwe_params_h10.csv (h = 10).
{
  echo "$HDR"
  $PIN ./bench_rlwe linear_chain 2048 16 36 1 20
  $PIN ./bench_rlwe two_state 4096 16 65 2 20
  $PIN ./bench_rlwe merge_heavy 4096 16 97 2 20
  $PIN ./bench_rlwe linear_chain 2048 10 36 1 20
  $PIN ./bench_rlwe two_state 4096 10 65 2 20
  $PIN ./bench_rlwe merge_heavy 4096 10 97 2 20
} > $OUT/results_rlwe_bench.csv

# Damgard-Jurik at h = 10 with the board proof, then at larger h.
{ echo "$HDR"; $PIN ./bench_dj 10 5 3 1 1; } > $OUT/results_dj_h10.csv
{
  echo "$HDR"
  $PIN ./bench_dj 20 1 1 0 1
  $PIN ./bench_dj 50 1 1 0 1
  $PIN ./bench_dj 100 1 1 0 1
} > $OUT/results_dj_large_h.csv

kill $SAMPLER
SAMPLER=
}

# The ElGamal board at full size: one entry timed three times on one thread,
# then the board proved and verified in chunks of 4096 entries.
board() {
  { echo "$HDR"; ./bench_eg_board 3 "$BOARD_M" "$THREADS" 4096 8 2> $OUT/results_eg_board_1e6_progress.txt; } \
    > $OUT/results_eg_board_1e6.csv
}

# The correctness tests of make tests; a failure stops the run.
tests() {
  make tests > $OUT/results_tests.txt 2>&1
}

case $MODE in
  suite) suite ;;
  board) board ;;
  all)   suite; board ;;
  everything)
    tests
    suite
    python3 policy_sim.py > $OUT/results_policy_sim.txt
    python3 powersum_check.py > $OUT/results_powersum_check.txt
    python3 paper_tables.py > $OUT/results_tables.txt
    board ;;
esac
