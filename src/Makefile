CC       := gcc
CFLAGS   := -std=c11 -Wall -Wextra -Wpedantic
LDFLAGS  :=
LIBS     := -lrelic -lgmp -lcrypto -lm

# Override if RELIC is not on the default search path:
#   make RELIC_INC=/opt/relic/include RELIC_LIB=/opt/relic/lib
RELIC_INC ?=
RELIC_LIB ?=

ifneq ($(RELIC_INC),)
    CFLAGS += -I$(RELIC_INC)
endif
ifneq ($(RELIC_LIB),)
    LDFLAGS += -L$(RELIC_LIB) -Wl,-rpath,$(RELIC_LIB)
endif

BUILD ?= release
ifeq ($(BUILD),debug)
    CFLAGS  += -O0 -g -fsanitize=address,undefined
    LDFLAGS += -fsanitize=address,undefined
else
    CFLAGS += -O2
endif

TARGET  := ect_demo
TESTBIN := test_mitm
BENCHBIN := bench_trace
MICROBIN := micro_bench
SOURCES := main.c elgamal.c utils.c ect.c nizk.c block.c schnorr.c sim.c scenario.c
OBJECTS := $(SOURCES:.c=.o)
HEADERS := params.h elgamal.h utils.h ect.h nizk.h block.h schnorr.h sim.h scenario.h \
           prg.h sig.h dj.h dj_board.h rlwe.h eg_board.h bench_util.h board_util.h

# Damgard-Jurik (GMP), Ring-LWE (own NTT) and the board proofs (OpenSSL for
# randomness, hashing and the Ed25519 / ML-DSA-44 issuance signatures).
DJOBJS   := dj.o dj_board.o prg.o sig.o board_util.o
RLWEOBJS := rlwe.o prg.o sig.o
NEWBINS  := test_dj test_dj_board test_rlwe test_eg_board bench_dj bench_rlwe bench_eg_board

.PHONY: all clean debug release test micro-bench tests benchmarks board-1e6 all-benchmarks everything powersum-check

all: $(TARGET)

debug:
	$(MAKE) BUILD=debug

release:
	$(MAKE) BUILD=release

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

# Correctness tests for Trace recovery (links every object but main.o).
TESTOBJS := $(filter-out main.o,$(OBJECTS))

$(TESTBIN): test_mitm.o $(TESTOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

test: $(TESTBIN)
	./$(TESTBIN)

$(BENCHBIN): bench_trace.o $(TESTOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

$(MICROBIN): micro_bench.o $(TESTOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

micro-bench: $(MICROBIN)
	./$(MICROBIN)

test_dj: test_dj.o $(DJOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lgmp -lcrypto -lm -lpthread

test_dj_board: test_dj_board.o $(DJOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lgmp -lcrypto -lm -lpthread

bench_dj: bench_dj.o $(DJOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lgmp -lcrypto -lm -lpthread

test_rlwe: test_rlwe.o $(RLWEOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lcrypto -lm

bench_rlwe: bench_rlwe.o $(RLWEOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lcrypto -lm

EGBOBJS  := eg_board.o board_util.o prg.o sig.o

test_eg_board: test_eg_board.o $(EGBOBJS) $(TESTOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS) -lpthread

bench_eg_board: bench_eg_board.o $(EGBOBJS) $(TESTOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS) -lpthread

# Correctness tests of the Damgard-Jurik and Ring-LWE instantiations and of
# the board proofs (test_dj at h = 10 takes several minutes).
tests: test_dj test_dj_board test_rlwe test_eg_board
	./test_eg_board
	./test_dj_board
	./test_rlwe 2048 16 36 1 1 32
	./test_rlwe 4096 16 65 2 2 32
	./test_rlwe 4096 16 97 2 4 32
	./test_rlwe 2048 10 36 1 1 32
	./test_rlwe 4096 10 65 2 2 32
	./test_rlwe 4096 10 97 2 4 32
	./test_dj 10 3072

# make benchmarks: the single-threaded suite; make board-1e6: the ElGamal
# board of 10^6 entries on 16 threads (about six hours); make
# all-benchmarks: both, one after the other; make everything: the tests, the
# suite, the policy simulation, the table rows and the board (about eight
# hours).
BENCHBINS := $(MICROBIN) $(BENCHBIN) bench_dj bench_rlwe bench_eg_board test_rlwe

benchmarks: $(BENCHBINS)
	./run_benchmarks.sh suite

board-1e6: bench_eg_board
	./run_benchmarks.sh board

all-benchmarks: $(BENCHBINS)
	./run_benchmarks.sh all

# The check of the power-sum encoding, in Python with GMP through ctypes
# (about fifteen minutes).
powersum-check:
	python3 powersum_check.py > ../results_powersum_check.txt

everything: $(BENCHBINS) test_dj test_dj_board test_rlwe test_eg_board
	./run_benchmarks.sh everything

clean:
	rm -f $(OBJECTS) test_mitm.o bench_trace.o micro_bench.o \
		$(TARGET) $(TESTBIN) $(BENCHBIN) $(MICROBIN) \
		$(NEWBINS) $(NEWBINS:=.o) dj.o dj_board.o prg.o sig.o rlwe.o eg_board.o board_util.o
