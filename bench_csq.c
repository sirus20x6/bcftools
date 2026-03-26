/*
 * Micro-benchmarks for bcftools csq hot paths.
 * Comparable to zig/bench/bench.zig
 *
 * Build: gcc -O2 -o bench_csq bench_csq.c csq.c gff.c regidx.c filter.c -I. -I../htslib ../htslib/libhts.a -lz -lpthread -lm -llzma -lbz2 -lcurl
 * Or simply: make bench_csq
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include "gff.h"

/* ---- Codon translation (from csq.c) ---- */
typedef struct { int id; const char *name, *code, *stop; } gencode_t;
static gencode_t gencode_std = {
    .id=0, .name="Standard simplified",
    .code="KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
    .stop="--------------M---------------------------------*-*-----*-------"
};

static const uint8_t nt4[] = {
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 3,4,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 3
};
static const uint8_t cnt4[] = {
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,3,4,2, 4,4,4,1, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 0,4,4,4, 4,4,4,4, 4,4,4,4,
    4,3,4,2, 4,4,4,1, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 0
};
#define _codon_idx(a,b,c) ((a)<<4 | (b)<<2 | (c))
#define _dna_idx(x) _codon_idx(nt4[(uint8_t)(x)[0]], nt4[(uint8_t)(x)[1]], nt4[(uint8_t)(x)[2]])
#define _cdna_idx(x) _codon_idx(cnt4[(uint8_t)(x)[2]], cnt4[(uint8_t)(x)[1]], cnt4[(uint8_t)(x)[0]])
#define dna2aa(x) (_dna_idx(x) > 63 ? 'X' : gencode_std.code[_dna_idx(x)])
#define cdna2aa(x) (_cdna_idx(x) > 63 ? 'X' : gencode_std.code[_cdna_idx(x)])

/* ---- regidx (from regidx.c) ---- */
#include "regidx.h"

/* ---- timing ---- */
static double elapsed_ms(struct timespec *start, struct timespec *end) {
    return (end->tv_sec - start->tv_sec) * 1000.0 + (end->tv_nsec - start->tv_nsec) / 1e6;
}

/* Simple PRNG (xorshift64) */
static uint64_t rng_state = 0x12345678ABCDEF01ULL;
static uint64_t rng_next(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

/* ---- Benchmark 1: Codon Translation ---- */
static void bench_codon_translation(void) {
    const int N = 1000000;
    const char bases[] = "ACGT";
    char codons[N][3];

    rng_state = 0x12345678ABCDEF01ULL;
    for (int i = 0; i < N; i++) {
        codons[i][0] = bases[rng_next() % 4];
        codons[i][1] = bases[rng_next() % 4];
        codons[i][2] = bases[rng_next() % 4];
    }

    struct timespec t0, t1;
    volatile uint8_t sink = 0;

    /* Forward */
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < N; i++) {
        sink += dna2aa(codons[i]);
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double ms = elapsed_ms(&t0, &t1);
    printf("  dna2aa (forward): %.1f ms  (%.0f ops/sec, n=%d)\n", ms, N / (ms / 1000.0), N);

    /* Reverse complement */
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < N; i++) {
        sink += cdna2aa(codons[i]);
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    ms = elapsed_ms(&t0, &t1);
    printf("  cdna2aa (revcomp): %.1f ms  (%.0f ops/sec, n=%d)\n", ms, N / (ms / 1000.0), N);
    (void)sink;
}

/* ---- Benchmark 2: regidx overlap queries ---- */
static void bench_regidx_overlap(void) {
    const int N_INTERVALS = 50000;
    const int N_QUERIES = 1000000;
    const int N_CHROMS = 25;

    /* Build regidx */
    regidx_t *idx = regidx_init(NULL, NULL, NULL, 0, NULL);
    char chr[16];
    for (int i = 0; i < N_INTERVALS; i++) {
        snprintf(chr, sizeof(chr), "chr%d", (i % N_CHROMS) + 1);
        uint32_t beg = (rng_next() % 100000000);
        uint32_t end = beg + (rng_next() % 1000) + 1;
        regidx_push(idx, chr, chr + strlen(chr) - 1, beg, end, NULL);
    }

    regitr_t *itr = regitr_init(idx);

    struct timespec t0, t1;
    int total_hits = 0;

    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < N_QUERIES; i++) {
        int ci = (rng_next() % N_CHROMS) + 1;
        snprintf(chr, sizeof(chr), "chr%d", ci);
        uint32_t qbeg = rng_next() % 100000000;
        uint32_t qend = qbeg + (rng_next() % 500);
        if (regidx_overlap(idx, chr, qbeg, qend, itr)) {
            while (regitr_overlap(itr)) {
                total_hits++;
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double ms = elapsed_ms(&t0, &t1);
    printf("  overlap queries: %.1f ms  (%.0f ops/sec, n=%d)\n", ms, N_QUERIES / (ms / 1000.0), N_QUERIES);
    printf("    %d intervals across %d chroms, %d total hits\n", N_INTERVALS, N_CHROMS, total_hits);

    regitr_destroy(itr);
    regidx_destroy(idx);
}

/* ---- Benchmark 3: malloc/free vs arena pattern ---- */
static void bench_malloc_pattern(void) {
    /* Simulate the C pattern: alloc per-transcript, free individually */
    const int N = 100000;
    struct timespec t0, t1;

    /* C pattern: individual malloc/free */
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < N; i++) {
        void *p1 = calloc(1, 256);  /* ref */
        void *p2 = calloc(1, 512);  /* sref */
        void *p3 = calloc(1, 64);   /* hap_node_t root */
        void *p4 = calloc(1, 128);  /* hap array */
        /* simulate some work */
        memset(p1, 'A', 256);
        memset(p2, 'C', 512);
        free(p4);
        free(p3);
        free(p2);
        free(p1);
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double ms = elapsed_ms(&t0, &t1);
    printf("  malloc/free (4 allocs per iter): %.1f ms  (%.0f ops/sec, n=%d)\n", ms, N / (ms / 1000.0), N);
}

int main(void) {
    printf("bcftools C micro-benchmarks\n");
    printf("===========================\n\n");

    printf("[1] Codon Translation (dna2aa)\n");
    bench_codon_translation();

    printf("\n[2] RegionIndex Overlap Queries\n");
    rng_state = 0xDEADBEEFCAFE1234ULL;
    bench_regidx_overlap();

    printf("\n[3] Malloc/Free Pattern (per-transcript)\n");
    bench_malloc_pattern();

    printf("\nDone.\n");
    return 0;
}
