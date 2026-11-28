/* Enumerate the cyclic flats of a binary matroid (column matroid of a k x n matrix over GF(2), n <= 32).
 * Input (stdin):  n k   then n columns, each as an unsigned integer (bit j = row j).
 * Output (stdout): one line per cyclic flat: "<mask as unsigned> <rank> <size>".
 * Method: circuits = supports of dual codewords w with rank(supp w) = |supp w| - 1;
 * BFS from cl(empty set): every cyclic flat Y > X is reached by joining one circuit at a time,
 * cl(X u C) for circuits C not contained in X (each intermediate join is a cyclic flat).
 * Build: cc -O2 -o cycflats cycflats.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int n, k;
static uint32_t col[32];

/* Gaussian basis: basis[b] has leading bit b (as row-vector in GF(2)^k, k <= 32) */
typedef struct { uint32_t v[32]; uint32_t have; int r; } Basis;

static void binit(Basis *B) { memset(B, 0, sizeof *B); }
static uint32_t breduce(const Basis *B, uint32_t x) {
    for (int b = 31; b >= 0; b--) if ((x >> b & 1) && (B->have >> b & 1)) x ^= B->v[b];
    return x;
}
static int badd(Basis *B, uint32_t x) {
    x = breduce(B, x);
    if (!x) return 0;
    int b = 31; while (!(x >> b & 1)) b--;
    B->v[b] = x; B->have |= 1u << b; B->r++; return 1;
}
static uint32_t closure_from(Basis *B) {
    uint32_t m = 0;
    for (int e = 0; e < n; e++) if (!breduce(B, col[e])) m |= 1u << e;
    return m;
}
static int rank_of(uint32_t mask) {
    Basis B; binit(&B);
    for (int e = 0; e < n; e++) if (mask >> e & 1) badd(&B, col[e]);
    return B.r;
}

/* hash set of uint32 masks (mask+1 stored, 0 = empty) */
static uint32_t *H; static uint64_t Hcap, Hcnt;
static int hins(uint32_t m) {
    uint64_t h = ((uint64_t)m * 0x9E3779B97F4A7C15ull) >> 20;
    for (uint64_t i = h & (Hcap - 1);; i = (i + 1) & (Hcap - 1)) {
        if (H[i] == 0) { H[i] = m + 1; Hcnt++; return 1; }
        if (H[i] == m + 1) return 0;
    }
}

int main(void) {
    if (scanf("%d %d", &n, &k) != 2 || n > 32 || k > 32) { fprintf(stderr, "bad input\n"); return 1; }
    for (int e = 0; e < n; e++) { unsigned x; if (scanf("%u", &x) != 1) return 1; col[e] = x; }
    /* dual code: kernel of the k x n matrix, via elimination on the n columns */
    /* express dependencies: track combination masks */
    uint32_t vec[32], comb[32]; int piv[32]; int nb = 0;
    uint32_t dualbasis[32]; int nd = 0;
    for (int e = 0; e < n; e++) {
        uint32_t x = col[e], c = 1u << e;
        for (int i = 0; i < nb; i++) if (x >> piv[i] & 1) { x ^= vec[i]; c ^= comb[i]; }
        if (x) { int b = 31; while (!(x >> b & 1)) b--; vec[nb] = x; comb[nb] = c; piv[nb] = b; nb++;
                 /* keep reduced: eliminate b from earlier vectors not needed for correctness */ }
        else dualbasis[nd++] = c;
    }
    int rk = nb;
    if (nd > 26) { fprintf(stderr, "dual dimension %d too large\n", nd); return 1; }
    /* circuits */
    uint32_t *circ = malloc(sizeof(uint32_t) << nd); int nc = 0;
    for (uint64_t s = 1; s < (1ull << nd); s++) {
        uint32_t w = 0; for (int i = 0; i < nd; i++) if (s >> i & 1) w ^= dualbasis[i];
        int sz = __builtin_popcount(w);
        if (rank_of(w) == sz - 1) circ[nc++] = w;
    }
    fprintf(stderr, "n=%d rank=%d dual_dim=%d circuits=%d\n", n, rk, nd, nc);
    Hcap = 1ull << 24; H = calloc(Hcap, sizeof(uint32_t));
    uint32_t *queue = malloc(sizeof(uint32_t) * (Hcap / 2)); uint64_t qh = 0, qt = 0;
    Basis B0; binit(&B0);
    uint32_t bot = closure_from(&B0);
    hins(bot); queue[qt++] = bot;
    while (qh < qt) {
        uint32_t X = queue[qh++];
        Basis BX; binit(&BX);
        for (int e = 0; e < n; e++) if (X >> e & 1) badd(&BX, col[e]);
        for (int c = 0; c < nc; c++) {
            uint32_t C = circ[c];
            if ((C & X) == C) continue;
            Basis B = BX;
            uint32_t rest = C & ~X;
            for (int e = 0; e < n; e++) if (rest >> e & 1) badd(&B, col[e]);
            uint32_t Y = closure_from(&B);
            if (hins(Y)) { if (qt >= Hcap / 2) { fprintf(stderr, "queue overflow\n"); return 1; } queue[qt++] = Y; }
        }
        if (Hcnt * 2 > Hcap) { fprintf(stderr, "hash too full\n"); return 1; }
    }
    for (uint64_t i = 0; i < qt; i++) printf("%u %d %d\n", queue[i], rank_of(queue[i]), __builtin_popcount(queue[i]));
    fprintf(stderr, "cyclic flats: %llu\n", (unsigned long long)qt);
    return 0;
}
