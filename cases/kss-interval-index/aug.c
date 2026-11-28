/* aug.c -- exhaustive augmentation search for small noncol.
   Reads connected graphs on m vertices (graph6, one per line) from stdin.
   Mode "aug": for each graph and each subset S of its vertices with |S|>=2 containing every
   degree-1 vertex, forms G+v (v adjacent to S) on N=m+1 vertices (so every delta>=2 connected
   graph on N vertices arises, since any connected graph has a non-cut vertex) and counts
   noncol exactly with early exit above CAP. Prints "g6 noncol diam" for noncol<=CAP.
   Mode "eval": prints "noncol diam" for each input graph (for cross-checking against Python).
   Usage: aug aug CAP < conn_m.g6     |   aug eval < graphs.g6                     */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXN 16
static int N;
static unsigned adj[MAXN];
static int D[MAXN][MAXN];

static int parse_g6(const char *s, unsigned *a) {
    int n = s[0] - 63; int k = 0;
    for (int i = 0; i < n; i++) a[i] = 0;
    const char *p = s + 1; int bit = 0; int val = 0; int nb = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
        if (nb == 0) { val = (*p++) - 63; nb = 6; }
        nb--; bit = (val >> nb) & 1;
        if (bit) { a[i] |= 1u << j; a[j] |= 1u << i; }
        k++;
    }
    return n;
}
static void to_g6(int n, unsigned *a, char *out) {
    int pos = 0; out[pos++] = (char)(n + 63);
    int val = 0, nb = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
        val = (val << 1) | ((a[i] >> j) & 1); nb++;
        if (nb == 6) { out[pos++] = (char)(val + 63); val = 0; nb = 0; }
    }
    if (nb) { val <<= (6 - nb); out[pos++] = (char)(val + 63); }
    out[pos] = 0;
}
static int apsp(int n, unsigned *a) { /* returns diameter, -1 if disconnected */
    int diam = 0;
    for (int s = 0; s < n; s++) {
        unsigned seen = 1u << s, front = 1u << s; int d = 0;
        for (int i = 0; i < n; i++) D[s][i] = -1;
        D[s][s] = 0;
        while (front) {
            unsigned nxt = 0;
            for (int x = 0; x < n; x++) if (front >> x & 1) nxt |= a[x];
            nxt &= ~seen; d++;
            for (int x = 0; x < n; x++) if (nxt >> x & 1) D[s][x] = d;
            seen |= nxt; front = nxt;
        }
        if (seen != ((n == 32) ? 0xffffffffu : ((1u << n) - 1))) return -1;
        if (d - 1 > diam) diam = d - 1;
    }
    return diam;
}
static int noncol(int n, int cap) {
    int c = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) for (int k = j + 1; k < n; k++) {
        int a = D[i][j], b = D[j][k], e = D[i][k];
        if (a + b == e || a + e == b || b + e == a) continue;
        if (++c > cap) return c;
    }
    return c;
}
int main(int argc, char **argv) {
    char line[256], out[64];
    unsigned a[MAXN];
    int mode_eval = (argc > 1 && strcmp(argv[1], "eval") == 0);
    int cap = (argc > 2) ? atoi(argv[2]) : 1000000;
    long long cand = 0, hits = 0;
    while (fgets(line, sizeof line, stdin)) {
        line[strcspn(line, "\r\n")] = 0;
        if (!line[0]) continue;
        int m = parse_g6(line, a);
        if (mode_eval) {
            int dm = apsp(m, a);
            printf("%d %d\n", dm < 0 ? -1 : noncol(m, 1000000), dm);
            continue;
        }
        int n = m + 1;
        unsigned must = 0;
        for (int i = 0; i < m; i++) if (__builtin_popcount(a[i]) == 1) must |= 1u << i;
        for (unsigned S = 0; S < (1u << m); S++) {
            if ((S & must) != must || __builtin_popcount(S) < 2) continue;
            unsigned b[MAXN];
            for (int i = 0; i < m; i++) b[i] = a[i] | (((S >> i) & 1) << m);
            b[m] = S;
            cand++;
            int dm = apsp(n, b);
            if (dm < 0) continue;
            int c = noncol(n, cap);
            if (c <= cap) { to_g6(n, b, out); printf("%s %d %d\n", out, c, dm); hits++; }
        }
    }
    fprintf(stderr, "candidates %lld hits %lld\n", cand, hits);
    return 0;
}
