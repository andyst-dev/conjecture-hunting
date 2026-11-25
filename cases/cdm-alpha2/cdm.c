/* CDM checker (Costa-Luu-Wood-Yip, arXiv:2512.17114, Conjecture 1.12).
 *
 * Input: graph6 lines on stdin, each a triangle-free graph H (e.g. from `geng -tc n`).
 * G = complement(H) has alpha(G) <= 2.
 *
 * Reduction: if G has a dominating edge xy (N(x) u N(y) = V), CDM holds for G.
 *   In H terms: xy is a non-edge of H whose ends have no common H-neighbour.
 *   So G has no dominating edge <=> every non-edge of H has a common neighbour
 *   <=> H is maximal triangle-free (MTF; diameter <= 2 for n >= 3).
 *   G is disconnected <=> H = K_{a,b} (CDM does not apply).
 * For each MTF H != K_{a,b}, we search G exhaustively for a connected dominating matching M
 * (|M| >= 2; |M| = 1 is impossible here):
 *   - edges of M are non-edges of H, pairwise disjoint;
 *   - connected: for e=ab, f=cd in M, not all four pairs {a,b}x{c,d} are H-edges;
 *   - dominating: for every e=ab in M, C(e) = N_H(a) & N_H(b) is contained in V(M).
 * Search is complete: we branch on an uncovered required vertex z, and every valid M through the
 * current partial matching must contain a G-edge at z.
 *
 * Output: counts per run; any graph without a CDM is printed as "COUNTEREXAMPLE <graph6>".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef uint32_t set;
static int n;
static set adj[32];       /* H adjacency */
static set common[32][32];
static long long n_in = 0, n_mtf = 0, n_kab = 0, n_checked = 0, n_counter = 0;
static long long size_hist[17];

static int parse_g6(const char *s) {
    int nn = s[0] - 63; int k = 0, bit = 5; const char *p = s + 1;
    n = nn;
    for (int i = 0; i < n; i++) adj[i] = 0;
    for (int j = 1; j < n; j++)
        for (int i = 0; i < j; i++) {
            int v = (*p - 63) >> bit & 1;
            if (v) { adj[i] |= 1u << j; adj[j] |= 1u << i; }
            if (--bit < 0) { bit = 5; p++; }
        }
    (void)k;
    return n;
}

static int is_mtf(void) {
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (!(adj[i] >> j & 1) && !(adj[i] & adj[j])) return 0;
    return 1;
}

static int is_complete_bipartite(void) {
    /* H connected triangle-free; complete bipartite iff every vertex's non-neighbourhood
       (incl. itself) is independent and equals the part. */
    set all = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    set A = all & ~adj[0];         /* part containing 0 */
    for (int v = 0; v < n; v++) {
        set part = (A >> v & 1) ? A : (all & ~A);
        if ((adj[v] | part) != all || (adj[v] & part)) return 0;
    }
    return 1;
}

/* M stored as list of edges */
static int ma[16], mb[16], msz;
static set used;
static int best_found;

static int g_edge(int a, int b) { return a != b && !(adj[a] >> b & 1); }

static int connected_to_all(int a, int b) {
    for (int t = 0; t < msz; t++) {
        int c = ma[t], d = mb[t];
        /* need some G-edge between {a,b} and {c,d}: i.e. not all 4 pairs are H-edges */
        if ((adj[a] >> c & 1) && (adj[a] >> d & 1) && (adj[b] >> c & 1) && (adj[b] >> d & 1)) return 0;
    }
    return 1;
}

static int search(set required) {
    set missing = required & ~used;
    if (!missing) {
        if (msz >= 1) { size_hist[msz]++; return 1; }
        return 0;
    }
    int z = __builtin_ctz(missing);
    /* branch: z must be matched with some partner w (G-edge zw), disjoint, connected, and its
       common set adds to required */
    for (int w = 0; w < n; w++) {
        if (w == z || (used >> w & 1) || !g_edge(z, w)) continue;
        if (!connected_to_all(z, w)) continue;
        ma[msz] = z; mb[msz] = w; msz++;
        used |= (1u << z) | (1u << w);
        int ok = search(required | common[z][w]);
        msz--; used &= ~((1u << z) | (1u << w));
        if (ok) return 1;
    }
    return 0;
}

static int has_cdm(void) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) common[i][j] = adj[i] & adj[j];
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++) {
            if (!g_edge(a, b)) continue;
            msz = 1; ma[0] = a; mb[0] = b; used = (1u << a) | (1u << b);
            if (search(common[a][b])) return 1;
        }
    return 0;
}

int main(int argc, char **argv) {
    char line[1024];
    int rawmode = (argc > 1 && strcmp(argv[1], "--raw") == 0);   /* test mode: print has_cdm for every input */
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == '>' || line[0] == '\n') continue;
        line[strcspn(line, "\r\n")] = 0;
        parse_g6(line);
        n_in++;
        if (rawmode) { printf("%s %d\n", line, has_cdm()); continue; }
        if (n < 3 || !is_mtf()) continue;
        n_mtf++;
        if (is_complete_bipartite()) { n_kab++; continue; }
        n_checked++;
        if (!has_cdm()) { n_counter++; printf("COUNTEREXAMPLE %s\n", line); fflush(stdout); }
    }
    printf("n=%d input=%lld mtf=%lld complete_bipartite=%lld checked=%lld counterexamples=%lld\n",
           n, n_in, n_mtf, n_kab, n_checked, n_counter);
    printf("size of first CDM found:");
    for (int s = 1; s < 17; s++) if (size_hist[s]) printf(" |M|=%d:%lld", s, size_hist[s]);
    printf("\n");
    return 0;
}
