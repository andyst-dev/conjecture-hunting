/* PRUNE hook for nauty's geng (compile geng.c with -DPRUNE=mtfprune), to be used with -t.
 *
 * geng builds each graph by adding vertices; the intermediate graphs are induced subgraphs
 * H' = H[0..k-1] of the final graph H on maxn vertices.
 * If H is maximal triangle-free (MTF), every "deficient" pair of H' (non-adjacent, with no common
 * neighbour inside H') must get a common neighbour among the rem = maxn - k vertices added later.
 * A later vertex w has N(w) n H' independent (H is triangle-free), so w covers only deficient
 * pairs that lie together inside one independent set of H'.
 * Two deficient pairs are "incompatible" if their union is not independent in H'; incompatible
 * pairs need different later vertices. Greedily collecting pairwise-incompatible deficient pairs
 * gives a lower bound on the number of later vertices needed, and we prune if it exceeds rem.
 * At k = maxn we keep only MTF graphs (rem = 0: no deficient pair allowed).
 * The bound is a necessary condition, so no MTF graph is lost. test_counts.sh checks the counts
 * against geng without pruning and against OEIS A216783.
 */
#include "gtools.h"

int mtfprune(graph *g, int k, int maxn) {
    int rem = maxn - k;
    int pa[1024], pb[1024]; int np = 0;
    setword un[1024];
    for (int a = 0; a < k; a++) {
        setword ra = g[a];
        for (int b = a + 1; b < k; b++) {
            if (ra & bit[b]) continue;
            if (ra & g[b]) continue;          /* has common neighbour */
            /* deficient pair */
            setword u = bit[a] | bit[b];
            int ok = 1;
            for (int t = 0; t < np; t++) {
                /* compatible if union {a,b,pa,pb} independent */
                setword w = u | un[t];
                int indep = 1;
                for (int x = 0; x < k && indep; x++)
                    if ((w & bit[x]) && (g[x] & w)) indep = 0;
                if (indep) { ok = 0; break; }
            }
            if (ok) {
                if (np >= rem) return 1;      /* need > rem later vertices */
                un[np] = u; pa[np] = a; pb[np] = b; np++;
            }
        }
    }
    return 0;
}
