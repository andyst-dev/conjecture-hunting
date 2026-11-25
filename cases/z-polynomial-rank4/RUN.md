# Running the long Z-polynomial search on your Mac

```
cd cases/z-polynomial-rank4
pip install sympy numpy
python3 test_zpoly.py            # must print known KL values and "True" everywhere
# rank 5 over GF(2), GF(3); rank 6 over GF(2); truncations of rank-6 configurations to rank 5
python3 search_z.py --rank 5 --p 2 --maxn 16 --minutes 240 --seed 1 &
python3 search_z.py --rank 5 --p 3 --maxn 14 --minutes 240 --seed 2 &
python3 search_z.py --rank 6 --p 2 --maxn 16 --minutes 240 --seed 3 &
python3 search_z.py --rank 6 --trunc 5 --p 3 --maxn 14 --minutes 240 --seed 4 &
```
- Progress goes to `data/z_search.jsonl` (best γ-root margin per restart; smaller means closer to a violation).
- Anything non-real-rooted goes to `data/z_candidates.jsonl` with `exact_real_rooted` from an exact Sturm count. Only `exact_real_rooted: false` matters. Then re-check it from scratch and contact Proudfoot before saying anything publicly.
- With Claude Code: "run RUN.md overnight in 4 parallel processes, then summarise data/z_search.jsonl and report any candidate".
