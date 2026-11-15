#!/usr/bin/env python3
"""G16. Converged gaps of the SU(2) strip: raise the electric-energy cut until the gap stops moving, then Aitken-extrapolate.

For each (N, g2) the cut is raised in steps until the gap changes by less than `tol` between successive cuts or the number of
states exceeds `max_states`.  The cut-limit is estimated by Aitken's delta-squared from the last three cuts (only if the successive
differences have the same sign and a ratio in (0, 1)); it is NOT reliable when the last step is still larger than about 1e-3.

Prints one JSON line per cut and one per (N, g2).  Cost at g2 = 0.4: N = 4, cut 80 (113,051 states) takes about 30 s; the matrix
assembly in Python and the Lanczos solve dominate (state enumeration is milliseconds).  Memory: about 2 GB at 475,000 states.

Run: python3 python/qg/strip_gap_scan.py N g2[,g2...] step start max_states   e.g.  4 0.8,0.6 8 32 250000
"""
from __future__ import annotations

import json
import sys
import time

try:
    from qg import spin_network_strip as sn
except ImportError:                                    # run as a script from python/qg
    import spin_network_strip as sn


def aitken(gaps):
    """Aitken estimate of the limit of the last three values, or None."""
    if len(gaps) < 3:
        return None
    d1, d0 = gaps[-1] - gaps[-2], gaps[-2] - gaps[-3]
    if d0 == 0 or not 0 < d1 / d0 < 1:
        return None
    r = d1 / d0
    return gaps[-1] + d1 * r / (1 - r)


def scan(N, g2s, step, start, max_states, tol=2e-5, out=print):
    results = {}
    for g2 in g2s:
        gaps, lams, lam = [], [], start
        while True:
            if len(sn.Net(N).states(lam)) > max_states:
                break
            t0 = time.time()
            w, n = sn.lowest(N, lam, g2, k=3)
            gaps.append(float(w[1] - w[0])); lams.append(lam)
            out(json.dumps(dict(kind="cut", N=N, g2=g2, lam=lam, states=n, gap=gaps[-1], E0=float(w[0]), sec=round(time.time() - t0, 1))))
            if len(gaps) >= 2 and abs(gaps[-1] - gaps[-2]) < tol:
                break
            lam += step
        fin = dict(kind="final", N=N, g2=g2, last_lam=lams[-1], last_gap=gaps[-1],
                   last_diff=(gaps[-1] - gaps[-2]) if len(gaps) > 1 else None, aitken=aitken(gaps))
        out(json.dumps(fin))
        results[g2] = fin
    return results


if __name__ == "__main__":
    N = int(sys.argv[1]); g2s = [float(x) for x in sys.argv[2].split(",")]
    scan(N, g2s, int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), out=lambda s: print(s, flush=True))
