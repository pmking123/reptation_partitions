"""
reptation_mc.py
===============
Correct Metropolis-Hastings sampler over integer partitions of L
with weight  w_tilde(lambda) = (k! / prod m_j!) * 2^k.

For each move type the proposal is UNIFORM over all valid moves of
that type, so q_fwd = 1/N_fwd and q_rev = 1/N_rev.
MH acceptance: min(1, exp(Delta_log_w) * N_fwd / N_rev).

Move types:
  Transfer:  enumerate all (a,b) pairs with a,b in parts, b>=2, a+1!=b,
             b-1!=a; propose one uniformly.
  Split:     enumerate all (a,c) with a in parts a>=2, 1<=c<=a//2;
             propose one uniformly.
  Merge:     enumerate all unordered pairs {c,d} of individual parts;
             propose one uniformly.

This ensures exact detailed balance and correct stationary distribution.
"""

import math
import random
import time
from collections import Counter


# ============================================================
# Partition
# ============================================================

class Partition:
    __slots__ = ('counts', 'k', 'L')

    def __init__(self, parts):
        self.counts = Counter(parts)
        self.k = len(parts)
        self.L = sum(parts)

    def copy(self):
        q = object.__new__(Partition)
        q.counts = Counter(self.counts)
        q.k = self.k
        q.L = self.L
        return q

    def to_tuple(self):
        out = []
        for v in sorted(self.counts, reverse=True):
            out.extend([v] * self.counts[v])
        return tuple(out)

    def log_w(self):
        r = math.lgamma(self.k + 1) + self.k * math.log(2)
        for m in self.counts.values():
            r -= math.lgamma(m + 1)
        return r

    def parts_list(self):
        """All individual part values as a list (with repetition)."""
        out = []
        for v, m in self.counts.items():
            out.extend([v] * m)
        return out


def _set(c, v, m):
    if m == 0:
        c.pop(v, None)
    else:
        c[v] = m


# ============================================================
# Enumerate valid moves
# ============================================================

def all_transfers(p):
    """
    All valid transfer moves (a->a+1, b->b-1).
    Conditions: a in parts, b in parts with b>=2.
    Excludes:
      - a+1==b: null move (a->a+1 then b=a+1->a cancels out)
      - a==b: requires m_a>=2; handled separately below
    Returns list of (a, b) tuples.
    """
    vals = list(p.counts.keys())
    dec_vals = [v for v in vals if v >= 2]
    moves = []
    for a in vals:
        for b in dec_vals:
            if a + 1 == b:
                continue  # null move
            if a == b:
                # a->a+1 and a->a-1: need two copies of a
                if p.counts[a] >= 2:
                    moves.append((a, b))
                continue
            moves.append((a, b))
    return moves


def all_splits(p):
    """
    All valid split moves: split part a into (c, a-c) with 1<=c<=a//2.
    Returns list of (a, c) tuples.
    """
    moves = []
    for a in p.counts:
        for c in range(1, a // 2 + 1):
            moves.append((a, c))
    return moves


def all_merges(p):
    """
    All valid merge moves: merge two individual parts.
    Returns list of unordered pairs as (c, d) with c <= d.
    """
    parts = p.parts_list()
    moves = set()
    k = len(parts)
    for i in range(k):
        for j in range(i+1, k):
            c, d = parts[i], parts[j]
            moves.add((min(c,d), max(c,d)))
    return list(moves)


# ============================================================
# Apply moves to get new partition
# ============================================================

def apply_transfer(p, a, b):
    """Apply transfer a->a+1, b->b-1. Returns new partition."""
    new_p = p.copy()
    c = new_p.counts
    if a == b:
        _set(c, a,   c.get(a, 0) - 2)
        _set(c, a+1, c.get(a+1, 0) + 1)
        _set(c, a-1, c.get(a-1, 0) + 1)
    else:
        _set(c, a,   c.get(a, 0) - 1)
        _set(c, a+1, c.get(a+1, 0) + 1)
        _set(c, b,   c.get(b, 0) - 1)
        _set(c, b-1, c.get(b-1, 0) + 1)
    return new_p


def apply_split(p, a, c):
    """Split part a into (c, a-c). Returns new partition."""
    d = a - c
    new_p = p.copy()
    ct = new_p.counts
    _set(ct, a, ct.get(a,0) - 1)
    _set(ct, c, ct.get(c,0) + 1)
    _set(ct, d, ct.get(d,0) + 1)
    new_p.k += 1
    return new_p


def apply_merge(p, c, d):
    """Merge parts c and d into c+d. Returns new partition."""
    new_p = p.copy()
    ct = new_p.counts
    _set(ct, c,   ct.get(c,0)   - 1)
    _set(ct, d,   ct.get(d,0)   - 1)
    _set(ct, c+d, ct.get(c+d,0) + 1)
    new_p.k -= 1
    return new_p


# ============================================================
# Metropolis-Hastings step with uniform proposals
# ============================================================

def metropolis_step(p, rng, p_transfer=0.5):
    """
    One MH step.
    With prob p_transfer: propose a transfer (uniform over valid transfers).
    Otherwise: with equal prob propose split or merge.

    Acceptance: min(1, exp(Delta_log_w) * N_fwd / N_rev)
    where N_fwd = number of valid moves of this type from p,
          N_rev = number of valid moves of this type from p'.
    """
    u = rng.random()

    if u < p_transfer:
        # Transfer move
        moves = all_transfers(p)
        if not moves:
            return p, False
        a, b   = rng.choice(moves)
        new_p  = apply_transfer(p, a, b)
        N_fwd  = len(moves)
        N_rev  = len(all_transfers(new_p))
        lw_old = p.log_w()
        lw_new = new_p.log_w()

    elif u < p_transfer + (1 - p_transfer) / 2:
        # Split move: reverse is a merge
        moves = all_splits(p)
        if not moves:
            return p, False
        a, c  = rng.choice(moves)
        new_p = apply_split(p, a, c)
        N_fwd = len(moves)
        N_rev = len(all_merges(new_p))   # reverse of split is merge
        lw_old = p.log_w()
        lw_new = new_p.log_w()

    else:
        # Merge move: reverse is a split
        moves = all_merges(p)
        if not moves:
            return p, False
        c, d  = rng.choice(moves)
        new_p = apply_merge(p, c, d)
        N_fwd = len(moves)
        N_rev = len(all_splits(new_p))   # reverse of merge is split
        lw_old = p.log_w()
        lw_new = new_p.log_w()

    # MH acceptance: exp(Delta_log_w + log(N_fwd/N_rev))
    log_alpha = (lw_new - lw_old) + math.log(N_fwd) - math.log(N_rev)
    if log_alpha >= 0 or rng.random() < math.exp(log_alpha):
        return new_p, True
    return p, False


# ============================================================
# Initialisation
# ============================================================

def init_partition(L, rng):
    k0   = max(1, int(1.28 * math.sqrt(L)))
    base = L // k0
    rem  = L % k0
    return Partition([base+1]*rem + [base]*(k0-rem))


# ============================================================
# Diagnostics
# ============================================================

def autocorrelation(series, max_lag=500):
    n    = len(series)
    ml   = min(max_lag, n//4)
    mean = sum(series)/n
    var  = sum((x-mean)**2 for x in series)/n
    if var < 1e-12:
        return [1.0]+[0.0]*ml
    return [sum((series[i]-mean)*(series[i+t]-mean) for i in range(n-t))
            / (var*(n-t)) for t in range(ml+1)]


def tau_int(acf):
    tau = 0.5
    for t in range(1, len(acf)):
        tau += acf[t]
        if t >= 5*tau:
            break
    return tau


# ============================================================
# Validation against exact w_tilde
# ============================================================

def validate(L_max=14, n_mc=300000):
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    try:
        from reptation_core import partitions as ep
    except ImportError:
        print("reptation_core.py not found; skipping.")
        return

    print("Validation: MC vs exact w_tilde")
    print()

    for L in range(4, L_max+1, 2):
        exact = {}
        Z = 0.0
        for p_t in ep(L):
            w = math.exp(Partition(list(p_t)).log_w())
            exact[p_t] = w; Z += w
        for p in exact: exact[p] /= Z
        km_ex = sum(len(p)*exact[p] for p in exact)
        k2_ex = sum(len(p)**2*exact[p] for p in exact)
        ks_ex = math.sqrt(max(k2_ex - km_ex**2, 0))

        rng  = random.Random(7)
        part = init_partition(L, rng)
        for _ in range(50000): part,_ = metropolis_step(part, rng)

        k_list = []; visits = Counter()
        for _ in range(n_mc):
            part, _ = metropolis_step(part, rng)
            k_list.append(part.k)
            visits[part.to_tuple()] += 1

        km_mc = sum(k_list)/n_mc
        ks_mc = math.sqrt(sum((k-km_mc)**2 for k in k_list)/n_mc)
        tv    = sum(abs(exact.get(p,0)-visits.get(p,0)/n_mc)
                    for p in set(exact)|set(visits))/2
        err   = abs(km_mc-km_ex)/km_ex*100

        print(f"  L={L:3d}  <k>_ex={km_ex:.4f}  <k>_MC={km_mc:.4f}  "
              f"err={err:.2f}%  TV={tv:.4f}  "
              f"sig_ex={ks_ex:.4f}  sig_MC={ks_mc:.4f}")
    print()


# ============================================================
# Run-length distribution
# ============================================================

def run_length_distribution(L, n_steps=500000, seed=42):
    PI     = math.pi
    pred_k = PI/math.sqrt(6)*math.sqrt(L)
    pred_r = math.sqrt(6)/PI*math.sqrt(L)

    rng  = random.Random(seed)
    part = init_partition(L, rng)
    for _ in range(n_steps//5): part,_ = metropolis_step(part, rng)

    dist=Counter(); k_sum=0
    for _ in range(n_steps):
        part,_ = metropolis_step(part, rng)
        k_sum += part.k
        for v,m in part.counts.items(): dist[v] += m

    total  = sum(dist.values())
    k_mean = k_sum/n_steps
    r_bar  = L/k_mean
    q      = max(0, 1-1/r_bar)

    print(f"  L={L}: <k>={k_mean:.2f} (pred {pred_k:.2f}), "
          f"r_bar={r_bar:.2f} (pred {pred_r:.2f})")
    print(f"  {'j':>4}  {'f_meas':>10}  {'f_pred':>10}  {'ratio':>7}")
    for j in range(1, min(max(dist)+1, int(6*r_bar)+2)):
        fm = dist.get(j,0)/total
        fp = (1/r_bar)*q**(j-1) if q>0 else 0
        if fm<5e-5 and fp<5e-5: break
        print(f"  {j:>4}  {fm:>10.5f}  {fp:>10.5f}  "
              f"  {fm/fp if fp>1e-9 else 0:>7.4f}")
    print()


# ============================================================
# Scaling study
# ============================================================

def scaling_study(L_vals, n_per_L=None, seed=42):
    PI       = math.pi
    print("Scaling study (weighted measure w_tilde = k! * 2^k / prod m_j!)")
    print("  NOTE: under this weighted measure E[k] ~ c * sqrt(L) * log(L),")
    print("  not sqrt(L). The ratio <k>/sqrt(L) grows as log(L) and does")
    print("  not converge. This is consistent with the canonical/GC")
    print("  inequivalence established in the paper (Section 6.3).")
    print("  The GC prediction pi/sqrt(6) applies only under the")
    print("  Boltzmann-weighted grand canonical ensemble, not here.")
    print()
    pred_k   = math.sqrt(6) / PI   # Erdos-Lehner coefficient for sqrt(L)*log(L)
    print(f"  Expected scaling: <k>/sqrt(L) ~ {pred_k:.4f} * log(L)  (growing)")
    print()
    hdr = (f"{'L':>6}  {'<k>/sqL':>9}  {'err%':>6}  "
           f"{'sig/L^.25':>10}  {'err%':>6}  "
           f"{'tau':>7}  {'acc':>6}  {'t(s)':>5}")
    print(hdr); print("-"*len(hdr))

    results = []
    for L in L_vals:
        n      = n_per_L or max(100000, 500*L)
        thin   = max(1, L//20)
        burnin = max(20000, 200*L)
        t0     = time.time()
        rng    = random.Random(seed)
        part   = init_partition(L, rng)

        for _ in range(burnin): part,_ = metropolis_step(part, rng)

        k_list=[]; acc=0
        for i in range(n*thin):
            part,a = metropolis_step(part, rng)
            acc += a
            if i%thin==0: k_list.append(part.k)
        t1 = time.time()

        sqL  = math.sqrt(L); L14=L**0.25
        km   = sum(k_list)/len(k_list)
        kvar = sum((k-km)**2 for k in k_list)/len(k_list)
        ksig = math.sqrt(max(kvar,0))
        ti   = tau_int(autocorrelation(k_list))
        ar   = acc/(n*thin)

        ratio_logL = km / (sqL * math.log(L)) if L > 1 else 0
        print(f"{L:>6}  {km/sqL:>9.4f}  {ratio_logL:>14.5f}  "
              f"{ksig/L14:>10.4f}  "
              f"{ti:>7.1f}  {ar:>6.3f}  {t1-t0:>5.1f}")
        results.append({'L':L,'km':km,'ksig':ksig,'tau':ti})
    return results


# ============================================================
# Entry point
# ============================================================

if __name__ == '__main__':
    print("\n" + "="*62)
    print("Reptation partition MC sampler")
    print("="*62 + "\n")

    validate(L_max=14, n_mc=300000)

    print("Run-length distributions")
    print()
    for L in [50, 100, 200]:
        run_length_distribution(L, n_steps=300000, seed=42)

    print("="*62)
    scaling_study([20,50,100,200,300,500], seed=42)
