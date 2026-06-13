"""
partition_numpy.py
==================
Numpy-accelerated partition representation and move operations,
shared by mc_fast.py (weighted sampler) and mc_uniform_fast.py
(uniform sampler).

State representation
--------------------
A partition of L into k parts is stored as two short numpy arrays,
always kept co-sorted in ascending order of part value:

    vals  : int32[n_d]   distinct part values  v_1 < v_2 < ... < v_{n_d}
    mults : int32[n_d]   multiplicities        m_1,  m_2, ...,  m_{n_d}

with  sum(vals * mults) == L  and  sum(mults) == k.

n_d = O(sqrt(L)) for typical partitions, so all operations below are
O(sqrt(L)) in both time and memory.

Move count formulas (O(n_d))
-----------------------------
n_transfers:
    For each (a, b) pair with a in vals, b in vals, b >= 2, a+1 != b,
    and (a != b or mult[a] >= 2).
    = n_v * n_d - n_adj - n_ei
    where n_v = n_d (all distinct values),
          n_d = number of decrementable values (val >= 2),
          n_adj = #{v in vals : v+1 also in vals}  (null-move pairs),
          n_ei  = #{v in dec  : mult[v] == 1}      (invalid same-value pairs).

n_splits:
    sum_a floor(a/2)  over all distinct a in vals.

n_merges:
    C(n_d, 2) + #{i : mults[i] >= 2}
    (distinct-value pairs + same-value pairs from multiplicity >= 2).

All move-count functions accept (vals, mults) numpy arrays and return
a Python int.
"""

import numpy as np
import math


# ============================================================
# Construction helpers
# ============================================================

def partition_to_arrays(parts_list):
    """
    Convert a list of part values (with repetition) to (vals, mults).
    Returns int32 numpy arrays sorted by ascending value.
    """
    from collections import Counter
    c = Counter(parts_list)
    vs = np.array(sorted(c.keys()), dtype=np.int32)
    ms = np.array([c[v] for v in vs], dtype=np.int32)
    return vs, ms


def arrays_to_parts_list(vals, mults):
    """Expand (vals, mults) back to a flat list of parts."""
    out = []
    for v, m in zip(vals.tolist(), mults.tolist()):
        out.extend([v] * m)
    return out


def arrays_to_tuple(vals, mults):
    """Canonical sorted-descending tuple, for exact comparisons."""
    out = []
    for v, m in zip(reversed(vals.tolist()), reversed(mults.tolist())):
        out.extend([v] * m)
    return tuple(out)


def init_arrays(L, rng_seed=42):
    """
    Initialise a partition near the expected part count.
    Returns (vals, mults, k, L).
    """
    import random
    k0   = max(1, int(1.28 * math.sqrt(L)))
    base = L // k0
    rem  = L % k0
    parts = [base + 1] * rem + [base] * (k0 - rem)
    vs, ms = partition_to_arrays(parts)
    return vs, ms


# ============================================================
# Insert / remove a single value in sorted arrays (O(n_d))
# ============================================================

def inc(vals, mults, v):
    return _inc(vals, mults, v)

def dec(vals, mults, v):
    return _dec(vals, mults, v)

def _inc(vals, mults, v):
    """
    Return new (vals, mults) with count of v incremented by 1.
    Inserts v if absent.
    """
    idx = np.searchsorted(vals, v)
    if idx < len(vals) and vals[idx] == v:
        new_m = mults.copy()
        new_m[idx] += 1
        return vals, new_m
    else:
        new_v = np.insert(vals,  idx, v).astype(np.int32)
        new_m = np.insert(mults, idx, 1).astype(np.int32)
        return new_v, new_m


def _dec(vals, mults, v):
    """
    Return new (vals, mults) with count of v decremented by 1.
    Removes v if count reaches 0.
    """
    idx = np.searchsorted(vals, v)
    # idx must be valid since v is present
    new_m = mults.copy()
    new_m[idx] -= 1
    if new_m[idx] == 0:
        new_v = np.delete(vals,  idx).astype(np.int32)
        new_m = np.delete(new_m, idx).astype(np.int32)
        return new_v, new_m
    return vals, new_m


# ============================================================
# Apply moves — all return (new_vals, new_mults, new_k)
# ============================================================

def apply_transfer_np(vals, mults, k, a, b):
    """
    Transfer: a -> a+1, b -> b-1.
    Preserves k and L.
    """
    v, m = _dec(vals,  mults,  a)
    v, m = _inc(v, m, a + 1)
    v, m = _dec(v, m, b)
    v, m = _inc(v, m, b - 1)
    return v, m, k


def apply_split_np(vals, mults, k, a, c):
    """
    Split: a -> (c, a-c).  k increases by 1.
    """
    d = a - c
    v, m = _dec(vals, mults, a)
    v, m = _inc(v, m, c)
    v, m = _inc(v, m, d)
    return v, m, k + 1


def apply_merge_np(vals, mults, k, c, d):
    """
    Merge: (c, d) -> c+d.  k decreases by 1.
    """
    v, m = _dec(vals,  mults,  c)
    v, m = _dec(v, m, d)
    v, m = _inc(v, m, c + d)
    return v, m, k - 1


# ============================================================
# Move count formulas (O(n_d), fully vectorised)
# ============================================================

def n_transfers_np(vals, mults):
    """
    Number of valid transfer moves.
    = n_v * n_dec - n_adj - n_ei
    """
    dec_mask = vals >= 2
    n_v   = len(vals)
    n_dec = int(dec_mask.sum())
    if n_dec == 0:
        return 0
    # n_adj: how many v satisfy v+1 also in vals
    vals_set = set(vals.tolist())
    n_adj = int(sum(1 for v in vals.tolist() if v + 1 in vals_set))
    # n_ei: decrementable values with multiplicity exactly 1
    n_ei  = int(((mults[dec_mask]) == 1).sum())
    return n_v * n_dec - n_adj - n_ei


def n_splits_np(vals):
    """Number of valid split moves = sum floor(a/2)."""
    return int((vals // 2).sum())


def n_merges_np(vals, mults):
    """
    Number of valid merge moves.
    = C(n_d, 2) + #{i : mults[i] >= 2}
    """
    n_d  = len(vals)
    n_rep = int((mults >= 2).sum())
    return n_d * (n_d - 1) // 2 + n_rep


# ============================================================
# All-moves enumeration (used for exact uniform transfer sampling
# and for the weighted sampler)
# ============================================================

def all_transfers_np(vals, mults):
    """
    Return list of (a, b) valid transfer pairs.
    Equivalent to mc.all_transfers but works on numpy arrays.
    """
    dec_vals = vals[vals >= 2].tolist()
    vals_set = set(vals.tolist())
    mults_d  = {int(v): int(m) for v, m in zip(vals, mults)}
    moves = []
    for a in vals.tolist():
        for b in dec_vals:
            if a + 1 == b:
                continue
            if a == b:
                if mults_d[a] >= 2:
                    moves.append((a, b))
                continue
            moves.append((a, b))
    return moves


def all_splits_np(vals):
    """Return list of (a, c) valid split pairs."""
    moves = []
    for a in vals.tolist():
        for c in range(1, a // 2 + 1):
            moves.append((a, c))
    return moves


def all_merges_np(vals, mults):
    """
    Return list of (c, d) valid merge pairs with c <= d.
    Includes same-value pairs when multiplicity >= 2.
    """
    vl = vals.tolist()
    ml = mults.tolist()
    moves = set()
    n_d = len(vl)
    for i in range(n_d):
        for j in range(i + 1, n_d):
            moves.add((vl[i], vl[j]))
        if ml[i] >= 2:
            moves.add((vl[i], vl[i]))
    return list(moves)


# ============================================================
# log_w for the weighted measure (mc_fast.py)
# ============================================================

def log_w_np(mults, k):
    """
    log w_tilde(lambda) = lgamma(k+1) + k*log(2) - sum lgamma(m_j+1)
    """
    r = math.lgamma(k + 1) + k * math.log(2)
    for m in mults.tolist():
        r -= math.lgamma(m + 1)
    return r


# ============================================================
# Autocorrelation and integrated autocorrelation time
# ============================================================

def autocorrelation_np(series, max_lag=500):
    """Vectorised autocorrelation using numpy."""
    x   = np.asarray(series, dtype=np.float64)
    n   = len(x)
    ml  = min(max_lag, n // 4)
    x  -= x.mean()
    var = np.dot(x, x) / n
    if var < 1e-12:
        return np.concatenate([[1.0], np.zeros(ml)])
    acf = np.empty(ml + 1)
    for t in range(ml + 1):
        acf[t] = np.dot(x[:n - t], x[t:]) / (var * (n - t))
    return acf


def tau_int_np(acf):
    tau = 0.5
    for t in range(1, len(acf)):
        tau += acf[t]
        if t >= 5 * tau:
            break
    return tau
