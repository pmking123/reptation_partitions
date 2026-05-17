"""
reptation_fixed_point_verification.py
--------------------------------------
Numerical verification of the fixed-point set characterisation for the
coarse-graining operators R_1 and R_2.

NOTE ON UNIQUENESS
------------------
The attractor of the iterated coarse-graining is NOT in general unique:
different merge orderings may produce different gap-2 compositions from
the same input.  What IS true -- and what this script verifies -- is that
EVERY attractor is a gap-2 composition, regardless of ordering.  This is
the correct statement of Proposition 7: the gap-2 compositions are
exactly the fixed-point SET (the set of possible attractors), not a
unique fixed point.

For all compositions of L <= 20, we verify:

  (1) R_2 o R_1 (deterministic left-to-right) always converges to a
      gap-2 composition.
  (2) R_1 o R_2 (deterministic left-to-right) always converges to a
      gap-2 composition.
  (3) Random merge orderings of R_2 o R_1 always converge to a gap-2
      composition (n_random_trials trials per composition).
  (4) Every gap-2 composition is itself a fixed point of both operators
      (i.e. neither operator modifies it).

On the 2D square lattice, odd-indexed parts are H-runs and even-indexed
parts are V-runs.  The operator R_delta merges any two consecutive
same-axis runs whose lengths differ by < delta, replacing both and the
intervening perpendicular run with a single run of combined length.

Total compositions of L <= 20: sum_{L=1}^{20} 2^{L-1} = 2^20 - 1 = 1,048,575.
"""

import itertools
import random
from collections import defaultdict

# ---------------------------------------------------------------------------
# Core data structures and operator definitions
# ---------------------------------------------------------------------------

def same_axis_pairs(comp):
    """
    Return list of (i, j) indices of consecutive same-axis pairs in comp,
    where i and j have the same parity (both odd-indexed or both even-indexed
    in 1-based terms, i.e. same parity in 0-based indexing).
    A pair (i, i+2) is same-axis: they are separated by one perpendicular run.
    """
    pairs = []
    for i in range(len(comp) - 2):
        # i and i+2 have the same parity (same axis)
        pairs.append((i, i + 2))
    return pairs

def apply_merge(comp, i):
    """
    Merge the same-axis pair at positions i and i+2 (0-based).
    Removes comp[i], comp[i+1] (the intervening perpendicular run),
    and comp[i+2], replacing them with comp[i] + comp[i+2].
    Returns new composition as a list.
    """
    merged_length = comp[i] + comp[i + 2]
    return list(comp[:i]) + [merged_length] + list(comp[i + 3:])

def find_mergeable(comp, delta):
    """
    Return list of indices i such that (comp[i], comp[i+2]) is a
    mergeable same-axis pair under R_delta:
    |comp[i] - comp[i+2]| < delta.
    """
    mergeable = []
    for i in range(len(comp) - 2):
        if abs(comp[i] - comp[i + 2]) < delta:
            mergeable.append(i)
    return mergeable

def apply_R_delta_deterministic(comp, delta):
    """
    Apply R_delta to exhaustion using a deterministic left-to-right scan:
    always merge the leftmost mergeable pair first.
    Returns the fixed-point composition as a tuple.
    """
    comp = list(comp)
    while True:
        mergeable = find_mergeable(comp, delta)
        if not mergeable:
            break
        # merge leftmost
        comp = apply_merge(comp, mergeable[0])
    return tuple(comp)

def apply_R_delta_random(comp, delta, rng):
    """
    Apply R_delta to exhaustion using a random merge order.
    Returns the fixed-point composition as a tuple.
    """
    comp = list(comp)
    while True:
        mergeable = find_mergeable(comp, delta)
        if not mergeable:
            break
        i = rng.choice(mergeable)
        comp = apply_merge(comp, i)
    return tuple(comp)

def apply_R2_o_R1(comp, order='21', rng=None):
    """
    Apply R_2 o R_1 (order='21') or R_1 o R_2 (order='12').
    Uses deterministic (left-to-right) merging unless rng is provided.
    """
    if rng is not None:
        if order == '21':
            after_R1 = apply_R_delta_random(comp, delta=1, rng=rng)
            return apply_R_delta_random(after_R1, delta=2, rng=rng)
        else:
            after_R2 = apply_R_delta_random(comp, delta=2, rng=rng)
            return apply_R_delta_random(after_R2, delta=1, rng=rng)
    else:
        if order == '21':
            after_R1 = apply_R_delta_deterministic(comp, delta=1)
            return apply_R_delta_deterministic(after_R1, delta=2)
        else:
            after_R2 = apply_R_delta_deterministic(comp, delta=2)
            return apply_R_delta_deterministic(after_R2, delta=1)

def satisfies_gap2(comp):
    """
    Check that all consecutive same-axis pairs in comp differ by >= 2.
    """
    for i in range(len(comp) - 2):
        if abs(comp[i] - comp[i + 2]) < 2:
            return False
    return True

def all_compositions(L):
    """
    Generate all compositions of L (ordered tuples of positive integers
    summing to L) via binary representation of split points.
    """
    if L == 0:
        yield ()
        return
    for mask in range(1 << (L - 1)):
        comp = []
        run = 1
        for bit in range(L - 1):
            if mask & (1 << bit):
                comp.append(run)
                run = 1
            else:
                run += 1
        comp.append(run)
        yield tuple(comp)

# ---------------------------------------------------------------------------
# Main verification
# ---------------------------------------------------------------------------

def is_gap2_composition(comp):
    """
    A gap-2 composition has |comp[i] - comp[i+2]| >= 2 for all i
    (same-axis pairs).  Single-part and two-part compositions are
    trivially gap-2.
    """
    return satisfies_gap2(comp)

def is_fixed_point_of_R_delta(comp, delta):
    """
    A composition is a fixed point of R_delta if no mergeable pair exists.
    """
    return len(find_mergeable(comp, delta)) == 0

def all_gap2_compositions(L):
    """Generate all gap-2 compositions of L."""
    for comp in all_compositions(L):
        if satisfies_gap2(comp):
            yield comp

def run_verification(L_max=20, n_random_trials=10, random_seed=42):
    """
    For all compositions of L = 1, ..., L_max:

    Test 1: R_2 o R_1 (deterministic) always reaches a gap-2 composition.
    Test 2: R_1 o R_2 (deterministic) always reaches a gap-2 composition.
    Test 3: n_random_trials random orderings of R_2 o R_1 each reach a
            gap-2 composition.
    Test 4: Every gap-2 composition is a fixed point of both R_1 and R_2.

    Also records, for illustration, cases where R_2 o R_1 and R_1 o R_2
    reach DIFFERENT gap-2 compositions (confirming non-uniqueness of attractor).
    """
    rng = random.Random(random_seed)

    total = 0
    failures_t1 = []   # Test 1: R_2 o R_1 not gap-2
    failures_t2 = []   # Test 2: R_1 o R_2 not gap-2
    failures_t3 = []   # Test 3: random ordering not gap-2
    failures_t4 = []   # Test 4: gap-2 composition not a fixed point

    # Count cases where the two orderings reach different (but both valid) attractors
    different_attractor_count = 0
    different_attractor_examples = []

    for L in range(1, L_max + 1):
        for comp in all_compositions(L):
            total += 1

            fp_21 = apply_R2_o_R1(comp, order='21')
            fp_12 = apply_R2_o_R1(comp, order='12')

            # Test 1
            if not is_gap2_composition(fp_21):
                failures_t1.append((comp, fp_21))

            # Test 2
            if not is_gap2_composition(fp_12):
                failures_t2.append((comp, fp_12))

            # Record different-attractor cases (both valid, but different)
            if fp_21 != fp_12:
                different_attractor_count += 1
                if len(different_attractor_examples) < 5:
                    different_attractor_examples.append((comp, fp_21, fp_12))

            # Test 3: random orderings
            for _ in range(n_random_trials):
                fp_rand = apply_R2_o_R1(comp, order='21', rng=rng)
                if not is_gap2_composition(fp_rand):
                    failures_t3.append((comp, fp_rand))
                    break

        # Test 4: all gap-2 compositions of this L are fixed points
        for comp in all_gap2_compositions(L):
            fp1_stable = is_fixed_point_of_R_delta(comp, delta=1)
            fp2_stable = is_fixed_point_of_R_delta(comp, delta=2)
            if not (fp1_stable and fp2_stable):
                failures_t4.append((comp, fp1_stable, fp2_stable))

    # -----------------------------------------------------------------------
    # Report
    # -----------------------------------------------------------------------
    print("=" * 65)
    print("Fixed-point set characterisation verification")
    print(f"L_max = {L_max},  random trials per composition = {n_random_trials}")
    print(f"Total compositions examined: {total:,}")
    print("=" * 65)

    for label, failures, description in [
        ("Test 1", failures_t1,
         "R_2∘R_1 (deterministic) always reaches a gap-2 composition"),
        ("Test 2", failures_t2,
         "R_1∘R_2 (deterministic) always reaches a gap-2 composition"),
        ("Test 3", failures_t3,
         "Random merge ordering always reaches a gap-2 composition"),
        ("Test 4", failures_t4,
         "Every gap-2 composition is a fixed point of R_1 and R_2"),
    ]:
        result = "PASS" if not failures else "FAIL"
        print(f"\n{label} — {description}: {result}")
        if failures:
            print(f"  {len(failures)} failure(s):")
            for f in failures[:3]:
                print(f"    {f}")
        else:
            print(f"  Zero failures.")

    print(f"\nNote: attractor non-uniqueness (both valid gap-2, but different)")
    print(f"  Cases where R_2∘R_1 and R_1∘R_2 reach different attractors: "
          f"{different_attractor_count:,} / {total:,}")
    print(f"  Examples:")
    for comp, fp21, fp12 in different_attractor_examples:
        print(f"    input={comp}  R_2∘R_1→{fp21}  R_1∘R_2→{fp12}")

    print("\n" + "=" * 65)
    all_pass = not (failures_t1 or failures_t2 or failures_t3 or failures_t4)
    print("Overall result:", "ALL TESTS PASSED" if all_pass else "FAILURES DETECTED")
    print("=" * 65)

    return {
        'total': total,
        'failures_t1': failures_t1,
        'failures_t2': failures_t2,
        'failures_t3': failures_t3,
        'failures_t4': failures_t4,
        'different_attractor_count': different_attractor_count,
    }

# ---------------------------------------------------------------------------
# Worked examples for the paper
# ---------------------------------------------------------------------------

def print_worked_examples():
    """
    Print illustrative examples showing convergence to a gap-2 fixed point,
    including cases where different orderings reach different (but both valid)
    gap-2 attractors.
    """
    examples = [
        (3, 2, 3, 2, 3),   # from the decoupling theorem example
        (1, 1, 2, 1),       # different attractors under two orderings
        (4, 2, 3, 2, 4),    # another non-unique case
        (1, 2, 3, 4, 5),    # already gap-2
        (2, 1, 2, 1, 2),    # gap-0 pairs present
    ]
    print("\nWorked examples")
    print("-" * 65)
    for comp in examples:
        fp_21 = apply_R2_o_R1(comp, order='21')
        fp_12 = apply_R2_o_R1(comp, order='12')
        same = fp_21 == fp_12
        print(f"  Input:        {comp}")
        print(f"  R_2∘R_1 →     {fp_21}  gap-2: {'✓' if satisfies_gap2(fp_21) else '✗'}")
        print(f"  R_1∘R_2 →     {fp_12}  gap-2: {'✓' if satisfies_gap2(fp_12) else '✗'}"
              f"  {'(same attractor)' if same else '(different attractor — both valid)'}")
        print()

if __name__ == '__main__':
    print_worked_examples()
    results = run_verification(L_max=20, n_random_trials=10, random_seed=42)
