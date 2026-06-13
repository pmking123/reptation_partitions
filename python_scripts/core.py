"""
core.py
=================
Core combinatorial utilities shared across reptation scripts.
"""

def partitions(L):
    """
    Generate all partitions of L as tuples in non-increasing order.
    Uses a standard recursive algorithm.
    """
    def _helper(n, max_part):
        if n == 0:
            yield ()
            return
        for p in range(min(n, max_part), 0, -1):
            for rest in _helper(n - p, p):
                yield (p,) + rest

    return _helper(L, L)