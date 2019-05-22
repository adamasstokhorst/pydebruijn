__all__ = ['powerset', 'hamming_weight']


def powerset(iterable, reverse=True):
    """Return a generator that yields all subsets in decreasing order.
    """
    import itertools
    s = list(iterable)
    if not reverse:
        return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))
    else:
        return itertools.chain.from_iterable(itertools.combinations(s, r) for r in reversed(range(len(s)+1)))


def hamming_weight(n):
    """Calculate the Hamming weight of a positive integer.
    """
    c = 0
    while n:
        c += 1
        n &= n - 1
    return c
