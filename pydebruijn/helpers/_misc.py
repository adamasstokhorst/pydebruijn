__all__ = ['powerset', 'hamming_weight']


def powerset(iterable, reverse=True):
    """
    Returns a generator that yields all subsets of an iterable.

    This recipe originated from Python 2's itertools documentation.

    Parameters
    ----------
    iterable : any iterable object

    reverse : bool (default=True)
        If True, yields subsets in decreasing order.

    Yields
    ------
    subset : tuple
        A subset of `iterable`.
    """
    import itertools
    s = list(iterable)
    if not reverse:
        return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))
    else:
        return itertools.chain.from_iterable(itertools.combinations(s, r) for r in reversed(range(len(s)+1)))


def hamming_weight(n):
    """
    Calculates the Hamming weight of a positive integer.

    The Hamming weight of an integer is defined as the number of
    1s in its binary representation.
    """
    c = 0
    while n:
        c += 1
        n &= n - 1
    return c
