__all__ = ['powerset', 'hamming_weight', 'lambda_left_shift', 'is_necklace', 'is_conecklace']


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


def lambda_left_shift(state, repeat=1):
    """
    Left shifts a state repeatedly such that the first 1 is at the end.
    """
    if repeat < 0:
        raise ValueError('repeat must be non-negative. (got {})'.format(repeat))
    elif repeat == 0:
        return state
    else:
        shift_index = state.index(1) + 1
        return lambda_left_shift(state[shift_index:] + state[:shift_index], repeat - 1)


def is_necklace(state):
    p = 1
    for j in range(2, len(state) + 1):
        if state[j-1] < state[j-1 - p]:
            return False
        if state[j-1] > state[j-1 - p]:
            p = j
    return len(state) % p == 0


def is_conecklace(state):
    return is_necklace(state + [1 - s for s in state])
