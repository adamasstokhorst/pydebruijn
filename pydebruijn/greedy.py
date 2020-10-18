import sympy as _sympy
from .helpers import (
    order_from_anf as _order_from_anf,
    anf_from_truth_table as _anf_from_truth_table,
    lambda_left_shift as _lambda,
    theta_left_shift as _theta,
    is_necklace as _is_necklace,
    is_conecklace as _is_conecklace
)
from .fsr import FeedbackShiftRegister as _FSR

__all__ = ['gpo', 'theorem_7']


def gpo(func, init_state, order=None):
    """
    Implementation of Generalized Prefer-Opposite algorithm.

    Parameters
    ----------
    func : algebraic normal form
        A SymPy expression describing the algebraic normal form of a
        feedback shift register.  Must be using integer symbols named
        `x_k`, where `k=0, 1, ...`.

    init_state : initial state
        A list containing only 0s and 1s as its elements.  Must be the
        same length as the order of `func`.

    order : integer, optional (default=None)
        The order of the algebraic normal form. If None, then it will
        be deduced from `func`.

    Returns
    -------
    fsr : FeedbackShiftRegister object
        This FSR will be using the given algebraic normal form
        at its core.

    Raises
    ------
    ValueError
        If the length of `init_state` is less than the order of the FSR.
    """
    if order is None:
        order = _order_from_anf(func)
    if len(init_state) < order:
        raise ValueError('init_state too short: len={}, expected at least {}'.format(len(init_state), order))
    order = max(order, len(init_state))
    syms = _sympy.symbols('x_:{}'.format(order), integer=True)
    table = {'0' * (order - 1): 1, '1' * (order - 1): 1}
    state = init_state[:]
    seq = ''.join(map(str, state))

    # DEBUG
    # step = 0
    while True:
        # step += 1

        args = zip(syms, state)
        next_bit = int(func.subs(args)) % 2
        last_bit = state[0]
        state = state[1:] + [1 - next_bit]
        if ''.join(map(str, state)) in seq:
            state[-1] = 1 - state[-1]
        seq += str(state[-1])
        table[''.join(map(str, state[:-1]))] = state[-1] ^ last_bit
        if len(table) == 2**(order - 1):
            # print step
            break

    return _FSR(_anf_from_truth_table(table, order), order, [0] * order)


def theorem_7(n, *k_values):
    """
    Implementation of Theorem 7 from "An Efficiently Generated Family of Binary de Bruijn Sequences".

    Parameters
    ----------
    n : integer
        The order of the resulting de Bruijn sequence. Must be at least 3.

    k_values : list of integers
        Distinct integers, as specified in the paper. Exclude k_1 and k_t.
        The largest of these must be strictly less than n - 1.

    Returns
    -------
    fsr : FeedbackShiftRegister object
        This FSR object represents the resulting de Bruijn sequence.

    Raises
    ------
    ValueError
        If any of the arguments does not satisfy the constraints.
    """
    # Validate values
    if n < 3:
        raise ValueError('n too small (got {})'.format(n))

    k_values = list(k_values) + [1, n]
    k_values.sort()
    if k_values[-2] >= n - 1:
        raise ValueError('k_{{t-1}} too large (got {})'.format(k_values[-2]))
    if k_values[1] <= 1:
        raise ValueError('k_{{1}} too small (got {})'.format(k_values[1]))

    for i in range(len(k_values) - 1):
        if k_values[i] >= k_values[i+1]:
            raise ValueError('k values not distinct')

    # Construct truth table
    table = {'0' * (n - 1): 1, '1' * (n - 1): 1}
    for i in range(1, 2**(n - 1) - 1):
        binary = bin(i)[2:]
        binary = '{{:0>{}}}'.format(n - 1).format(binary)

        state = [int(s) for s in binary]
        next_state = state[0] + state[-1]
        if state[0] == 0 and _is_conecklace(state):
            next_state += 1
        elif state[0] == 1:
            state_wt = sum(state)
            k_comparison = [k <= state_wt for k in k_values]
            index = k_comparison.index(False) - 1
            if index >= 0 and _is_necklace(_lambda(state, k_values[index])):
                next_state += 1

        table[binary] = next_state % 2

    return _FSR(_anf_from_truth_table(table, n), n, [0] * n)


def theorem_9(n, k):
    """
    Implementation of Theorem 9 from "An Efficiently Generated Family of Binary de Bruijn Sequences".

    Parameters
    ----------
    n : integer
        The order of the resulting de Bruijn sequence. Must be at least 3.

    k : integer
        A positive integer strictly less than lcm(1, 2, 3, ..., n-2).

    Returns
    -------
    fsr : FeedbackShiftRegister object
        This FSR object represents the resulting de Bruijn sequence.

    Raises
    ------
    ValueError
        If any of the arguments does not satisfy the constraints.
    """
    # Validate values
    if n < 3:
        raise ValueError('n too small (got {})'.format(n))

    delta = _sympy.lcm_list(range(1, n - 1))
    if not 1 <= k <= delta:
        raise ValueError('k out of bounds (got {})'.format(k))

    # Construct truth table
    table = {'0' * (n - 1): 1, '1' * (n - 1): 1}
    for i in range(1, 2**(n - 1) - 1):
        binary = bin(i)[2:]
        binary = '{{:0>{}}}'.format(n - 1).format(binary)

        state = [int(s) for s in binary]
        next_state = state[0] + state[-1]
        if state[0] == 0 and _is_conecklace(state):
            next_state += 1
        elif state[0] == 1 and _is_necklace(_lambda(state, k)):
            next_state += 1

        table[binary] = next_state % 2

    return _FSR(_anf_from_truth_table(table, n), n, [0] * n)


def theorem_11(n, *k_values):
    """
    Implementation of Theorem 11 from "An Efficiently Generated Family of Binary de Bruijn Sequences".

    Parameters
    ----------
    n : integer
        The order of the resulting de Bruijn sequence. Must be at least 3.

    k_values : list of integers
        Distinct integers, as specified in the paper. Exclude k_1 and k_t.
        The largest of these must be strictly less than n - 1.

    Returns
    -------
    fsr : FeedbackShiftRegister object
        This FSR object represents the resulting de Bruijn sequence.

    Raises
    ------
    ValueError
        If any of the arguments does not satisfy the constraints.
    """
    # Validate values
    if n < 3:
        raise ValueError('n too small (got {})'.format(n))

    k_values = list(k_values) + [1, n]
    k_values.sort()
    if k_values[-2] >= n - 1:
        raise ValueError('k_{{t-1}} too large (got {})'.format(k_values[-2]))
    if k_values[1] <= 1:
        raise ValueError('k_2 too small (got {})'.format(k_values[1]))

    for i in range(len(k_values) - 1):
        if k_values[i] >= k_values[i+1]:
            raise ValueError('k values not distinct')

    # Construct truth table
    table = {'0' * (n - 1): 1, '1' * (n - 1): 1}
    for i in range(1, 2**(n - 1) - 1):
        binary = bin(i)[2:]
        binary = '{{:0>{}}}'.format(n - 1).format(binary)

        state = [int(s) for s in binary]
        u_state = [1 - s for s in state[:-1]] + [0]
        w_state = [0] + state[:-1]
        w_inv_state = [1 - s for s in w_state]
        next_state = state[0] + state[-1]
        if state[-1] == 1 and _is_conecklace(u_state):
            next_state += 1
        elif state[-1] == 0:
            state_wt = sum(w_inv_state)
            k_comparison = [k <= state_wt for k in k_values]
            index = k_comparison.index(False) - 1
            if index >= 0 and _is_necklace(_theta(w_state, k_values[index] - 1)):
                next_state += 1

        table[binary] = next_state % 2

    return _FSR(_anf_from_truth_table(table, n), n, [0] * n)


def theorem_12(n, k):
    """
    Implementation of Theorem 12 from "An Efficiently Generated Family of Binary de Bruijn Sequences".

    Parameters
    ----------
    n : integer
        The order of the resulting de Bruijn sequence. Must be at least 2.

    k : integer
        A nonnegative integer.

    Returns
    -------
    fsr : FeedbackShiftRegister object
        This FSR object represents the resulting de Bruijn sequence.

    Raises
    ------
    ValueError
        If any of the arguments does not satisfy the constraints.
    """
    # Validate values
    if n < 2:
        raise ValueError('n too small (got {})'.format(n))
    if k < 0:
        raise ValueError('k must be nonnegative (got {})'.format(k))

    # Construct truth table
    table = {'0' * (n - 1): 1, '1' * (n - 1): 1}
    for i in range(1, 2**(n - 1) - 1):
        binary = bin(i)[2:]
        binary = '{{:0>{}}}'.format(n - 1).format(binary)

        state = [int(s) for s in binary]
        u_state = [1 - s for s in state[:-1]] + [0]
        w_state = [0] + state[:-1]
        next_state = state[0] + state[-1]
        if state[-1] == 1 and _is_conecklace(u_state):
            next_state += 1
        elif state[-1] == 0 and _is_necklace(_theta(w_state, k)):
            next_state += 1

        table[binary] = next_state % 2

    return _FSR(_anf_from_truth_table(table, n), n, [0] * n)
