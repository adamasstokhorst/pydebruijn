import sympy as _sympy
from .helpers import order_from_anf as _order_from_anf, anf_from_truth_table as _anf_from_truth_table
from .fsr import FeedbackShiftRegister as _FSR

__all__ = ['gpo']


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
