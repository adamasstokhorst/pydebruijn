import sympy as _sympy
from .helpers import order_from_anf as _order_from_anf, anf_from_truth_table as _anf_from_truth_table
from .fsr import FeedbackShiftRegister as _FSR

__all__ = ['gpo']


def gpo(func, init_state, order=None):
    """Implementation of Generalized Prefer-Opposite algorithm.
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
    step = 0
    while True:
        step += 1

        args = zip(syms, state)
        next_bit = int(func.subs(args)) % 2
        last_bit = state[0]
        state = state[1:] + [1 - next_bit]
        if ''.join(map(str, state)) in seq:
            state[-1] = 1 - state[-1]
        seq += str(state[-1])
        table[''.join(map(str, state[:-1]))] = state[-1] ^ last_bit
        if len(table) == 2**(order - 1):
            print step
            break

    return _FSR(_anf_from_truth_table(table, order), order, [0] * order)
