import sympy
from ._misc import powerset

__all__ = ['toggle_anf', 'order_from_anf']


def toggle_anf(anf, order=None):
    """changes an ANF to its modified version and vice versa"""
    new_anf = anf.func(*anf.args)
    if order is None:
        order = order_from_anf(new_anf)
    syms = sympy.symbols('x_:{}'.format(order), integer=True)
    for subset in powerset(range(1, order), reverse=False):
        new_anf += reduce(lambda a, b: a * b, map(lambda a: syms[a], list(subset)), 1)
    return sympy.Poly(new_anf, modulus=2).as_expr()


def order_from_anf(anf):
    """counts the number of variables in an ANF"""
    n = 0
    for arg in anf.args:
        if isinstance(arg, sympy.Symbol):
            m = int(str(arg).split('_')[-1]) + 1
            n = m if m > n else n
        else:
            m = order_from_anf(arg)
            n = m if m > n else n
    return n
