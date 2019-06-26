import sympy
from ._misc import powerset

__all__ = ['toggle_anf', 'order_from_anf']


def toggle_anf(anf, order=None):
    """
    Removes or adds the zero cycle from a given feedback shift register.

    Parameters
    ----------
    anf : algebraic normal form
        A SymPy expression describing the algebraic normal form of a
        feedback shift register.  Must be using integer symbols named
        `x_k`, where `k=0, 1, ...`.

    order : integer, optional (default=None)
        The order of `anf`.  If None, then it will be deduced from `anf`.

    Returns
    -------
    new_anf : algebraic normal form
        The algebraic normal form of the feedback shift register
        corresponding to `anf` but with the zero cycle added or removed.

    Examples
    --------
    >>> x0, x1 = sympy.symbols('x_0 x_1', integer=True)
    >>> toggle_anf(x0 + x1)
    x_0 + 1
    >>> toggle_anf(x0 + x1, 3)
    x_0 + x_1*x_2 + x_2 + 1
    """
    new_anf = anf.func(*anf.args)
    if order is None:
        order = order_from_anf(new_anf)
    syms = sympy.symbols('x_:{}'.format(order), integer=True)
    for subset in powerset(range(1, order), reverse=False):
        new_anf += reduce(lambda a, b: a * b, map(lambda a: syms[a], list(subset)), 1)
    return sympy.Poly(new_anf, modulus=2).as_expr()


def order_from_anf(anf):
    """
    Counts the minimum order of a feedback shift register.

    This is accomplished by stringifying the symbols in the expression
    and finding the greatest index.

    Parameters
    ----------
    anf : algebraic normal form
        A SymPy expression describing the algebraic normal form of a
        feedback shift register.  Must be using integer symbols named
        `x_k`, where `k=0, 1, ...`.

    Returns
    -------
    n : integer
        The minimum order of `anf`.

    Examples
    --------
    >>> x0, x1, x2 = sympy.symbols('x_0 x_1 x_2', integer=True)
    >>> order_from_anf(x0 + x1)
    2
    >>> order_from_anf(x0 + x2)
    3
    """
    n = 0
    for arg in anf.args:
        if isinstance(arg, sympy.Symbol):
            m = int(str(arg).split('_')[-1]) + 1
            n = m if m > n else n
        else:
            m = order_from_anf(arg)
            n = m if m > n else n
    return n
