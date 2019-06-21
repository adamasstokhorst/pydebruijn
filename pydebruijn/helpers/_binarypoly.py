import sympy
from sympy.abc import x
from ._misc import hamming_weight
from ._zechlogs import get_representatives

__all__ = ['is_primitive', 'generate_primitives', 'get_associate_poly',
           'lfsr_from_poly', 'seq_decimation', 'poly_decimation', 'get_special_state']


def is_primitive(poly):
    """
    Checks whether a binary polynomial is primitive over GF(2).

    Parameters
    ----------
    poly : SymPy polynomial
        The polynomial must be using `x` as its generator and has
        `modulus` set to 2.

    Returns
    -------
    b : boolean
        True if `poly` is primitive over GF(2), False otherwise.

    Examples
    --------
    >>> is_primitive(sympy.Poly(x**3 + x + 1, x, modulus=2))
    True
    >>> is_primitive(sympy.Poly(x**4 + 1, x, modulus=2))  # reducible
    False
    >>> is_primitive(sympy.Poly(x**4 + x**3 + x**2 + x + 1, x, modulus=2))  # irreducible non-primitive
    False
    """
    if not poly.is_irreducible:
        return False

    degree = int(sympy.degree(poly))
    if sympy.isprime(2**degree - 1):
        return True

    for k in (d for d in sympy.divisors(2**degree-1) if degree < d < (2**degree - 1)):
        q = sympy.Poly(x**k + 1, x, modulus=2)
        if sympy.rem(q, poly, modulus=2).is_zero:
            return False
    return True


def generate_primitives(degree, limit=None):
    """
    Returns a list of primitive polynomials of a given degree.

    This function searches for the first primitive polynomial using
    brute force, and then uses decimations to obtain the rest.

    Parameters
    ----------
    degree : integer
        The degree of polynomials to be returned.

    limit : integer, optional (default=None)
        If None, returns all primitive polynomials of the given degree.
        Otherwise, returns at most `limit` polynomials.

    Returns
    -------
    l : list
        List of primitive polynomials of degree `degree`.
    """
    if degree == 1:
        return [sympy.Poly(x+1, x, modulus=2)]

    poly = None
    for k in xrange(1, 2**(degree-1)):
        if hamming_weight(k) % 2 == 1:
            if int(bin(k)[:1:-1] + '0' * (degree - int(sympy.log(k, 2)) - 2), 2) < k:
                continue
            proto_poly = 0
            power = 0
            while k:
                if k % 2:
                    proto_poly += x**power
                k >>= 1
                power += 1
            proto_poly *= x
            proto_poly += x**degree + 1
            poly = sympy.Poly(proto_poly, x, modulus=2)
            if is_primitive(poly):
                break

    decimations = get_representatives(degree)[:limit]
    temp_out = [poly_decimation(poly, t) for t in decimations if sympy.gcd(t, 2**degree - 1) == 1]

    return sorted(temp_out, key=lambda a: a.all_coeffs())


def get_associate_poly(poly_list):
    """
    Computes the associates of polynomials in a list.

    The associate polynomial of a polynomial `p` is defined as the
    primitive polynomial which yields `p` when decimated by a certain
    value `t`, also called its order.

    Parameters
    ----------
    poly_list : list of SymPy polynomials
        The polynomials must be using `x` as its generator and has
        `modulus` set to 2.

    Returns
    -------
    associates : list of dict
        The dictionary will have four keys:

        - 'poly' : the original polynomial
        - 'associate' : the associate polynomial
        - 'order' : the `t` value such that when t-decimated, the
          associate polynomial yields the input polynomial
        - 'period' : the period of the LFSR corresponding to the
          input polynomial

        The elements of this list will correspond to the elements of
        `poly_list` in the same order.
    """
    associates = []
    indexes = []
    for d in set(map(sympy.degree, poly_list)):
        poly_subset = [poly for poly in poly_list if sympy.degree(poly) == d]
        primitive_polys = generate_primitives(d)
        for poly in poly_subset:
            indexes.append(poly_list.index(poly))
            if poly in primitive_polys:
                associates.append({'poly': poly,
                                   'order': 1,
                                   'period': 2**d-1,
                                   'associate': poly})
                continue

            # routine to find sequence's period
            sequence = [1] + [0]*(d-1)
            e = 1
            divisor = sympy.divisors(2**d-1)
            for e in divisor:
                if e < d:
                    continue
                while len(sequence) < 2*e:
                    sequence.append(lfsr_from_poly(poly, sequence[-d:])[-1])

                is_periodic = True
                for i in range(e):
                    is_periodic = is_periodic and sequence[i] == sequence[i+e]
                    if not is_periodic:
                        is_periodic = False
                        break
                if is_periodic:
                    break

            # routine to find associate polynomial
            for primitive in primitive_polys:
                if poly_decimation(primitive, (2**d-1)/e) == poly:
                    associates.append({'poly': poly,
                                       'order': (2**d-1)/e,
                                       'period': e,
                                       'associate': primitive})
                    break
    return [p for _, p in sorted(zip(indexes, associates))]


def lfsr_from_poly(poly, state):
    r"""
    Computes the next state of a LFSR with a given feedback polynomial.

    If the polynomial given as :math:`x^n + a_{n-1}x_^{n-1} + \ldots + a_0`
    and state is given as :math:`(s_0, s_1, \ldots, s_{n-1})`, then
    the next bit is computed as

    .. math::

        \textup{next\_bit} = a_{n-1}s_0 + a_{n-2}s_1 + \dots + a_0s_{n-1}.

    Parameters
    ----------
    poly : SymPy polynomial
        The binary polynomial corresponding to the LFSR. Note that
        this module takes the leading term as a dummy term, and

    state : list
        This list is expected to only contain 0 or 1, and its length
        must match the degree of the polynomial.

    Returns
    -------
    next_state : list
        The state succeeding `state` under the given LFSR.
    """
    lfsr = list(reversed(poly.all_coeffs()))[:-1]
    next_state = state + [sum(map(lambda a, b: a & int(b), state, lfsr)) % 2]
    return next_state[1:]


def seq_decimation(p, t, offset=0, c_state=None):
    r"""
    Decimates the sequence corresponding to the given polynomial.

    If the sequence is :math:'\{a_n\}', then a k-decimation is defined
    as the subsequence :math:'\{a_kn\}'.

    Parameters
    ----------
    p : SymPy polynomial
        A binary polynomial corresponding to a LFSR.

    t : integer
        The decimation value.

    offset : integer, optional (default=0)
        Shift the sequence by this amount before decimating.

    c_state : list, optional (default=None)
        The initial state of the LFSR.  It must be the same length as
        the degree of `p` and contains only 0 or 1.  If None, the
        initial state will all 1s.

    Returns
    -------
    l : list
        The decimated sequence.
    """
    deg = sympy.degree(p)
    if c_state is None:
        c_state = [1] * deg
    ret = [0] * (2*deg)
    for _ in range(offset):
        c_state = lfsr_from_poly(p, c_state)
    indexes = map(lambda a: (a*t) % (2**deg - 1), range(2*deg))
    ctr = 0
    i = 0
    while ctr < 2*deg:
        if i in indexes:
            pos = indexes.index(i)
            ret[pos] = c_state[0]
            indexes[pos] = -1
            ctr += 1
        c_state = lfsr_from_poly(p, c_state)
        i += 1
        if i >= 2**deg - 1:
            i -= 2**deg - 1
    return ret


def poly_decimation(p, t):
    """
    Decimates the given polynomial and returns another polynomial.

    This function decimates the calculated sequence and reconstructs
    the polynomial corresponding to that sequence using the
    Berlekamp-Massey algorithm.

    Parameters
    ----------
    p : SymPy polynomial
        A binary polynomial corresponding to a LFSR.

    t : integer
        The decimation value.

    Returns
    -------
    q : SymPy polynomial
        The polynomial corresponding to the decimated sequence.
        Note that if the degree of `q` does not equal the degree of `p`,
        this function returns None.
    """
    from operator import mul
    n = sympy.degree(p)
    s = seq_decimation(p, t)
    while len(s) < 2*n:
        s += s

    cd = sympy.Poly(1, x, modulus=2)
    l, m, bd = 0, -1, 1
    for i in range(2*n):
        sub_cd = list(reversed(cd.all_coeffs()))[1:l+1]
        sub_s = list(reversed(s[i-l:i]))
        sub_cd += [0] * (len(sub_s) - len(sub_cd))
        disc = s[i] + sum(map(mul, sub_cd, sub_s))
        if disc % 2 == 1:
            td = cd
            cd += bd * sympy.Poly(x**(i - m), x, modulus=2)
            if l <= i/2:
                l = i + 1 - l
                m = i
                bd = td
    if sympy.degree(cd) == n:
        cd = sympy.Poly(reversed(cd.all_coeffs()), x, modulus=2)
        return cd
    else:
        return None


def get_special_state(p, t):
    """
    Finds the special state corresponding to a decimated polynomial.

    The special state is such that when the sequence corresponding
    to `p` is decimated with the special state as the initial state,
    the state yielded with zero offset is the conjugate state to
    the zero state.

    Parameters
    ----------
    p : SymPy polynomial
        A primitive binary polynomial.

    t : integer
        The decimation value.

    Returns
    -------
    l : list
        The special state.
    """
    # TODO: rewrite this function to be more clear
    from operator import add
    deg = sympy.degree(p)

    base = [1] + [0] * (deg - 1)
    base_state = seq_decimation(p, t, c_state=base)[:deg]
    ones = []
    init = base[:]
    for i in range(1, deg):
        init[i] = 1
        init[i-1] = 0
        if i % t != 0:
            state = seq_decimation(p, t, c_state=init)[:deg]
            ones.append((init[:], state))

    i = 0
    for i in range(2**len(ones)):
        cstate = base_state[:]
        for j in range(len(ones)):
            if i & 2**j != 0:
                cstate = map(add, cstate, ones[j][1])
        cstate = map(lambda x: x % 2, cstate)
        if cstate == base:
            break

    for j in range(len(ones)):
        if i & 2**j != 0:
            base = map(add, base, ones[j][0])
        base = map(lambda x: x % 2, base)

    return base
