import sympy
from sympy.abc import x
from ._misc import hamming_weight

__all__ = ['is_primitive', 'generate_primitives', 'get_associate_poly',
           'lfsr_from_poly', 'seq_decimation', 'poly_decimation', 'get_special_state']


def is_primitive(poly):
    """Check whether a binary polynomial is primitive over GF(2).
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
    """Return a list of (possibly all) primitive polynomials over GF(2) of a given degree.
    """
    if degree == 1:
        return [sympy.Poly(x+1, x, modulus=2)]

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

    # TODO: use get representatives from _zechlogs.py?
    decimations = []
    for k in range(1, 2**(degree-1)):
        if sympy.gcd(k, 2**degree - 1) != 1:
            continue
        is_new_coset = True
        for i in range(1, degree):
            if (k * 2**i) % (2**degree - 1) in decimations:
                is_new_coset = False
                break
        if is_new_coset:
            decimations.append(k)
        if len(decimations) >= limit:
            break

    return sorted(map(lambda a: poly_decimation(poly, a), decimations), key=lambda a: a.all_coeffs())


def get_associate_poly(poly_list):
    """Compute associate primitive polynomials for all polynomials in a given list.
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
    """Compute the next state of a LFSR with a given feedback polynomial.
    """
    lfsr = list(reversed(poly.all_coeffs()))[:-1]
    next_state = state + [sum(map(lambda a, b: a & int(b), state, lfsr)) % 2]
    return next_state[1:]


def seq_decimation(p, t, offset=0, c_state=None):
    """Return the t-decimated sequence of a polynomial.
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
    """Return the t-decimated polynomial of a polynomial.
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
