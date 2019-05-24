__all__ = ['ZechContainer', 'get_representatives', 'get_decimation_value', 'retrieve_zech_log', 'fetch_and_save']


class ZechContainer(dict):
    def __init__(self, order):
        self._ord = order
        self._modulo = 2 ** order - 1
        super(ZechContainer, self).__init__()

    def __getitem__(self, key):
        if key in self:
            return super(ZechContainer, self).__getitem__(key)
        else:
            for i in xrange(1, self._ord):
                key *= 2
                key %= self._modulo
                if key in self:
                    return 2 ** (self._ord - i) * self[key] % self._modulo

        # if it doesn't find anything, then it's not covered by the representatives
        raise KeyError(str(2 * key % self._modulo))

    @property
    def order(self):
        return self._ord


def get_representatives(n):
    import collections

    modulus = 2**n - 1

    t_flag = collections.defaultdict(lambda: -1)
    for t in xrange(1, modulus):
        if t_flag[t] == -1:
            t_flag[t] = 1
            for i in xrange(1, n):
                t *= 2
                t %= modulus
                if t_flag[t] == -1:
                    t_flag[t] = 0
                else:
                    break

    return [i for i in xrange(1, modulus) if t_flag[i]]


def get_decimation_value(p, q):
    # returns t such that p decimated by t is q
    import sympy
    from ._binarypoly import poly_decimation

    # quick checks
    if p == q:
        return 1
    if sympy.degree(p) != sympy.degree(q):
        return None

    degree = sympy.degree(p)
    modulus = 2 ** degree - 1
    t_values = [i for i in get_representatives(degree) if sympy.gcd(i, modulus) == 1]

    for t in t_values:
        if t == 1:
            continue
        if poly_decimation(p, t) == q:
            return t

    # if all else fails...
    return None


def retrieve_zech_log(p):
    import os
    import pickle
    import bz2
    import sympy

    degree = sympy.degree(p)
    modulus = 2**degree - 1

    file_path = os.path.abspath(os.path.dirname(__file__))
    file_path = os.path.join(file_path, 'zechdata')
    file_path = os.path.join(file_path, 'zechdata_{}'.format(degree))

    try:
        with bz2.BZ2File(file_path, 'r') as f:
            s = pickle.load(f)
            base_poly = s['poly']
            base_data = s['data']
    except IOError:
        raise IOError('data non-existent, call fetch_and_save({})'.format(degree))

    if p == base_poly:
        return base_data
    else:
        new_data = ZechContainer(degree)

        t_val = get_decimation_value(base_poly, p)
        t_inv = sympy.mod_inverse(t_val, modulus)

        for r in base_data:
            new_data[r] = t_inv * base_data[t_val * r % modulus] % modulus

        return new_data


def fetch_and_save(n):
    import os
    import pickle
    import itertools
    import bz2
    from ._binarypoly import generate_primitives

    import requests
    from lxml import etree

    p = generate_primitives(n, limit=1)[0]
    data = ZechContainer(n)
    reps = get_representatives(n)

    poly_string = str(p.as_expr()).replace('**', '^')
    uri = "http://magma.maths.usyd.edu.au/xml/calculator.xml"

    for sublist in itertools.izip_longest(*[iter(reps)]*200):
        m_input = ("P<x>:=PolynomialRing(GF(2));\n"
                   "f := {};\n".format(poly_string) +
                   "K<w> := ExtensionField<GF(2), x|f>;\n"
                   "Init := {};\n".format([a for a in sublist if a is not None]) +
                   "for m := 1 to #Init do\n"
                   "    z := ZechLog(K, Init[m]);\n"
                   "    printf \"%o,%o,\", Init[m], z;\n"
                   "end for;")
        result = None
        while not result:
            r = requests.post(uri, data={'input': m_input})
            result = etree.fromstring(r.text).xpath('//line')
            if result:
                result = result[0]
            else:
                continue
            result = zip(*[iter(result.text.split(','))]*2)
            for key, value in result:
                data[int(key)] = int(value)

    file_path = os.path.abspath(os.path.dirname(__file__))
    file_path = os.path.join(file_path, 'zechdata')
    file_path = os.path.join(file_path, 'zechdata_{}'.format(degree))

    # this shouldn't fail
    with bz2.BZ2File(file_path, 'w') as f:
        pickle.dump({'poly': p, 'data': data}, f, protocol=1)
