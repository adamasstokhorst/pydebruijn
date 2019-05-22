import sympy
from ._misc import powerset

__all__ = ['padded_bin', 'seq_from_truth_table', 'truth_table_from_seq', 'anf_from_truth_table']


def padded_bin(k, n):
    """Returns a number in binary, padded with zeroes to a given length.
    """
    s = bin(k)
    return '0' * (n - len(bin(k)) + 2) + s[2:]


# TODO: anything dealing with sequences should assume it to be a list, not a string.
def seq_from_truth_table(func_dict, n):
    """Computes the corresponding binary sequence string from a given truth table.
    """
    seq = '0' * n
    while len(seq) < 2**n:
        seq += str(func_dict[seq[-(n-1):]] ^ int(seq[-n]))
    return seq


def truth_table_from_seq(seq, degree):
    """Converts de Bruijn sequence string to its underlying truth table.
    """
    table = {}
    s = seq + seq[:degree+1]
    for i in xrange(2**degree):
        if seq[i+1:i+degree] not in table:
            table[s[i+1:i+degree]] = int(s[i]) ^ int(s[i+degree])

        if len(table) == 2**(degree - 1):
            return table


def anf_from_truth_table(table, degree):
    """Computes the algebraic normal form corresponding to a given truth table.
    """
    sym = [sympy.Symbol('x_{}'.format(a), integer=True) for a in range(degree)]
    anf = sym[0]
    for s in table:
        if table[s]:
            zero_index = [a+1 for a in range(degree-1) if s[a] == '0']
            one_index = [a+1 for a in range(degree-1) if s[a] == '1']
            init = reduce(lambda a, b: a * b, map(lambda a: sym[a], one_index), 1)

            for subset in powerset(zero_index, reverse=False):
                anf += init * reduce(lambda a, b: a * b, map(lambda a: sym[a], list(subset)), 1)
    return sympy.poly(anf, modulus=2).as_expr()
