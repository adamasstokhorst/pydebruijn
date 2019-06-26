import sympy
from ._misc import powerset

__all__ = ['padded_bin', 'seq_from_truth_table', 'truth_table_from_seq', 'anf_from_truth_table']


def padded_bin(k, n):
    """
    Converts a number to binary form, left-padded to a given length.

    This function might not work correctly with negative integers.

    Parameters
    ----------
    k : integer
        The number to be converted to binary.

    n : integer
        The target length of the binary string. The returned string
        will be at least of this length.

    Returns
    -------
    s : string
        The binary string, which is at least `n` characters long.
    """
    s = bin(k)
    return '0' * (n - len(bin(k)) + 2) + s[2:]


# TODO: anything dealing with sequences should assume it to be a list, not a string.
def seq_from_truth_table(table, order):
    """
    Computes the binary sequence string from a given truth table.

    Currently, this function assumes that the truth table describes
    a de Bruijn sequence.

    Parameters
    ----------
    table : dict
        A mapping from binary strings of length `n` to 0 or 1.

    order : integer
        The order of the boolean function.

    Returns
    -------
    seq : string
        The binary sequence represented by the truth table.
    """
    seq = '0' * order
    while len(seq) < 2**order:
        seq += str(table[seq[-(order-1):]] ^ int(seq[-order]))
    return seq


def truth_table_from_seq(seq, order):
    """
    Converts a de Bruijn sequence string to truth table form.

    Parameters
    ----------
    seq : string
        A de Bruijn sequence as a binary string.

    order : integer
        The order of the de Bruijn sequence.

    Returns
    -------
    table : string
        The truth table form of the underlying boolean function.
    """
    table = {}
    s = seq + seq[:order+1]
    for i in xrange(2**order):
        if seq[i+1:i+order] not in table:
            table[s[i+1:i+order]] = int(s[i]) ^ int(s[i+order])

        if len(table) == 2**(order - 1):
            return table


def anf_from_truth_table(table, order):
    """
    Computes the algebraic normal form from a given truth table.

    Parameters
    ----------
    table : dict
        A mapping from binary strings of length `n` to 0 or 1.

    order : integer
        The order of the boolean function.

    Returns
    -------
    expr : SymPy expression
        The algebraic normal form of the boolean function.
    """
    sym = [sympy.Symbol('x_{}'.format(a), integer=True) for a in range(order)]
    anf = sym[0]
    for s in table:
        if table[s]:
            zero_index = [a+1 for a in range(order-1) if s[a] == '0']
            one_index = [a+1 for a in range(order-1) if s[a] == '1']
            init = reduce(lambda a, b: a * b, map(lambda a: sym[a], one_index), 1)

            for subset in powerset(zero_index, reverse=False):
                anf += init * reduce(lambda a, b: a * b, map(lambda a: sym[a], list(subset)), 1)
    return sympy.poly(anf, modulus=2).as_expr()
