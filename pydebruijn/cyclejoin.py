import sympy as _sympy
import networkx as _nx
import itertools as _iters
import collections as _collections
from sympy.abc import x as _x
from .helpers import (poly_decimation as _poly_decimation,
                      is_primitive as _is_primitive,
                      get_associate_poly as _get_associate_poly,
                      seq_decimation as _seq_decimation,
                      lfsr_from_poly as _lfsr_from_poly,
                      spanning_trees as _spanning_trees,
                      retrieve_zech_log as _retrieve_zech_log,
                      get_special_state as _get_special_state)
from .fsr import FeedbackShiftRegister as _FSR

__all__ = ['DeBruijnPoly', 'DeBruijnZech']


class DeBruijnPoly(object):
    def __init__(self, *args):
        """
        Initializes a de Bruijn sequence generator.

        This class makes use of irreducible binary polynomials to generate
        de Bruijn sequences.

        Parameters
        ----------
        args : binary string(s) (required)
            Polynomials to be used to generate de Bruijn sequences.
            Input polynomials must be irreducible; reducible polynomials
            are silently ignored.  Coefficients are given in decreasing
            power -- e.g. `x**3 + x + 1` is written as `1011`.
        """
        if not args:
            raise ValueError('no arguments passed (at least 1 expected)')

        # properties modifiable by user
        self._state = None

        # properties that are read-only
        self._polys = []
        self._states = []
        self._associates = []
        self._graph = _nx.MultiGraph()
        self._poly = None
        self._order = 0
        self._p_matrix = None
        self._adjacency_matrix = None
        self._param_generator = None
        self._fsr = None

        for binary_string in args:
            binary_seq = map(lambda a: 0 if a == '0' else 1, binary_string)
            proto_poly = reduce(lambda a, b: a * _x + b, binary_seq, 0)
            poly = _sympy.Poly(proto_poly, _x, modulus=2)
            if poly.is_irreducible and poly not in self._polys:
                self._polys.append(poly)
        if not self._polys:
            raise ValueError('no irreducible polynomial supplied.')

        self._polys.sort(key=_sympy.degree)
        self._poly = reduce(lambda a, b: a * b, self._polys)
        self._order = self._poly.degree()

        self._state = [0] * self._order
        self._sym = _sympy.symbols('x_:{}'.format(self._order), integer=True)

        self.__initialize()

    def __initialize(self):
        """
        Method for actually initializing the de Bruijn sequence generator.

        Users need not to run this method.
        """
        # populate states
        self._associates = _get_associate_poly(self._polys)
        for entry in self._associates:
            entry_state = []
            degree = _sympy.degree(entry['associate'])
            init_state = [1] * degree
            for i in xrange(entry['order']):
                entry_state.append(_seq_decimation(entry['associate'], entry['order'], i, init_state)[:degree])
            entry_state.append([0] * degree)
            self._states.append(entry_state)

        # find special state
        p_matrix = []
        for poly in self._polys:
            degree = _sympy.degree(poly)
            for i in xrange(degree):
                state = [0] * degree
                state[i] = 1
                for j in xrange(self._order - degree):
                    state.append(_lfsr_from_poly(poly, state[-degree:])[-1])
                p_matrix += state
        p_matrix = _sympy.Matrix(self._order, self._order, p_matrix)
        self._p_matrix = p_matrix
        special_state = map(int, _sympy.Matrix(1, self._order, [1] + [0] * (self._order - 1)) * p_matrix.inv_mod(2))
        special_states = []
        i = 0
        for poly in self._polys:
            special_states.append(special_state[i:i + _sympy.degree(poly)])
            i += _sympy.degree(poly)

        # find viable pairs
        # notes: all_pairs = list of dictionary for each polynomial, where the keys are pairs of states,
        #                    and the entries are the corresponding shifts for the states.
        # this code looks really messy, is there a better way to do this?
        all_pairs = []
        for i, entry in enumerate(self._associates):
            cur_pairs = _collections.defaultdict(list)
            special_list = [special_states[i]]
            for _ in xrange(entry['period'] - 1):
                special_list.append(_lfsr_from_poly(entry['poly'], special_list[-1]))
            for state_1 in xrange(entry['order'] + 1):
                added_state = map(lambda a: map(lambda b, c: b ^ c, self._states[i][state_1], a), special_list)
                for state_2 in xrange(state_1, entry['order'] + 1):
                    cur_state = self._states[i][state_2][:]
                    for shift_2 in xrange(entry['period'] if state_2 != entry['order'] else 1):
                        if cur_state in added_state:
                            shift_1 = added_state.index(cur_state)
                            if state_2 == entry['order']:
                                shift_2 = shift_1
                            cur_pairs[(state_1, state_2)].append((-shift_1 % entry['period'],
                                                                  (shift_2 - shift_1) % entry['period']))
                            if state_1 != state_2:
                                cur_pairs[(state_2, state_1)].append(((shift_2 - shift_1) % entry['period'],
                                                                      -shift_1 % entry['period']))
                        cur_state = _lfsr_from_poly(entry['poly'], cur_state)
            all_pairs.append(cur_pairs)

        # find conjugate pairs and construct adjacency graph
        graph = self._graph
        for param_1, param_2, shifts_1, shifts_2 in self.__conjugate_pair_generator(all_pairs):
            if param_1 < param_2:  # was != !!
                graph.add_edge(param_1, param_2, shift={param_1: shifts_1, param_2: shifts_2})

        # clean up
        del all_pairs

        self._adjacency_matrix = -_sympy.Matrix(_nx.to_numpy_matrix(self._graph)).applyfunc(int)
        for i in range(self._adjacency_matrix.rows):
            self._adjacency_matrix[i, i] -= sum(self._adjacency_matrix[i, :])

        simple_graph = _nx.Graph(self._graph)
        self._param_generator = self.__param_generator(_spanning_trees(simple_graph))

        anf = self.__get_algebraic_normal_form(*self._param_generator.next())
        self._fsr = _FSR(anf, order=self._order, init_state=self._state)

    # this could be done better?
    def __conjugate_pair_generator(self, all_pairs):
        """
        Method for generating conjugate pairs.

        Users need not to run this method.

        See also
        --------
        __initialize
        """
        # notes: pairs_product = cartesian product of all the pairs of states for each polynomials,
        #                        i.e. it's a tuple of size s, where s is the number of polynomials
        #        pairs = a tuple of size s, entry at index i is a state-pair for i-th polynomial.
        pairs_product = _iters.product(*map(lambda a: a.keys(), all_pairs))
        for pairs in pairs_product:
            for shifts in _iters.product(*[all_pairs[i][pair] for i, pair in enumerate(pairs)]):
                # there's probably a more efficient way of doing this
                periods_1 = [entry['period'] if pairs[i][0] != entry['order'] else 1
                             for i, entry in enumerate(self._associates)]
                periods_2 = [entry['period'] if pairs[i][1] != entry['order'] else 1
                             for i, entry in enumerate(self._associates)]
                found = True
                for l1 in _iters.product(*[range(a) for a in periods_1]):
                    for l2 in _iters.product(*[range(a) for a in periods_2]):
                        found = True
                        for (i, j) in _iters.combinations(range(len(self._polys)), 2):
                            if _sympy.gcd(periods_1[i], periods_2[j]) != 1:
                                diff_1 = shifts[i][0] - shifts[j][0]
                                diff_2 = shifts[i][1] - shifts[j][1]
                                found &= (diff_1 - l1[i] + l1[j]) % _sympy.gcd(periods_1[i], periods_2[j]) == 0
                                found &= (diff_2 - l2[i] + l2[j]) % _sympy.gcd(periods_1[i], periods_2[j]) == 0
                                if not found:
                                    break
                        if found:
                            # we can reduce number of pairs by messing around in here, possibly
                            x_list_1 = []
                            x_list_2 = []
                            for i in range(len(self._polys)):
                                x_list_1.append(_sympy.gcd(periods_1[i], reduce(_sympy.lcm, periods_1[:i], 1)))
                                x_list_2.append(_sympy.gcd(periods_2[i], reduce(_sympy.lcm, periods_2[:i], 1)))
                            eff_shift_1 = [(shifts[i][0] - shifts[0][0]) % x_list_1[i]
                                           for i in range(len(self._polys))]
                            eff_shift_2 = [(shifts[i][1] - shifts[0][1]) % x_list_2[i]
                                           for i in range(len(self._polys))]
                            param_1 = tuple([a[0] for a in pairs] + eff_shift_1)
                            param_2 = tuple([a[1] for a in pairs] + eff_shift_2)
                            yield param_1, param_2, \
                                  [shifts[i][0] for i in range(len(self._polys))], \
                                  [shifts[i][1] for i in range(len(self._polys))]
                            break
                    if found:
                        break

    def __param_generator(self, trees):
        """
        Method for generating sequence parameters from spanning trees.

        Users need not to run this method.

        See also
        --------
        __initialize
        """
        for tree in trees:
            param_list = map(lambda a: [self._graph.get_edge_data(*a)[k]['shift'][a[0]]
                                        for k in self._graph.get_edge_data(*a)], tree)
            for param in _iters.product(*param_list):
                yield tree, param

    def __get_algebraic_normal_form(self, tree, param):
        """
        Method for generating the algebraic normal form of the feedback shift register.

        Users need not to run this method.

        See also
        --------
        __initialize
        """
        terms = [a[0][0] for a in self._poly.terms() if a[0][0] != self._order]
        anf = sum([self._sym[a] for a in terms])

        # pass the state_set into bits, probably by storing as bound variable?
        for w, (p1, p2) in enumerate(tree):
            state = []
            cur_state = list(p1)[:len(self._polys)]
            for k in range(len(self._polys)):
                sub_state = self._states[k][cur_state[k]]
                for l in range(param[w][k]):
                    sub_state = _lfsr_from_poly(self._polys[k], sub_state)
                state += sub_state
            state = (_sympy.Matrix(1, self._order, state) * self._p_matrix).applyfunc(lambda a: a % 2)[:]

            anf += reduce(lambda a, b: a * b, [self._sym[a] + state[a] + 1 for a in range(1, self._order)])

        return _sympy.Poly(anf, modulus=2).as_expr()

    @property
    def polys(self):
        """
        Returns the list of polynomials used in constructing the sequence
        generator.
        """
        return self._polys

    @property
    def poly(self):
        """
        Returns the sequence's generating polynomial.

        Equivalent to multiplying the elements of `.polys` together.
        """
        return self._poly

    @property
    def order(self):
        """
        Returns the degree of the sequence's generating polynomial.

        Equivalent to `.poly.degree()`.
        """
        return self._order

    @property
    def adjacency_matrix(self):
        """
        Returns the adjacency matrix of the connectivity graph.
        """
        return self._adjacency_matrix

    @property
    def fsr(self):
        """
        Returns the `FeedbackShiftRegister` object corresponding to the
        current sequence.

        See also
        --------
        FeedbackShiftRegister
        """
        return self._fsr

    @property
    def state(self):
        """
        Returns the current state of the feedback shift register.
        """
        return self.fsr.state

    @state.setter
    def state(self, iterable):
        """
        Sets the current state of the feedback shift register.

        Parameters
        ----------
        iterable : any iterable object
            Replace the current state with this iterable.  The resulting
            state may not be the same as the given iterable.

        See also
        --------
        FeedbackShiftRegister
        """
        self.fsr.state = iterable

    def next_sequence(self):
        """
        Changes the parameters of the generator to generate a different
        sequence.

        Raises
        ------
        StopIteration
            If the connectivity graph has yielded all possible
            spanning trees.
        """
        # Method will raise StopIteration when sequences are exhausted, don't forget to handle it.
        anf = self.__get_algebraic_normal_form(*self._param_generator.next())
        self._fsr = _FSR(anf, order=self._order, init_state=self._state)


class DeBruijnZech(object):
    def __init__(self, p, t):
        """
        Initializes a de Bruijn sequence generator.

        This class makes use of Zech logarithms to generate
        de Bruijn sequences.

        Parameters
        ----------
        p : Sympy.Poly object (required)
            A primitive binary polynomial.

        t : integer (required)
            The decimation value to be applied to `p`.  If `p` is of
            degree `n`, this value must necessarily divide `2**n - 1`,
            but this is not a sufficient requirement.

        Raises
        ------
        ValueError
            If `p` is not primitive, or if `t` does not fulfill
            a certain criterion.
        """
        # properties modifiable by user
        self._state = None

        # properties that are read-only
        self._states = []
        self._graph = _nx.MultiGraph()
        self._associate = p
        self._poly = None
        self._t_value = t
        self._order = 0
        self._adjacency_matrix = None
        self._param_generator = None
        self._fsr = None

        self._order = self._associate.degree()
        if not _is_primitive(p):
            raise ValueError('polynomial not primitive: {}'.format(p.as_expr()))
        if (2 ** self._order - 1) % t != 0 or _poly_decimation(p, t) is None:
            raise ValueError('inappropriate t-value: {}'.format(t))

        self._state = [0] * self._order
        self._sym = _sympy.symbols('x_:{}'.format(self._order), integer=True)
        self._poly = _poly_decimation(p, t)

        self.__initialize()

    def __initialize(self):
        """
        Method for actually initializing the de Bruijn sequence generator.

        Users need not to run this method.
        """
        # populate states
        degree = self._order
        init_state = _get_special_state(self._associate, self._t_value)
        for i in xrange(self._t_value):
            self._states.append(_seq_decimation(self._associate, self._t_value, i, init_state)[:degree])
        self._states.append([0] * degree)

        zech_log = _retrieve_zech_log(self._associate)

        # find conjugate pairs and construct adjacency graph
        graph = self._graph
        for z1, z2 in ((a, zech_log[a]) for a in xrange(1, 2**degree - 1)):
            if z1 < z2:
                param_1, param_2 = (z1 % self._t_value,), (z2 % self._t_value,)
                shifts_1, shifts_2 = (z1 / self._t_value,), (z2 / self._t_value,)
                graph.add_edge(param_1, param_2, shift={param_1: shifts_1, param_2: shifts_2})

        self._adjacency_matrix = -_sympy.Matrix(_nx.to_numpy_matrix(self._graph)).applyfunc(int)
        for i in range(self._adjacency_matrix.rows):
            self._adjacency_matrix[i, i] -= sum(self._adjacency_matrix[i, :])

        simple_graph = _nx.Graph(self._graph)
        self._param_generator = self.__param_generator(_spanning_trees(simple_graph))

        # if auto_arm:
        #     anf = self.__get_algebraic_normal_form(*self._param_generator.next())
        #     self._fsr = _FSR(anf, order=self._order, init_state=self._state)

    def __param_generator(self, trees):
        """
        Method for generating sequence parameters from spanning trees.

        Users need not to run this method.

        See also
        --------
        __initialize
        """
        for tree in trees:
            param_list = map(lambda a: [self._graph.get_edge_data(*a)[k]['shift'][a[0]]
                                        for k in self._graph.get_edge_data(*a)], tree)
            for param in _iters.product(*param_list):
                yield tree, param

    def __get_algebraic_normal_form(self, tree, param):
        """
        Method for generating the algebraic normal form of the feedback shift register.

        Users need not to run this method.

        See also
        --------
        __initialize
        """
        from .helpers import powerset

        terms = [a[0][0] for a in self._poly.terms() if a[0][0] != self._order]
        anf = sum([self._sym[a] for a in terms])

        for w, (p1, p2) in enumerate(tree):
            cur_state = list(p1)
            state = self._states[cur_state[0]][:]
            for l in range(param[w][0]):
                state = _lfsr_from_poly(self._poly, state)
            anf += reduce(lambda a, b: a * b, [self._sym[a] + state[a] + 1 for a in range(1, self._order)])

        for subset in powerset(range(1, self._order), reverse=False):
            anf += reduce(lambda a, b: a * b, map(lambda a: self._sym[a], list(subset)), 1)

        return _sympy.Poly(anf, modulus=2).as_expr()

    @property
    def poly(self):
        """
        Returns the sequence's generating polynomial.

        Equivalent to applying `t`-decimation to `p`.
        """
        return self._poly

    @property
    def order(self):
        """
        Returns the degree of the sequence's generating polynomial.

        Equivalent to `.poly.degree()`.
        """
        return self._order

    @property
    def adjacency_matrix(self):
        """
        Returns the adjacency matrix of the connectivity graph.
        """
        return self._adjacency_matrix

    @property
    def fsr(self):
        """
        Returns the `FeedbackShiftRegister` object corresponding to the
        current sequence.

        See also
        --------
        FeedbackShiftRegister
        """
        return self._fsr

    @property
    def state(self):
        """
        Returns the current state of the feedback shift register.
        """
        return self.fsr.state

    @state.setter
    def state(self, iterable):
        """
        Sets the current state of the feedback shift register.

        Parameters
        ----------
        iterable : any iterable object
            Replace the current state with this iterable.  The resulting
            state may not be the same as the given iterable.

        See also
        --------
        FeedbackShiftRegister
        """
        self.fsr.state = iterable

    def next_sequence(self):
        """
        Changes the parameters of the generator to generate a different
        sequence.

        Raises
        ------
        StopIteration
            If the connectivity graph has yielded all possible
            spanning trees.
        """
        # Method will raise StopIteration when sequences are exhausted, don't forget to handle it.
        anf = self.__get_algebraic_normal_form(*self._param_generator.next())
        self._fsr = _FSR(anf, order=self._order, init_state=self._state)

