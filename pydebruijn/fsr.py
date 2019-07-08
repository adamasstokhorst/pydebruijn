import sympy as _sympy
from .helpers import toggle_anf as _toggle_anf, order_from_anf as _order_from_anf

__all__ = ['FeedbackShiftRegister']


class FeedbackShiftRegister(object):
    """
    Class for implementing feedback shift registers.

    This class does not check for validity of its inputs.  Users may
    invoke the constructor but must take care to provide valid inputs.
    """
    def __init__(self, anf, order=None, init_state=None):
        """
        Initializes a feedback shift register.

        Parameters
        ----------
        anf : algebraic normal form
            A SymPy expression describing the algebraic normal form of a
            feedback shift register.  Must be using integer symbols named
            `x_k`, where `k=0, 1, ...`.

        order : integer, optional (default=None)
            The order of `anf`.  If None, then it will be deduced from `anf`.

        init_state : list, optional (default=None)
            A list containing only 0s and 1s as its elements.  Must be the
            same length as the order of `func`.
        """
        # this will come from the greedy function so there's no need to check for validity
        self._anf = _sympy.Poly(anf, modulus=2).as_expr()
        self._mod_anf = _toggle_anf(anf, order)
        if order is None:
            self._order = _order_from_anf(anf)
        else:
            self._order = order
        self._state = None
        if init_state is None:
            self.state = [0] * (self.order - 1) + [1]
        else:
            self.state = init_state[:]

    @property
    def modified_anf(self):
        """
        Returns the algebraic normal form of the feedback shift register
        with the zero cycle removed, if present, and vice-versa.

        Returns
        -------
        mod_anf : algebraic normal form
            The modified algebraic normal form.

        See also
        --------
        anf
        """
        return self._mod_anf

    @property
    def modified_anf_degree(self):
        """
        Returns the degree of the modified algebraic normal form.

        Returns
        -------
        degree : integer
            The degree of the modified algebraic normal form.

        See also
        --------
        modified_anf
        """
        return _sympy.Poly(self.modified_anf, modulus=2).total_degree()

    @property
    def anf(self):
        """
        Returns the algebraic normal form of the feedback shift register.

        Returns
        -------
        anf : algebraic normal form
            The algebraic normal form as a SymPy expression.

        See also
        --------
        modified_anf
        """
        return self._anf

    @property
    def anf_degree(self):
        """
        Returns the degree of the algebraic normal form.

        Returns
        -------
        degree : integer
            The degree of the algebraic normal form.

        See also
        --------
        anf
        """
        return _sympy.Poly(self.anf, modulus=2).total_degree()

    @property
    def order(self):
        """
        Returns the order of the algebraic normal form.

        The order of a feedback shift register is also
        the size of its state.

        Returns
        -------
        order : integer
            The order of the algebraic normal form.
        """
        return self._order

    @property
    def state(self):
        """
        Returns the current state of the algebraic normal form.

        Returns
        -------
        state : list
            The state of the algebraic normal form. Its elements
            are always only 0s and 1s.
        """
        return self._state

    @state.setter
    def state(self, iterable):
        """
        Sets the current state of the algebraic normal form.

        Parameters
        ----------
        iterable : iterable object
            The iterable is passed through a filter that converts all
            truthy values to 1 and falsy values to 0.  If its length is
            greater than the order of the FSR, it will be trimmed.  If
            its length is less than the order of the FSR, it will be
            padded with 0s at the end.
        """
        s = map(lambda x: 1 if x else 0, iterable)
        if len(s) < self._order:
            s += [0] * (self._order - len(s))
        elif len(s) > self._order:
            s = s[:self._order]
        self._state = s

    def advance(self):
        """
        Advances the current state of feedback shift register by one step.

        See also
        --------
        next_state
        """
        # this is really expensive for some reason
        args = zip(_sympy.symbols('x_:{}'.format(self.order), integer=True), self.state)
        self.state = self.state[1:] + [self.anf.subs(args) % 2]

    def next_state(self):
        """
        Returns the state after the current state of the
        feedback shift register.

        This method modifies the current state of the feedback shift
        register as a side-effect.

        Returns
        -------
        next_state : list
            The next state of the feedback shift register.

        See also
        --------
        advance
        """
        self.advance()
        return self.state

    def bit(self):
        """
        Returns the next bit of feedback shift register.

        This method modifies the current state of the feedback shift
        register as a side-effect.

        This method is synonymous with `fsr.next_state()[-1]`.

        Returns
        -------
        next_bit : int
            The next bit of the feedback shift register.

        See also
        --------
        next_state
        """
        self.advance()
        return self.state[-1]

    def sequence(self):
        """
        Returns the full sequence of the feedback shift register.

        Do keep the order of the feedback shift register in mind,
        as this method will try to return the entire sequence,
        no matter how big it is.

        Returns
        -------
        sequence : list
            The full sequence of the feedback shift register.
        """
        cur_state = self.state[:]
        seq = cur_state[:]
        self.advance()
        while self.state != cur_state:
            seq.append(self.state[-1])
            self.advance()
        return seq[:-(self.order-1)]
