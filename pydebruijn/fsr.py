import sympy as _sympy
from .helpers import toggle_anf as _toggle_anf, order_from_anf as _order_from_anf

__all__ = ['FeedbackShiftRegister']


class FeedbackShiftRegister(object):
    def __init__(self, anf, order=None, init_state=None):
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
        return self._mod_anf

    @property
    def modified_anf_degree(self):
        return _sympy.Poly(self.modified_anf, modulus=2).total_degree()

    @property
    def anf(self):
        return self._anf

    @property
    def anf_degree(self):
        return _sympy.Poly(self.anf, modulus=2).total_degree()

    @property
    def order(self):
        return self._order

    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, iterable):
        s = map(lambda x: 1 if x else 0, iterable)
        if len(s) < self._order:
            s += [0] * (self._order - len(s))
        elif len(s) > self._order:
            s = s[:self._order]
        self._state = s

    def advance(self):
        # this is really expensive for some reason
        args = zip(_sympy.symbols('x_:{}'.format(self.order), integer=True), self.state)
        self.state = self.state[1:] + [self.anf.subs(args) % 2]

    def next_state(self):
        self.advance()
        return self.state

    def bit(self):
        self.advance()
        return self.state[-1]

    def sequence(self):
        cur_state = self.state[:]
        seq = cur_state[:]
        self.advance()
        while self.state != cur_state:
            seq.append(self.state[-1])
            self.advance()
        return seq[:-(self.order-1)]
