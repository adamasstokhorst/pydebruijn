import pydebruijn as pdb
import sympy as s

from sympy.abc import x

db = pdb.cyclejoin.DeBruijnZech(s.Poly(x**3+x+1, x, modulus=2), 1, s.Poly(x**6+x+1, x, modulus=2), 7)
print ''.join(map(str, db.fsr.sequence()))
print len((map(str, db.fsr.sequence())))

db = pdb.cyclejoin.DeBruijnZech(s.Poly(x**3+x+1, x, modulus=2), 1, s.Poly(x**3+x**2+1, x, modulus=2), 1)
print ''.join(map(str, db.fsr.sequence()))
print len((map(str, db.fsr.sequence())))

db = pdb.cyclejoin.DeBruijnZech(s.Poly(x**2+x+1, x, modulus=2), 1, s.Poly(x**4 + x**3+x**2+x+1, x, modulus=2), 5)
print ''.join(map(str, db.fsr.sequence()))
print len((map(str, db.fsr.sequence())))
