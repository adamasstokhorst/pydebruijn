import pydebruijn as pdb
from pprint import pprint

polys = '111',  #'11111', '1001001', '1010111', '1110101', '100011011'

for poly in polys:
    print 'Polynomial: 1011, 11111', poly

    print '-- Zechlog based --'
    db = pdb.cyclejoin.DeBruijnZech('1011', '11111', poly)
    pprint(db._p_matrix)
    pprint(db.adjacency_matrix)
    print db.adjacency_matrix.minor(0, 0)
    for _ in range(0):
        s = db.fsr.sequence()
        print len(s), ''.join(map(str, s))
        db.next_sequence()
    print
    
    print '-- ProdIrreds based --'
    db = pdb.cyclejoin.DeBruijnPoly('1011', '11111', poly)
    pprint(db._p_matrix)
    pprint(db.adjacency_matrix)
    print db.adjacency_matrix.minor(0, 0)
    for _ in range(0):
        s = db.fsr.sequence()
        print len(s), ''.join(map(str, s))
        db.next_sequence()
    print
    
