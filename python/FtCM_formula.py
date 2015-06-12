
def multisets( n, k ):
    if k == 0: return [(0,)*n]
    if n == 0: return []
    if n == 1: return [(k,)]
    return [(0,)+val for val in multisets(n-1,k)] + \
            [(val[0]+1,)+val[1:] for val in multisets(n,k-1)]

def submultisets( m, k ):        
    msets = multisets( len(m),k)
    smsets = []
    for mset in msets: 
        if (all((e_s <= e_m) for e_s, e_m in zip(mset, m))):
            smsets.append(mset)
    return smsets        

print multisets(1,2)
print multisets(2,3)
print multisets(3,3)
print submultisets((1,1,1),3)
print submultisets((3,3,3),3)


        
