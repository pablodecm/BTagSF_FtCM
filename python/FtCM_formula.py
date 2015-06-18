try:
    from symengine import Symbol, sympify
    from symengine.lib.symengine_wrapper import binomial
except ImportError:
    print "symengine (C++ CAS) not found, using sympy"
    from sympy import Symbol, sympify, binomial

def multiset( n, k ):
    if k == 0: return [(0,)*n]
    if n == 0: return []
    if n == 1: return [(k,)]
    return [(0,)+val for val in multiset(n-1,k)] + \
            [(val[0]+1,)+val[1:] for val in multiset(n,k-1)]

def submultiset( m, k ):        
    mset = multiset( len(m),k)
    smset = []
    for comp in mset: 
        if (all((e_c <= e_m) for e_c, e_m in zip(comp, m))):
            smset.append(comp)
    return smset        

def prob_component( comp, m_comp):
    n_t = len(comp) # number of types
    effs = [Symbol("eff_{}".format(i)) for i in range(n_t)]
    product = sympify(1)
    for t in range(n_t): 
        c = comp[t]
        m = m_comp[t]
        product *= binomial(m,c)*effs[t]**c*(1-effs[t])**(m-c)
    return product

def prob_submultiset( m, k):
    smset = submultiset(m, k) # get possible permutations
    comp_list = [prob_component(comp, m) for comp in smset] 
    return sum(comp_list) # sum the prob of all components 

def formula_categories(cat_fractions, k):    
    f_cats = sympify(0)
    for cat, frac in cat_fractions.items():
        f_cats += frac*prob_submultiset(cat, k)
    return f_cats.expand()

    
    



        
