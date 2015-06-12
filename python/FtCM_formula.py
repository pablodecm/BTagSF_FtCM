from sympy import Symbol, IndexedBase, Product
from sympy import binomial, sympify
from ast import literal_eval
import json

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
    t = Symbol('t') # type of jet
    C = IndexedBase('C') # componenent number vector 
    N = IndexedBase('N') # max componenent number vector
    E = IndexedBase('E') # efficiency vector 
    product =  Product(binomial(N[t],C[t])*E[t]**(C[t])*(1-E[t])**(N[t]-C[t]), (t,0,len(comp)-1))
    C_dict = { C[i] : comp[i] for i in range(n_t) }
    N_dict = { N[i] : m_comp[i] for i in range(n_t) }
    subs_dict = dict(C_dict.items() + N_dict.items())
    return product.doit().subs(subs_dict)

def prob_submultiset( m, k):
    smset = submultiset(m, k) # get possible permutations
    comp_list = [prob_component(comp, m) for comp in smset] 
    return sum(comp_list) # sum the prob of all components 

def load_cat_counts(filename):
    with open(filename) as f:
        return { literal_eval(k) : v for k,v in json.load(f).items()} 

def formula_categories(filename, k):    
    cat_counts = load_cat_counts(filename)
    f_cats = sympify(0)
    for cat, count in cat_counts.items():
        f_cats += count[0]*prob_submultiset(cat, k)
    return f_cats.expand()

    
    



        
