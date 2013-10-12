#! /usr/bin/env ipython

import sys
import pickle
import sympy

ip = get_ipython()

for order in [int(o) for o in sys.argv[1:]]:
    print("Calculating with polynomials to order {0}".format(order))
    Num = sympy.var('Num:{0}'.format(order), real=True)
    Den = sympy.var('Den:{0}'.format(order), real=True)
    PolynomialVariable = sympy.Symbol('PolynomialVariable', real=True)
    p_Num = sum(Num[i]*PolynomialVariable**i for i in range(order))
    p_Den = sum(Den[i]*PolynomialVariable**i for i in range(order))
    ip.magic("time p_Ratio = sympy.series(p_Num/p_Den,x=PolynomialVariable,x0=0,n=order)")
    pickle.dump(p_Ratio, file('PolynomialRatioSeries_Order{0}.dat'.format(order),'w'))

# Load later, in a different python session, with:
# p_Ratio = pickle.load(file('PolynomialRatioSeries_Order{0}.dat'.format(order)))
