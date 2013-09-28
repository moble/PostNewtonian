# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Set up python

# <codecell>

# Make sure division of integers does not round to the nearest integer
from __future__ import division

# Make everything in python's symbolic math package available
from sympy import *

# Print symbolic expressions nicely
init_printing()

# <headingcell level=1>

# Define the basic variables

# <markdowncell>

# Only the most basic variables should be defined here.

# <codecell>

m1, m2, v = symbols('m1, m2, v', real=True)
chi1, chi2, ell = symbols('chi1, chi2, ell', seq=True, real=True)

# <headingcell level=1>

# Derived variables

# <markdowncell>

# Any variable that can be derived from the variables above should be put in this section.

# <codecell>

# These are just combinations of the two masses
m = symbols('m', real=True) # = m1+m2
delta = symbols('delta', real=True) # = (m1-m2)/m
nu = symbols('nu', real=True) # = m1*m2/m
q = symbols('q', real=True) # = m1/m2

# These are just components of the spin vectors
chi1_ell, chi1_n, chi1_la = symbols('chi1_ell, chi1_n, chi1_la', real=True)
chi2_ell, chi2_n, chi2_la = symbols('chi2_ell, chi2_n, chi2_la', real=True)

# These are combinations of the individual spin vectors
S = symbol('S', seq=True, real=True) #  =  S_1 + S_2  =  chi1*(1+delta)/2 + chi2*(1-delta)/2
Sigma = symbol('Sigma', seq=True, real=True) #  =  Sigma_1 + S_2  =  chi1*(1+delta)/2 + chi2*(1-delta)/2
S_ell, S_n, S_la = symbols('S_ell, S_n, S_la', real=True)
Sigma_ell, Sigma_n, Sigma_la = symbols('Sigma_ell, Sigma_n, Sigma_la', real=True)

