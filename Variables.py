# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Set up python

# <codecell>

# Make sure division of integers does not round to the nearest integer
from __future__ import division

# Make everything in python's symbolic math package available
from sympy import * # Make sure sympy functions are used in preference to numpy
import sympy # Make sympy. constructions available
from sympy import Rational as frac # Rename for similarity to latex
from sympy import log as ln

# Print symbolic expressions nicely
init_printing()

# We'll use the numpy `array` object for vectors
from numpy import array, cross, dot

# <headingcell level=1>

# Define the basic variables

# <markdowncell>

# Only the most basic variables should be defined here.

# <codecell>

# Dimensionful quantities, just in case anybody uses them...
G,c = symbols('G,c', real=True)

# Unit orbital angular velocity vector ("Newtonian" angular momentum)
Lhat_N, Lhat_Nx, Lhat_Ny, Lhat_Nz = symbols('Lhat_N, Lhat_Nx, Lhat_Ny, Lhat_Nz', real=True)

# Unit separation vector between the compact objects
nhat, nhat_x, nhat_y, nhat_z = symbols('nhat, nhat_x, nhat_y, nhat_z', real=True)

# Masses
m1, m2 = symbols('m1, m2', real=True)

# Spin vectors
chi1, chi1_l, chi1_n, chi1_la = symbols('chi1, chi1_l, chi1_n, chi1_lambda', real=True)
chi2, chi2_l, chi2_n, chi2_la = symbols('chi2, chi2_l, chi2_n, chi2_lambda', real=True)

# Tidal deformabilities
lambda1, lambda2 = symbols('lambda1, lambda2', real=True)

# v = x**(1/2) = Omega_orb**(1/3)
v = symbols('v', real=True)

# <headingcell level=1>

# Derived variables

# <markdowncell>

# Any variable that can be derived from the variables above should be put in this section.
# 
# These variables should probably be left in arbitrary form, unless a particular simplification is desired.  The `BasicSubstitutions` dictionary should map between the general names and their definitions in terms of basic variables.  In numerical codes, their values can be calculated once and then stored, so that the values do not have to be re-calculated every time they appear in an expression.

# <codecell>

BasicSubstitutions = {} # For now, just initialize the dictionary

# <markdowncell>

# Various common combinations of the two masses:

# <codecell>

m = symbols('m', real=True);
BasicSubstitutions[m] = m1+m2

delta = symbols('delta', real=True);
BasicSubstitutions[delta] = (m1-m2)/(m1+m2)

nu = symbols('nu', real=True);
BasicSubstitutions[nu] = m1*m2/(m1+m2)

q = symbols('q', real=True);
BasicSubstitutions[q] = m1/m2

# <markdowncell>

# The system's vector basis is given by $(\hat{L}_{\text{N}}, \hat{n}, \hat{\lambda})$.  Here, we define the remaining element, and give the relevant substitutions in terms of Cartesian basis elements.

# <codecell>

lambdahat, lambdahat_x, lambdahat_y, lambdahat_z = symbols('lambdahat, lambdahat_x, lambdahat_y, lambdahat_z', real=True)

BasicSubstitutions[Lhat_N] = array([Lhat_Nx, Lhat_Ny, Lhat_Nz])
BasicSubstitutions[nhat] = array([nhat_x, nhat_y, nhat_z])
BasicSubstitutions[lambdahat] = cross(BasicSubstitutions[Lhat_N], BasicSubstitutions[nhat])

BasicSubstitutions[lambdahat_x] = BasicSubstitutions[lambdahat][0]
BasicSubstitutions[lambdahat_y] = BasicSubstitutions[lambdahat][1]
BasicSubstitutions[lambdahat_z] = BasicSubstitutions[lambdahat][2]

# <markdowncell>

# Various spin components and combinations:

# <codecell>

BasicSubstitutions[chi1] = array([chi1_l, chi1_n, chi1_la])
BasicSubstitutions[chi2] = array([chi2_l, chi2_n, chi2_la])

chi1chi1, chi1chi2, chi2chi2 = symbols('chi1chi1, chi1chi2, chi2chi2', real=True)
BasicSubstitutions[chi1chi1] = dot(BasicSubstitutions[chi1], BasicSubstitutions[chi1])
BasicSubstitutions[chi1chi2] = dot(BasicSubstitutions[chi1], BasicSubstitutions[chi2])
BasicSubstitutions[chi2chi2] = dot(BasicSubstitutions[chi2], BasicSubstitutions[chi2])

S, S_l, S_n, S_la = symbols('S, S_l, S_n, S_lambda', real=True)
BasicSubstitutions[S] = ((m1+m2)**2)*chi1*((1+delta)**2/4) + chi2*((1-delta)**2/4)
BasicSubstitutions[S_l] = ((m1+m2)**2)*chi1_l*((1+delta)**2/4) + chi2_l*((1-delta)**2/4)
BasicSubstitutions[S_n] = ((m1+m2)**2)*chi1_n*((1+delta)**2/4) + chi2_n*((1-delta)**2/4)
BasicSubstitutions[S_la] = ((m1+m2)**2)*chi1_la*((1+delta)**2/4) + chi2_la*((1-delta)**2/4)

Sigma, Sigma_l, Sigma_n, Sigma_la = symbols('Sigma, Sigma_l, Sigma_n, Sigma_lambda', real=True)
BasicSubstitutions[Sigma] = ((m1+m2)**2)*(chi2*(1-delta)/2 - chi1*(1+delta)/2)
BasicSubstitutions[Sigma_l] = ((m1+m2)**2)*(chi2_l*(1-delta)/2 - chi1_l*(1+delta)/2)
BasicSubstitutions[Sigma_n] = ((m1+m2)**2)*(chi2_n*(1-delta)/2 - chi1_n*(1+delta)/2)
BasicSubstitutions[Sigma_la] = ((m1+m2)**2)*(chi2_la*(1-delta)/2 - chi1_la*(1+delta)/2)

# <markdowncell>

# Other powers of the angular velocity that find frequent use:

# <codecell>

x, Omega_orb = symbols('x, Omega_orb', real=True)
BasicSubstitutions[x] = v**2
BasicSubstitutions[Omega_orb] = v**3

