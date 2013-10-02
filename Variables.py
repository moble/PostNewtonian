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

# We'll use an orderd dictionary to satisfy dependencies
from collections import OrderedDict

# <headingcell level=1>

# Fundamental variables

# <markdowncell>

# Only the most basic variables should be defined here.

# <codecell>

# Dimensionful quantities, just in case anybody uses them...
G,c = symbols('G,c', real=True)

# Unit orbital angular velocity vector ("Newtonian" angular momentum)
Lhat_Nx, Lhat_Ny, Lhat_Nz = symbols('Lhat_Nx, Lhat_Ny, Lhat_Nz', real=True)

# Unit separation vector between the compact objects
nhat_x, nhat_y, nhat_z = symbols('nhat_x, nhat_y, nhat_z', real=True)

# Mass of object 1.  Note that the total mass is assumed to be 1, so m2 is derived
m1 = symbols('m1', real=True)

# Spin vectors
chi1_x, chi1_y, chi1_z = symbols('chi1_x, chi1_y, chi1_z', real=True)
chi2_x, chi2_y, chi2_z = symbols('chi2_x, chi2_x, chi2_z', real=True)

# Tidal deformabilities, in units where the total mass is 1
lambda1, lambda2 = symbols('lambda1, lambda2', real=True)

# v = x**(1/2) = Omega_orb**(1/3)
v = symbols('v', real=True)

# <headingcell level=1>

# Derived variables

# <markdowncell>

# Any variable that can be derived from the variables above should be put in this section.
# 
# These variables should probably be left in arbitrary form, unless a particular simplification is desired.  The `BasicSubstitutions` dictionary should map from the general names and their definitions in terms of basic variables.  In numerical codes, their values can be calculated once per time step and then stored, so that the values do not have to be re-calculated every time they appear in an expression.

# <codecell>

BasicSubstitutions = OrderedDict() # For now, just initialize the dictionary

# <markdowncell>

# Of course, some variables will only need to be computed once for the entire system, and stored.  Just so that these quantities are not continually recalculated, we make a list of them for future reference.

# <codecell>

VariableConstants = []

# <markdowncell>

# Various common combinations of the two masses:

# <codecell>

m2 = symbols('m2', real=True)
BasicSubstitutions[m2] = 1-m1
VariableConstants += [m2]

m = symbols('m', real=True);
BasicSubstitutions[m] = m1+m2
VariableConstants += [m]

delta = symbols('delta', real=True);
BasicSubstitutions[delta] = m1-m2
VariableConstants += [delta]

nu, nu__2, nu__3 = symbols('nu, nu__2, nu__3', real=True);
BasicSubstitutions[nu] = m1*m2
BasicSubstitutions[nu__2] = (m1*m2)**2
BasicSubstitutions[nu__3] = (m1*m2)**3
VariableConstants += [nu, nu__2, nu__3]

q = symbols('q', real=True);
BasicSubstitutions[q] = m1/m2
VariableConstants += [q]

# <markdowncell>

# The system's vector basis is given by $(\hat{L}_{\text{N}}, \hat{n}, \hat{\lambda})$.  Here, we define the remaining element, and give the relevant substitutions in terms of Cartesian basis elements.

# <codecell>

Lhat_N, nhat = symbols('Lhat_N, nhat', real=True)

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

chi1, chi2 = symbols('chi1, chi2', real=True)
BasicSubstitutions[chi1] = array([chi1_x, chi1_y, chi1_z])
BasicSubstitutions[chi2] = array([chi2_x, chi2_y, chi2_z])

chi1chi1, chi1chi2, chi2chi2 = symbols('chi1chi1, chi1chi2, chi2chi2', real=True)
BasicSubstitutions[chi1chi1] = dot(BasicSubstitutions[chi1], BasicSubstitutions[chi1])
BasicSubstitutions[chi1chi2] = dot(BasicSubstitutions[chi1], BasicSubstitutions[chi2])
BasicSubstitutions[chi2chi2] = dot(BasicSubstitutions[chi2], BasicSubstitutions[chi2])

chi1_l, chi1_n, chi1_la = symbols('chi1_l, chi1_n, chi1_lambda', real=True)
chi2_l, chi2_n, chi2_la = symbols('chi2_l, chi2_n, chi2_lambda', real=True)
BasicSubstitutions[chi1_l]  = dot(BasicSubstitutions[chi1], BasicSubstitutions[Lhat_N])
BasicSubstitutions[chi1_n]  = dot(BasicSubstitutions[chi1], BasicSubstitutions[nhat])
BasicSubstitutions[chi1_la] = dot(BasicSubstitutions[chi1], BasicSubstitutions[lambdahat])
BasicSubstitutions[chi2_l]  = dot(BasicSubstitutions[chi2], BasicSubstitutions[Lhat_N])
BasicSubstitutions[chi2_n]  = dot(BasicSubstitutions[chi2], BasicSubstitutions[nhat])
BasicSubstitutions[chi2_la] = dot(BasicSubstitutions[chi2], BasicSubstitutions[lambdahat])

sqrt1Mchi1chi1, sqrt1Mchi2chi2 = symbols('sqrt1Mchi1chi1, sqrt1Mchi2chi2', real=True)
BasicSubstitutions[sqrt1Mchi1chi1] = sqrt(1-chi1chi1)
BasicSubstitutions[sqrt1Mchi2chi2] = sqrt(1-chi2chi2)

S, S_l, S_n, S_la = symbols('S, S_l, S_n, S_lambda', real=True)
BasicSubstitutions[S] = chi1*m1**2 + chi2*m2**2
BasicSubstitutions[S_l] = chi1_l*m1**2 + chi2_l*m2**2
BasicSubstitutions[S_n] = chi1_n*m1**2 + chi2_l*m2**2
BasicSubstitutions[S_la] = chi1_la*m1**2 + chi2_la*m2**2

Sigma, Sigma_l, Sigma_n, Sigma_la = symbols('Sigma, Sigma_l, Sigma_n, Sigma_lambda', real=True)
BasicSubstitutions[Sigma] = chi2*m2 - chi1*m1
BasicSubstitutions[Sigma_l] = chi2_l*m2 - chi1_l*m1
BasicSubstitutions[Sigma_n] = chi2_n*m2 - chi1_n*m1
BasicSubstitutions[Sigma_la] = chi2_la*m2 - chi1_la*m1

# <markdowncell>

# Other powers of the angular velocity that find frequent use:

# <codecell>

x, Omega_orb, logv = symbols('x, Omega_orb, logv', real=True)
BasicSubstitutions[x] = v**2
BasicSubstitutions[Omega_orb] = v**3
BasicSubstitutions[logv] = log(v)

