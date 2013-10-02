# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# Always run this first
# NOTE: Do not define new variables in this notebook or in Variables.py; 
#       define them in Variables.ipynb.  Variables.py will be automatically
#       recreated when Variables.ipynb is saved.

from __future__ import division # This needs to be here, even though it's in Variables.py
execfile('ExecNotebook.ipy')
execnotebook('Variables.ipynb')
execnotebook('CodeOutput.ipynb')

# <markdowncell>

# This collection of absorption terms is possibly incomplete; only [Alvi (2001)](http://link.aps.org/doi/10.1103/PhysRevD.64.104020) is used.
# 
# See also (the last three referenced by Marsat et al.)
# 
#   - [Tagoshi, Sasaki (1994)](http://arxiv.org/abs/gr-qc/9405062)
#   - E. Poisson and M. Sasaki, Phys. Rev. D 51, 5753 (1995)
#   - H. Tagoshi, S. Mano, and E. Takasugi, Prog. Theor. Phys. 98, 829 (1997), arXiv:gr-qc/9711072
#   - K. Chatziioannou, E. Poisson, and N. Yunes, (2012), arXiv:1211.1686

# <codecell>

AbsorptionCoefficient = frac(32,5)*nu**2*v**10

# <codecell>

AlviTerm =  (v**5/4)*(m1**3)*(1+3*chi1chi1)*(-chi1_l + 2*(1+sqrt(1-chi1chi1))*m1*v**3) \
           +(v**5/4)*(m2**3)*(1+3*chi2chi2)*(-chi2_l + 2*(1+sqrt(1-chi2chi2))*m2*v**3)

# <codecell>

# Absorption = AbsorptionCoefficient*AlviTerm

