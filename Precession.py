# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# Always run this first
# NOTE: Do not define new variables in this notebook or in Variables.py; 
#       define them in Variables.ipynb.  Variables.py will be automatically
#       recreated when Variables.ipynb is saved.

%run -i Variables.py
from __future__ import division # This needs to be here, even though it's in Variables.py

# <headingcell level=1>

# Precession of orbital angular velocity $\vec{\Omega}_{\vec{\Omega}_{\text{orb}}}$

# <markdowncell>

# Eqs. (4.3) and (4.4) of [Bohé et al. (2013)](http://arxiv.org/abs/1212.5520v2) say that the precession of the orbital angular velocity is along $\hat{n}$, with magnitude (in their notation) $a_{\ell}/r\omega = \gamma\, a_{\ell} / v^3$.

# <codecell>

gamma_Omega_orbFactor = v**2
gamma_Omega_orbTerms = {}
gamma_Omega_orbTerms[0] = 1
gamma_Omega_orbTerms[1] = 0
gamma_Omega_orbTerms[2] = -nu/3 + 1
gamma_Omega_orbTerms[3] = 5*S_l/3 + Sigma_l*delta
gamma_Omega_orbTerms[4] = -65*nu/12 + 1
gamma_Omega_orbTerms[5] = 8*S_l*nu/9 + 10*S_l/3 + 2*Sigma_l*delta
gamma_Omega_orbTerms[6] = nu**3/81 + 229*nu**2/36 - 41*pi**2*nu/192 - 2203*nu/2520 + 1
gamma_Omega_orbTerms[7] = -6*S_l*nu**2 - 127*S_l*nu/12 + 5*S_l - 8*Sigma_l*delta*nu**2/3 - 61*Sigma_l*delta*nu/6 + 3*Sigma_l*delta

gamma_Omega_orb = gamma_Omega_orbFactor * sum([gamma_Omega_orbTerms[key]*v**key for key in sorted(gamma_Omega_orbTerms)]).simplify()

# <codecell>

a_ell_Omega_orbFactor = v**7
a_ell_Omega_orbTerms = {}
a_ell_Omega_orbTerms[0] = 7*S_n + 3*Sigma_n*delta
a_ell_Omega_orbTerms[1] = 0
a_ell_Omega_orbTerms[2] = -29*S_n*nu/3 - 10*S_n - 9*Sigma_n*delta*nu/2 - 6*Sigma_n*delta
a_ell_Omega_orbTerms[3] = 0
a_ell_Omega_orbTerms[4] = 52*S_n*nu**2/9 + 59*S_n*nu/4 + 3*S_n/2 + 17*Sigma_n*delta*nu**2/6 + 73*Sigma_n*delta*nu/8 + 3*Sigma_n*delta/2
a_ell_Omega_orbTerms[5] = 0
a_ell_Omega_orbTerms[6] = 0
a_ell_Omega_orbTerms[7] = 0

a_ell_Omega_orb = a_ell_Omega_orbFactor * sum([a_ell_Omega_orbTerms[key]*v**key for key in sorted(a_ell_Omega_orbTerms)])

# <codecell>

Omega_Omegavec_orb = gamma_Omega_orb*a_ell_Omega_orb/v**3

# <codecell>

%run -i /tmp/test.py

# <codecell>

print CCodeOutput([['return ', 'Omega_Omegavec_orb']])

# <codecell>


