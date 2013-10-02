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

# <headingcell level=1>

# Precession of orbital angular velocity $\vec{\Omega}_{\hat{L}_{\text{N}}}$

# <markdowncell>

# [Bohé et al. (2013)](http://arxiv.org/abs/1212.5520v2) say that the precession of the orbital angular velocity is along $\hat{n}$, with magnitude (in their notation) $a_{\ell}/r\omega = \gamma\, a_{\ell} / v^3$.
# 
# *NOTE:* There is a 3pN gauge term in $\gamma$ that I have simply dropped here.  It is $\ln(r/r_0')$.
# 
# The following are Eqs. (4.3) and (4.4) of [Bohé et al. (2013)](http://arxiv.org/abs/1212.5520v2).

# <codecell>

gamma_Lhat_NFactor = v**2
gamma_Lhat_NTerms = {}
gamma_Lhat_NTerms[0] = 1
gamma_Lhat_NTerms[1] = 0
gamma_Lhat_NTerms[2] = -nu/3 + 1
gamma_Lhat_NTerms[3] = 5*S_l/3 + Sigma_l*delta
gamma_Lhat_NTerms[4] = -65*nu/12 + 1
gamma_Lhat_NTerms[5] = (frac(8,9)*nu + frac(10,3))*S_l + 2*Sigma_l*delta
gamma_Lhat_NTerms[6] = nu**3/81 + 229*nu**2/36 - 41*pi**2*nu/192 - 2203*nu/2520 + 1
gamma_Lhat_NTerms[7] = (-6*nu**2 - 127*nu/12 + 5)*S_l - 8*Sigma_l*delta*nu**2/3 + (-61*nu/6 + 3)*Sigma_l*delta

gamma_Lhat_N = gamma_Lhat_NFactor * sum([gamma_Lhat_NTerms[key]*v**key
                                         for key in sorted(gamma_Lhat_NTerms)]).simplify()

# <codecell>

a_ell_Lhat_NFactor = v**7
a_ell_Lhat_NTerms = {}
a_ell_Lhat_NTerms[0] = 7*S_n + 3*Sigma_n*delta
a_ell_Lhat_NTerms[1] = 0
a_ell_Lhat_NTerms[2] = (-29*nu/3-10)*S_n + (-9*nu/2-6)*delta*Sigma_n
a_ell_Lhat_NTerms[3] = 0
a_ell_Lhat_NTerms[4] = (frac(52,9)*nu**2 + frac(59,4)*nu + frac(3,2))*S_n \
                          + (frac(17,6)*nu**2 + frac(73,8)*nu + frac(3,2))*delta*Sigma_n
a_ell_Lhat_NTerms[5] = 0
a_ell_Lhat_NTerms[6] = 0
a_ell_Lhat_NTerms[7] = 0

a_ell_Lhat_N = a_ell_Lhat_NFactor * sum([a_ell_Lhat_NTerms[key]*v**key for key in sorted(a_ell_Lhat_NTerms)])
a_ell_overvcubed = a_ell_Lhat_N/v**3

# <codecell>

# Omega_Lhat_N = gamma_Lhat_N*a_ell_Lhat_N/v**3

# <headingcell level=1>

# Precession of spins $\vec{\Omega}_{1,2}$

# <markdowncell>

# In this formulation, the spins precess entirely about the orbital angular velocity $\hat{L}_{\text{N}}$.
# 
# The frequency of that precession is given by Eq. (4.5) of [Bohé et al. (2013)](http://arxiv.org/abs/1212.5520v2):

# <codecell>

Omega1_Factor = v**5
Omega1_Terms = {}
Omega1_Terms[0] = frac(3,4) + frac(1,2)*nu - frac(3,4)*delta
Omega1_Terms[1] = 0
Omega1_Terms[2] = frac(9,16) + frac(5,4)*nu - frac(1,24)*nu**2 + delta*(-frac(9,16) + frac(5,8)*nu)
Omega1_Terms[3] = 0
Omega1_Terms[4] = frac(27,32) + frac(3,16)*nu - frac(105,32)*nu**2 - frac(1,48)*nu**3 \
                  + delta*(-frac(27,32) + frac(39,8)*nu - frac(5,32)*nu**2)

Omega1 = Omega1_Factor * sum([Omega1_Terms[key]*v**key for key in sorted(Omega1_Terms)])
Omega2 = Omega1.subs(delta, -delta)

