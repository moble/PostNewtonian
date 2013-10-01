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

# Individual energy terms

# <markdowncell>

# In this notebook, every term will be multiplied by the following coefficient.  [Note that fractions need to be entered as, e.g., `frac(32,5)` so that they are not converted to finite-precision decimals.]

# <codecell>

EnergyCoefficient = -(nu*v**2)/2

# <markdowncell>

# We will define a dictionary of dictionaries to contain all the energy terms.  The top-level keys will be different types of terms: nonspinning, spin-spin, etc.  The dictionaries they contain will give those terms at different PN orders. 

# <codecell>

EnergyTerms = {'Nonspinning':{},
               'ExtraEMRI':{},
               'SpinSpin':{},
               'SpinOrbit':{},
               'NSTidal':{},
               }

# <markdowncell>

# These terms come from Eq. (194) of [Blanchet (2006)](http://www.livingreviews.org/lrr-2006-4):

# <codecell>

EnergyTerms['Nonspinning'][0] = 1
EnergyTerms['Nonspinning'][1] = 0
EnergyTerms['Nonspinning'][2] = -frac(3,4) - frac(1,12)*nu
EnergyTerms['Nonspinning'][3] = 0
EnergyTerms['Nonspinning'][4] = -frac(27,8) + frac(19,8)*nu - frac(1,24)*nu**2
EnergyTerms['Nonspinning'][5] = 0
EnergyTerms['Nonspinning'][6] = -frac(675,64) + (frac(34445,576) - frac(205,96)*pi**2)*nu - frac(155,96)*nu**2 \
                                - frac(35,5184)*nu**3
EnergyTerms['Nonspinning'][7] = 0

# <markdowncell>

# The EMRI terms...

# <codecell>

EnergyTerms['ExtraEMRI'][8] = 

# <markdowncell>

# The spin-spin terms in the energy are known to...

# <codecell>

EnergyTerms['SpinSpin'][0] = 
EnergyTerms['SpinSpin'][1] = 
EnergyTerms['SpinSpin'][2] = 
EnergyTerms['SpinSpin'][3] = 
EnergyTerms['SpinSpin'][4] = 
EnergyTerms['SpinSpin'][5] = 
EnergyTerms['SpinSpin'][6] = 
EnergyTerms['SpinSpin'][7] = 
EnergyTerms['SpinSpin'][8] = 

# <markdowncell>

# The spin-orbit terms in the energy are now complete to 4.0pN (the last term is zero).  These terms come from Eq. (4.6) of [Bohé et al. (2012)](http://arxiv.org/abs/1212.5520v2):

# <codecell>

EnergyTerms['SpinOrbit'][0] = 0
EnergyTerms['SpinOrbit'][1] = 0
EnergyTerms['SpinOrbit'][2] = 0
EnergyTerms['SpinOrbit'][3] = frac(14,3)*S_l + 2*delta*Sigma_l
EnergyTerms['SpinOrbit'][4] = 0
EnergyTerms['SpinOrbit'][5] = (11-61*nu/9)*S_l + (3-10*nu/3)*delta*Sigma_l
EnergyTerms['SpinOrbit'][6] = 0
EnergyTerms['SpinOrbit'][7] = (frac(135,4)-frac(367,4)*nu+frac(29,12)*nu**2)*S_l \
                              + (frac(27,4)-39*nu+frac(5,4)*nu**2)*delta*Sigma_l
EnergyTerms['SpinOrbit'][8] = 0

# <markdowncell>

# The tidal-coupling terms come in to the energy at relative 5pN order, and are known to 6pN order.
# 
# These terms come from Eq. (2.11) of [Vines et al. (2011)](http://prd.aps.org/abstract/PRD/v83/i8/e084051).  Note their unusual convention for mass ratios, where $\chi_1 = m_1/m$ in their notation; in particular, $\chi$ is not a spin parameter.  Also note that $\hat{\lambda} = \lambda_2 v^{10}/(m_1+m_2)^5$, and we need to add the coupling terms again with $1 \leftrightarrow 2$.  Finally, note the normalization difference, where the overall factor is different by $-2$.

# <codecell>

EnergyTerms['NSTidal'][0] = 0
EnergyTerms['NSTidal'][1] = 0
EnergyTerms['NSTidal'][2] = 0
EnergyTerms['NSTidal'][3] = 0
EnergyTerms['NSTidal'][4] = 0
EnergyTerms['NSTidal'][5] = 0
EnergyTerms['NSTidal'][6] = 0
EnergyTerms['NSTidal'][7] = 0
EnergyTerms['NSTidal'][8] = 0
EnergyTerms['NSTidal'][9] = 0
EnergyTerms['NSTidal'][10] = -9*(m1/m2)*lambda2 - 9*(m2/m1)*lambda1
EnergyTerms['NSTidal'][11] = 0
EnergyTerms['NSTidal'][12] = -frac(11,2)*(m1/m2)*(3+2*m2+3*m2**2)*lambda2 - frac(11,2)*(m2/m1)*(3+2*m1+3*m1**2)*lambda1

# Note that the above terms should be divided by (m1+m2)**5, except that here we use units with m1+m2=1

