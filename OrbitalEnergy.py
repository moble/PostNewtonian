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

# Individual energy terms

# <markdowncell>

# In this notebook, every term will be multiplied by the following coefficient.  [Note that fractions need to be entered as, e.g., `frac(32,5)` so that they are not converted to finite-precision decimals.]

# <codecell>

EnergyCoefficient = -(nu*v**2)/2

# <markdowncell>

# We will define a dictionary of dictionaries to contain all the energy terms.  The top-level keys will be different types of terms: nonspinning, spin-spin, etc.  The dictionaries they contain will give those terms at different PN orders. 

# <codecell>

EnergyTerms = {'Nonspinning':{},
               'IncompleteNonspinning':{},
               'SpinOrbit':{},
               'SpinSquared':{},
               'NSTidal':{},
               }

# <markdowncell>

# The nonspinning orbital binding energy is known through 3.5pN come from Eq. (194) of [Blanchet (2006)](http://www.livingreviews.org/lrr-2006-4).

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

# The 4pN term from Eq. (5.2d) of [Jaranowski and Schäfer](http://arxiv.org/abs/1303.3225v1) is almost known exactly, except that the $\nu$-linear piece of the 4pN term is known only numerically, to 7 decimal places.  The remaining terms are best summarized as Eq. (3.1) of [Barausse et al.](http://arxiv.org/abs/1111.5610v2)  They are only known exactly at lowest order in $\nu$, and numerically to first order in $\nu$, leaving remaining unkown terms at higher orders in $\nu$.

# <codecell>

EnergyTerms['IncompleteNonspinning'][8] = \
    -frac(3969,128) + (153.8803)*nu + (-frac(498449,3456) + frac(3157,576)*pi**2)*nu**2 \
    + frac(301,1728)*nu**3 + frac(77,31104)*nu**4 + frac(896,15)*nu*ln(v)
EnergyTerms['IncompleteNonspinning'][9] = 0
EnergyTerms['IncompleteNonspinning'][10] = -frac(45927,512) + (-55.13)*nu + (-frac(9976,35)-frac(3808,15)*nu)*nu*ln(v)
EnergyTerms['IncompleteNonspinning'][11] = 0
EnergyTerms['IncompleteNonspinning'][12] = -frac(264627,1024) + (588.)*nu - (2288.)*nu*ln(v)

# <markdowncell>

# ***(Is the following true?  What about [this paper](http://arxiv.org/abs/1302.6723v2)?)***
# The spin-squared terms (by which I mean both spin-spin and spin-orbit squared terms) in the energy are known only at 2pN order (from [Kidder (1995)](http://link.aps.org/doi/10.1103/PhysRevD.52.821) and [Will and Wiseman (1996)](http://link.aps.org/doi/10.1103/PhysRevD.54.4813)).  They are most conveniently given in Eq. (C4) of [Arun et al.](http://arxiv.org/abs/0810.5336v3)  We first need to convert from Arun et al.'s slightly inconvenient spin definitions:

# <codecell>

chis = array([chi1_l+chi2_l,chi1_n+chi2_n,chi1_la+chi2_la])/2
chia = array([chi1_l-chi2_l,chi1_n-chi2_n,chi1_la-chi2_la])/2
chis_l = (chi1_l+chi2_l)/2
chia_l = (chi1_l-chi2_l)/2
SSTerm = nu*((dot(chis,chis) - dot(chia,chia)) - 3*((chis_l)**2 - (chia_l)**2)) \
    + (frac(1,2) - nu)*(dot(chis,chis) + dot(chia,chia) - 3*((chis_l)**2 + (chia_l)**2)) \
    + delta*(dot(chis,chia) - 3*(chis_l)*(chia_l))

# <codecell>

EnergyTerms['SpinSquared'][4] = SSTerm.expand().simplify()

# <markdowncell>

# The spin-orbit terms in the energy are now complete to 4.0pN (the last term is zero).  These terms come from Eq. (4.6) of [Bohé et al. (2012)](http://arxiv.org/abs/1212.5520v2):

# <codecell>

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

EnergyTerms['NSTidal'][10] = -9*(m1/m2)*lambda2 - 9*(m2/m1)*lambda1
EnergyTerms['NSTidal'][11] = 0
EnergyTerms['NSTidal'][12] = -frac(11,2)*(m1/m2)*(3+2*m2+3*m2**2)*lambda2 - frac(11,2)*(m2/m1)*(3+2*m1+3*m1**2)*lambda1

# Note that the above terms should be divided by (m1+m2)**5, and each occurence of m1 or m2 should be divided by (m1+m2),
# except that here we use units with m1+m2=1

# <headingcell level=1>

# Collected terms

# <codecell>

def EnergySum(SpinTerms=True, IncompleteNonspinningTerms=False, NSTidalTerms=False) :
    """
    Return an expression for the orbital binding energy with the given options.
    
    """
    E = sum([EnergyTerms['Nonspinning'][i]*v**i for i in sorted(EnergyTerms['Nonspinning'])])
    if (SpinTerms) :
        E += sum([EnergyTerms['SpinOrbit'][i]*v**i for i in sorted(EnergyTerms['SpinOrbit'])])
        E += sum([EnergyTerms['SpinSquared'][i]*v**i for i in sorted(EnergyTerms['SpinSquared'])])
    if (IncompleteNonspinningTerms) :
        E += sum([EnergyTerms['IncompleteNonspinning'][i]*v**i for i in sorted(EnergyTerms['IncompleteNonspinning'])])
    if (NSTidalTerms) :
        E += sum([EnergyTerms['NSTidal'][i]*v**i for i in sorted(EnergyTerms['NSTidal'])])
    return EnergyCoefficient*E

