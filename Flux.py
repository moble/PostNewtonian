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

# Individual flux terms

# <markdowncell>

# In this notebook, every term will be multiplied by the following coefficient.  [Note that fractions need to be entered as, e.g., `frac(32,5)` so that they are not converted to finite-precision decimals.]

# <codecell>

FluxCoefficient = frac(32,5)*nu**2*v**10

# <markdowncell>

# We will define a dictionary of dictionaries to contain all the flux terms.  The top-level keys will be different types of terms: nonspinning, spin-spin, etc.  The dictionaries they contain will give those terms at different PN orders. 

# <codecell>

FluxTerms = {'Nonspinning':{},
             'IncompleteNonspinning':{},
             'SpinSquared':{},
             'SpinOrbit':{},
             'NSTidal':{},
             }

# <markdowncell>

# The nonspinning flux terms are complete to 3.5pN order.  These terms are given by Eq. (231) of [Blanchet (2006)](http://www.livingreviews.org/lrr-2006-4):

# <codecell>

FluxTerms['Nonspinning'][0] = 1
FluxTerms['Nonspinning'][1] = 0
FluxTerms['Nonspinning'][2] = -frac(1247,336) - frac(35,12)*nu
FluxTerms['Nonspinning'][3] = 4*pi
FluxTerms['Nonspinning'][4] = -frac(44711,9072) + frac(9271,504)*nu + frac(65,18)*nu**2
FluxTerms['Nonspinning'][5] = (-frac(8191,672) - frac(583,24)*nu)*pi
FluxTerms['Nonspinning'][6] = \
    frac(6643739519,69854400) + frac(16,3)*pi**2 - EulerGamma*frac(1712,105) - frac(1712,105)*ln(4) \
     - frac(1712,105)*ln(v)+ (-frac(134543,7776) + frac(41,48)*pi**2)*nu - frac(94403,3024)*nu**2 - frac(775,324)*nu**3
FluxTerms['Nonspinning'][7] = (-frac(16285,504) + frac(214745,1728)*nu + frac(193385,3024)*nu**2)*pi

# <markdowncell>

# This term is the 4.0pN nonspinning term in the limit as the mass ratio goes to zero.  This is given by Eq. (43) of [Tagoshi and Sasaki (1994)](http://arxiv.org/abs/gr-qc/9405062v1):

# <codecell>

FluxTerms['IncompleteNonspinning'][8] = \
    -frac(323105549467,3178375200) - frac(1369,126)*pi**2 + frac(232597,4410)*EulerGamma \
    + frac(39931,294)*ln(2) - frac(47385,1568)*ln(3) + frac(232597,4410)*ln(v)

# <markdowncell>

# The following are EMRI terms from Appendix A of [Fujita (2012)](http://arxiv.org/abs/1211.5535v1).  He computed them up to 22pN.  That seems like overkill, so we'll just go up to 6pN.

# <codecell>

FluxTerms['IncompleteNonspinning'][8] = \
    232597*log(v)/4410 - 1369*pi**2/126 - 323105549467/3178375200 - 47385*log(3)/1568 + 232597*EulerGamma/4410 \
    + 39931*log(2)/294
FluxTerms['IncompleteNonspinning'][9] = \
    -6848*pi*log(v)/105 - 13696*pi*log(2)/105 - 6848*EulerGamma*pi/105 + 265978667519*pi/745113600
FluxTerms['IncompleteNonspinning'][10] = \
    916628467*log(v)/7858620 - 2500861660823683/2831932303200 - 424223*pi**2/6804 - 83217611*log(2)/1122660 \
    + 916628467*EulerGamma/7858620 + 47385*log(3)/196
FluxTerms['IncompleteNonspinning'][11] = \
    177293*pi*log(v)/1176 - 142155*pi*log(3)/784 + 8399309750401*pi/101708006400 + 177293*EulerGamma*pi/1176 \
    + 8521283*pi*log(2)/17640
FluxTerms['IncompleteNonspinning'][12] = \
    1465472*log(v)**2/11025 - 246137536815857*log(v)/157329572400 - 27392*pi**2*log(v)/315 \
    + 2930944*EulerGamma*log(v)/11025 + 5861888*log(2)*log(v)/11025 - 271272899815409*log(2)/157329572400 \
    - 54784*pi**2*log(2)/315 - 246137536815857*EulerGamma/157329572400 - 437114506833*log(3)/789268480 - 256*pi**4/45 \
    - 27392*EulerGamma*pi**2/315 - 27392*zeta(3)/105 - 37744140625*log(5)/260941824 + 1465472*EulerGamma**2/11025 \
    + 5861888*EulerGamma*log(2)/11025 + 5861888*log(2)**2/11025 + 2067586193789233570693/602387400044430000 \
    + 3803225263*pi**2/10478160

# <markdowncell>

# The spin-squared terms (by which I mean both spin-spin and spin-orbit squared terms) in the flux are known only at 2pN order (from [Kidder (1995)](http://link.aps.org/doi/10.1103/PhysRevD.52.821) and [Will and Wiseman (1996)](http://link.aps.org/doi/10.1103/PhysRevD.54.4813)).  They are most conveniently given in Eq. (C10) of [Arun et al.](http://arxiv.org/abs/0810.5336v3)  We first need to convert from Arun et al.'s slightly inconvenient spin definitions:

# <codecell>

chis = array([chi1_l+chi2_l,chi1_n+chi2_n,chi1_la+chi2_la])/2
chia = array([chi1_l-chi2_l,chi1_n-chi2_n,chi1_la-chi2_la])/2
chis_l = (chi1_l+chi2_l)/2
chia_l = (chi1_l-chi2_l)/2
SSTerm = (frac(287,96) + nu/24)*(chis_l)**2 - (frac(89,96) + frac(7,24)*nu)*dot(chis,chis) \
    + (frac(287,96) - 12*nu)*(chia_l)**2 + (-(frac(89,96)) + 4*nu)*dot(chia,chia) \
    + frac(287,48)*delta*chis_l*chia_l - frac(89,48)*delta*dot(chis,chia)

# <codecell>

FluxTerms['SpinSquared'][4] = horner(SSTerm.expand().simplify())

# <markdowncell>

# The spin-orbit terms in the flux are now known to 4.0pN.  These terms come from Eq. (4.9) of [Marsat et al. (2013)](http://arxiv.org/abs/1307.6793v1):

# <codecell>

FluxTerms['SpinOrbit'][3] = 4*S_l - frac(5,4)*delta*Sigma_l
FluxTerms['SpinOrbit'][4] = 0
FluxTerms['SpinOrbit'][5] = (-frac(9,2)+frac(272,9)*nu)*S_l + (-frac(13,16)+frac(43,4)*nu)*delta*Sigma_l
FluxTerms['SpinOrbit'][6] = (-16*pi)*S_l + (-frac(31,6)*pi)*delta*Sigma_l
FluxTerms['SpinOrbit'][7] = (frac(476645,6804)+frac(6172,189)*nu-frac(2810,27)*nu**2)*S_l \
                            + (frac(9535,336)+frac(1849,126)*nu-frac(1501,36)*nu**2)*delta*Sigma_l
FluxTerms['SpinOrbit'][8] = (-frac(3485,96)+frac(13879,72)*nu)*pi*S_l \
                            + (-frac(7163,672)+frac(130583,2016)*nu)*pi*delta*Sigma_l

# <markdowncell>

# The tidal-coupling terms come in to the energy at relative 5pN order, and are known partially at 6pN order.
# 
# These terms come from Eq. (3.6) of [Vines et al. (2011)](http://prd.aps.org/abstract/PRD/v83/i8/e084051).  Note their unusual convention for mass ratios, where $\chi_1 = m_1/m$ in their notation; in particular, $\chi$ is not a spin parameter.  Also note that $\hat{\lambda} = \lambda_2 v^{10}/(m_1+m_2)^5$, and we need to add the coupling terms again with $1 \leftrightarrow 2$.  Finally, note the normalization difference, where a different overall factor is used, leading to a sign difference.

# <codecell>

FluxTerms['NSTidal'][10] = (12-18/m2)*lambda2 + (12-18/m1)*lambda1
FluxTerms['NSTidal'][11] = 0
FluxTerms['NSTidal'][12] = (704+1803*m2-4501*m2**2+2170*m2**3)*lambda2/(28*m2) \
                           + (704+1803*m1-4501*m1**2+2170*m1**3)*lambda1/(28*m1)

# Note that the above terms should be divided by (m1+m2)**5, and m1 and m2 should be divided by (m1+m2),
# except that here we use units with m1+m2=1

# <codecell>

#print FluxTerms

# <headingcell level=1>

# Collected flux terms

# <codecell>

def FluxSum(SpinTerms=True, IncompleteNonspinningTerms=False, NSTidalTerms=False) :
    """
    Return an expression for the GW flux with the given options.
    
    """
    F = sum([FluxTerms['Nonspinning'][i]*v**i for i in sorted(FluxTerms['Nonspinning'])])
    if (SpinTerms) :
        F += sum([FluxTerms['SpinOrbit'][i]*v**i for i in sorted(FluxTerms['SpinOrbit'])])
        F += sum([FluxTerms['SpinSquared'][i]*v**i for i in sorted(FluxTerms['SpinSquared'])])
    if (IncompleteNonspinningTerms) :
        F += sum([FluxTerms['IncompleteNonspinning'][i]*v**i for i in sorted(FluxTerms['IncompleteNonspinning'])])
    if (NSTidalTerms) :
        F += sum([FluxTerms['NSTidal'][i]*v**i for i in sorted(FluxTerms['NSTidal'])])
    return FluxCoefficient*F

# <codecell>

#Flux = FluxSum(IncompleteNonspinningTerms=True).subs(log(v), logv).subs(Pow(nu,3), nu__3).subs(Pow(nu,2), nu__2)
#print CCodeOutput(['Flux'])

