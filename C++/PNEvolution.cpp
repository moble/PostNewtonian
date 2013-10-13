#include "PNEvolution.hpp"
#include <cmath>

using Quaternions::Quaternion;
using Quaternions::exp;
using Quaternions::conjugate;

void AngularVelocityIntegrand(const double r0, const double r1, const double r2,
                              const double omega0, const double omega1, const double omega2,
                              double& rdot0, double& rdot1, double& rdot2)
{
  const double rmag = sqrt(r0*r0+r1*r1+r2*r2);

  // If rmag is basically zero, return an answer assuming exactly zero
  if(rmag<Quaternion_Epsilon) {
    rdot0 = omega0/2.0;
    rdot1 = omega1/2.0;
    rdot2 = omega2/2.0;
    return;
  }

  // Otherwise, actually do the calculation
  const double dotTerm = (r0*omega0+r1*omega1+r2*omega2)/(rmag*rmag);
  const double cotTerm = rmag/(2*tan(rmag));
  rdot0 = (omega0 - r0*dotTerm)*cotTerm + r0*dotTerm + 0.5*(omega1*r2 - omega2*r1);
  rdot1 = (omega1 - r1*dotTerm)*cotTerm + r1*dotTerm + 0.5*(omega2*r0 - omega0*r2);
  rdot2 = (omega2 - r2*dotTerm)*cotTerm + r2*dotTerm + 0.5*(omega0*r1 - omega1*r0);
  return;
}

void AngularVelocityIntegrand(const double r0, const double r1,
                              const double omega0, const double omega1, double omega2,
                              double& rdot0, double& rdot1)
{
  const double rmag = sqrt(r0*r0+r1*r1);

  // If rmag is basically zero, return an answer assuming exactly zero
  if(rmag<Quaternion_Epsilon) {
    rdot0 = omega0/2.0;
    rdot1 = omega1/2.0;
    return;
  }

  // Otherwise, actually do the calculation
  const double dotTerm = (r0*omega0+r1*omega1)/(rmag*rmag);
  const double cotTerm = rmag/(2*tan(rmag));
  rdot0 = (omega0 - r0*dotTerm)*cotTerm + r0*dotTerm - 0.5*omega2*r1;
  rdot1 = (omega1 - r1*dotTerm)*cotTerm + r1*dotTerm + 0.5*omega2*r0;
  return;
}

const Quaternions::Quaternion xHat(0,1.0,0,0);
const Quaternions::Quaternion yHat(0,0,1.0,0);
const Quaternions::Quaternion zHat(0,0,0,1.0);

#include "PNApproximants.ipp"
