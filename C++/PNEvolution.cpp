#include "PNEvolution.hpp"
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_odeiv2.h>

using Quaternions::Quaternion;
using Quaternions::exp;
using Quaternions::conjugate;
using Quaternions::xHat;
using Quaternions::yHat;
using Quaternions::zHat;
using Quaternions::FrameFromAngularVelocity_Integrand;
using Quaternions::FrameFromAngularVelocity_2D_Integrand;


#include "PNApproximants.ipp"


void PostNewtonian::EvolvePN(const std::string& Approximant,
			     const double v_i, const double m1,
			     const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
			     std::vector<double>& v,
			     std::vector<std::vector<double> >& chi1,
			     std::vector<std::vector<double> >& chi2,
			     std::vector<Quaternions::Quaternion>& R_frame,
			     std::vector<double>& Phi
			     )
{
  
  return;
}

void PostNewtonian::EvolvePN(const std::string& Approximant, const double PNOrder,
			     const double v0, const double v_i,
			     const double m1,
			     const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
			     const Quaternions::Quaternion& R_frame_i,
			     std::vector<double>& v,
			     std::vector<std::vector<double> >& chi1,
			     std::vector<std::vector<double> >& chi2,
			     std::vector<Quaternions::Quaternion>& R_frame,
			     std::vector<double>& Phi
			     )
{
  
  return;
}

std::vector<std::vector<double> > ellHat(const std::vector<Quaternions::Quaternion>& R_frame) {
  std::vector<std::vector<double> > ellhat(R_frame.size());
  for(unsigned int i=0; i<ellhat.size(); ++i) {
    ellhat[i] = (R_frame[i]*zHat*R_frame[i].conjugate()).vec();
  }
  return ellhat;
}
std::vector<std::vector<double> > nHat(const std::vector<Quaternions::Quaternion>& R_frame) {
  std::vector<std::vector<double> > nhat(R_frame.size());
  for(unsigned int i=0; i<nhat.size(); ++i) {
    nhat[i] = (R_frame[i]*xHat*R_frame[i].conjugate()).vec();
  }
  return nhat;
}
std::vector<std::vector<double> > lambdaHat(const std::vector<Quaternions::Quaternion>& R_frame) {
  std::vector<std::vector<double> > lambdahat(R_frame.size());
  for(unsigned int i=0; i<lambdahat.size(); ++i) {
    lambdahat[i] = (R_frame[i]*yHat*R_frame[i].conjugate()).vec();
  }
  return lambdahat;
}
