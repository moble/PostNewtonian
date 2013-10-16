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

#define QuatLogDiscontinuity 1.4142135623730951
inline Quaternion Unflipped(const Quaternion& R0, const Quaternion& R1) {
  if((R1-R0).abs() > QuatLogDiscontinuity) {
    return -R1;
  } else {
    return R1;
  }
}



#include "PNApproximants.ipp"


void PostNewtonian::EvolvePN(const std::string& Approximant,
			     const double v_i, const double m1,
			     const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
			     std::vector<double>& t, std::vector<double>& v,
			     std::vector<std::vector<double> >& chi1, std::vector<std::vector<double> >& chi2,
			     std::vector<Quaternions::Quaternion>& R_frame,
			     std::vector<double>& Phi
			     )
{
  EvolvePN(Approximant, 3.5,
	   v_i, v_i,
	   m1,
	   chi1_i, chi2_i,
	   Quaternion(1,0,0,0),
	   t, v,
	   chi1, chi2,
	   R_frame,
	   Phi
	   );
  return;
}

// This will be the right-hand side for the ODE integration; params
// will point to a TaylorTn object.
int funcT1 (double t, const double y[], double dydt[], void* params) {
  TaylorTn* Tn = (TaylorTn*) params;
  return Tn->TaylorT1_3p5PN(t, y, dydt);
}
int funcT4 (double t, const double y[], double dydt[], void* params) {
  TaylorTn* Tn = (TaylorTn*) params;
  return Tn->TaylorT4_3p5PN(t, y, dydt);
}
int funcT5 (double t, const double y[], double dydt[], void* params) {
  TaylorTn* Tn = (TaylorTn*) params;
  return Tn->TaylorT5_3p5PN(t, y, dydt);
}

void PostNewtonian::EvolvePN(const std::string& Approximant, const double PNOrder,
			     const double v0, const double v_i,
			     const double m1,
			     const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
			     const Quaternions::Quaternion& R_frame_i,
			     std::vector<double>& t, std::vector<double>& v,
			     std::vector<std::vector<double> >& chi1, std::vector<std::vector<double> >& chi2,
			     std::vector<Quaternions::Quaternion>& R_frame,
			     std::vector<double>& Phi
			     )
{
  std::cerr << __FILE__ << ":" << __LINE__ << ": Add option to evolve in reverse, and option to run both ways." << std::endl;

  // Transform the input into the forms we will actually use
  const double chi1Mag_i = std::sqrt(chi1_i[0]*chi1_i[0] + chi1_i[1]*chi1_i[1] + chi1_i[2]*chi1_i[2]);
  const double chi2Mag_i = std::sqrt(chi2_i[0]*chi2_i[0] + chi2_i[1]*chi2_i[1] + chi2_i[2]*chi2_i[2]);
  const std::vector<double> rfrak_chi1_i = sqrtOfRotor(-Quaternion(chi1_i).normalized()*zHat).log().vec();
  const std::vector<double> rfrak_chi2_i = sqrtOfRotor(-Quaternion(chi2_i).normalized()*zHat).log().vec();
  const std::vector<double> rfrak_ell_i = R_frame_i.log().vec();

  // These are the basic variables to be evolved
  std::vector<double> y(9);
  y[0] = v_i;
  y[1] = rfrak_chi1_i[0];
  y[2] = rfrak_chi1_i[1];
  y[3] = rfrak_chi2_i[0];
  y[4] = rfrak_chi2_i[1];
  y[5] = rfrak_ell_i[0];
  y[6] = rfrak_ell_i[1];
  y[7] = rfrak_ell_i[2];
  y[8] = 0.0;

  // Tn encapsulates all the actual PN calculations -- especially the
  // right-hand sides of the evolution system
  TaylorTn Tn(xHat, yHat, zHat, m1, v_i,
	      chi1Mag_i, chi2Mag_i,
	      rfrak_chi1_i[0], rfrak_chi1_i[1], rfrak_chi2_i[0], rfrak_chi2_i[1],
	      rfrak_ell_i[0], rfrak_ell_i[1], rfrak_ell_i[2]);

  // Here are the parameters for the evolution
  const double nu = m1*(1-m1);
  double time = -5.0/(256.0*nu*std::pow(y[0],8)); // This is the lowest-order pN time-to-merger
  double endtime = -3*time; // Give ourselves a large margin of error in case inspiral runs longer than expected
  const unsigned int MinSteps = 100000; // This is only a very rough lower limit
  const unsigned int MaxSteps = 10000000; // This is a hard upper limit
  double h = 1.0;
  const double eps_abs = 1.e-13;
  const double eps_rel = 1.e-13;
  const double hmin = 1.0e-7;
  const double hmax = (endtime-time) / (2.0*MinSteps);

  // We will be using `push_back`, so we first reserve the rough lower
  // limit we will need (after clearing out any content the input
  // vectors had)
  t.clear(); t.reserve(MinSteps);
  v.clear(); v.reserve(MinSteps);
  chi1.clear(); chi1.reserve(MinSteps);
  chi2.clear(); chi2.reserve(MinSteps);
  R_frame.clear(); R_frame.reserve(MinSteps);
  Phi.clear(); Phi.reserve(MinSteps);

  // Declare and initialize the GSL ODE integrator
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_step* s = gsl_odeiv2_step_alloc(T, y.size());
  gsl_odeiv2_control* c = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
  gsl_odeiv2_evolve* e = gsl_odeiv2_evolve_alloc(y.size());
  gsl_odeiv2_system sysT1 = {funcT1, NULL, y.size(), (void *) &Tn};
  gsl_odeiv2_system sysT4 = {funcT4, NULL, y.size(), (void *) &Tn};
  gsl_odeiv2_system sysT5 = {funcT5, NULL, y.size(), (void *) &Tn};
  gsl_odeiv2_system* sys;
  std::cerr << __FILE__ << ":" << __LINE__ << ": Add more options for PN orders here; possibly use static member function." << std::endl;
  if(Approximant.compare("TaylorT1")==0) {
    sys = &sysT1;
  } else if(Approximant.compare("TaylorT4")==0) {
    sys = &sysT4;
  } else if(Approximant.compare("TaylorT5")==0) {
    sys = &sysT5;
  } else {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Approximant '" << Approximant << "' is not yet implemented." << std::endl;
    throw(-1);
  }

  // Store the data at the first step
  {
    Tn.Recalculate(time, &y[0]);
    t.push_back(time);
    v.push_back(y[0]);
    const Quaternion R_chi1_i = exp(Quaternion(0.0, y[1], y[2], 0.0));
    chi1.push_back((chi1Mag_i*R_chi1_i*zHat*R_chi1_i.conjugate()).vec());
    const Quaternion R_chi2_i = exp(Quaternion(0.0, y[3], y[4], 0.0));
    chi2.push_back((chi2Mag_i*R_chi2_i*zHat*R_chi2_i.conjugate()).vec());
    const Quaternion R_frame_i = exp(Quaternion(0.0, y[5], y[6], y[7]));
    // const Quaternion R_frame_i = Unflipped(R_frame.back(), exp(Quaternion(0.0, y[5], y[6], y[7])));
    R_frame.push_back(R_frame_i);
    Phi.push_back(y[8]);
  }

  // Run the integration
  unsigned int NSteps = 0;
  unsigned int nSteps = 0;
  while (time < endtime) {
    // Take a step
    int status = gsl_odeiv2_evolve_apply(e, c, s, sys, &time, time+hmax, &h, &y[0]);
    ++NSteps;
    ++nSteps;

    // Check if it worked and the system is still reasonable
    if(status == GSL_EDOM) {
      // std::cout << "Velocity v has become greater than 1.0" << std::endl;
      break;
    } else if(status == GSL_EDIVERGE) {
      // std::cout << "Velocity is no longer increasing" << std::endl;
      break;
    } else if(status != GSL_SUCCESS) {
      std::cerr << "GSL odeiv2 error.  Return value=" << status << "\n" << std::endl;
      break;
    }

    // If it worked, store the data
    {
      Tn.Recalculate(time, &y[0]);
      t.push_back(time);
      v.push_back(y[0]);
      const Quaternion R_chi1_i = exp(Quaternion(0.0, y[1], y[2], 0.0));
      chi1.push_back((chi1Mag_i*R_chi1_i*zHat*R_chi1_i.conjugate()).vec());
      const Quaternion R_chi2_i = exp(Quaternion(0.0, y[3], y[4], 0.0));
      chi2.push_back((chi2Mag_i*R_chi2_i*zHat*R_chi2_i.conjugate()).vec());
      const Quaternion R_frame_i = Unflipped(R_frame.back(), exp(Quaternion(0.0, y[5], y[6], y[7])));
      R_frame.push_back(R_frame_i);
      Phi.push_back(y[8]);
    }

    // Check if we should stop because this has gone on suspiciously long
    if(time>=endtime) {
      std::cerr << "Time has gone on four times as long as expected.  This seems strange, so we'll stop."
		<< "\nNote that this is unusual.  You may have a short waveform that stops before merger." << std::endl;
      break;
    }

    // Check if we should stop because there have been too many steps
    if(NSteps>MaxSteps) {
      std::cerr << "\n\nThe integration has taken " << NSteps << ".  This seems excessive, so we'll stop."
		<< "\nNote that this is unusual.  You may have a short waveform that stops before merger." << std::endl;
      break;
    }

    // Check if we should stop because the step has gotten too small,
    // but make sure we at least take 500 steps since the last
    // disruption.  [This is the condition that we expect to stop us
    // near merger.]
    if(nSteps>500 && h<hmin) { break; }

    // Reset values of quaternion logarithms to smaller sizes, if
    // necessary.  If this resets, we reset nSteps to zero, because
    // this may make the time stepper take smaller steps
    const double rfrakMag_chi1 = std::sqrt(y[1]*y[1]+y[2]*y[2]);
    if(rfrakMag_chi1>M_PI/2.) {
      y[1] = (rfrakMag_chi1-M_PI)*y[1]/rfrakMag_chi1;
      y[2] = (rfrakMag_chi1-M_PI)*y[2]/rfrakMag_chi1;
      nSteps=0; // This may make the integrator take a few small steps at first
    }
    const double rfrakMag_chi2 = std::sqrt(y[3]*y[3]+y[4]*y[4]);
    if(rfrakMag_chi2>M_PI/2.) {
      y[3] = (rfrakMag_chi2-M_PI)*y[3]/rfrakMag_chi2;
      y[4] = (rfrakMag_chi2-M_PI)*y[4]/rfrakMag_chi2;
      nSteps=0; // This may make the integrator take a few small steps at first
    }
    const double rfrakMag_ellHat = std::sqrt(y[5]*y[5]+y[6]*y[6]+y[7]*y[7]);
    if(rfrakMag_ellHat>M_PI/2.) {
      y[5] = (rfrakMag_ellHat-M_PI)*y[5]/rfrakMag_ellHat;
      y[6] = (rfrakMag_ellHat-M_PI)*y[6]/rfrakMag_ellHat;
      y[7] = (rfrakMag_ellHat-M_PI)*y[7]/rfrakMag_ellHat;
      nSteps=0; // This may make the integrator take a few small steps at first
    }
  }

  // Free the gsl storage
  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);

  // Make the last time = 0.0
  const double tback = t.back();
  for(unsigned int i=0; i<t.size(); ++i) {
    t[i] -= tback;
  }

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
