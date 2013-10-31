#include "PNEvolution.hpp"
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_odeiv2.h>

using Quaternions::Quaternion;
using Quaternions::exp;
using Quaternions::log;
using Quaternions::conjugate;
using Quaternions::xHat;
using Quaternions::yHat;
using Quaternions::zHat;
using Quaternions::FrameFromAngularVelocity_Integrand;
using Quaternions::FrameFromAngularVelocity_2D_Integrand;
using std::vector;

/////////////////////////////////
//// Local utility functions ////
/////////////////////////////////

#define QuatLogDiscontinuity 1.4142135623730951
inline Quaternion Unflipped(const Quaternion& R0, const Quaternion& R1) {
  if((R1-R0).abs() > QuatLogDiscontinuity) {
    return -R1;
  } else {
    return R1;
  }
}

std::vector<double> operator-(const std::vector<double>& b) {
  std::vector<double> c(b.size());
  for(unsigned int i=0; i<b.size(); ++i) {
    c[i] = -b[i];
  }
  return c;
}

std::vector<double> operator*(const double a, const std::vector<double>& b) {
  std::vector<double> c(b.size());
  for(unsigned int i=0; i<b.size(); ++i) {
    c[i] = a*b[i];
  }
  return c;
}

std::vector<double> operator*(const std::vector<double>& b, const double a) {
  std::vector<double> c(b.size());
  for(unsigned int i=0; i<b.size(); ++i) {
    c[i] = a*b[i];
  }
  return c;
}

std::vector<double> operator/(const std::vector<double>& b, const double a) {
  std::vector<double> c(b.size());
  for(unsigned int i=0; i<b.size(); ++i) {
    c[i] = b[i]/a;
  }
  return c;
}

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
  std::vector<double> c(b.size());
  for(unsigned int i=0; i<b.size(); ++i) {
    c[i] = a[i]+b[i];
  }
  return c;
}

std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b) {
  std::vector<double> c(b.size());
  for(unsigned int i=0; i<b.size(); ++i) {
    c[i] = a[i]-b[i];
  }
  return c;
}

void cross(const double a[3], const double b[3], double c[3]) {
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = -a[0]*b[2]+a[2]*b[0];
  c[2] = a[0]*b[1]-a[1]*b[0];
  return;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}


///////////////////////////////////
//// Defining the approximants ////
///////////////////////////////////
class TaylorTn {
public:
  virtual int TaylorT1(double t, const double* y, double* dydt) { return GSL_FAILURE; }
  virtual int TaylorT4(double t, const double* y, double* dydt) { return GSL_FAILURE; }
  virtual int TaylorT5(double t, const double* y, double* dydt) { return GSL_FAILURE; }
  virtual void Recalculate(double t, const double* y) { }
};
#include "PNApproximants.ipp"

/////////////////////////////
//// Evolution functions ////
/////////////////////////////
void PostNewtonian::EvolvePN(const std::string& Approximant,
			     const double v_i, const double m1, const double m2,
			     const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
			     std::vector<double>& t, std::vector<double>& v,
			     std::vector<std::vector<double> >& chi1, std::vector<std::vector<double> >& chi2,
			     std::vector<Quaternions::Quaternion>& R_frame,
			     std::vector<double>& Phi
			     )
{
  EvolvePN(Approximant, 3.5,
	   v_i, v_i,
	   m1, m2,
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
  return Tn->TaylorT1(t, y, dydt);
}
int funcT4 (double t, const double y[], double dydt[], void* params) {
  TaylorTn* Tn = (TaylorTn*) params;
  return Tn->TaylorT4(t, y, dydt);
}
int funcT5 (double t, const double y[], double dydt[], void* params) {
  TaylorTn* Tn = (TaylorTn*) params;
  return Tn->TaylorT5(t, y, dydt);
}

void PostNewtonian::EvolvePN(const std::string& Approximant, const double PNOrder,
			     const double v0, const double v_i,
			     const double m1, const double m2,
			     const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
			     const Quaternions::Quaternion& R_frame_i,
			     std::vector<double>& t, std::vector<double>& v,
			     std::vector<std::vector<double> >& chi1, std::vector<std::vector<double> >& chi2,
			     std::vector<Quaternions::Quaternion>& R_frame,
			     std::vector<double>& Phi
			     )
{
  std::cerr << __FILE__ << ":" << __LINE__ << ": Add easier interface for running both ways." << std::endl;

  // Transform the input into the forms we will actually use
  const std::vector<double> ellHat_i = (R_frame_i*zHat*R_frame_i.conjugate()).vec();
  const std::vector<double> nHat_i = (R_frame_i*xHat*R_frame_i.conjugate()).vec();
  const Quaternion Rax =
    sqrtOfRotor(-normalized(Quaternion(0., ellHat_i[0], ellHat_i[1], ellHat_i[2]))*zHat);
  const double gamma_i = (log(Rax.conjugate() * R_frame_i)).dot(zHat) * 2.0;

  // These are the basic variables to be evolved
  std::vector<double> y(12);
  y[0] = v_i;
  y[1] = chi1_i[0];
  y[2] = chi1_i[1];
  y[3] = chi1_i[2];
  y[4] = chi2_i[0];
  y[5] = chi2_i[1];
  y[6] = chi2_i[2];
  y[7] = ellHat_i[0];
  y[8] = ellHat_i[1];
  y[9] = ellHat_i[2];
  y[10] = 0.0; // Phi
  y[11] = gamma_i;

  {
    const Quaternion R = Rax * exp(((gamma_i)/2.)*zHat);
    std::cout << "d nHat_i: \t" << R_frame_i*xHat*R_frame_i.conjugate() - R*xHat*R.conjugate() << std::endl;
    std::cout << "d lambdaHat_i: \t" << R_frame_i*yHat*R_frame_i.conjugate() - R*yHat*R.conjugate() << std::endl;
    std::cout << "d ellHat_i: \t" << R_frame_i*zHat*R_frame_i.conjugate() - R*zHat*R.conjugate() << std::endl << std::endl;
  }

  // Tn encapsulates all the actual PN calculations -- especially the
  // right-hand sides of the evolution system
  TaylorTn* Tn;
  switch(int(2*PNOrder)) {
  case 0:
    Tn = new TaylorTn_0PN(m1, m2, v_i,
			  chi1_i[0], chi1_i[1], chi1_i[2],
			  chi2_i[0], chi2_i[1], chi2_i[2],
			  ellHat_i[0], ellHat_i[1], ellHat_i[2],
			  nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  case 1:
    Tn = new TaylorTn_0p50PN(m1, m2, v_i,
			     chi1_i[0], chi1_i[1], chi1_i[2],
			     chi2_i[0], chi2_i[1], chi2_i[2],
			     ellHat_i[0], ellHat_i[1], ellHat_i[2],
			     nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  case 2:
    Tn = new TaylorTn_1p0PN(m1, m2, v_i,
			    chi1_i[0], chi1_i[1], chi1_i[2],
			    chi2_i[0], chi2_i[1], chi2_i[2],
			    ellHat_i[0], ellHat_i[1], ellHat_i[2],
			    nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  case 3:
    Tn = new TaylorTn_1p5PN(m1, m2, v_i,
			    chi1_i[0], chi1_i[1], chi1_i[2],
			    chi2_i[0], chi2_i[1], chi2_i[2],
			    ellHat_i[0], ellHat_i[1], ellHat_i[2],
			    nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  case 4:
    Tn = new TaylorTn_2p0PN(m1, m2, v_i,
			    chi1_i[0], chi1_i[1], chi1_i[2],
			    chi2_i[0], chi2_i[1], chi2_i[2],
			    ellHat_i[0], ellHat_i[1], ellHat_i[2],
			    nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  case 5:
    Tn = new TaylorTn_2p5PN(m1, m2, v_i,
			    chi1_i[0], chi1_i[1], chi1_i[2],
			    chi2_i[0], chi2_i[1], chi2_i[2],
			    ellHat_i[0], ellHat_i[1], ellHat_i[2],
			    nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  case 6:
    Tn = new TaylorTn_3p0PN(m1, m2, v_i,
			    chi1_i[0], chi1_i[1], chi1_i[2],
			    chi2_i[0], chi2_i[1], chi2_i[2],
			    ellHat_i[0], ellHat_i[1], ellHat_i[2],
			    nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  case 7:
    Tn = new TaylorTn_3p5PN(m1, m2, v_i,
			    chi1_i[0], chi1_i[1], chi1_i[2],
			    chi2_i[0], chi2_i[1], chi2_i[2],
			    ellHat_i[0], ellHat_i[1], ellHat_i[2],
			    nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  case 8:
    Tn = new TaylorTn_4p0PN(m1, m2, v_i,
			    chi1_i[0], chi1_i[1], chi1_i[2],
			    chi2_i[0], chi2_i[1], chi2_i[2],
			    ellHat_i[0], ellHat_i[1], ellHat_i[2],
			    nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  case 9:
    Tn = new TaylorTn_4p5PN(m1, m2, v_i,
			    chi1_i[0], chi1_i[1], chi1_i[2],
			    chi2_i[0], chi2_i[1], chi2_i[2],
			    ellHat_i[0], ellHat_i[1], ellHat_i[2],
			    nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  case 10:
    Tn = new TaylorTn_5p0PN(m1, m2, v_i,
			    chi1_i[0], chi1_i[1], chi1_i[2],
			    chi2_i[0], chi2_i[1], chi2_i[2],
			    ellHat_i[0], ellHat_i[1], ellHat_i[2],
			    nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  case 11:
    Tn = new TaylorTn_5p5PN(m1, m2, v_i,
			    chi1_i[0], chi1_i[1], chi1_i[2],
			    chi2_i[0], chi2_i[1], chi2_i[2],
			    ellHat_i[0], ellHat_i[1], ellHat_i[2],
			    nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  case 12:
    Tn = new TaylorTn_6p0PN(m1, m2, v_i,
			    chi1_i[0], chi1_i[1], chi1_i[2],
			    chi2_i[0], chi2_i[1], chi2_i[2],
			    ellHat_i[0], ellHat_i[1], ellHat_i[2],
			    nHat_i[0], nHat_i[1], nHat_i[2]);
    break;
  default:
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": PN order " << PNOrder << " is not yet implemented." << std::endl;
    throw(-1);
  }

  // Here are the parameters for the evolution
  const double nu = m1*m2/((m1+m2)*(m1+m2));
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
  gsl_odeiv2_system sysT1 = {funcT1, NULL, y.size(), (void *) Tn};
  gsl_odeiv2_system sysT4 = {funcT4, NULL, y.size(), (void *) Tn};
  gsl_odeiv2_system sysT5 = {funcT5, NULL, y.size(), (void *) Tn};
  gsl_odeiv2_system* sys;
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
    Tn->Recalculate(time, &y[0]);
    vector<double> chi1_i(3), chi2_i(3);
    chi1_i[0] = y[1];
    chi1_i[1] = y[2];
    chi1_i[2] = y[3];
    chi2_i[0] = y[4];
    chi2_i[1] = y[5];
    chi2_i[2] = y[6];
    const Quaternion Rax =
      sqrtOfRotor(-normalized(Quaternion(0., y[7], y[8], y[9]))*zHat);
    const Quaternion R_frame_i = Rax * exp(((y[10]+y[11])/2.)*zHat);
    t.push_back(time);
    v.push_back(y[0]);
    chi1.push_back(chi1_i);
    chi2.push_back(chi2_i);
    R_frame.push_back(R_frame_i);
    Phi.push_back(y[10]);
    std::cout << time << " " << y[11] << std::endl;
  }

  // Run the integration
  unsigned int NSteps = 0;
  while (time < endtime) {
    // Take a step
    int status = gsl_odeiv2_evolve_apply(e, c, s, sys, &time, time+hmax, &h, &y[0]);
    ++NSteps;

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
      Tn->Recalculate(time, &y[0]);
      vector<double> chi1_i(3), chi2_i(3);
      chi1_i[0] = y[1];
      chi1_i[1] = y[2];
      chi1_i[2] = y[3];
      chi2_i[0] = y[4];
      chi2_i[1] = y[5];
      chi2_i[2] = y[6];
      const Quaternion Rax =
	sqrtOfRotor(-normalized(Quaternion(0., y[7], y[8], y[9]))*zHat);
      const Quaternion R_frame_i = Rax * exp(((y[10]+y[11])/2.)*zHat);
      t.push_back(time);
      v.push_back(y[0]);
      chi1.push_back(chi1_i);
      chi2.push_back(chi2_i);
      R_frame.push_back(R_frame_i);
      Phi.push_back(y[10]);
      std::cout << time << " " << y[11] << std::endl;
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
    // but make sure we at least take 100 steps since the beginning.
    // [This is the condition that we expect to stop us near merger.]
    if(NSteps>100 && h<hmin) { break; }
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

std::vector<std::vector<double> > PostNewtonian::ellHat(const std::vector<Quaternions::Quaternion>& R_frame) {
  std::vector<std::vector<double> > ellhat(R_frame.size());
  for(unsigned int i=0; i<ellhat.size(); ++i) {
    ellhat[i] = (R_frame[i]*zHat*R_frame[i].conjugate()).vec();
  }
  return ellhat;
}
std::vector<std::vector<double> > PostNewtonian::nHat(const std::vector<Quaternions::Quaternion>& R_frame) {
  std::vector<std::vector<double> > nhat(R_frame.size());
  for(unsigned int i=0; i<nhat.size(); ++i) {
    nhat[i] = (R_frame[i]*xHat*R_frame[i].conjugate()).vec();
  }
  return nhat;
}
std::vector<std::vector<double> > PostNewtonian::lambdaHat(const std::vector<Quaternions::Quaternion>& R_frame) {
  std::vector<std::vector<double> > lambdahat(R_frame.size());
  for(unsigned int i=0; i<lambdahat.size(); ++i) {
    lambdahat[i] = (R_frame[i]*yHat*R_frame[i].conjugate()).vec();
  }
  return lambdahat;
}
