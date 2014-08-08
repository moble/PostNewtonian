#include "PNEvolution.hpp"
#include <cmath>
#include <iomanip>

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

// This macro is useful for debugging
#define INFOTOCERR std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << ": "
#define INFOTOCOUT std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << ": "

#define QuatLogDiscontinuity 1.4142135623730951
inline Quaternion Unflipped(const Quaternion& R0, const Quaternion& R1) {
  if((R1-R0).abs() > QuatLogDiscontinuity) {
    return -R1;
  } else {
    return R1;
  }
}


///////////////////////////////////
//// Defining the approximants ////
///////////////////////////////////
class TaylorTn_Q {
public:
  virtual Quaternion OrbitalAngularMomentum() = 0;
  virtual int TaylorT1(double t, const double* y, double* dydt) { return GSL_FAILURE; }
  virtual int TaylorT4(double t, const double* y, double* dydt) { return GSL_FAILURE; }
  virtual int TaylorT5(double t, const double* y, double* dydt) { return GSL_FAILURE; }
  virtual void Recalculate(double t, const double* y) { }
};
#include "PNApproximants_Q.ipp"


void PostNewtonian::EvolvePN_Q(const std::string& Approximant,
                               const double v_i, const double m1, const double m2,
                               const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
                               std::vector<double>& t, std::vector<double>& v,
                               std::vector<std::vector<double> >& chi1, std::vector<std::vector<double> >& chi2,
                               std::vector<Quaternions::Quaternion>& R_frame,
                               std::vector<double>& Phi, std::vector<std::vector<double> >& L
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
           Phi,
           L
           );
  return;
}

// These will be the right-hand sides for the ODE integration; params
// will point to a TaylorTn_Q object.
namespace {
  int funcT1 (double t, const double y[], double dydt[], void* params) {
    TaylorTn_Q* Tn = (TaylorTn_Q*) params;
    return Tn->TaylorT1(t, y, dydt);
  }
  int funcT4 (double t, const double y[], double dydt[], void* params) {
    TaylorTn_Q* Tn = (TaylorTn_Q*) params;
    return Tn->TaylorT4(t, y, dydt);
  }
  int funcT5 (double t, const double y[], double dydt[], void* params) {
    TaylorTn_Q* Tn = (TaylorTn_Q*) params;
    return Tn->TaylorT5(t, y, dydt);
  }
};

// This function combines the forward and backward vectors we need
template <class T>
void CombineForwardAndBackward(std::vector<T>& Forward, const std::vector<T>& Backward) {
  std::vector<T> Combined(Forward.size()+Backward.size());
  typename std::vector<T>::iterator Cit = Combined.begin();
  for(typename std::vector<T>::const_reverse_iterator Bit = Backward.rbegin(); Bit != Backward.rend(); ++Bit,++Cit) {
    *Cit=*Bit;
  }
  for(typename std::vector<T>::iterator Fit = Forward.begin(); Fit != Forward.end(); ++Fit,++Cit) {
    *Cit=*Fit;
  }
  Forward.swap(Combined);
  return;
}

void PostNewtonian::EvolvePN_Q(const std::string& Approximant, const double PNOrbitalEvolutionOrder,
                               const double v_0, const double v_i,
                               const double m1, const double m2,
                               const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
                               const Quaternions::Quaternion& R_frame_i,
                               std::vector<double>& t, std::vector<double>& v,
                               std::vector<std::vector<double> >& chi1, std::vector<std::vector<double> >& chi2,
                               std::vector<Quaternions::Quaternion>& R_frame,
                               std::vector<double>& Phi, std::vector<std::vector<double> >& L,
                               const unsigned int MinStepsPerOrbit, const bool ForwardInTime)
{
  // Transform the input into the forms we will actually use
  const double chi1Mag = Quaternions::Quaternion(chi1_i).abs();
  const double chi2Mag = Quaternions::Quaternion(chi2_i).abs();
  const Quaternion S_chi1_i = ( chi1Mag>1e-12
                                ? std::sqrt(chi1Mag) * Quaternions::sqrtOfRotor(-Quaternions::Quaternion(chi1_i).normalized()*Quaternions::zHat)
                                : Quaternions::Zero);
  const Quaternion S_chi2_i = ( chi2Mag>1e-12
                                ? std::sqrt(chi2Mag) * Quaternions::sqrtOfRotor(-Quaternions::Quaternion(chi2_i).normalized()*Quaternions::zHat)
                                : Quaternions::Zero);
  const std::vector<double> rfrak_frame_i = R_frame_i.log().vec();

  // These are the basic variables to be evolved
  std::vector<double> y(9);
  y[0] = v_i; // PN expansion parameter v [ = sqrt(x) = Omega_orb^{1/3} ]
  y[1] = 0.0; // x component of logarithm of R_chi1 dynamics rotor
  y[2] = 0.0; // y component of logarithm of R_chi1 dynamics rotor
  y[3] = 0.0; // x component of logarithm of R_chi2 dynamics rotor
  y[4] = 0.0; // y component of logarithm of R_chi2 dynamics rotor
  y[5] = rfrak_frame_i[0]; // x component of logarithm of R_frame rotor
  y[6] = rfrak_frame_i[1]; // y component of logarithm of R_frame rotor
  y[7] = rfrak_frame_i[2]; // z component of logarithm of R_frame rotor
  y[8] = 0.0; // Orbital phase Phi [integrated as a convenient diagnostic; not actually needed]

  // Tn encapsulates all the actual PN calculations -- especially the
  // right-hand sides of the evolution system
  TaylorTn_Q* Tn = 0;
  switch(int(2*PNOrbitalEvolutionOrder)) {
  case 0:
    Tn = new TaylorTn_0PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                            S_chi1_i, S_chi2_i,
                            0.0, 0.0, 0.0, 0.0,
                            rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  case 1:
    Tn = new TaylorTn_0p50PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                               S_chi1_i, S_chi2_i,
                               0.0, 0.0, 0.0, 0.0,
                               rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  case 2:
    Tn = new TaylorTn_1p0PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                              S_chi1_i, S_chi2_i,
                              0.0, 0.0, 0.0, 0.0,
                              rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  case 3:
    Tn = new TaylorTn_1p5PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                              S_chi1_i, S_chi2_i,
                              0.0, 0.0, 0.0, 0.0,
                              rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  case 4:
    Tn = new TaylorTn_2p0PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                              S_chi1_i, S_chi2_i,
                              0.0, 0.0, 0.0, 0.0,
                              rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  case 5:
    Tn = new TaylorTn_2p5PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                              S_chi1_i, S_chi2_i,
                              0.0, 0.0, 0.0, 0.0,
                              rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  case 6:
    Tn = new TaylorTn_3p0PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                              S_chi1_i, S_chi2_i,
                              0.0, 0.0, 0.0, 0.0,
                              rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  case 7:
    Tn = new TaylorTn_3p5PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                              S_chi1_i, S_chi2_i,
                              0.0, 0.0, 0.0, 0.0,
                              rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  case 8:
    Tn = new TaylorTn_4p0PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                              S_chi1_i, S_chi2_i,
                              0.0, 0.0, 0.0, 0.0,
                              rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  case 9:
    Tn = new TaylorTn_4p5PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                              S_chi1_i, S_chi2_i,
                              0.0, 0.0, 0.0, 0.0,
                              rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  case 10:
    Tn = new TaylorTn_5p0PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                              S_chi1_i, S_chi2_i,
                              0.0, 0.0, 0.0, 0.0,
                              rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  case 11:
    Tn = new TaylorTn_5p5PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                              S_chi1_i, S_chi2_i,
                              0.0, 0.0, 0.0, 0.0,
                              rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  case 12:
    Tn = new TaylorTn_6p0PN_Q(xHat, yHat, zHat, m1, m2, v_i,
                              S_chi1_i, S_chi2_i,
                              0.0, 0.0, 0.0, 0.0,
                              rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2]);
    break;
  default:
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": PN order " << PNOrbitalEvolutionOrder << " is not yet implemented." << std::endl;
    throw(-1);
  }
  if(Tn==0) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Tn did not get assigned with PNOrbitalEvolutionOrder=" << PNOrbitalEvolutionOrder << "." << std::endl;
    throw(-1);
  }

  // Here are the parameters for the evolution
  const double nu = m1*m2/((m1+m2)*(m1+m2));
  const double time_to_merger = 5.0/(256.0*nu*std::pow(y[0],8)); // This is the lowest-order pN time-to-merger
  double time = 0.0;
  double endtime =
    ForwardInTime
    ? 4*time_to_merger // Happily run far into positive times just for a comfy margin of error
    : 3*(-5.0/(256.0*nu*std::pow(v_0,8))); // This should be (3x) a pretty good estimate, considering that it should be very early...
  const unsigned int MinSteps = 100000; // This is only a very rough lower limit
  const unsigned int MaxSteps = 200000; // This is a hard upper limit; much bigger than this and Mike's laptop has been known to crash
  double h = ForwardInTime ? 1.0 : -1.0;
  const double eps_abs = 1.e-6;
  const double eps_rel = 1.e-6;
  const double hmin = ForwardInTime ? 1.0e-7 : -1.0e-7;
  const double hmin_storage = ForwardInTime ? 1.0e-5 : -1.0e-5;
  const double hmax = (endtime-time) / (2.0*MinSteps); // Time-direction is taken care of
  double hnext = hmax;

  // We will be using `push_back`, so we first reserve the rough lower
  // limit we will need (after clearing out any content the input
  // vectors had)
  t.clear(); t.reserve(MinSteps);
  v.clear(); v.reserve(MinSteps);
  chi1.clear(); chi1.reserve(MinSteps);
  chi2.clear(); chi2.reserve(MinSteps);
  R_frame.clear(); R_frame.reserve(MinSteps);
  Phi.clear(); Phi.reserve(MinSteps);
  L.clear(); L.reserve(MinSteps);

  // Declare and initialize the GSL ODE integrator
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_step* s = gsl_odeiv2_step_alloc(T, y.size());
  gsl_odeiv2_control* c = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
  gsl_odeiv2_evolve* e = gsl_odeiv2_evolve_alloc(y.size());
  gsl_odeiv2_system sysT1 = {::funcT1, NULL, y.size(), (void *) Tn};
  gsl_odeiv2_system sysT4 = {::funcT4, NULL, y.size(), (void *) Tn};
  gsl_odeiv2_system sysT5 = {::funcT5, NULL, y.size(), (void *) Tn};
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
    t.push_back(time);
    v.push_back(y[0]);
    const Quaternion R_chi1_i = exp(Quaternion(0.0, y[1], y[2], 0.0));
    chi1.push_back((S_chi1_i*R_chi1_i*zHat*R_chi1_i.conjugate()*S_chi1_i.conjugate()).vec());
    const Quaternion R_chi2_i = exp(Quaternion(0.0, y[3], y[4], 0.0));
    chi2.push_back((S_chi2_i*R_chi2_i*zHat*R_chi2_i.conjugate()*S_chi2_i.conjugate()).vec());
    const Quaternion R_frame_i = exp(Quaternion(0.0, y[5], y[6], y[7]));
    // const Quaternion R_frame_i = Unflipped(R_frame.back(), exp(Quaternion(0.0, y[5], y[6], y[7])));
    R_frame.push_back(R_frame_i);
    Phi.push_back(y[8]);
    L.push_back((Tn->OrbitalAngularMomentum()).vec());
  }

  // Run the integration
  unsigned int NSteps = 0;
  unsigned int nSteps = 0;
  while ((ForwardInTime && time < endtime) ||
         (!ForwardInTime && y[0]>v_0)) {
    // Take a step
    int status = gsl_odeiv2_evolve_apply(e, c, s, sys, &time, time+hnext, &h, &y[0]);
    ++NSteps;
    ++nSteps;

    // Check if it worked and the system is still reasonable
    if(status == GSL_EDOM) {
      INFOTOCOUT << "Velocity v has become greater than 1.0.  This is a nice way for PN to stop." << std::endl;
      break;
    } else if(status == GSL_EDIVERGE) {
      INFOTOCOUT << "Velocity is no longer increasing.  This is not an uncommon way for PN to stop, and probably okay." << std::endl;
      break;
    } else if(status != GSL_SUCCESS) {
      INFOTOCERR << "GSL odeiv2 error.  Return value=" << status << "\nThis is potentially a very bad way for PN to stop." << std::endl;
      break;
    }

    // If the time step was large enough, store the data
    Tn->Recalculate(time, &y[0]);
    if(std::abs(t.back()-time)>=std::abs(hmin_storage)) {
      t.push_back(time);
      v.push_back(y[0]);
      const Quaternion R_chi1_i = exp(Quaternion(0.0, y[1], y[2], 0.0));
      chi1.push_back((S_chi1_i*R_chi1_i*zHat*R_chi1_i.conjugate()*S_chi1_i.conjugate()).vec());
      const Quaternion R_chi2_i = exp(Quaternion(0.0, y[3], y[4], 0.0));
      chi2.push_back((S_chi2_i*R_chi2_i*zHat*R_chi2_i.conjugate()*S_chi2_i.conjugate()).vec());
      const Quaternion R_frame_i = Unflipped(R_frame.back(), exp(Quaternion(0.0, y[5], y[6], y[7])));
      R_frame.push_back(R_frame_i);
      Phi.push_back(y[8]);
      L.push_back((Tn->OrbitalAngularMomentum()).vec());
    }

    // Check if we should stop because this has gone on suspiciously long
    if((ForwardInTime && time>=endtime)
       || (!ForwardInTime && time<endtime)) {
      INFOTOCERR << "Time has gone on four times as long as expected.  This seems strange, so we'll stop."
                 << "\nNote that this is unusual.  You may have a short waveform that stops before merger,"
                 << "\nor one of the stopping criteria may have gotten fooled." << std::endl;
      break;
    }

    // Check if we should stop because there have been too many steps
    if(NSteps>MaxSteps) {
      INFOTOCERR << "\n\nThe integration has taken " << NSteps << ".  This seems excessive, so we'll stop."
                 << "\nNote that this is unusual.  You may have a short waveform that stops before merger." << std::endl;
      break;
    }

    // Check if we should stop because the step has gotten too small,
    // but make sure we at least take 500 steps since the last
    // disruption.  [This is the condition that we expect to stop us
    // near merger.]
    if(nSteps>500 && std::abs(h)<std::abs(hmin)) {
      INFOTOCOUT << "The step size " << h << " has become smaller than the lower limit of " << hmin
                 << ".\nThis is probably fine, as it indicates the evolution has naturally reached its end." << std::endl;
      break;
    }

    // Set the next time-step size (actually just an upper limit)
    if(MinStepsPerOrbit!=0) {
      const double OrbitalPeriod = 2*M_PI/(y[0]*y[0]*y[0]);
      hnext = (ForwardInTime
               ? std::min(h, std::min(hmax,  OrbitalPeriod/(MinStepsPerOrbit+1.)))
               : std::max(h, std::max(hmax, -OrbitalPeriod/(MinStepsPerOrbit+1.))) );
      // INFOTOCERR << time << " " << h << " " << hnext << " " << y[0] << " " << NSteps << " " << nSteps << std::endl;
    }

    // Reset values of quaternion logarithms to smaller sizes, if
    // necessary.  If this resets, we reset nSteps to zero, because
    // this may make the time stepper take smaller steps
    const double rfrakMag_chi1 = std::sqrt(y[1]*y[1]+y[2]*y[2]);
    if(rfrakMag_chi1>M_PI/2.) {
      y[1] = (rfrakMag_chi1-M_PI)*y[1]/rfrakMag_chi1;
      y[2] = (rfrakMag_chi1-M_PI)*y[2]/rfrakMag_chi1;
      nSteps=0; // This may make the integrator take a few small steps at first
      gsl_odeiv2_evolve_reset(e); // This should be called whenever the next use of `e` will not be a continuation of the previous step
      // INFOTOCERR << time << std::endl;
    }
    const double rfrakMag_chi2 = std::sqrt(y[3]*y[3]+y[4]*y[4]);
    if(rfrakMag_chi2>M_PI/2.) {
      y[3] = (rfrakMag_chi2-M_PI)*y[3]/rfrakMag_chi2;
      y[4] = (rfrakMag_chi2-M_PI)*y[4]/rfrakMag_chi2;
      nSteps=0; // This may make the integrator take a few small steps at first
      gsl_odeiv2_evolve_reset(e); // This should be called whenever the next use of `e` will not be a continuation of the previous step
      // INFOTOCERR << time << std::endl;
    }
    const double rfrakMag_ellHat = std::sqrt(y[5]*y[5]+y[6]*y[6]+y[7]*y[7]);
    if(rfrakMag_ellHat>M_PI/2.) {
      y[5] = (rfrakMag_ellHat-M_PI)*y[5]/rfrakMag_ellHat;
      y[6] = (rfrakMag_ellHat-M_PI)*y[6]/rfrakMag_ellHat;
      y[7] = (rfrakMag_ellHat-M_PI)*y[7]/rfrakMag_ellHat;
      nSteps=0; // This may make the integrator take a few small steps at first
      gsl_odeiv2_evolve_reset(e); // This should be called whenever the next use of `e` will not be a continuation of the previous step
      // INFOTOCERR << time << std::endl;
    }
  }

  if( ! ((ForwardInTime && time < endtime) || (!ForwardInTime && y[0]<v_0)) ) {
    if(ForwardInTime) {
      INFOTOCERR << "Stepping ended because the time " << time << " has exceeded the endtime of " << endtime << ".\n"
                 << "This could possibly indicate that the equations are not well behaved." << std::endl;
    } else {
      INFOTOCERR << std::setprecision(15)
                 << "Stepping stopped even though backwards evolution has only reached v=" << y[0] << ", which is greater than the target v_0=" << v_0
                 << ".\nThis is short of the expected stopping criterion, and may indicate a problem." << std::endl;
    }
  }

  // Free the gsl storage
  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);

  // Combine two segments if we need to also run backwards
  if(v_0<v_i && ForwardInTime) {
    // Run the reverse evolution
    std::vector<double> tBackward, vBackward, PhiBackward;
    std::vector<std::vector<double> > chi1Backward, chi2Backward, LBackward;
    std::vector<Quaternions::Quaternion> R_frameBackward;
    EvolvePN_Q(Approximant, PNOrbitalEvolutionOrder, v_0, v_i, m1, m2, chi1_i, chi2_i, R_frame_i,
               tBackward, vBackward, chi1Backward, chi2Backward, R_frameBackward, PhiBackward, LBackward,
               MinStepsPerOrbit, false);

    // Remove the first element of each of the backward data vectors
    // so that point is not duplicated
    tBackward.erase(tBackward.begin());
    vBackward.erase(vBackward.begin());
    PhiBackward.erase(PhiBackward.begin());
    chi1Backward.erase(chi1Backward.begin());
    chi2Backward.erase(chi2Backward.begin());
    LBackward.erase(LBackward.begin());
    R_frameBackward.erase(R_frameBackward.begin());

    // Combine the data
    CombineForwardAndBackward(t, tBackward);
    CombineForwardAndBackward(v, vBackward);
    CombineForwardAndBackward(Phi, PhiBackward);
    CombineForwardAndBackward(chi1, chi1Backward);
    CombineForwardAndBackward(chi2, chi2Backward);
    CombineForwardAndBackward(L, LBackward);
    CombineForwardAndBackward(R_frame, R_frameBackward);
  }

  // // Make the last time = 0.0
  // if(ForwardInTime) {
  //   const double tback = t.back();
  //   for(unsigned int i=0; i<t.size(); ++i) {
  //     t[i] -= tback;
  //   }
  // }

  delete Tn;

  return;
}
