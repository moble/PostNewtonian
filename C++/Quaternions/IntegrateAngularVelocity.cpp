// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#include <cmath>
#include <iostream>
#include <iomanip>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_odeiv2.h>
// Are these necessary?
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

#include "Quaternions.hpp"
#include "IntegrateAngularVelocity.hpp"
using Quaternions::Quaternion;
using Quaternions::QuaternionArray;

// Note: Don't do 'using namespace std' because we don't want to
// confuse which log, exp, etc., is being used in any instance.
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::vector;


// Define some error codes, which will be used in python
#define FailedGSLCall 6


#ifndef DOXYGEN
typedef struct {
  gsl_interp_accel* accX;
  gsl_spline* splineX;
  gsl_interp_accel* accY;
  gsl_spline* splineY;
  gsl_interp_accel* accZ;
  gsl_spline* splineZ;
} ParameterSet;
int FrameFromAngularVelocity_RHS(double t, const double ri[], double drdt[], void* mu) {
  // Interpolate and unpack some values
  ParameterSet* params = (ParameterSet*) mu;
  vector<double> rfrak(3);
  rfrak[0] = ri[0];
  rfrak[1] = ri[1];
  rfrak[2] = ri[2];
  vector<double> Omega(3);
  Omega[0] = gsl_spline_eval(params->splineX, t, params->accX);
  Omega[1] = gsl_spline_eval(params->splineY, t, params->accY);
  Omega[2] = gsl_spline_eval(params->splineZ, t, params->accZ);
  // Evaluate the RHS and unpack
  const std::vector<double> rfrakDot = Quaternions::FrameFromAngularVelocity_Integrand(rfrak, Omega);
  drdt[0] = rfrakDot[0];
  drdt[1] = rfrakDot[1];
  drdt[2] = rfrakDot[2];
  // GSL wants to hear that everything went okay
  return GSL_SUCCESS;
}
#endif // DOXYGEN
std::vector<Quaternion> Quaternions::FrameFromAngularVelocity(const std::vector<Quaternion>& Omega, const std::vector<double>& T) {
  ///
  /// \param Omega Vector of Quaternions.
  /// \param T Vector of corresponding times.
  ///
  /// Note that each element of Omega should be a pure-vector
  /// Quaternion, corresponding to the angular-velocity vector at the
  /// instant of time.
  const vector<double> OmegaX = Quaternions::Component1(Omega);
  const vector<double> OmegaY = Quaternions::Component2(Omega);
  const vector<double> OmegaZ = Quaternions::Component3(Omega);

  // Set up the spline to interpolate Omega
  gsl_interp_accel* accX = gsl_interp_accel_alloc();
  gsl_spline* splineX = gsl_spline_alloc(gsl_interp_cspline, OmegaX.size());
  gsl_spline_init(splineX, &T[0], &OmegaX[0], OmegaX.size());
  gsl_interp_accel* accY = gsl_interp_accel_alloc();
  gsl_spline* splineY = gsl_spline_alloc(gsl_interp_cspline, OmegaY.size());
  gsl_spline_init(splineY, &T[0], &OmegaY[0], OmegaY.size());
  gsl_interp_accel* accZ = gsl_interp_accel_alloc();
  gsl_spline* splineZ = gsl_spline_alloc(gsl_interp_cspline, OmegaZ.size());
  gsl_spline_init(splineZ, &T[0], &OmegaZ[0], OmegaZ.size());

  // Set up the integrator
  const double hstart = (T[1]-T[0])/10.;
  const double epsabs = 1.e-9;
  const double epsrel = 1.e-9;
  ParameterSet params;
  params.accX = accX;
  params.splineX = splineX;
  params.accY = accY;
  params.splineY = splineY;
  params.accZ = accZ;
  params.splineZ = splineZ;
  gsl_odeiv2_system sys = {FrameFromAngularVelocity_RHS, NULL, 3, (void *) &params};
  // gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, hstart, epsabs, epsrel);
  double t = T[0];
  double r[3] = {0.0, 0.0, 0.0};

  // Run the integration
  vector<vector<double> > rs;
  rs.push_back(vector<double>(r, r+3));
  for(unsigned int i=1; i<T.size(); ++i) {
    int status = gsl_odeiv2_driver_apply(d, &t, T[i], r);
    if(status != GSL_SUCCESS) {
      gsl_spline_free(splineX);
      gsl_interp_accel_free(accX);
      gsl_spline_free(splineY);
      gsl_interp_accel_free(accY);
      gsl_spline_free(splineZ);
      gsl_interp_accel_free(accZ);
      gsl_odeiv2_driver_free(d);
      // cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": gsl_odeiv2_driver_apply returned an error; return value=" << status << "\n\n";
      throw(FailedGSLCall);
    }
    rs.push_back(vector<double>(r, r+3));
    // Reset the value of r if necessary
    const double absquatlogR = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    if(absquatlogR>M_PI/2.) {
      r[0] = (absquatlogR-M_PI)*r[0]/absquatlogR;
      r[1] = (absquatlogR-M_PI)*r[1]/absquatlogR;
      r[2] = (absquatlogR-M_PI)*r[2]/absquatlogR;
      // cerr << "Flipping r at t=" << t << endl;
    }
  }

  // Free the gsl storage
  gsl_spline_free(splineX);
  gsl_interp_accel_free(accX);
  gsl_spline_free(splineY);
  gsl_interp_accel_free(accY);
  gsl_spline_free(splineZ);
  gsl_interp_accel_free(accZ);
  gsl_odeiv2_driver_free(d);

  return UnflipRotors(exp(QuaternionArray(rs)));
}
