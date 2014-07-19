#ifndef PNWAVEFORMMODES_HPP
#define PNWAVEFORMMODES_HPP

#include <vector>
#include <complex>
#include "Quaternions.hpp"

namespace PostNewtonian {

  /// These functions return all known modes of a post-Newtonian
  /// waveform, given the present speed and spins of the black holes.
  /// The BHs are assumed to be along the x axis, with individual
  /// velocities along the y axis.  (That is, the modes are given in
  /// the corotating frame.)  This waveform will need to be
  /// transformed back into an inertial frame for many applications.
  ///
  /// The returned data are given as vectors of complex numbers,
  /// representing mode (2,-2), followed by (2,-1), etc., all the way
  /// up to the highest known mode -- currently (8,8).

  std::vector<std::vector<std::complex<double> > > WaveformModes
  (const double m1, const double m2, const std::vector<double>& v,
   const std::vector<std::vector<double> >& chi1, const std::vector<std::vector<double> >& chi2,
   const double PNWaveformModeOrder=3.5);

};



#endif // PNWAVEFORMMODES_HPP
