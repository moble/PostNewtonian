#include "PNWaveformModes.hpp"
#include "PNWaveformModes.ipp"

using PostNewtonian::WaveformModes;

std::vector<std::vector<std::complex<double> > > WaveformModes
(const double m1, const double m2, const std::vector<double>& v,
 const std::vector<std::vector<double> >& chi1, const std::vector<std::vector<double> >& chi2) {
  WaveformModes_3p5PN WM(m1, m2, v[0],
			 chi1[0][0], chi1[0][1], chi1[0][2],
			 chi2[0][0], chi2[0][1], chi2[0][2],
			 0.0, 0.0, 0.0);
  std::vector<std::vector<std::complex<double> > > Modes(v.size(), std::vector<std::complex<double> >(ellMax*(ellMax+2) - 3));
  for(unsigned int i=0; i<v.size(); ++i) {
    Modes[i] = WM(v[i], chi1[i], chi2[i]);
  }
  return Modes;
}
