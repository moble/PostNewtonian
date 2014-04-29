#include "PNWaveformModes.hpp"
#include "PNWaveformModes.ipp"

std::vector<std::vector<std::complex<double> > > PostNewtonian::WaveformModes
(const double m1, const double m2, const std::vector<double>& v,
 const std::vector<std::vector<double> >& chi1, const std::vector<std::vector<double> >& chi2) {
  WaveformModes_3p5PN WM(m1, m2, v[0],
                         chi1[0][0], chi1[0][1], chi1[0][2],
                         chi2[0][0], chi2[0][1], chi2[0][2],
                         0., 0., 1., 1., 0., 0.);
  const unsigned int N_lm = ellMax*(ellMax+2) - 3;
  const unsigned int N_t = v.size();
  std::vector<std::vector<std::complex<double> > > Modes(N_lm, std::vector<std::complex<double> >(N_t));
  std::vector<std::complex<double> > modes(N_lm);
  for(unsigned int i_t=0; i_t<N_t; ++i_t) {
    modes = WM(v[i_t], chi1[i_t], chi2[i_t]);
    for(unsigned int i_lm=0; i_lm<N_lm; ++i_lm) {
      Modes[i_lm][i_t] = modes[i_lm];
    }
  }
  return Modes;
}
