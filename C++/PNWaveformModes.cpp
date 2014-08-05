#include "PNWaveformModes.hpp"
#include "PNWaveformModes.ipp"


std::vector<std::vector<std::complex<double> > > PostNewtonian::WaveformModes
(const double m1, const double m2, const std::vector<double>& v,
 const std::vector<std::vector<double> >& chi1, const std::vector<std::vector<double> >& chi2, const double PNWaveformModeOrder) {
  WaveformModes_Base* WM;
  switch(int(2*PNWaveformModeOrder)) {
  case 0:
    WM = new WaveformModes_0PN(m1, m2, v[0]);
    break;
  case 1:
    WM = new WaveformModes_0p50PN(m1, m2, v[0]);
    break;
  case 2:
    WM = new WaveformModes_1p0PN(m1, m2, v[0],
                                 Quaternions::Quaternion(0.0, chi1[0][0], chi1[0][1], chi1[0][2]),
                                 Quaternions::Quaternion(0.0, chi2[0][0], chi2[0][1], chi2[0][2]));
    break;
  case 3:
    WM = new WaveformModes_1p5PN(m1, m2, v[0],
                                 Quaternions::Quaternion(0.0, chi1[0][0], chi1[0][1], chi1[0][2]),
                                 Quaternions::Quaternion(0.0, chi2[0][0], chi2[0][1], chi2[0][2]));
    break;
  case 4:
    WM = new WaveformModes_2p0PN(m1, m2, v[0],
                                 Quaternions::Quaternion(0.0, chi1[0][0], chi1[0][1], chi1[0][2]),
                                 Quaternions::Quaternion(0.0, chi2[0][0], chi2[0][1], chi2[0][2]));
    break;
  case 5:
    WM = new WaveformModes_2p5PN(m1, m2, v[0],
                                 Quaternions::Quaternion(0.0, chi1[0][0], chi1[0][1], chi1[0][2]),
                                 Quaternions::Quaternion(0.0, chi2[0][0], chi2[0][1], chi2[0][2]));
    break;
  case 6:
    WM = new WaveformModes_3p0PN(m1, m2, v[0],
                                 Quaternions::Quaternion(0.0, chi1[0][0], chi1[0][1], chi1[0][2]),
                                 Quaternions::Quaternion(0.0, chi2[0][0], chi2[0][1], chi2[0][2]));
    break;
  case 7:
    WM = new WaveformModes_3p5PN(m1, m2, v[0],
                                 Quaternions::Quaternion(0.0, chi1[0][0], chi1[0][1], chi1[0][2]),
                                 Quaternions::Quaternion(0.0, chi2[0][0], chi2[0][1], chi2[0][2]));
    break;
  default:
    std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << ": PNWaveformModeOrder " << PNWaveformModeOrder << " not yet implemented." << std::endl;
    throw(0);
  }

  const unsigned int N_lm = ellMax*(ellMax+2) - 3;
  const unsigned int N_t = v.size();
  std::vector<std::vector<std::complex<double> > > Modes(N_lm, std::vector<std::complex<double> >(N_t));
  std::vector<std::complex<double> > modes(N_lm);
  for(unsigned int i_t=0; i_t<N_t; ++i_t) {
    modes = (*WM)(v[i_t], chi1[i_t], chi2[i_t]);
    for(unsigned int i_lm=0; i_lm<N_lm; ++i_lm) {
      Modes[i_lm][i_t] = modes[i_lm];
    }
  }
  delete WM;
  return Modes;
}
