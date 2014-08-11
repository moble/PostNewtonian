// File produced automatically by WaveformModeCodeGen.ipynb

class WaveformModes_Base {
public:
  virtual std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k) = 0;
  virtual std::vector<std::complex<double> > operator()(
    const double v_k, const std::vector<double>& chi1, const std::vector<double>& chi2)
  {
    return this->operator()(v_k, chi1[0], chi1[1], chi1[2], chi2[0], chi2[1], chi2[2]);
  }
};

const unsigned int ellMax = 8;
const std::complex<double> I(0,1.0);
inline std::complex<double> conjugate(const std::complex<double>& a) { return std::conj(a); }

class WaveformModes_0PN : public WaveformModes_Base {
private:
  const double M1, M2;
  double v;
  const double M, nu;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_2_0, hHat_4_0_0;

public:
  WaveformModes_0PN(const double M1_i, const double M2_i, const double v_i) :
    M1(M1_i), M2(M2_i), v(v_i), M(M1 + M2), nu(M1*M2/pow(M, 2)), rhOverM_coeff(6.34132367616962*nu*pow(v, 2)),
    hHat_2_0_0(-0.145802960879951), hHat_2_2_0(1.00000000000000), hHat_4_0_0(-0.00140298964521140)
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;

    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = hHat_2_0_0*rhOverM_coeff;
    Asymm = 0;
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = hHat_2_2_0*rhOverM_coeff;
    Asymm = 0;
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = 0;
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = 0;
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_0PN : public WaveformModes_Base


class WaveformModes_0p50PN : public WaveformModes_Base {
private:
  const double M1, M2;
  double v;
  const double M, delta, nu;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_1_1, hHat_2_2_0, hHat_3_1_1, hHat_3_3_1, hHat_4_0_0;

public:
  WaveformModes_0p50PN(const double M1_i, const double M2_i, const double v_i) :
    M1(M1_i), M2(M2_i), v(v_i), M(M1 + M2), delta((M1 - M2)/M), nu(M1*M2/pow(M, 2)),
    rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), hHat_2_0_0(-0.145802960879951), hHat_2_1_1(0.333333333333333*I*delta),
    hHat_2_2_0(1.00000000000000), hHat_3_1_1(0.0222717701593687*I*delta), hHat_3_3_1(-0.776323754260148*I*delta),
    hHat_4_0_0(-0.00140298964521140)
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;

    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = hHat_2_0_0*rhOverM_coeff;
    Asymm = 0;
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = hHat_2_1_1*rhOverM_coeff*v;
    Asymm = 0;
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = hHat_2_2_0*rhOverM_coeff;
    Asymm = 0;
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = hHat_3_1_1*rhOverM_coeff*v;
    Asymm = 0;
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = hHat_3_3_1*rhOverM_coeff*v;
    Asymm = 0;
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = 0;
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = 0;
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_0p50PN : public WaveformModes_Base


class WaveformModes_1p0PN : public WaveformModes_Base {
private:
  const double M1, M2;
  double v;
  const double M, delta, nu;
  Quaternions::Quaternion chiVec1, chiVec2;
  double chi1_n, chi1_lambda, chi1_ell, chi2_n, chi2_lambda, chi2_ell, Sigma_ell, Sigma_n, Sigma_lambda;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_1_1, hHat_2_2_0, hHat_2_2_2, hHat_3_1_1, hHat_3_2_2, hHat_3_3_1,
                             hHat_4_0_0, hHat_4_2_2, hHat_4_4_2;
  std::complex<double> hHat_spin_Symm_2_1_2, hHat_spin_Asymm_2_2_2, hHat_spin_Asymm_2_0_2;

public:
  WaveformModes_1p0PN(const double M1_i, const double M2_i, const double v_i, const Quaternions::Quaternion chiVec1_i,
                      const Quaternions::Quaternion chiVec2_i) :
    M1(M1_i), M2(M2_i), v(v_i), M(M1 + M2), delta((M1 - M2)/M), nu(M1*M2/pow(M, 2)), chiVec1(chiVec1_i),
    chiVec2(chiVec2_i), chi1_n(chiVec1[1]), chi1_lambda(chiVec1[2]), chi1_ell(chiVec1[3]), chi2_n(chiVec2[1]),
    chi2_lambda(chiVec2[2]), chi2_ell(chiVec2[3]), Sigma_ell(M*(-M1*chi1_ell + M2*chi2_ell)), Sigma_n(M*(-M1*chi1_n +
    M2*chi2_n)), Sigma_lambda(M*(-M1*chi1_lambda + M2*chi2_lambda)), rhOverM_coeff(6.34132367616962*nu*pow(v, 2)),
    hHat_2_0_0(-0.145802960879951), hHat_2_1_1(0.333333333333333*I*delta), hHat_2_2_0(1.00000000000000),
    hHat_2_2_2(1.30952380952381*nu - 2.54761904761905), hHat_3_1_1(0.0222717701593687*I*delta),
    hHat_3_2_2(-0.845154254728516*nu + 0.281718084909506), hHat_3_3_1(-0.776323754260148*I*delta),
    hHat_4_0_0(-0.00140298964521140), hHat_4_2_2(-0.10647942749999*nu + 0.0354931424999967),
    hHat_4_4_2(2.25374467927604*nu - 0.751248226425348), hHat_spin_Symm_2_1_2(0.5*I*Sigma_ell/pow(M, 2)),
    hHat_spin_Asymm_2_2_2(0.5*(-Sigma_lambda - I*Sigma_n)/pow(M, 2)),
    hHat_spin_Asymm_2_0_2(0.408248290463863*I*Sigma_n/pow(M, 2))
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;
    chiVec1 = Quaternions::Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
    chiVec2 = Quaternions::Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

    chi1_n = chiVec1[1];
    chi1_lambda = chiVec1[2];
    chi1_ell = chiVec1[3];
    chi2_n = chiVec2[1];
    chi2_lambda = chiVec2[2];
    chi2_ell = chiVec2[3];
    Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell);
    Sigma_n = M*(-M1*chi1_n + M2*chi2_n);
    Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda);
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    hHat_spin_Symm_2_1_2 = 0.5*I*Sigma_ell/pow(M, 2);
    hHat_spin_Asymm_2_2_2 = 0.5*(-Sigma_lambda - I*Sigma_n)/pow(M, 2);
    hHat_spin_Asymm_2_0_2 = 0.408248290463863*I*Sigma_n/pow(M, 2);

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = hHat_2_0_0*rhOverM_coeff;
    Asymm = hHat_spin_Asymm_2_0_2*rhOverM_coeff*pow(v, 2);
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_2_1_1 + hHat_spin_Symm_2_1_2*v);
    Asymm = 0;
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = rhOverM_coeff*(hHat_2_2_0 + hHat_2_2_2*pow(v, 2));
    Asymm = hHat_spin_Asymm_2_2_2*rhOverM_coeff*pow(v, 2);
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = hHat_3_1_1*rhOverM_coeff*v;
    Asymm = 0;
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = hHat_3_2_2*rhOverM_coeff*pow(v, 2);
    Asymm = 0;
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = hHat_3_3_1*rhOverM_coeff*v;
    Asymm = 0;
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = 0;
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = hHat_4_2_2*rhOverM_coeff*pow(v, 2);
    Asymm = 0;
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = hHat_4_4_2*rhOverM_coeff*pow(v, 2);
    Asymm = 0;
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = 0;
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_1p0PN : public WaveformModes_Base


class WaveformModes_1p5PN : public WaveformModes_Base {
private:
  const double M1, M2;
  double v;
  const double M, delta, nu;
  Quaternions::Quaternion chiVec1, chiVec2;
  double chi1_n, chi1_lambda, chi1_ell, chi2_n, chi2_lambda, chi2_ell, S_ell, S_n, S_lambda, Sigma_ell, Sigma_n,
         Sigma_lambda;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_1_1, hHat_2_1_3, hHat_2_2_0, hHat_2_2_2, hHat_2_2_3, hHat_3_1_1,
                             hHat_3_1_3, hHat_3_2_2, hHat_3_3_1, hHat_3_3_3, hHat_4_0_0, hHat_4_1_3, hHat_4_2_2,
                             hHat_4_3_3, hHat_4_4_2, hHat_5_1_3, hHat_5_3_3, hHat_5_5_3;
  std::complex<double> hHat_spin_Symm_2_2_3, hHat_spin_Symm_2_1_2, hHat_spin_Symm_2_0_3, hHat_spin_Symm_3_2_3,
                       hHat_spin_Asymm_2_2_2, hHat_spin_Asymm_2_1_3, hHat_spin_Asymm_2_0_2, hHat_spin_Asymm_3_3_3,
                       hHat_spin_Asymm_3_1_3;

public:
  WaveformModes_1p5PN(const double M1_i, const double M2_i, const double v_i, const Quaternions::Quaternion chiVec1_i,
                      const Quaternions::Quaternion chiVec2_i) :
    M1(M1_i), M2(M2_i), v(v_i), M(M1 + M2), delta((M1 - M2)/M), nu(M1*M2/pow(M, 2)), chiVec1(chiVec1_i),
    chiVec2(chiVec2_i), chi1_n(chiVec1[1]), chi1_lambda(chiVec1[2]), chi1_ell(chiVec1[3]), chi2_n(chiVec2[1]),
    chi2_lambda(chiVec2[2]), chi2_ell(chiVec2[3]), S_ell(pow(M1, 2)*chi1_ell + pow(M2, 2)*chi2_ell), S_n(pow(M1,
    2)*chi1_n + pow(M2, 2)*chi2_n), S_lambda(pow(M1, 2)*chi1_lambda + pow(M2, 2)*chi2_lambda), Sigma_ell(M*(-M1*chi1_ell
    + M2*chi2_ell)), Sigma_n(M*(-M1*chi1_n + M2*chi2_n)), Sigma_lambda(M*(-M1*chi1_lambda + M2*chi2_lambda)),
    rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), hHat_2_0_0(-0.145802960879951), hHat_2_1_1(0.333333333333333*I*delta),
    hHat_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)), hHat_2_2_0(1.00000000000000),
    hHat_2_2_2(1.30952380952381*nu - 2.54761904761905), hHat_2_2_3(6.28318530717959),
    hHat_3_1_1(0.0222717701593687*I*delta), hHat_3_1_3(-0.0148478467729125*I*delta*(nu + 4.0)),
    hHat_3_2_2(-0.845154254728516*nu + 0.281718084909506), hHat_3_3_1(-0.776323754260148*I*delta),
    hHat_3_3_3(-1.5526475085203*I*delta*(nu - 2.0)), hHat_4_0_0(-0.00140298964521140),
    hHat_4_1_3(0.00376461626210521*I*delta*(-2.0*nu + 1.0)), hHat_4_2_2(-0.10647942749999*nu + 0.0354931424999967),
    hHat_4_3_3(0.268926437100239*I*delta*(2.0*nu - 1.0)), hHat_4_4_2(2.25374467927604*nu - 0.751248226425348),
    hHat_5_1_3(0.000176960830360287*I*delta*(-2.0*nu + 1.0)), hHat_5_3_3(0.0464469088412683*I*delta*(2.0*nu - 1.0)),
    hHat_5_5_3(-0.801376894396698*I*delta*(2.0*nu - 1.0)), hHat_spin_Symm_2_2_3(0.166666666666667*(3.0*S_ell +
    5.0*Sigma_ell*delta)/pow(M, 2)), hHat_spin_Symm_2_1_2(0.5*I*Sigma_ell/pow(M, 2)),
    hHat_spin_Symm_2_0_3(0.408248290463863*(5.0*S_ell + 3.0*Sigma_ell*delta)/pow(M, 2)),
    hHat_spin_Symm_3_2_3(0.563436169819011*(S_ell + Sigma_ell*delta)/pow(M, 2)),
    hHat_spin_Asymm_2_2_2(0.5*(-Sigma_lambda - I*Sigma_n)/pow(M, 2)),
    hHat_spin_Asymm_2_1_3(0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
    13.0*Sigma_n*delta)/pow(M, 2)), hHat_spin_Asymm_2_0_2(0.408248290463863*I*Sigma_n/pow(M, 2)),
    hHat_spin_Asymm_3_3_3(0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(M, 2)),
    hHat_spin_Asymm_3_1_3(0.17817416127495*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/pow(M, 2))
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;
    chiVec1 = Quaternions::Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
    chiVec2 = Quaternions::Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

    chi1_n = chiVec1[1];
    chi1_lambda = chiVec1[2];
    chi1_ell = chiVec1[3];
    chi2_n = chiVec2[1];
    chi2_lambda = chiVec2[2];
    chi2_ell = chiVec2[3];
    S_ell = pow(M1, 2)*chi1_ell + pow(M2, 2)*chi2_ell;
    S_n = pow(M1, 2)*chi1_n + pow(M2, 2)*chi2_n;
    S_lambda = pow(M1, 2)*chi1_lambda + pow(M2, 2)*chi2_lambda;
    Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell);
    Sigma_n = M*(-M1*chi1_n + M2*chi2_n);
    Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda);
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    hHat_spin_Symm_2_2_3 = 0.166666666666667*(3.0*S_ell + 5.0*Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_2_1_2 = 0.5*I*Sigma_ell/pow(M, 2);
    hHat_spin_Symm_2_0_3 = 0.408248290463863*(5.0*S_ell + 3.0*Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_3_2_3 = 0.563436169819011*(S_ell + Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Asymm_2_2_2 = 0.5*(-Sigma_lambda - I*Sigma_n)/pow(M, 2);
    hHat_spin_Asymm_2_1_3 = 0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
      13.0*Sigma_n*delta)/pow(M, 2);
    hHat_spin_Asymm_2_0_2 = 0.408248290463863*I*Sigma_n/pow(M, 2);
    hHat_spin_Asymm_3_3_3 = 0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(M, 2);
    hHat_spin_Asymm_3_1_3 = 0.17817416127495*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/pow(M, 2);

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_3*pow(v, 3));
    Asymm = hHat_spin_Asymm_2_0_2*rhOverM_coeff*pow(v, 2);
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_2_1_3*v + hHat_spin_Symm_2_1_2));
    Asymm = hHat_spin_Asymm_2_1_3*rhOverM_coeff*pow(v, 3);
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3)));
    Asymm = hHat_spin_Asymm_2_2_2*rhOverM_coeff*pow(v, 2);
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_3_1_1 + hHat_3_1_3*pow(v, 2));
    Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3);
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 + hHat_spin_Symm_3_2_3*v);
    Asymm = 0;
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = rhOverM_coeff*v*(hHat_3_3_1 + hHat_3_3_3*pow(v, 2));
    Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3);
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = 0;
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = hHat_4_1_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = hHat_4_2_2*rhOverM_coeff*pow(v, 2);
    Asymm = 0;
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = hHat_4_3_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = hHat_4_4_2*rhOverM_coeff*pow(v, 2);
    Asymm = 0;
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = hHat_5_1_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = hHat_5_3_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = hHat_5_5_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = 0;
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_1p5PN : public WaveformModes_Base


class WaveformModes_2p0PN : public WaveformModes_Base {
private:
  const double M1, M2;
  double v;
  const double M, delta, nu;
  Quaternions::Quaternion chiVec1, chiVec2;
  double chi1_n, chi1_lambda, chi1_ell, chi2_n, chi2_lambda, chi2_ell, S_ell, S_n, S_lambda, Sigma_ell, Sigma_n,
         Sigma_lambda, S1_ell, S1_n, S1_lambda, S2_ell, S2_n, S2_lambda;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_1_1, hHat_2_1_3, hHat_2_1_4, hHat_2_2_0, hHat_2_2_2, hHat_2_2_3,
                             hHat_2_2_4, hHat_3_1_1, hHat_3_1_3, hHat_3_1_4, hHat_3_2_2, hHat_3_2_4, hHat_3_3_1,
                             hHat_3_3_3, hHat_3_3_4, hHat_4_0_0, hHat_4_1_3, hHat_4_2_2, hHat_4_2_4, hHat_4_3_3,
                             hHat_4_4_2, hHat_4_4_4, hHat_5_1_3, hHat_5_2_4, hHat_5_3_3, hHat_5_4_4, hHat_5_5_3,
                             hHat_6_2_4, hHat_6_4_4, hHat_6_6_4;
  std::complex<double> hHat_spin_Symm_2_2_3, hHat_spin_Symm_2_2_4, hHat_spin_Symm_2_1_2, hHat_spin_Symm_2_1_4,
                       hHat_spin_Symm_2_0_3, hHat_spin_Symm_2_0_4, hHat_spin_Symm_3_3_4, hHat_spin_Symm_3_2_3,
                       hHat_spin_Symm_3_1_4, hHat_spin_Symm_4_3_4, hHat_spin_Symm_4_1_4, hHat_spin_Asymm_2_2_2,
                       hHat_spin_Asymm_2_2_4, hHat_spin_Asymm_2_1_3, hHat_spin_Asymm_2_1_4, hHat_spin_Asymm_2_0_2,
                       hHat_spin_Asymm_2_0_4, hHat_spin_Asymm_3_3_3, hHat_spin_Asymm_3_2_4, hHat_spin_Asymm_3_1_3,
                       hHat_spin_Asymm_3_0_4, hHat_spin_Asymm_4_4_4, hHat_spin_Asymm_4_2_4, hHat_spin_Asymm_4_0_4;

public:
  WaveformModes_2p0PN(const double M1_i, const double M2_i, const double v_i, const Quaternions::Quaternion chiVec1_i,
                      const Quaternions::Quaternion chiVec2_i) :
    M1(M1_i), M2(M2_i), v(v_i), M(M1 + M2), delta((M1 - M2)/M), nu(M1*M2/pow(M, 2)), chiVec1(chiVec1_i),
    chiVec2(chiVec2_i), chi1_n(chiVec1[1]), chi1_lambda(chiVec1[2]), chi1_ell(chiVec1[3]), chi2_n(chiVec2[1]),
    chi2_lambda(chiVec2[2]), chi2_ell(chiVec2[3]), S_ell(pow(M1, 2)*chi1_ell + pow(M2, 2)*chi2_ell), S_n(pow(M1,
    2)*chi1_n + pow(M2, 2)*chi2_n), S_lambda(pow(M1, 2)*chi1_lambda + pow(M2, 2)*chi2_lambda), Sigma_ell(M*(-M1*chi1_ell
    + M2*chi2_ell)), Sigma_n(M*(-M1*chi1_n + M2*chi2_n)), Sigma_lambda(M*(-M1*chi1_lambda + M2*chi2_lambda)),
    S1_ell(pow(M1, 2)*chi1_ell), S1_n(pow(M1, 2)*chi1_n), S1_lambda(pow(M1, 2)*chi1_lambda), S2_ell(pow(M2,
    2)*chi2_ell), S2_n(pow(M2, 2)*chi2_n), S2_lambda(pow(M2, 2)*chi2_lambda), rhOverM_coeff(6.34132367616962*nu*pow(v,
    2)), hHat_2_0_0(-0.145802960879951), hHat_2_1_1(0.333333333333333*I*delta),
    hHat_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)), hHat_2_1_4(delta*(0.628764787039964 + 1.0471975511966*I)),
    hHat_2_2_0(1.00000000000000), hHat_2_2_2(1.30952380952381*nu - 2.54761904761905), hHat_2_2_3(6.28318530717959),
    hHat_2_2_4(0.000661375661375661*nu*(2047.0*nu - 7483.0) - 1.43716931216931), hHat_3_1_1(0.0222717701593687*I*delta),
    hHat_3_1_3(-0.0148478467729125*I*delta*(nu + 4.0)), hHat_3_1_4(delta*(0.0620557076072073 + 0.0699688295151131*I)),
    hHat_3_2_2(-0.845154254728516*nu + 0.281718084909506), hHat_3_2_4(0.00313020094343895*nu*(-365.0*nu + 725.0) -
    0.604128782083717), hHat_3_3_1(-0.776323754260148*I*delta), hHat_3_3_3(-1.5526475085203*I*delta*(nu - 2.0)),
    hHat_3_3_4(delta*(-1.37192659820446 - 7.31667900957279*I)), hHat_4_0_0(-0.00140298964521140),
    hHat_4_1_3(0.00376461626210521*I*delta*(-2.0*nu + 1.0)), hHat_4_2_2(-0.10647942749999*nu + 0.0354931424999967),
    hHat_4_2_4(0.000107554977272717*nu*(-285.0*nu + 4025.0) - 0.141004575204532),
    hHat_4_3_3(0.268926437100239*I*delta*(2.0*nu - 1.0)), hHat_4_4_2(2.25374467927604*nu - 0.751248226425348),
    hHat_4_4_4(0.0113825488852325*nu*(525.0*nu - 1273.0) + 4.04991089336574),
    hHat_5_1_3(0.000176960830360287*I*delta*(-2.0*nu + 1.0)), hHat_5_2_4(0.00998814611056655*nu*(5.0*nu - 5.0) +
    0.00998814611056655), hHat_5_3_3(0.0464469088412683*I*delta*(2.0*nu - 1.0)),
    hHat_5_4_4(-0.276799624590764*nu*(5.0*nu - 5.0) - 0.276799624590764), hHat_5_5_3(-0.801376894396698*I*delta*(2.0*nu
    - 1.0)), hHat_6_2_4(0.000835250737974468*nu*(5.0*nu - 5.0) + 0.000835250737974468),
    hHat_6_4_4(-0.058558165806266*nu*(5.0*nu - 5.0) - 0.058558165806266), hHat_6_6_4(0.903141370807658*nu*(5.0*nu - 5.0)
    + 0.903141370807658), hHat_spin_Symm_2_2_3(0.166666666666667*(3.0*S_ell + 5.0*Sigma_ell*delta)/pow(M, 2)),
    hHat_spin_Symm_2_2_4(0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
    15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(M, 4)*nu)), hHat_spin_Symm_2_1_2(0.5*I*Sigma_ell/pow(M, 2)),
    hHat_spin_Symm_2_1_4(0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(M, 2)),
    hHat_spin_Symm_2_0_3(0.408248290463863*(5.0*S_ell + 3.0*Sigma_ell*delta)/pow(M, 2)),
    hHat_spin_Symm_2_0_4(0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(M, 4)*nu)),
    hHat_spin_Symm_3_3_4(0.388161877130074*I*(7.0*S_ell*delta - 3.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2)),
    hHat_spin_Symm_3_2_3(0.563436169819011*(S_ell + Sigma_ell*delta)/pow(M, 2)),
    hHat_spin_Symm_3_1_4(0.0111358850796843*I*(S_ell*delta - 5.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2)),
    hHat_spin_Symm_4_3_4(0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(M, 2)),
    hHat_spin_Symm_4_1_4(0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(M, 2)),
    hHat_spin_Asymm_2_2_2(0.5*(-Sigma_lambda - I*Sigma_n)/pow(M, 2)),
    hHat_spin_Asymm_2_2_4(0.0119047619047619*(19.0*S_lambda*delta + 182.0*I*S_n*delta - 43.0*Sigma_lambda*nu +
    5.0*Sigma_lambda - 280.0*I*Sigma_n*nu + 98.0*I*Sigma_n)/pow(M, 2)),
    hHat_spin_Asymm_2_1_3(0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
    13.0*Sigma_n*delta)/pow(M, 2)), hHat_spin_Asymm_2_1_4(0.5*(-3.0*S1_ell*S2_n - 3.0*S1_n*S2_ell)/(pow(M, 4)*nu)),
    hHat_spin_Asymm_2_0_2(0.408248290463863*I*Sigma_n/pow(M, 2)),
    hHat_spin_Asymm_2_0_4(0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(M, 2)),
    hHat_spin_Asymm_3_3_3(0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(M, 2)),
    hHat_spin_Asymm_3_2_4(0.0352147606136882*(-Sigma_lambda*(83.0*nu - 17.0) + 4.0*I*Sigma_n*(55.0*nu - 13.0) +
    25.0*delta*(S_lambda - 4.0*I*S_n))/pow(M, 2)), hHat_spin_Asymm_3_1_3(0.17817416127495*(I*S_lambda + S_n +
    delta*(I*Sigma_lambda + Sigma_n))/pow(M, 2)), hHat_spin_Asymm_3_0_4(0.038575837490523*(-17.0*S_lambda*delta +
    Sigma_lambda*(35.0*nu - 9.0))/pow(M, 2)), hHat_spin_Asymm_4_4_4(0.950798536569581*(-3.0*Sigma_lambda*nu +
    Sigma_lambda - I*Sigma_n*(3.0*nu - 1.0) + delta*(S_lambda + I*S_n))/pow(M, 2)),
    hHat_spin_Asymm_4_2_4(0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) + 14.0*I*Sigma_n*(3.0*nu - 1.0) +
    delta*(13.0*S_lambda - 14.0*I*S_n))/pow(M, 2)), hHat_spin_Asymm_4_0_4(0.00841793787126842*I*(S_n*delta -
    3.0*Sigma_n*nu + Sigma_n)/pow(M, 2))
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;
    chiVec1 = Quaternions::Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
    chiVec2 = Quaternions::Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

    chi1_n = chiVec1[1];
    chi1_lambda = chiVec1[2];
    chi1_ell = chiVec1[3];
    chi2_n = chiVec2[1];
    chi2_lambda = chiVec2[2];
    chi2_ell = chiVec2[3];
    S_ell = pow(M1, 2)*chi1_ell + pow(M2, 2)*chi2_ell;
    S_n = pow(M1, 2)*chi1_n + pow(M2, 2)*chi2_n;
    S_lambda = pow(M1, 2)*chi1_lambda + pow(M2, 2)*chi2_lambda;
    Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell);
    Sigma_n = M*(-M1*chi1_n + M2*chi2_n);
    Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda);
    S1_ell = pow(M1, 2)*chi1_ell;
    S1_n = pow(M1, 2)*chi1_n;
    S1_lambda = pow(M1, 2)*chi1_lambda;
    S2_ell = pow(M2, 2)*chi2_ell;
    S2_n = pow(M2, 2)*chi2_n;
    S2_lambda = pow(M2, 2)*chi2_lambda;
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    hHat_spin_Symm_2_2_3 = 0.166666666666667*(3.0*S_ell + 5.0*Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_2_2_4 = 0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
      15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(M, 4)*nu);
    hHat_spin_Symm_2_1_2 = 0.5*I*Sigma_ell/pow(M, 2);
    hHat_spin_Symm_2_1_4 = 0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(M, 2);
    hHat_spin_Symm_2_0_3 = 0.408248290463863*(5.0*S_ell + 3.0*Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_2_0_4 = 0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(M, 4)*nu);
    hHat_spin_Symm_3_3_4 = 0.388161877130074*I*(7.0*S_ell*delta - 3.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2);
    hHat_spin_Symm_3_2_3 = 0.563436169819011*(S_ell + Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_3_1_4 = 0.0111358850796843*I*(S_ell*delta - 5.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2);
    hHat_spin_Symm_4_3_4 = 0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(M, 2);
    hHat_spin_Symm_4_1_4 = 0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(M, 2);
    hHat_spin_Asymm_2_2_2 = 0.5*(-Sigma_lambda - I*Sigma_n)/pow(M, 2);
    hHat_spin_Asymm_2_2_4 = 0.0119047619047619*(19.0*S_lambda*delta + 182.0*I*S_n*delta - 43.0*Sigma_lambda*nu +
      5.0*Sigma_lambda - 280.0*I*Sigma_n*nu + 98.0*I*Sigma_n)/pow(M, 2);
    hHat_spin_Asymm_2_1_3 = 0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
      13.0*Sigma_n*delta)/pow(M, 2);
    hHat_spin_Asymm_2_1_4 = 0.5*(-3.0*S1_ell*S2_n - 3.0*S1_n*S2_ell)/(pow(M, 4)*nu);
    hHat_spin_Asymm_2_0_2 = 0.408248290463863*I*Sigma_n/pow(M, 2);
    hHat_spin_Asymm_2_0_4 = 0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(M, 2);
    hHat_spin_Asymm_3_3_3 = 0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(M, 2);
    hHat_spin_Asymm_3_2_4 = 0.0352147606136882*(-Sigma_lambda*(83.0*nu - 17.0) + 4.0*I*Sigma_n*(55.0*nu - 13.0) +
      25.0*delta*(S_lambda - 4.0*I*S_n))/pow(M, 2);
    hHat_spin_Asymm_3_1_3 = 0.17817416127495*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/pow(M, 2);
    hHat_spin_Asymm_3_0_4 = 0.038575837490523*(-17.0*S_lambda*delta + Sigma_lambda*(35.0*nu - 9.0))/pow(M, 2);
    hHat_spin_Asymm_4_4_4 = 0.950798536569581*(-3.0*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3.0*nu - 1.0) +
      delta*(S_lambda + I*S_n))/pow(M, 2);
    hHat_spin_Asymm_4_2_4 = 0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) + 14.0*I*Sigma_n*(3.0*nu - 1.0) +
      delta*(13.0*S_lambda - 14.0*I*S_n))/pow(M, 2);
    hHat_spin_Asymm_4_0_4 = 0.00841793787126842*I*(S_n*delta - 3.0*Sigma_n*nu + Sigma_n)/pow(M, 2);

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = rhOverM_coeff*(hHat_2_0_0 + pow(v, 3)*(hHat_spin_Symm_2_0_3 + hHat_spin_Symm_2_0_4*v));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*pow(v, 2));
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 +
      hHat_spin_Symm_2_1_4))));
    Asymm = rhOverM_coeff*pow(v, 3)*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v);
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 +
      hHat_spin_Symm_2_2_4))));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*pow(v, 2));
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = 0;
    Asymm = hHat_spin_Asymm_3_0_4*rhOverM_coeff*pow(v, 4);
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_3_1_1 + pow(v, 2)*(hHat_3_1_3 + v*(hHat_3_1_4 + hHat_spin_Symm_3_1_4)));
    Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3);
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 + v*(hHat_3_2_4*v + hHat_spin_Symm_3_2_3));
    Asymm = hHat_spin_Asymm_3_2_4*rhOverM_coeff*pow(v, 4);
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = rhOverM_coeff*v*(hHat_3_3_1 + pow(v, 2)*(hHat_3_3_3 + v*(hHat_3_3_4 + hHat_spin_Symm_3_3_4)));
    Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3);
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = hHat_spin_Asymm_4_0_4*rhOverM_coeff*pow(v, 4);
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_1_3 + hHat_spin_Symm_4_1_4*v);
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_2_2 + hHat_4_2_4*pow(v, 2));
    Asymm = hHat_spin_Asymm_4_2_4*rhOverM_coeff*pow(v, 4);
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_3_3 + hHat_spin_Symm_4_3_4*v);
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_4_2 + hHat_4_4_4*pow(v, 2));
    Asymm = hHat_spin_Asymm_4_4_4*rhOverM_coeff*pow(v, 4);
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = hHat_5_1_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = hHat_5_2_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = hHat_5_3_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = hHat_5_4_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = hHat_5_5_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = hHat_6_2_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = hHat_6_4_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = hHat_6_6_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = 0;
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_2p0PN : public WaveformModes_Base


class WaveformModes_2p5PN : public WaveformModes_Base {
private:
  const double M1, M2;
  double v;
  const double M, delta, nu;
  Quaternions::Quaternion chiVec1, chiVec2;
  double chi1_n, chi1_lambda, chi1_ell, chi2_n, chi2_lambda, chi2_ell, S_ell, S_n, S_lambda, Sigma_ell, Sigma_n,
         Sigma_lambda, S1_ell, S1_n, S1_lambda, S2_ell, S2_n, S2_lambda;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_1_1, hHat_2_1_3, hHat_2_1_4, hHat_2_1_5, hHat_2_2_0, hHat_2_2_2,
                             hHat_2_2_3, hHat_2_2_4, hHat_2_2_5, hHat_3_0_5, hHat_3_1_1, hHat_3_1_3, hHat_3_1_4,
                             hHat_3_1_5, hHat_3_2_2, hHat_3_2_4, hHat_3_2_5, hHat_3_3_1, hHat_3_3_3, hHat_3_3_4,
                             hHat_3_3_5, hHat_4_0_0, hHat_4_1_3, hHat_4_1_5, hHat_4_2_2, hHat_4_2_4, hHat_4_2_5,
                             hHat_4_3_3, hHat_4_3_5, hHat_4_4_2, hHat_4_4_4, hHat_4_4_5, hHat_5_1_3, hHat_5_1_5,
                             hHat_5_2_4, hHat_5_3_3, hHat_5_3_5, hHat_5_4_4, hHat_5_5_3, hHat_5_5_5, hHat_6_1_5,
                             hHat_6_2_4, hHat_6_3_5, hHat_6_4_4, hHat_6_5_5, hHat_6_6_4, hHat_7_1_5, hHat_7_3_5,
                             hHat_7_5_5, hHat_7_7_5;
  std::complex<double> hHat_spin_Symm_2_2_3, hHat_spin_Symm_2_2_4, hHat_spin_Symm_2_1_2, hHat_spin_Symm_2_1_4,
                       hHat_spin_Symm_2_0_3, hHat_spin_Symm_2_0_4, hHat_spin_Symm_3_3_4, hHat_spin_Symm_3_2_3,
                       hHat_spin_Symm_3_1_4, hHat_spin_Symm_4_3_4, hHat_spin_Symm_4_1_4, hHat_spin_Asymm_2_2_2,
                       hHat_spin_Asymm_2_2_4, hHat_spin_Asymm_2_1_3, hHat_spin_Asymm_2_1_4, hHat_spin_Asymm_2_0_2,
                       hHat_spin_Asymm_2_0_4, hHat_spin_Asymm_3_3_3, hHat_spin_Asymm_3_2_4, hHat_spin_Asymm_3_1_3,
                       hHat_spin_Asymm_3_0_4, hHat_spin_Asymm_4_4_4, hHat_spin_Asymm_4_2_4, hHat_spin_Asymm_4_0_4;

public:
  WaveformModes_2p5PN(const double M1_i, const double M2_i, const double v_i, const Quaternions::Quaternion chiVec1_i,
                      const Quaternions::Quaternion chiVec2_i) :
    M1(M1_i), M2(M2_i), v(v_i), M(M1 + M2), delta((M1 - M2)/M), nu(M1*M2/pow(M, 2)), chiVec1(chiVec1_i),
    chiVec2(chiVec2_i), chi1_n(chiVec1[1]), chi1_lambda(chiVec1[2]), chi1_ell(chiVec1[3]), chi2_n(chiVec2[1]),
    chi2_lambda(chiVec2[2]), chi2_ell(chiVec2[3]), S_ell(pow(M1, 2)*chi1_ell + pow(M2, 2)*chi2_ell), S_n(pow(M1,
    2)*chi1_n + pow(M2, 2)*chi2_n), S_lambda(pow(M1, 2)*chi1_lambda + pow(M2, 2)*chi2_lambda), Sigma_ell(M*(-M1*chi1_ell
    + M2*chi2_ell)), Sigma_n(M*(-M1*chi1_n + M2*chi2_n)), Sigma_lambda(M*(-M1*chi1_lambda + M2*chi2_lambda)),
    S1_ell(pow(M1, 2)*chi1_ell), S1_n(pow(M1, 2)*chi1_n), S1_lambda(pow(M1, 2)*chi1_lambda), S2_ell(pow(M2,
    2)*chi2_ell), S2_n(pow(M2, 2)*chi2_n), S2_lambda(pow(M2, 2)*chi2_lambda), rhOverM_coeff(6.34132367616962*nu*pow(v,
    2)), hHat_2_0_0(-0.145802960879951), hHat_2_1_1(0.333333333333333*I*delta),
    hHat_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)), hHat_2_1_4(delta*(0.628764787039964 + 1.0471975511966*I)),
    hHat_2_1_5(0.000661375661375661*I*delta*(nu*(237.0*nu - 2036.0) - 172.0)), hHat_2_2_0(1.00000000000000),
    hHat_2_2_2(1.30952380952381*nu - 2.54761904761905), hHat_2_2_3(6.28318530717959),
    hHat_2_2_4(0.000661375661375661*nu*(2047.0*nu - 7483.0) - 1.43716931216931), hHat_2_2_5(5.08638810581205*nu -
    24.0*I*nu - 16.0071625682908), hHat_3_0_5(-0.370328039909021*I*nu), hHat_3_1_1(0.0222717701593687*I*delta),
    hHat_3_1_3(-0.0148478467729125*I*delta*(nu + 4.0)), hHat_3_1_4(delta*(0.0620557076072073 + 0.0699688295151131*I)),
    hHat_3_1_5(-0.000112483687673579*I*delta*(nu*(247.0*nu + 272.0) - 607.0)), hHat_3_2_2(-0.845154254728516*nu +
    0.281718084909506), hHat_3_2_4(0.00313020094343895*nu*(-365.0*nu + 725.0) - 0.604128782083717),
    hHat_3_2_5(-5.31026079561053*nu + 3.71867872080547*I*nu + 1.77008693187018 - 0.845154254728517*I),
    hHat_3_3_1(-0.776323754260148*I*delta), hHat_3_3_3(-1.5526475085203*I*delta*(nu - 2.0)),
    hHat_3_3_4(delta*(-1.37192659820446 - 7.31667900957279*I)), hHat_3_3_5(-0.00235249622503075*I*delta*(nu*(887.0*nu -
    3676.0) + 369.0)), hHat_4_0_0(-0.00140298964521140), hHat_4_1_3(0.00376461626210521*I*delta*(-2.0*nu + 1.0)),
    hHat_4_1_5(-2.85198201674637e-5*I*delta*(nu*(332.0*nu - 1011.0) + 404.0)), hHat_4_2_2(-0.10647942749999*nu +
    0.0354931424999967), hHat_4_2_4(0.000107554977272717*nu*(-285.0*nu + 4025.0) - 0.141004575204532),
    hHat_4_2_5(0.00709862849999933*nu*(-94.2477796076938 + 84.0*I) + 0.22300999146161 - 0.149071198499986*I),
    hHat_4_3_3(0.268926437100239*I*delta*(2.0*nu - 1.0)), hHat_4_3_5(0.00203732149318363*I*delta*(nu*(524.0*nu - 1267.0)
    + 468.0)), hHat_4_4_2(2.25374467927604*nu - 0.751248226425348), hHat_4_4_4(0.0113825488852325*nu*(525.0*nu - 1273.0)
    + 4.04991089336574), hHat_4_4_5(28.3213909099228*nu + 0.0187812056606337*I*(-527.578706662453*nu + 114.192902220818)
    - 9.44046363664094), hHat_5_1_3(0.000176960830360287*I*delta*(-2.0*nu + 1.0)),
    hHat_5_1_5(-4.5374571887253e-6*I*delta*(nu*(4.0*nu - 352.0) + 179.0)), hHat_5_2_4(0.00998814611056655*nu*(5.0*nu -
    5.0) + 0.00998814611056655), hHat_5_3_3(0.0464469088412683*I*delta*(2.0*nu - 1.0)),
    hHat_5_3_5(0.00119094638054534*I*delta*(8.0*nu*(11.0*nu - 58.0) + 207.0)), hHat_5_4_4(-0.276799624590764*nu*(5.0*nu
    - 5.0) - 0.276799624590764), hHat_5_5_3(-0.801376894396698*I*delta*(2.0*nu - 1.0)),
    hHat_5_5_5(-0.0205481254973512*I*delta*(16.0*nu*(16.0*nu - 43.0) + 263.0)),
    hHat_6_1_5(2.35829888333555e-5*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), hHat_6_2_4(0.000835250737974468*nu*(5.0*nu - 5.0)
    + 0.000835250737974468), hHat_6_3_5(-0.0163097621781264*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    hHat_6_4_4(-0.058558165806266*nu*(5.0*nu - 5.0) - 0.058558165806266), hHat_6_5_5(0.299357979653564*I*delta*(nu -
    1.0)*(3.0*nu - 1.0)), hHat_6_6_4(0.903141370807658*nu*(5.0*nu - 5.0) + 0.903141370807658),
    hHat_7_1_5(8.17593033339979e-7*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), hHat_7_3_5(-0.0018582230503756*I*delta*(nu -
    1.0)*(3.0*nu - 1.0)), hHat_7_5_5(0.0733861624905401*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    hHat_7_7_5(-1.05422444934392*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), hHat_spin_Symm_2_2_3(0.166666666666667*(3.0*S_ell +
    5.0*Sigma_ell*delta)/pow(M, 2)), hHat_spin_Symm_2_2_4(0.166666666666667*(12.0*S1_ell*S2_ell +
    10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n - 15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(M, 4)*nu)),
    hHat_spin_Symm_2_1_2(0.5*I*Sigma_ell/pow(M, 2)), hHat_spin_Symm_2_1_4(0.0238095238095238*I*(-86.0*S_ell*delta +
    Sigma_ell*(139.0*nu - 79.0))/pow(M, 2)), hHat_spin_Symm_2_0_3(0.408248290463863*(5.0*S_ell +
    3.0*Sigma_ell*delta)/pow(M, 2)), hHat_spin_Symm_2_0_4(0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(M,
    4)*nu)), hHat_spin_Symm_3_3_4(0.388161877130074*I*(7.0*S_ell*delta - 3.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2)),
    hHat_spin_Symm_3_2_3(0.563436169819011*(S_ell + Sigma_ell*delta)/pow(M, 2)),
    hHat_spin_Symm_3_1_4(0.0111358850796843*I*(S_ell*delta - 5.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2)),
    hHat_spin_Symm_4_3_4(0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(M, 2)),
    hHat_spin_Symm_4_1_4(0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(M, 2)),
    hHat_spin_Asymm_2_2_2(0.5*(-Sigma_lambda - I*Sigma_n)/pow(M, 2)),
    hHat_spin_Asymm_2_2_4(0.0119047619047619*(19.0*S_lambda*delta + 182.0*I*S_n*delta - 43.0*Sigma_lambda*nu +
    5.0*Sigma_lambda - 280.0*I*Sigma_n*nu + 98.0*I*Sigma_n)/pow(M, 2)),
    hHat_spin_Asymm_2_1_3(0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
    13.0*Sigma_n*delta)/pow(M, 2)), hHat_spin_Asymm_2_1_4(0.5*(-3.0*S1_ell*S2_n - 3.0*S1_n*S2_ell)/(pow(M, 4)*nu)),
    hHat_spin_Asymm_2_0_2(0.408248290463863*I*Sigma_n/pow(M, 2)),
    hHat_spin_Asymm_2_0_4(0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(M, 2)),
    hHat_spin_Asymm_3_3_3(0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(M, 2)),
    hHat_spin_Asymm_3_2_4(0.0352147606136882*(-Sigma_lambda*(83.0*nu - 17.0) + 4.0*I*Sigma_n*(55.0*nu - 13.0) +
    25.0*delta*(S_lambda - 4.0*I*S_n))/pow(M, 2)), hHat_spin_Asymm_3_1_3(0.17817416127495*(I*S_lambda + S_n +
    delta*(I*Sigma_lambda + Sigma_n))/pow(M, 2)), hHat_spin_Asymm_3_0_4(0.038575837490523*(-17.0*S_lambda*delta +
    Sigma_lambda*(35.0*nu - 9.0))/pow(M, 2)), hHat_spin_Asymm_4_4_4(0.950798536569581*(-3.0*Sigma_lambda*nu +
    Sigma_lambda - I*Sigma_n*(3.0*nu - 1.0) + delta*(S_lambda + I*S_n))/pow(M, 2)),
    hHat_spin_Asymm_4_2_4(0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) + 14.0*I*Sigma_n*(3.0*nu - 1.0) +
    delta*(13.0*S_lambda - 14.0*I*S_n))/pow(M, 2)), hHat_spin_Asymm_4_0_4(0.00841793787126842*I*(S_n*delta -
    3.0*Sigma_n*nu + Sigma_n)/pow(M, 2))
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;
    chiVec1 = Quaternions::Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
    chiVec2 = Quaternions::Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

    chi1_n = chiVec1[1];
    chi1_lambda = chiVec1[2];
    chi1_ell = chiVec1[3];
    chi2_n = chiVec2[1];
    chi2_lambda = chiVec2[2];
    chi2_ell = chiVec2[3];
    S_ell = pow(M1, 2)*chi1_ell + pow(M2, 2)*chi2_ell;
    S_n = pow(M1, 2)*chi1_n + pow(M2, 2)*chi2_n;
    S_lambda = pow(M1, 2)*chi1_lambda + pow(M2, 2)*chi2_lambda;
    Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell);
    Sigma_n = M*(-M1*chi1_n + M2*chi2_n);
    Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda);
    S1_ell = pow(M1, 2)*chi1_ell;
    S1_n = pow(M1, 2)*chi1_n;
    S1_lambda = pow(M1, 2)*chi1_lambda;
    S2_ell = pow(M2, 2)*chi2_ell;
    S2_n = pow(M2, 2)*chi2_n;
    S2_lambda = pow(M2, 2)*chi2_lambda;
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    hHat_spin_Symm_2_2_3 = 0.166666666666667*(3.0*S_ell + 5.0*Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_2_2_4 = 0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
      15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(M, 4)*nu);
    hHat_spin_Symm_2_1_2 = 0.5*I*Sigma_ell/pow(M, 2);
    hHat_spin_Symm_2_1_4 = 0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(M, 2);
    hHat_spin_Symm_2_0_3 = 0.408248290463863*(5.0*S_ell + 3.0*Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_2_0_4 = 0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(M, 4)*nu);
    hHat_spin_Symm_3_3_4 = 0.388161877130074*I*(7.0*S_ell*delta - 3.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2);
    hHat_spin_Symm_3_2_3 = 0.563436169819011*(S_ell + Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_3_1_4 = 0.0111358850796843*I*(S_ell*delta - 5.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2);
    hHat_spin_Symm_4_3_4 = 0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(M, 2);
    hHat_spin_Symm_4_1_4 = 0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(M, 2);
    hHat_spin_Asymm_2_2_2 = 0.5*(-Sigma_lambda - I*Sigma_n)/pow(M, 2);
    hHat_spin_Asymm_2_2_4 = 0.0119047619047619*(19.0*S_lambda*delta + 182.0*I*S_n*delta - 43.0*Sigma_lambda*nu +
      5.0*Sigma_lambda - 280.0*I*Sigma_n*nu + 98.0*I*Sigma_n)/pow(M, 2);
    hHat_spin_Asymm_2_1_3 = 0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
      13.0*Sigma_n*delta)/pow(M, 2);
    hHat_spin_Asymm_2_1_4 = 0.5*(-3.0*S1_ell*S2_n - 3.0*S1_n*S2_ell)/(pow(M, 4)*nu);
    hHat_spin_Asymm_2_0_2 = 0.408248290463863*I*Sigma_n/pow(M, 2);
    hHat_spin_Asymm_2_0_4 = 0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(M, 2);
    hHat_spin_Asymm_3_3_3 = 0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(M, 2);
    hHat_spin_Asymm_3_2_4 = 0.0352147606136882*(-Sigma_lambda*(83.0*nu - 17.0) + 4.0*I*Sigma_n*(55.0*nu - 13.0) +
      25.0*delta*(S_lambda - 4.0*I*S_n))/pow(M, 2);
    hHat_spin_Asymm_3_1_3 = 0.17817416127495*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/pow(M, 2);
    hHat_spin_Asymm_3_0_4 = 0.038575837490523*(-17.0*S_lambda*delta + Sigma_lambda*(35.0*nu - 9.0))/pow(M, 2);
    hHat_spin_Asymm_4_4_4 = 0.950798536569581*(-3.0*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3.0*nu - 1.0) +
      delta*(S_lambda + I*S_n))/pow(M, 2);
    hHat_spin_Asymm_4_2_4 = 0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) + 14.0*I*Sigma_n*(3.0*nu - 1.0) +
      delta*(13.0*S_lambda - 14.0*I*S_n))/pow(M, 2);
    hHat_spin_Asymm_4_0_4 = 0.00841793787126842*I*(S_n*delta - 3.0*Sigma_n*nu + Sigma_n)/pow(M, 2);

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = rhOverM_coeff*(hHat_2_0_0 + pow(v, 3)*(hHat_spin_Symm_2_0_3 + hHat_spin_Symm_2_0_4*v));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*pow(v, 2));
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_2_1_5*v +
      hHat_spin_Symm_2_1_4))));
    Asymm = rhOverM_coeff*pow(v, 3)*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v);
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 +
      hHat_2_2_5*v + hHat_spin_Symm_2_2_4))));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*pow(v, 2));
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = hHat_3_0_5*rhOverM_coeff*pow(v, 5);
    Asymm = hHat_spin_Asymm_3_0_4*rhOverM_coeff*pow(v, 4);
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_3_1_1 + pow(v, 2)*(hHat_3_1_3 + v*(hHat_3_1_4 + hHat_3_1_5*v + hHat_spin_Symm_3_1_4)));
    Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3);
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 + v*(hHat_spin_Symm_3_2_3 + v*(hHat_3_2_4 + hHat_3_2_5*v)));
    Asymm = hHat_spin_Asymm_3_2_4*rhOverM_coeff*pow(v, 4);
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = rhOverM_coeff*v*(hHat_3_3_1 + pow(v, 2)*(hHat_3_3_3 + v*(hHat_3_3_4 + hHat_3_3_5*v + hHat_spin_Symm_3_3_4)));
    Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3);
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = hHat_spin_Asymm_4_0_4*rhOverM_coeff*pow(v, 4);
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_1_3 + v*(hHat_4_1_5*v + hHat_spin_Symm_4_1_4));
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_2_2 + pow(v, 2)*(hHat_4_2_4 + hHat_4_2_5*v));
    Asymm = hHat_spin_Asymm_4_2_4*rhOverM_coeff*pow(v, 4);
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_3_3 + v*(hHat_4_3_5*v + hHat_spin_Symm_4_3_4));
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_4_2 + pow(v, 2)*(hHat_4_4_4 + hHat_4_4_5*v));
    Asymm = hHat_spin_Asymm_4_4_4*rhOverM_coeff*pow(v, 4);
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_1_3 + hHat_5_1_5*pow(v, 2));
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = hHat_5_2_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_3_3 + hHat_5_3_5*pow(v, 2));
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = hHat_5_4_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_5_3 + hHat_5_5_5*pow(v, 2));
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = hHat_6_1_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = hHat_6_2_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = hHat_6_3_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = hHat_6_4_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = hHat_6_5_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = hHat_6_6_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = hHat_7_1_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = hHat_7_3_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = hHat_7_5_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = hHat_7_7_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = 0;
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_2p5PN : public WaveformModes_Base


class WaveformModes_3p0PN : public WaveformModes_Base {
private:
  const double M1, M2;
  double v;
  const double M, delta, nu;
  Quaternions::Quaternion chiVec1, chiVec2;
  double chi1_n, chi1_lambda, chi1_ell, chi2_n, chi2_lambda, chi2_ell, S_ell, S_n, S_lambda, Sigma_ell, Sigma_n,
         Sigma_lambda, S1_ell, S1_n, S1_lambda, S2_ell, S2_n, S2_lambda, logv;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_1_1, hHat_2_1_3, hHat_2_1_4, hHat_2_1_5, hHat_2_1_6, hHat_2_2_0,
                             hHat_2_2_2, hHat_2_2_3, hHat_2_2_4, hHat_2_2_5, hHat_2_2_6, hHat_2_2_lnv_6, hHat_3_0_5,
                             hHat_3_1_1, hHat_3_1_3, hHat_3_1_4, hHat_3_1_5, hHat_3_1_6, hHat_3_2_2, hHat_3_2_4,
                             hHat_3_2_5, hHat_3_2_6, hHat_3_3_1, hHat_3_3_3, hHat_3_3_4, hHat_3_3_5, hHat_3_3_6,
                             hHat_4_0_0, hHat_4_1_3, hHat_4_1_5, hHat_4_1_6, hHat_4_2_2, hHat_4_2_4, hHat_4_2_5,
                             hHat_4_2_6, hHat_4_3_3, hHat_4_3_5, hHat_4_3_6, hHat_4_4_2, hHat_4_4_4, hHat_4_4_5,
                             hHat_4_4_6, hHat_5_1_3, hHat_5_1_5, hHat_5_1_6, hHat_5_2_4, hHat_5_2_6, hHat_5_3_3,
                             hHat_5_3_5, hHat_5_3_6, hHat_5_4_4, hHat_5_4_6, hHat_5_5_3, hHat_5_5_5, hHat_5_5_6,
                             hHat_6_1_5, hHat_6_2_4, hHat_6_2_6, hHat_6_3_5, hHat_6_4_4, hHat_6_4_6, hHat_6_5_5,
                             hHat_6_6_4, hHat_6_6_6, hHat_7_1_5, hHat_7_2_6, hHat_7_3_5, hHat_7_4_6, hHat_7_5_5,
                             hHat_7_6_6, hHat_7_7_5, hHat_8_2_6, hHat_8_4_6, hHat_8_6_6, hHat_8_8_6;
  std::complex<double> hHat_spin_Symm_2_2_3, hHat_spin_Symm_2_2_4, hHat_spin_Symm_2_1_2, hHat_spin_Symm_2_1_4,
                       hHat_spin_Symm_2_0_3, hHat_spin_Symm_2_0_4, hHat_spin_Symm_3_3_4, hHat_spin_Symm_3_2_3,
                       hHat_spin_Symm_3_1_4, hHat_spin_Symm_4_3_4, hHat_spin_Symm_4_1_4, hHat_spin_Asymm_2_2_2,
                       hHat_spin_Asymm_2_2_4, hHat_spin_Asymm_2_1_3, hHat_spin_Asymm_2_1_4, hHat_spin_Asymm_2_0_2,
                       hHat_spin_Asymm_2_0_4, hHat_spin_Asymm_3_3_3, hHat_spin_Asymm_3_2_4, hHat_spin_Asymm_3_1_3,
                       hHat_spin_Asymm_3_0_4, hHat_spin_Asymm_4_4_4, hHat_spin_Asymm_4_2_4, hHat_spin_Asymm_4_0_4;

public:
  WaveformModes_3p0PN(const double M1_i, const double M2_i, const double v_i, const Quaternions::Quaternion chiVec1_i,
                      const Quaternions::Quaternion chiVec2_i) :
    M1(M1_i), M2(M2_i), v(v_i), M(M1 + M2), delta((M1 - M2)/M), nu(M1*M2/pow(M, 2)), chiVec1(chiVec1_i),
    chiVec2(chiVec2_i), chi1_n(chiVec1[1]), chi1_lambda(chiVec1[2]), chi1_ell(chiVec1[3]), chi2_n(chiVec2[1]),
    chi2_lambda(chiVec2[2]), chi2_ell(chiVec2[3]), S_ell(pow(M1, 2)*chi1_ell + pow(M2, 2)*chi2_ell), S_n(pow(M1,
    2)*chi1_n + pow(M2, 2)*chi2_n), S_lambda(pow(M1, 2)*chi1_lambda + pow(M2, 2)*chi2_lambda), Sigma_ell(M*(-M1*chi1_ell
    + M2*chi2_ell)), Sigma_n(M*(-M1*chi1_n + M2*chi2_n)), Sigma_lambda(M*(-M1*chi1_lambda + M2*chi2_lambda)),
    S1_ell(pow(M1, 2)*chi1_ell), S1_n(pow(M1, 2)*chi1_n), S1_lambda(pow(M1, 2)*chi1_lambda), S2_ell(pow(M2,
    2)*chi2_ell), S2_n(pow(M2, 2)*chi2_n), S2_lambda(pow(M2, 2)*chi2_lambda), logv(log(v)),
    rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), hHat_2_0_0(-0.145802960879951), hHat_2_1_1(0.333333333333333*I*delta),
    hHat_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)), hHat_2_1_4(delta*(0.628764787039964 + 1.0471975511966*I)),
    hHat_2_1_5(0.000661375661375661*I*delta*(nu*(237.0*nu - 2036.0) - 172.0)),
    hHat_2_1_6(0.00595238095238095*delta*(nu*(722.635532333439 + 37.6991118430775*I) - 64.1340082780763 -
    106.814150222053*I)), hHat_2_2_0(1.00000000000000), hHat_2_2_2(1.30952380952381*nu - 2.54761904761905),
    hHat_2_2_3(6.28318530717959), hHat_2_2_4(0.000661375661375661*nu*(2047.0*nu - 7483.0) - 1.43716931216931),
    hHat_2_2_5(5.08638810581205*nu - 24.0*I*nu - 16.0071625682908), hHat_2_2_6(1.00208433541767e-5*nu*(nu*(114635.0*nu -
    729396.0) - 834555.0) + 4.21514354629858*nu + 32.3588011610077 + 12.8057300546327*I),
    hHat_2_2_lnv_6(-8.15238095238095), hHat_3_0_5(-0.370328039909021*I*nu), hHat_3_1_1(0.0222717701593687*I*delta),
    hHat_3_1_3(-0.0148478467729125*I*delta*(nu + 4.0)), hHat_3_1_4(delta*(0.0620557076072073 + 0.0699688295151131*I)),
    hHat_3_1_5(-0.000112483687673579*I*delta*(nu*(247.0*nu + 272.0) - 607.0)),
    hHat_3_1_6(0.000742392338645623*delta*(-46.5203026391962*nu - 15.707963267949*I*(7.0*nu + 16.0) -
    222.903548889591)), hHat_3_2_2(-0.845154254728516*nu + 0.281718084909506),
    hHat_3_2_4(0.00313020094343895*nu*(-365.0*nu + 725.0) - 0.604128782083717), hHat_3_2_5(-5.31026079561053*nu +
    3.71867872080547*I*nu + 1.77008693187018 - 0.845154254728517*I), hHat_3_2_6(7.11409305327034e-5*nu*(nu*(-16023.0*nu
    + 100026.0) - 17387.0) - 0.103225490202953), hHat_3_3_1(-0.776323754260148*I*delta),
    hHat_3_3_3(-1.5526475085203*I*delta*(nu - 2.0)), hHat_3_3_4(delta*(-1.37192659820446 - 7.31667900957279*I)),
    hHat_3_3_5(-0.00235249622503075*I*delta*(nu*(887.0*nu - 3676.0) + 369.0)),
    hHat_3_3_6(0.000319474795991831*delta*(-87338.4780856744*nu - 11451.1052223348*I*(3.0*nu - 8.0) +
    17177.2748951318)), hHat_4_0_0(-0.00140298964521140), hHat_4_1_3(0.00376461626210521*I*delta*(-2.0*nu + 1.0)),
    hHat_4_1_5(-2.85198201674637e-5*I*delta*(nu*(332.0*nu - 1011.0) + 404.0)),
    hHat_4_1_6(0.00012548720873684*delta*(-1744.17766166719*nu - 94.2477796076938*I*(2.0*nu - 1.0) + 105.588830833597)),
    hHat_4_2_2(-0.10647942749999*nu + 0.0354931424999967), hHat_4_2_4(0.000107554977272717*nu*(-285.0*nu + 4025.0) -
    0.141004575204532), hHat_4_2_5(0.00709862849999933*nu*(-94.2477796076938 + 84.0*I) + 0.22300999146161 -
    0.149071198499986*I), hHat_4_2_6(1.37890996503484e-7*nu*(115.0*nu*(3363.0*nu + 34822.0) - 5460759.0) +
    0.184032298439331), hHat_4_3_3(0.268926437100239*I*delta*(2.0*nu - 1.0)),
    hHat_4_3_5(0.00203732149318363*I*delta*(nu*(524.0*nu - 1267.0) + 468.0)),
    hHat_4_3_6(0.000332007947037332*delta*(12359.8791491886*nu + 7634.0701482232*I*(2.0*nu - 1.0) - 3213.43957459432)),
    hHat_4_4_2(2.25374467927604*nu - 0.751248226425348), hHat_4_4_4(0.0113825488852325*nu*(525.0*nu - 1273.0) +
    4.04991089336574), hHat_4_4_5(28.3213909099228*nu + 0.0187812056606337*I*(-527.578706662453*nu + 114.192902220818) -
    9.44046363664094), hHat_4_4_6(2.91860227826476e-6*nu*(5.0*nu*(678291.0*nu - 3231338.0) + 9793071.0) -
    4.0101757911199), hHat_5_1_3(0.000176960830360287*I*delta*(-2.0*nu + 1.0)),
    hHat_5_1_5(-4.5374571887253e-6*I*delta*(nu*(4.0*nu - 352.0) + 179.0)),
    hHat_5_1_6(2.52801186228981e-6*delta*(-8958.08121055678*nu - 219.911485751286*I*(2.0*nu - 1.0) + 278.040605278392)),
    hHat_5_2_4(0.00998814611056655*nu*(5.0*nu - 5.0) + 0.00998814611056655),
    hHat_5_2_6(7.68318931582042e-5*nu*(35.0*nu*(33.0*nu - 118.0) + 3079.0) - 0.0429270763059624),
    hHat_5_3_3(0.0464469088412683*I*delta*(2.0*nu - 1.0)), hHat_5_3_5(0.00119094638054534*I*delta*(8.0*nu*(11.0*nu -
    58.0) + 207.0)), hHat_5_3_6(9.10188297888856e-7*delta*(923537.386398884*nu + 480946.419338061*I*(2.0*nu - 1.0) -
    271701.693199442)), hHat_5_4_4(-0.276799624590764*nu*(5.0*nu - 5.0) - 0.276799624590764),
    hHat_5_4_6(-0.00212922788146741*nu*(5.0*nu*(339.0*nu - 1042.0) + 3619.0) + 1.35388475720164),
    hHat_5_5_3(-0.801376894396698*I*delta*(2.0*nu - 1.0)), hHat_5_5_5(-0.0205481254973512*I*delta*(16.0*nu*(16.0*nu -
    43.0) + 263.0)), hHat_5_5_6(1.83171861576388e-5*delta*(-679921.609610114*nu - 687223.392972767*I*(2.0*nu - 1.0) +
    164747.804805057)), hHat_6_1_5(2.35829888333555e-5*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    hHat_6_2_4(0.000835250737974468*nu*(5.0*nu - 5.0) + 0.000835250737974468),
    hHat_6_2_6(0.000417625368987234*nu*(nu*(7.0*nu - 64.0) + 59.0) - 0.00483252212685228),
    hHat_6_3_5(-0.0163097621781264*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), hHat_6_4_4(-0.058558165806266*nu*(5.0*nu - 5.0) -
    0.058558165806266), hHat_6_4_6(-0.029279082903133*nu*(nu*(19.0*nu - 88.0) + 71.0) + 0.388993529998767),
    hHat_6_5_5(0.299357979653564*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), hHat_6_6_4(0.903141370807658*nu*(5.0*nu - 5.0) +
    0.903141370807658), hHat_6_6_6(0.451570685403829*nu*(nu*(39.0*nu - 128.0) + 91.0) - 7.2896410643761),
    hHat_7_1_5(8.17593033339979e-7*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), hHat_7_2_6(-0.00134580482328584*nu*pow(nu - 1.0,
    2) + 0.000192257831897977), hHat_7_3_5(-0.0018582230503756*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    hHat_7_4_6(0.161597034311329*nu*pow(nu - 1.0, 2) - 0.0230852906159042), hHat_7_5_5(0.0733861624905401*I*delta*(nu -
    1.0)*(3.0*nu - 1.0)), hHat_7_6_6(-2.3464301109844*nu*pow(nu - 1.0, 2) + 0.3352043015692),
    hHat_7_7_5(-1.05422444934392*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), hHat_8_2_6(-8.42775671401151e-5*nu*pow(nu - 1.0, 2)
    + 1.20396524485879e-5), hHat_8_4_6(0.0226281108784145*nu*pow(nu - 1.0, 2) - 0.00323258726834493),
    hHat_8_6_6(-0.645290686836342*nu*pow(nu - 1.0, 2) + 0.0921843838337631), hHat_8_8_6(8.82604070589592*nu*pow(nu -
    1.0, 2) - 1.26086295798513), hHat_spin_Symm_2_2_3(0.166666666666667*(3.0*S_ell + 5.0*Sigma_ell*delta)/pow(M, 2)),
    hHat_spin_Symm_2_2_4(0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
    15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(M, 4)*nu)), hHat_spin_Symm_2_1_2(0.5*I*Sigma_ell/pow(M, 2)),
    hHat_spin_Symm_2_1_4(0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(M, 2)),
    hHat_spin_Symm_2_0_3(0.408248290463863*(5.0*S_ell + 3.0*Sigma_ell*delta)/pow(M, 2)),
    hHat_spin_Symm_2_0_4(0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(M, 4)*nu)),
    hHat_spin_Symm_3_3_4(0.388161877130074*I*(7.0*S_ell*delta - 3.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2)),
    hHat_spin_Symm_3_2_3(0.563436169819011*(S_ell + Sigma_ell*delta)/pow(M, 2)),
    hHat_spin_Symm_3_1_4(0.0111358850796843*I*(S_ell*delta - 5.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2)),
    hHat_spin_Symm_4_3_4(0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(M, 2)),
    hHat_spin_Symm_4_1_4(0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(M, 2)),
    hHat_spin_Asymm_2_2_2(0.5*(-Sigma_lambda - I*Sigma_n)/pow(M, 2)),
    hHat_spin_Asymm_2_2_4(0.0119047619047619*(19.0*S_lambda*delta + 182.0*I*S_n*delta - 43.0*Sigma_lambda*nu +
    5.0*Sigma_lambda - 280.0*I*Sigma_n*nu + 98.0*I*Sigma_n)/pow(M, 2)),
    hHat_spin_Asymm_2_1_3(0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
    13.0*Sigma_n*delta)/pow(M, 2)), hHat_spin_Asymm_2_1_4(0.5*(-3.0*S1_ell*S2_n - 3.0*S1_n*S2_ell)/(pow(M, 4)*nu)),
    hHat_spin_Asymm_2_0_2(0.408248290463863*I*Sigma_n/pow(M, 2)),
    hHat_spin_Asymm_2_0_4(0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(M, 2)),
    hHat_spin_Asymm_3_3_3(0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(M, 2)),
    hHat_spin_Asymm_3_2_4(0.0352147606136882*(-Sigma_lambda*(83.0*nu - 17.0) + 4.0*I*Sigma_n*(55.0*nu - 13.0) +
    25.0*delta*(S_lambda - 4.0*I*S_n))/pow(M, 2)), hHat_spin_Asymm_3_1_3(0.17817416127495*(I*S_lambda + S_n +
    delta*(I*Sigma_lambda + Sigma_n))/pow(M, 2)), hHat_spin_Asymm_3_0_4(0.038575837490523*(-17.0*S_lambda*delta +
    Sigma_lambda*(35.0*nu - 9.0))/pow(M, 2)), hHat_spin_Asymm_4_4_4(0.950798536569581*(-3.0*Sigma_lambda*nu +
    Sigma_lambda - I*Sigma_n*(3.0*nu - 1.0) + delta*(S_lambda + I*S_n))/pow(M, 2)),
    hHat_spin_Asymm_4_2_4(0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) + 14.0*I*Sigma_n*(3.0*nu - 1.0) +
    delta*(13.0*S_lambda - 14.0*I*S_n))/pow(M, 2)), hHat_spin_Asymm_4_0_4(0.00841793787126842*I*(S_n*delta -
    3.0*Sigma_n*nu + Sigma_n)/pow(M, 2))
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;
    chiVec1 = Quaternions::Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
    chiVec2 = Quaternions::Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

    chi1_n = chiVec1[1];
    chi1_lambda = chiVec1[2];
    chi1_ell = chiVec1[3];
    chi2_n = chiVec2[1];
    chi2_lambda = chiVec2[2];
    chi2_ell = chiVec2[3];
    S_ell = pow(M1, 2)*chi1_ell + pow(M2, 2)*chi2_ell;
    S_n = pow(M1, 2)*chi1_n + pow(M2, 2)*chi2_n;
    S_lambda = pow(M1, 2)*chi1_lambda + pow(M2, 2)*chi2_lambda;
    Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell);
    Sigma_n = M*(-M1*chi1_n + M2*chi2_n);
    Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda);
    S1_ell = pow(M1, 2)*chi1_ell;
    S1_n = pow(M1, 2)*chi1_n;
    S1_lambda = pow(M1, 2)*chi1_lambda;
    S2_ell = pow(M2, 2)*chi2_ell;
    S2_n = pow(M2, 2)*chi2_n;
    S2_lambda = pow(M2, 2)*chi2_lambda;
    logv = log(v);
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    hHat_spin_Symm_2_2_3 = 0.166666666666667*(3.0*S_ell + 5.0*Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_2_2_4 = 0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
      15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(M, 4)*nu);
    hHat_spin_Symm_2_1_2 = 0.5*I*Sigma_ell/pow(M, 2);
    hHat_spin_Symm_2_1_4 = 0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(M, 2);
    hHat_spin_Symm_2_0_3 = 0.408248290463863*(5.0*S_ell + 3.0*Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_2_0_4 = 0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(M, 4)*nu);
    hHat_spin_Symm_3_3_4 = 0.388161877130074*I*(7.0*S_ell*delta - 3.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2);
    hHat_spin_Symm_3_2_3 = 0.563436169819011*(S_ell + Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_3_1_4 = 0.0111358850796843*I*(S_ell*delta - 5.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2);
    hHat_spin_Symm_4_3_4 = 0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(M, 2);
    hHat_spin_Symm_4_1_4 = 0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(M, 2);
    hHat_spin_Asymm_2_2_2 = 0.5*(-Sigma_lambda - I*Sigma_n)/pow(M, 2);
    hHat_spin_Asymm_2_2_4 = 0.0119047619047619*(19.0*S_lambda*delta + 182.0*I*S_n*delta - 43.0*Sigma_lambda*nu +
      5.0*Sigma_lambda - 280.0*I*Sigma_n*nu + 98.0*I*Sigma_n)/pow(M, 2);
    hHat_spin_Asymm_2_1_3 = 0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
      13.0*Sigma_n*delta)/pow(M, 2);
    hHat_spin_Asymm_2_1_4 = 0.5*(-3.0*S1_ell*S2_n - 3.0*S1_n*S2_ell)/(pow(M, 4)*nu);
    hHat_spin_Asymm_2_0_2 = 0.408248290463863*I*Sigma_n/pow(M, 2);
    hHat_spin_Asymm_2_0_4 = 0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(M, 2);
    hHat_spin_Asymm_3_3_3 = 0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(M, 2);
    hHat_spin_Asymm_3_2_4 = 0.0352147606136882*(-Sigma_lambda*(83.0*nu - 17.0) + 4.0*I*Sigma_n*(55.0*nu - 13.0) +
      25.0*delta*(S_lambda - 4.0*I*S_n))/pow(M, 2);
    hHat_spin_Asymm_3_1_3 = 0.17817416127495*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/pow(M, 2);
    hHat_spin_Asymm_3_0_4 = 0.038575837490523*(-17.0*S_lambda*delta + Sigma_lambda*(35.0*nu - 9.0))/pow(M, 2);
    hHat_spin_Asymm_4_4_4 = 0.950798536569581*(-3.0*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3.0*nu - 1.0) +
      delta*(S_lambda + I*S_n))/pow(M, 2);
    hHat_spin_Asymm_4_2_4 = 0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) + 14.0*I*Sigma_n*(3.0*nu - 1.0) +
      delta*(13.0*S_lambda - 14.0*I*S_n))/pow(M, 2);
    hHat_spin_Asymm_4_0_4 = 0.00841793787126842*I*(S_n*delta - 3.0*Sigma_n*nu + Sigma_n)/pow(M, 2);

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = rhOverM_coeff*(hHat_2_0_0 + pow(v, 3)*(hHat_spin_Symm_2_0_3 + hHat_spin_Symm_2_0_4*v));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*pow(v, 2));
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_spin_Symm_2_1_4
      + v*(hHat_2_1_5 + hHat_2_1_6*v)))));
    Asymm = rhOverM_coeff*pow(v, 3)*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v);
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 +
      hHat_spin_Symm_2_2_4 + v*(hHat_2_2_5 + v*(hHat_2_2_6 + hHat_2_2_lnv_6*logv))))));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*pow(v, 2));
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = hHat_3_0_5*rhOverM_coeff*pow(v, 5);
    Asymm = hHat_spin_Asymm_3_0_4*rhOverM_coeff*pow(v, 4);
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_3_1_1 + pow(v, 2)*(hHat_3_1_3 + v*(hHat_3_1_4 + hHat_spin_Symm_3_1_4 + v*(hHat_3_1_5 +
      hHat_3_1_6*v))));
    Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3);
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 + v*(hHat_spin_Symm_3_2_3 + v*(hHat_3_2_4 + v*(hHat_3_2_5 +
      hHat_3_2_6*v))));
    Asymm = hHat_spin_Asymm_3_2_4*rhOverM_coeff*pow(v, 4);
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = rhOverM_coeff*v*(hHat_3_3_1 + pow(v, 2)*(hHat_3_3_3 + v*(hHat_3_3_4 + hHat_spin_Symm_3_3_4 + v*(hHat_3_3_5 +
      hHat_3_3_6*v))));
    Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3);
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = hHat_spin_Asymm_4_0_4*rhOverM_coeff*pow(v, 4);
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_1_3 + v*(hHat_spin_Symm_4_1_4 + v*(hHat_4_1_5 + hHat_4_1_6*v)));
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_2_2 + pow(v, 2)*(hHat_4_2_4 + v*(hHat_4_2_5 + hHat_4_2_6*v)));
    Asymm = hHat_spin_Asymm_4_2_4*rhOverM_coeff*pow(v, 4);
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_3_3 + v*(hHat_spin_Symm_4_3_4 + v*(hHat_4_3_5 + hHat_4_3_6*v)));
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_4_2 + pow(v, 2)*(hHat_4_4_4 + v*(hHat_4_4_5 + hHat_4_4_6*v)));
    Asymm = hHat_spin_Asymm_4_4_4*rhOverM_coeff*pow(v, 4);
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_1_3 + pow(v, 2)*(hHat_5_1_5 + hHat_5_1_6*v));
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_5_2_4 + hHat_5_2_6*pow(v, 2));
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_3_3 + pow(v, 2)*(hHat_5_3_5 + hHat_5_3_6*v));
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_5_4_4 + hHat_5_4_6*pow(v, 2));
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_5_3 + pow(v, 2)*(hHat_5_5_5 + hHat_5_5_6*v));
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = hHat_6_1_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_6_2_4 + hHat_6_2_6*pow(v, 2));
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = hHat_6_3_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_6_4_4 + hHat_6_4_6*pow(v, 2));
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = hHat_6_5_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_6_6_4 + hHat_6_6_6*pow(v, 2));
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = hHat_7_1_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = hHat_7_2_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = hHat_7_3_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = hHat_7_4_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = hHat_7_5_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = hHat_7_6_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = hHat_7_7_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = hHat_8_2_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = hHat_8_4_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = hHat_8_6_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = hHat_8_8_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_3p0PN : public WaveformModes_Base


class WaveformModes_3p5PN : public WaveformModes_Base {
private:
  const double M1, M2;
  double v;
  const double M, delta, nu;
  Quaternions::Quaternion chiVec1, chiVec2;
  double chi1_n, chi1_lambda, chi1_ell, chi2_n, chi2_lambda, chi2_ell, S_ell, S_n, S_lambda, Sigma_ell, Sigma_n,
         Sigma_lambda, S1_ell, S1_n, S1_lambda, S2_ell, S2_n, S2_lambda, logv;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_1_1, hHat_2_1_3, hHat_2_1_4, hHat_2_1_5, hHat_2_1_6, hHat_2_2_0,
                             hHat_2_2_2, hHat_2_2_3, hHat_2_2_4, hHat_2_2_5, hHat_2_2_6, hHat_2_2_lnv_6, hHat_2_2_7,
                             hHat_3_0_5, hHat_3_1_1, hHat_3_1_3, hHat_3_1_4, hHat_3_1_5, hHat_3_1_6, hHat_3_2_2,
                             hHat_3_2_4, hHat_3_2_5, hHat_3_2_6, hHat_3_3_1, hHat_3_3_3, hHat_3_3_4, hHat_3_3_5,
                             hHat_3_3_6, hHat_4_0_0, hHat_4_1_3, hHat_4_1_5, hHat_4_1_6, hHat_4_2_2, hHat_4_2_4,
                             hHat_4_2_5, hHat_4_2_6, hHat_4_3_3, hHat_4_3_5, hHat_4_3_6, hHat_4_4_2, hHat_4_4_4,
                             hHat_4_4_5, hHat_4_4_6, hHat_5_1_3, hHat_5_1_5, hHat_5_1_6, hHat_5_2_4, hHat_5_2_6,
                             hHat_5_3_3, hHat_5_3_5, hHat_5_3_6, hHat_5_4_4, hHat_5_4_6, hHat_5_5_3, hHat_5_5_5,
                             hHat_5_5_6, hHat_6_1_5, hHat_6_2_4, hHat_6_2_6, hHat_6_3_5, hHat_6_4_4, hHat_6_4_6,
                             hHat_6_5_5, hHat_6_6_4, hHat_6_6_6, hHat_7_1_5, hHat_7_2_6, hHat_7_3_5, hHat_7_4_6,
                             hHat_7_5_5, hHat_7_6_6, hHat_7_7_5, hHat_8_2_6, hHat_8_4_6, hHat_8_6_6, hHat_8_8_6;
  std::complex<double> hHat_spin_Symm_2_2_3, hHat_spin_Symm_2_2_4, hHat_spin_Symm_2_1_2, hHat_spin_Symm_2_1_4,
                       hHat_spin_Symm_2_0_3, hHat_spin_Symm_2_0_4, hHat_spin_Symm_3_3_4, hHat_spin_Symm_3_2_3,
                       hHat_spin_Symm_3_1_4, hHat_spin_Symm_4_3_4, hHat_spin_Symm_4_1_4, hHat_spin_Asymm_2_2_2,
                       hHat_spin_Asymm_2_2_4, hHat_spin_Asymm_2_1_3, hHat_spin_Asymm_2_1_4, hHat_spin_Asymm_2_0_2,
                       hHat_spin_Asymm_2_0_4, hHat_spin_Asymm_3_3_3, hHat_spin_Asymm_3_2_4, hHat_spin_Asymm_3_1_3,
                       hHat_spin_Asymm_3_0_4, hHat_spin_Asymm_4_4_4, hHat_spin_Asymm_4_2_4, hHat_spin_Asymm_4_0_4;

public:
  WaveformModes_3p5PN(const double M1_i, const double M2_i, const double v_i, const Quaternions::Quaternion chiVec1_i,
                      const Quaternions::Quaternion chiVec2_i) :
    M1(M1_i), M2(M2_i), v(v_i), M(M1 + M2), delta((M1 - M2)/M), nu(M1*M2/pow(M, 2)), chiVec1(chiVec1_i),
    chiVec2(chiVec2_i), chi1_n(chiVec1[1]), chi1_lambda(chiVec1[2]), chi1_ell(chiVec1[3]), chi2_n(chiVec2[1]),
    chi2_lambda(chiVec2[2]), chi2_ell(chiVec2[3]), S_ell(pow(M1, 2)*chi1_ell + pow(M2, 2)*chi2_ell), S_n(pow(M1,
    2)*chi1_n + pow(M2, 2)*chi2_n), S_lambda(pow(M1, 2)*chi1_lambda + pow(M2, 2)*chi2_lambda), Sigma_ell(M*(-M1*chi1_ell
    + M2*chi2_ell)), Sigma_n(M*(-M1*chi1_n + M2*chi2_n)), Sigma_lambda(M*(-M1*chi1_lambda + M2*chi2_lambda)),
    S1_ell(pow(M1, 2)*chi1_ell), S1_n(pow(M1, 2)*chi1_n), S1_lambda(pow(M1, 2)*chi1_lambda), S2_ell(pow(M2,
    2)*chi2_ell), S2_n(pow(M2, 2)*chi2_n), S2_lambda(pow(M2, 2)*chi2_lambda), logv(log(v)),
    rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), hHat_2_0_0(-0.145802960879951), hHat_2_1_1(0.333333333333333*I*delta),
    hHat_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)), hHat_2_1_4(delta*(0.628764787039964 + 1.0471975511966*I)),
    hHat_2_1_5(0.000661375661375661*I*delta*(nu*(237.0*nu - 2036.0) - 172.0)),
    hHat_2_1_6(0.00595238095238095*delta*(nu*(722.635532333439 + 37.6991118430775*I) - 64.1340082780763 -
    106.814150222053*I)), hHat_2_2_0(1.00000000000000), hHat_2_2_2(1.30952380952381*nu - 2.54761904761905),
    hHat_2_2_3(6.28318530717959), hHat_2_2_4(0.000661375661375661*nu*(2047.0*nu - 7483.0) - 1.43716931216931),
    hHat_2_2_5(5.08638810581205*nu - 24.0*I*nu - 16.0071625682908), hHat_2_2_6(1.00208433541767e-5*nu*(nu*(114635.0*nu -
    729396.0) - 834555.0) + 4.21514354629858*nu + 32.3588011610077 + 12.8057300546327*I),
    hHat_2_2_lnv_6(-8.15238095238095), hHat_2_2_7(0.00831109167616348*nu*(560.0*nu - 2459.0) -
    0.00017636684303351*I*nu*(24396.0*nu - 501655.0) - 9.03000110615162), hHat_3_0_5(-0.370328039909021*I*nu),
    hHat_3_1_1(0.0222717701593687*I*delta), hHat_3_1_3(-0.0148478467729125*I*delta*(nu + 4.0)),
    hHat_3_1_4(delta*(0.0620557076072073 + 0.0699688295151131*I)),
    hHat_3_1_5(-0.000112483687673579*I*delta*(nu*(247.0*nu + 272.0) - 607.0)),
    hHat_3_1_6(0.000742392338645623*delta*(-46.5203026391962*nu - 15.707963267949*I*(7.0*nu + 16.0) -
    222.903548889591)), hHat_3_2_2(-0.845154254728516*nu + 0.281718084909506),
    hHat_3_2_4(0.00313020094343895*nu*(-365.0*nu + 725.0) - 0.604128782083717), hHat_3_2_5(-5.31026079561053*nu +
    3.71867872080547*I*nu + 1.77008693187018 - 0.845154254728517*I), hHat_3_2_6(7.11409305327034e-5*nu*(nu*(-16023.0*nu
    + 100026.0) - 17387.0) - 0.103225490202953), hHat_3_3_1(-0.776323754260148*I*delta),
    hHat_3_3_3(-1.5526475085203*I*delta*(nu - 2.0)), hHat_3_3_4(delta*(-1.37192659820446 - 7.31667900957279*I)),
    hHat_3_3_5(-0.00235249622503075*I*delta*(nu*(887.0*nu - 3676.0) + 369.0)),
    hHat_3_3_6(0.000319474795991831*delta*(-87338.4780856744*nu - 11451.1052223348*I*(3.0*nu - 8.0) +
    17177.2748951318)), hHat_4_0_0(-0.00140298964521140), hHat_4_1_3(0.00376461626210521*I*delta*(-2.0*nu + 1.0)),
    hHat_4_1_5(-2.85198201674637e-5*I*delta*(nu*(332.0*nu - 1011.0) + 404.0)),
    hHat_4_1_6(0.00012548720873684*delta*(-1744.17766166719*nu - 94.2477796076938*I*(2.0*nu - 1.0) + 105.588830833597)),
    hHat_4_2_2(-0.10647942749999*nu + 0.0354931424999967), hHat_4_2_4(0.000107554977272717*nu*(-285.0*nu + 4025.0) -
    0.141004575204532), hHat_4_2_5(0.00709862849999933*nu*(-94.2477796076938 + 84.0*I) + 0.22300999146161 -
    0.149071198499986*I), hHat_4_2_6(1.37890996503484e-7*nu*(115.0*nu*(3363.0*nu + 34822.0) - 5460759.0) +
    0.184032298439331), hHat_4_3_3(0.268926437100239*I*delta*(2.0*nu - 1.0)),
    hHat_4_3_5(0.00203732149318363*I*delta*(nu*(524.0*nu - 1267.0) + 468.0)),
    hHat_4_3_6(0.000332007947037332*delta*(12359.8791491886*nu + 7634.0701482232*I*(2.0*nu - 1.0) - 3213.43957459432)),
    hHat_4_4_2(2.25374467927604*nu - 0.751248226425348), hHat_4_4_4(0.0113825488852325*nu*(525.0*nu - 1273.0) +
    4.04991089336574), hHat_4_4_5(28.3213909099228*nu + 0.0187812056606337*I*(-527.578706662453*nu + 114.192902220818) -
    9.44046363664094), hHat_4_4_6(2.91860227826476e-6*nu*(5.0*nu*(678291.0*nu - 3231338.0) + 9793071.0) -
    4.0101757911199), hHat_5_1_3(0.000176960830360287*I*delta*(-2.0*nu + 1.0)),
    hHat_5_1_5(-4.5374571887253e-6*I*delta*(nu*(4.0*nu - 352.0) + 179.0)),
    hHat_5_1_6(2.52801186228981e-6*delta*(-8958.08121055678*nu - 219.911485751286*I*(2.0*nu - 1.0) + 278.040605278392)),
    hHat_5_2_4(0.00998814611056655*nu*(5.0*nu - 5.0) + 0.00998814611056655),
    hHat_5_2_6(7.68318931582042e-5*nu*(35.0*nu*(33.0*nu - 118.0) + 3079.0) - 0.0429270763059624),
    hHat_5_3_3(0.0464469088412683*I*delta*(2.0*nu - 1.0)), hHat_5_3_5(0.00119094638054534*I*delta*(8.0*nu*(11.0*nu -
    58.0) + 207.0)), hHat_5_3_6(9.10188297888856e-7*delta*(923537.386398884*nu + 480946.419338061*I*(2.0*nu - 1.0) -
    271701.693199442)), hHat_5_4_4(-0.276799624590764*nu*(5.0*nu - 5.0) - 0.276799624590764),
    hHat_5_4_6(-0.00212922788146741*nu*(5.0*nu*(339.0*nu - 1042.0) + 3619.0) + 1.35388475720164),
    hHat_5_5_3(-0.801376894396698*I*delta*(2.0*nu - 1.0)), hHat_5_5_5(-0.0205481254973512*I*delta*(16.0*nu*(16.0*nu -
    43.0) + 263.0)), hHat_5_5_6(1.83171861576388e-5*delta*(-679921.609610114*nu - 687223.392972767*I*(2.0*nu - 1.0) +
    164747.804805057)), hHat_6_1_5(2.35829888333555e-5*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    hHat_6_2_4(0.000835250737974468*nu*(5.0*nu - 5.0) + 0.000835250737974468),
    hHat_6_2_6(0.000417625368987234*nu*(nu*(7.0*nu - 64.0) + 59.0) - 0.00483252212685228),
    hHat_6_3_5(-0.0163097621781264*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), hHat_6_4_4(-0.058558165806266*nu*(5.0*nu - 5.0) -
    0.058558165806266), hHat_6_4_6(-0.029279082903133*nu*(nu*(19.0*nu - 88.0) + 71.0) + 0.388993529998767),
    hHat_6_5_5(0.299357979653564*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), hHat_6_6_4(0.903141370807658*nu*(5.0*nu - 5.0) +
    0.903141370807658), hHat_6_6_6(0.451570685403829*nu*(nu*(39.0*nu - 128.0) + 91.0) - 7.2896410643761),
    hHat_7_1_5(8.17593033339979e-7*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), hHat_7_2_6(-0.00134580482328584*nu*pow(nu - 1.0,
    2) + 0.000192257831897977), hHat_7_3_5(-0.0018582230503756*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    hHat_7_4_6(0.161597034311329*nu*pow(nu - 1.0, 2) - 0.0230852906159042), hHat_7_5_5(0.0733861624905401*I*delta*(nu -
    1.0)*(3.0*nu - 1.0)), hHat_7_6_6(-2.3464301109844*nu*pow(nu - 1.0, 2) + 0.3352043015692),
    hHat_7_7_5(-1.05422444934392*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), hHat_8_2_6(-8.42775671401151e-5*nu*pow(nu - 1.0, 2)
    + 1.20396524485879e-5), hHat_8_4_6(0.0226281108784145*nu*pow(nu - 1.0, 2) - 0.00323258726834493),
    hHat_8_6_6(-0.645290686836342*nu*pow(nu - 1.0, 2) + 0.0921843838337631), hHat_8_8_6(8.82604070589592*nu*pow(nu -
    1.0, 2) - 1.26086295798513), hHat_spin_Symm_2_2_3(0.166666666666667*(3.0*S_ell + 5.0*Sigma_ell*delta)/pow(M, 2)),
    hHat_spin_Symm_2_2_4(0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
    15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(M, 4)*nu)), hHat_spin_Symm_2_1_2(0.5*I*Sigma_ell/pow(M, 2)),
    hHat_spin_Symm_2_1_4(0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(M, 2)),
    hHat_spin_Symm_2_0_3(0.408248290463863*(5.0*S_ell + 3.0*Sigma_ell*delta)/pow(M, 2)),
    hHat_spin_Symm_2_0_4(0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(M, 4)*nu)),
    hHat_spin_Symm_3_3_4(0.388161877130074*I*(7.0*S_ell*delta - 3.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2)),
    hHat_spin_Symm_3_2_3(0.563436169819011*(S_ell + Sigma_ell*delta)/pow(M, 2)),
    hHat_spin_Symm_3_1_4(0.0111358850796843*I*(S_ell*delta - 5.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2)),
    hHat_spin_Symm_4_3_4(0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(M, 2)),
    hHat_spin_Symm_4_1_4(0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(M, 2)),
    hHat_spin_Asymm_2_2_2(0.5*(-Sigma_lambda - I*Sigma_n)/pow(M, 2)),
    hHat_spin_Asymm_2_2_4(0.0119047619047619*(19.0*S_lambda*delta + 182.0*I*S_n*delta - 43.0*Sigma_lambda*nu +
    5.0*Sigma_lambda - 280.0*I*Sigma_n*nu + 98.0*I*Sigma_n)/pow(M, 2)),
    hHat_spin_Asymm_2_1_3(0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
    13.0*Sigma_n*delta)/pow(M, 2)), hHat_spin_Asymm_2_1_4(0.5*(-3.0*S1_ell*S2_n - 3.0*S1_n*S2_ell)/(pow(M, 4)*nu)),
    hHat_spin_Asymm_2_0_2(0.408248290463863*I*Sigma_n/pow(M, 2)),
    hHat_spin_Asymm_2_0_4(0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(M, 2)),
    hHat_spin_Asymm_3_3_3(0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(M, 2)),
    hHat_spin_Asymm_3_2_4(0.0352147606136882*(-Sigma_lambda*(83.0*nu - 17.0) + 4.0*I*Sigma_n*(55.0*nu - 13.0) +
    25.0*delta*(S_lambda - 4.0*I*S_n))/pow(M, 2)), hHat_spin_Asymm_3_1_3(0.17817416127495*(I*S_lambda + S_n +
    delta*(I*Sigma_lambda + Sigma_n))/pow(M, 2)), hHat_spin_Asymm_3_0_4(0.038575837490523*(-17.0*S_lambda*delta +
    Sigma_lambda*(35.0*nu - 9.0))/pow(M, 2)), hHat_spin_Asymm_4_4_4(0.950798536569581*(-3.0*Sigma_lambda*nu +
    Sigma_lambda - I*Sigma_n*(3.0*nu - 1.0) + delta*(S_lambda + I*S_n))/pow(M, 2)),
    hHat_spin_Asymm_4_2_4(0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) + 14.0*I*Sigma_n*(3.0*nu - 1.0) +
    delta*(13.0*S_lambda - 14.0*I*S_n))/pow(M, 2)), hHat_spin_Asymm_4_0_4(0.00841793787126842*I*(S_n*delta -
    3.0*Sigma_n*nu + Sigma_n)/pow(M, 2))
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;
    chiVec1 = Quaternions::Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
    chiVec2 = Quaternions::Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

    chi1_n = chiVec1[1];
    chi1_lambda = chiVec1[2];
    chi1_ell = chiVec1[3];
    chi2_n = chiVec2[1];
    chi2_lambda = chiVec2[2];
    chi2_ell = chiVec2[3];
    S_ell = pow(M1, 2)*chi1_ell + pow(M2, 2)*chi2_ell;
    S_n = pow(M1, 2)*chi1_n + pow(M2, 2)*chi2_n;
    S_lambda = pow(M1, 2)*chi1_lambda + pow(M2, 2)*chi2_lambda;
    Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell);
    Sigma_n = M*(-M1*chi1_n + M2*chi2_n);
    Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda);
    S1_ell = pow(M1, 2)*chi1_ell;
    S1_n = pow(M1, 2)*chi1_n;
    S1_lambda = pow(M1, 2)*chi1_lambda;
    S2_ell = pow(M2, 2)*chi2_ell;
    S2_n = pow(M2, 2)*chi2_n;
    S2_lambda = pow(M2, 2)*chi2_lambda;
    logv = log(v);
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    hHat_spin_Symm_2_2_3 = 0.166666666666667*(3.0*S_ell + 5.0*Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_2_2_4 = 0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
      15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(M, 4)*nu);
    hHat_spin_Symm_2_1_2 = 0.5*I*Sigma_ell/pow(M, 2);
    hHat_spin_Symm_2_1_4 = 0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(M, 2);
    hHat_spin_Symm_2_0_3 = 0.408248290463863*(5.0*S_ell + 3.0*Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_2_0_4 = 0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(M, 4)*nu);
    hHat_spin_Symm_3_3_4 = 0.388161877130074*I*(7.0*S_ell*delta - 3.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2);
    hHat_spin_Symm_3_2_3 = 0.563436169819011*(S_ell + Sigma_ell*delta)/pow(M, 2);
    hHat_spin_Symm_3_1_4 = 0.0111358850796843*I*(S_ell*delta - 5.0*Sigma_ell*(3.0*nu - 1.0))/pow(M, 2);
    hHat_spin_Symm_4_3_4 = 0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(M, 2);
    hHat_spin_Symm_4_1_4 = 0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(M, 2);
    hHat_spin_Asymm_2_2_2 = 0.5*(-Sigma_lambda - I*Sigma_n)/pow(M, 2);
    hHat_spin_Asymm_2_2_4 = 0.0119047619047619*(19.0*S_lambda*delta + 182.0*I*S_n*delta - 43.0*Sigma_lambda*nu +
      5.0*Sigma_lambda - 280.0*I*Sigma_n*nu + 98.0*I*Sigma_n)/pow(M, 2);
    hHat_spin_Asymm_2_1_3 = 0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
      13.0*Sigma_n*delta)/pow(M, 2);
    hHat_spin_Asymm_2_1_4 = 0.5*(-3.0*S1_ell*S2_n - 3.0*S1_n*S2_ell)/(pow(M, 4)*nu);
    hHat_spin_Asymm_2_0_2 = 0.408248290463863*I*Sigma_n/pow(M, 2);
    hHat_spin_Asymm_2_0_4 = 0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(M, 2);
    hHat_spin_Asymm_3_3_3 = 0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(M, 2);
    hHat_spin_Asymm_3_2_4 = 0.0352147606136882*(-Sigma_lambda*(83.0*nu - 17.0) + 4.0*I*Sigma_n*(55.0*nu - 13.0) +
      25.0*delta*(S_lambda - 4.0*I*S_n))/pow(M, 2);
    hHat_spin_Asymm_3_1_3 = 0.17817416127495*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/pow(M, 2);
    hHat_spin_Asymm_3_0_4 = 0.038575837490523*(-17.0*S_lambda*delta + Sigma_lambda*(35.0*nu - 9.0))/pow(M, 2);
    hHat_spin_Asymm_4_4_4 = 0.950798536569581*(-3.0*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3.0*nu - 1.0) +
      delta*(S_lambda + I*S_n))/pow(M, 2);
    hHat_spin_Asymm_4_2_4 = 0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) + 14.0*I*Sigma_n*(3.0*nu - 1.0) +
      delta*(13.0*S_lambda - 14.0*I*S_n))/pow(M, 2);
    hHat_spin_Asymm_4_0_4 = 0.00841793787126842*I*(S_n*delta - 3.0*Sigma_n*nu + Sigma_n)/pow(M, 2);

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = rhOverM_coeff*(hHat_2_0_0 + pow(v, 3)*(hHat_spin_Symm_2_0_3 + hHat_spin_Symm_2_0_4*v));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*pow(v, 2));
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_spin_Symm_2_1_4
      + v*(hHat_2_1_5 + hHat_2_1_6*v)))));
    Asymm = rhOverM_coeff*pow(v, 3)*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v);
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 +
      hHat_spin_Symm_2_2_4 + v*(hHat_2_2_5 + v*(hHat_2_2_6 + hHat_2_2_7*v + hHat_2_2_lnv_6*logv))))));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*pow(v, 2));
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = hHat_3_0_5*rhOverM_coeff*pow(v, 5);
    Asymm = hHat_spin_Asymm_3_0_4*rhOverM_coeff*pow(v, 4);
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_3_1_1 + pow(v, 2)*(hHat_3_1_3 + v*(hHat_3_1_4 + hHat_spin_Symm_3_1_4 + v*(hHat_3_1_5 +
      hHat_3_1_6*v))));
    Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3);
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 + v*(hHat_spin_Symm_3_2_3 + v*(hHat_3_2_4 + v*(hHat_3_2_5 +
      hHat_3_2_6*v))));
    Asymm = hHat_spin_Asymm_3_2_4*rhOverM_coeff*pow(v, 4);
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = rhOverM_coeff*v*(hHat_3_3_1 + pow(v, 2)*(hHat_3_3_3 + v*(hHat_3_3_4 + hHat_spin_Symm_3_3_4 + v*(hHat_3_3_5 +
      hHat_3_3_6*v))));
    Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3);
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = hHat_spin_Asymm_4_0_4*rhOverM_coeff*pow(v, 4);
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_1_3 + v*(hHat_spin_Symm_4_1_4 + v*(hHat_4_1_5 + hHat_4_1_6*v)));
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_2_2 + pow(v, 2)*(hHat_4_2_4 + v*(hHat_4_2_5 + hHat_4_2_6*v)));
    Asymm = hHat_spin_Asymm_4_2_4*rhOverM_coeff*pow(v, 4);
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_3_3 + v*(hHat_spin_Symm_4_3_4 + v*(hHat_4_3_5 + hHat_4_3_6*v)));
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_4_2 + pow(v, 2)*(hHat_4_4_4 + v*(hHat_4_4_5 + hHat_4_4_6*v)));
    Asymm = hHat_spin_Asymm_4_4_4*rhOverM_coeff*pow(v, 4);
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_1_3 + pow(v, 2)*(hHat_5_1_5 + hHat_5_1_6*v));
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_5_2_4 + hHat_5_2_6*pow(v, 2));
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_3_3 + pow(v, 2)*(hHat_5_3_5 + hHat_5_3_6*v));
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_5_4_4 + hHat_5_4_6*pow(v, 2));
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_5_3 + pow(v, 2)*(hHat_5_5_5 + hHat_5_5_6*v));
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = hHat_6_1_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_6_2_4 + hHat_6_2_6*pow(v, 2));
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = hHat_6_3_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_6_4_4 + hHat_6_4_6*pow(v, 2));
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = hHat_6_5_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_6_6_4 + hHat_6_6_6*pow(v, 2));
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = hHat_7_1_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = hHat_7_2_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = hHat_7_3_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = hHat_7_4_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = hHat_7_5_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = hHat_7_6_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = hHat_7_7_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = hHat_8_2_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = hHat_8_4_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = hHat_8_6_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = hHat_8_8_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_3p5PN : public WaveformModes_Base
