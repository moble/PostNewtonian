// File produced automatically by WaveformModeCodeGen.ipynb

class WaveformModes_Base {
public:
  virtual std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_x_k, const double chi1_y_k, const double chi1_z_k,
    const double chi2_x_k, const double chi2_y_k, const double chi2_z_k,
    const double ellHat_x_k, const double ellHat_y_k, const double ellHat_z_k) = 0;
  virtual std::vector<std::complex<double> > operator()(
    const double v_k, const std::vector<double>& chi1, const std::vector<double>& chi2)
  {
    return this->operator()(v_k, chi1[0], chi1[1], chi1[2], chi2[0], chi2[1], chi2[2], 0.0, 0.0, 1.0);
  }
};

const unsigned int ellMax = 8;
const std::complex<double> I(0,1.0);
inline std::complex<double> conjugate(const std::complex<double>& a) { return std::conj(a); }

class WaveformModes_0PN : public WaveformModes_Base {
private:
  const double m1, m2;
  double v;
  const double m, nu;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_2_0, hHat_4_0_0;

public:
  WaveformModes_0PN(const double m1_i, const double m2_i, const double v_i) :
    m1(m1_i), m2(m2_i), v(v_i), m(m1 + m2), nu(m1*m2/pow(m, 2)), rhOverM_coeff(6.34132367616962*nu*pow(v, 2)),
    hHat_2_0_0(-0.145802960879951), hHat_2_2_0(1.00000000000000), hHat_4_0_0(-0.00140298964521140)
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_x_k, const double chi1_y_k, const double chi1_z_k,
    const double chi2_x_k, const double chi2_y_k, const double chi2_z_k,
    const double ellHat_x_k, const double ellHat_y_k, const double ellHat_z_k)
  {
    v = v_k;

    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = conjugate(hHat_2_2_0)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = hHat_2_0_0*rhOverM_coeff;
    Modes[i++] = 0;
    Modes[i++] = hHat_2_2_0*rhOverM_coeff;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = hHat_4_0_0*rhOverM_coeff;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;

    return Modes;
  }

}; // class WaveformModes_0PN : public WaveformModes_Base


class WaveformModes_0p50PN : public WaveformModes_Base {
private:
  const double m1, m2;
  double v;
  const double m, delta, nu;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_1_1, hHat_2_2_0, hHat_3_1_1, hHat_3_3_1, hHat_4_0_0;

public:
  WaveformModes_0p50PN(const double m1_i, const double m2_i, const double v_i) :
    m1(m1_i), m2(m2_i), v(v_i), m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)),
    rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), hHat_2_0_0(-0.145802960879951), hHat_2_1_1(0.333333333333333*I*delta),
    hHat_2_2_0(1.00000000000000), hHat_3_1_1(0.0222717701593687*I*delta), hHat_3_3_1(-0.776323754260148*I*delta),
    hHat_4_0_0(-0.00140298964521140)
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_x_k, const double chi1_y_k, const double chi1_z_k,
    const double chi2_x_k, const double chi2_y_k, const double chi2_z_k,
    const double ellHat_x_k, const double ellHat_y_k, const double ellHat_z_k)
  {
    v = v_k;

    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = conjugate(hHat_2_2_0)*conjugate(rhOverM_coeff);
    Modes[i++] = conjugate(hHat_2_1_1)*conjugate(rhOverM_coeff)*conjugate(v);
    Modes[i++] = hHat_2_0_0*rhOverM_coeff;
    Modes[i++] = hHat_2_1_1*rhOverM_coeff*v;
    Modes[i++] = hHat_2_2_0*rhOverM_coeff;
    Modes[i++] = -conjugate(hHat_3_3_1)*conjugate(rhOverM_coeff)*conjugate(v);
    Modes[i++] = 0;
    Modes[i++] = -conjugate(hHat_3_1_1)*conjugate(rhOverM_coeff)*conjugate(v);
    Modes[i++] = 0;
    Modes[i++] = hHat_3_1_1*rhOverM_coeff*v;
    Modes[i++] = 0;
    Modes[i++] = hHat_3_3_1*rhOverM_coeff*v;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = hHat_4_0_0*rhOverM_coeff;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;

    return Modes;
  }

}; // class WaveformModes_0p50PN : public WaveformModes_Base


class WaveformModes_1p0PN : public WaveformModes_Base {
private:
  const Quaternions::Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const Quaternions::Quaternion S_chi1, S_chi2;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z;
  const double m, delta, nu;
  Quaternions::Quaternion R, nHat, lambdaHat, ellHat;
  double nHat_x, nHat_y, nHat_z, lambdaHat_x, lambdaHat_y, lambdaHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternions::Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, chi1_ell, chi1_n, chi1_lambda, chi2_ell, chi2_n, chi2_lambda,
         Sigma_ell, Sigma_n, Sigma_lambda;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_1_1, hHat_2_2_0, hHat_2_2_2, hHat_3_1_1, hHat_3_2_2, hHat_3_3_1,
                             hHat_4_0_0, hHat_4_2_2, hHat_4_4_2;
  std::complex<double> hHat_spin_Symm_2_1_2, hHat_spin_Asymm_2_2_2, hHat_spin_Asymm_2_0_2;

public:
  WaveformModes_1p0PN(const Quaternions::Quaternion xHat_i, const Quaternions::Quaternion yHat_i, const
                      Quaternions::Quaternion zHat_i, const double m1_i, const double m2_i, const double v_i, const
                      Quaternions::Quaternion S_chi1_i, const Quaternions::Quaternion S_chi2_i, const double
                      rfrak_chi1_x_i, const double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double
                      rfrak_chi2_y_i, const double rfrak_frame_x_i, const double rfrak_frame_y_i, const double
                      rfrak_frame_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), S_chi1(S_chi1_i), S_chi2(S_chi2_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_frame_x(rfrak_frame_x_i), rfrak_frame_y(rfrak_frame_y_i),
    rfrak_frame_z(rfrak_frame_z_i), m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_frame_x*xHat +
    rfrak_frame_y*yHat + rfrak_frame_z*zHat)), nHat(R*xHat*conjugate(R)), lambdaHat(R*yHat*conjugate(R)),
    ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]), lambdaHat_x(lambdaHat[1]),
    lambdaHat_y(lambdaHat[2]), lambdaHat_z(lambdaHat[3]), ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]),
    R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)), R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)),
    chiVec1(S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)),
    chiVec2(S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]),
    chi1_z(chiVec1[3]), chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1_ell(chi1_x*ellHat_x +
    chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z),
    chi1_lambda(chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z), chi2_ell(chi2_x*ellHat_x +
    chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z),
    chi2_lambda(chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z), Sigma_ell(m*(-chi1_ell*m1 +
    chi2_ell*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), Sigma_lambda(m*(-chi1_lambda*m1 + chi2_lambda*m2)),
    rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), hHat_2_0_0(-0.145802960879951), hHat_2_1_1(0.333333333333333*I*delta),
    hHat_2_2_0(1.00000000000000), hHat_2_2_2(1.30952380952381*nu - 2.54761904761905),
    hHat_3_1_1(0.0222717701593687*I*delta), hHat_3_2_2(-0.845154254728516*nu + 0.281718084909506),
    hHat_3_3_1(-0.776323754260148*I*delta), hHat_4_0_0(-0.00140298964521140), hHat_4_2_2(-0.10647942749999*nu +
    0.0354931424999967), hHat_4_4_2(2.25374467927604*nu - 0.751248226425348),
    hHat_spin_Symm_2_1_2(0.5*I*Sigma_ell/pow(m, 2)), hHat_spin_Asymm_2_2_2(0.5*(-Sigma_lambda - I*Sigma_n)/pow(m, 2)),
    hHat_spin_Asymm_2_0_2(0.408248290463863*I*Sigma_n/pow(m, 2))
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_x_k, const double chi1_y_k, const double chi1_z_k,
    const double chi2_x_k, const double chi2_y_k, const double chi2_z_k,
    const double ellHat_x_k, const double ellHat_y_k, const double ellHat_z_k)
  {
    v = v_k;
    chi1_x = chi1_x_k;
    chi1_y = chi1_y_k;
    chi1_z = chi1_z_k;
    chi2_x = chi2_x_k;
    chi2_y = chi2_y_k;
    chi2_z = chi2_z_k;
    ellHat_x = ellHat_x_k;
    ellHat_y = ellHat_y_k;
    ellHat_z = ellHat_z_k;

    R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat);
    nHat = R*xHat*conjugate(R);
    lambdaHat = R*yHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    lambdaHat_x = lambdaHat[1];
    lambdaHat_y = lambdaHat[2];
    lambdaHat_z = lambdaHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1);
    chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1_ell = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi1_lambda = chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z;
    chi2_ell = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    chi2_lambda = chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z;
    Sigma_ell = m*(-chi1_ell*m1 + chi2_ell*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    Sigma_lambda = m*(-chi1_lambda*m1 + chi2_lambda*m2);
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    hHat_spin_Symm_2_1_2 = 0.5*I*Sigma_ell/pow(m, 2);
    hHat_spin_Asymm_2_2_2 = 0.5*(-Sigma_lambda - I*Sigma_n)/pow(m, 2);
    hHat_spin_Asymm_2_0_2 = 0.408248290463863*I*Sigma_n/pow(m, 2);

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = (conjugate(hHat_2_2_0) + conjugate(hHat_2_2_2)*pow(conjugate(v), 2))*conjugate(rhOverM_coeff) -
      conjugate(hHat_spin_Asymm_2_2_2)*conjugate(rhOverM_coeff)*pow(conjugate(v), 2);
    Modes[i++] = (conjugate(hHat_2_1_1) +
      conjugate(hHat_spin_Symm_2_1_2)*conjugate(v))*conjugate(rhOverM_coeff)*conjugate(v);
    Modes[i++] = hHat_2_0_0*rhOverM_coeff + hHat_spin_Asymm_2_0_2*rhOverM_coeff*pow(v, 2);
    Modes[i++] = rhOverM_coeff*v*(hHat_2_1_1 + hHat_spin_Symm_2_1_2*v);
    Modes[i++] = hHat_spin_Asymm_2_2_2*rhOverM_coeff*pow(v, 2) + rhOverM_coeff*(hHat_2_2_0 + hHat_2_2_2*pow(v, 2));
    Modes[i++] = -conjugate(hHat_3_3_1)*conjugate(rhOverM_coeff)*conjugate(v);
    Modes[i++] = -conjugate(hHat_3_2_2)*conjugate(rhOverM_coeff)*pow(conjugate(v), 2);
    Modes[i++] = -conjugate(hHat_3_1_1)*conjugate(rhOverM_coeff)*conjugate(v);
    Modes[i++] = 0;
    Modes[i++] = hHat_3_1_1*rhOverM_coeff*v;
    Modes[i++] = hHat_3_2_2*rhOverM_coeff*pow(v, 2);
    Modes[i++] = hHat_3_3_1*rhOverM_coeff*v;
    Modes[i++] = conjugate(hHat_4_4_2)*conjugate(rhOverM_coeff)*pow(conjugate(v), 2);
    Modes[i++] = 0;
    Modes[i++] = conjugate(hHat_4_2_2)*conjugate(rhOverM_coeff)*pow(conjugate(v), 2);
    Modes[i++] = 0;
    Modes[i++] = hHat_4_0_0*rhOverM_coeff;
    Modes[i++] = 0;
    Modes[i++] = hHat_4_2_2*rhOverM_coeff*pow(v, 2);
    Modes[i++] = 0;
    Modes[i++] = hHat_4_4_2*rhOverM_coeff*pow(v, 2);
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;

    return Modes;
  }

}; // class WaveformModes_1p0PN : public WaveformModes_Base


class WaveformModes_1p5PN : public WaveformModes_Base {
private:
  const Quaternions::Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const Quaternions::Quaternion S_chi1, S_chi2;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z;
  const double m, delta, nu;
  Quaternions::Quaternion R, nHat, lambdaHat, ellHat;
  double nHat_x, nHat_y, nHat_z, lambdaHat_x, lambdaHat_y, lambdaHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternions::Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, chi1_ell, chi1_n, chi1_lambda, chi2_ell, chi2_n, chi2_lambda,
         S_ell, S_n, S_lambda, Sigma_ell, Sigma_n, Sigma_lambda;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_1_1, hHat_2_1_3, hHat_2_2_0, hHat_2_2_2, hHat_2_2_3, hHat_3_1_1,
                             hHat_3_1_3, hHat_3_2_2, hHat_3_3_1, hHat_3_3_3, hHat_4_0_0, hHat_4_1_3, hHat_4_2_2,
                             hHat_4_3_3, hHat_4_4_2, hHat_5_1_3, hHat_5_3_3, hHat_5_5_3;
  std::complex<double> hHat_spin_Symm_2_2_3, hHat_spin_Symm_2_1_2, hHat_spin_Symm_3_2_3, hHat_spin_Asymm_2_2_2,
                       hHat_spin_Asymm_2_1_3, hHat_spin_Asymm_2_0_2, hHat_spin_Asymm_3_3_3, hHat_spin_Asymm_3_1_3;

public:
  WaveformModes_1p5PN(const Quaternions::Quaternion xHat_i, const Quaternions::Quaternion yHat_i, const
                      Quaternions::Quaternion zHat_i, const double m1_i, const double m2_i, const double v_i, const
                      Quaternions::Quaternion S_chi1_i, const Quaternions::Quaternion S_chi2_i, const double
                      rfrak_chi1_x_i, const double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double
                      rfrak_chi2_y_i, const double rfrak_frame_x_i, const double rfrak_frame_y_i, const double
                      rfrak_frame_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), S_chi1(S_chi1_i), S_chi2(S_chi2_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_frame_x(rfrak_frame_x_i), rfrak_frame_y(rfrak_frame_y_i),
    rfrak_frame_z(rfrak_frame_z_i), m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_frame_x*xHat +
    rfrak_frame_y*yHat + rfrak_frame_z*zHat)), nHat(R*xHat*conjugate(R)), lambdaHat(R*yHat*conjugate(R)),
    ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]), lambdaHat_x(lambdaHat[1]),
    lambdaHat_y(lambdaHat[2]), lambdaHat_z(lambdaHat[3]), ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]),
    R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)), R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)),
    chiVec1(S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)),
    chiVec2(S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]),
    chi1_z(chiVec1[3]), chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1_ell(chi1_x*ellHat_x +
    chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z),
    chi1_lambda(chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z), chi2_ell(chi2_x*ellHat_x +
    chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z),
    chi2_lambda(chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z), S_ell(chi1_ell*pow(m1, 2) +
    chi2_ell*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)), S_lambda(chi1_lambda*pow(m1, 2) +
    chi2_lambda*pow(m2, 2)), Sigma_ell(m*(-chi1_ell*m1 + chi2_ell*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)),
    Sigma_lambda(m*(-chi1_lambda*m1 + chi2_lambda*m2)), rhOverM_coeff(6.34132367616962*nu*pow(v, 2)),
    hHat_2_0_0(-0.145802960879951), hHat_2_1_1(0.333333333333333*I*delta),
    hHat_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)), hHat_2_2_0(1.00000000000000),
    hHat_2_2_2(1.30952380952381*nu - 2.54761904761905), hHat_2_2_3(6.28318530717959),
    hHat_3_1_1(0.0222717701593687*I*delta), hHat_3_1_3(-0.0148478467729125*I*delta*(nu + 4.0)),
    hHat_3_2_2(-0.845154254728516*nu + 0.281718084909506), hHat_3_3_1(-0.776323754260148*I*delta),
    hHat_3_3_3(-1.5526475085203*I*delta*(nu - 2.0)), hHat_4_0_0(-0.00140298964521140),
    hHat_4_1_3(0.00376461626210521*I*delta*(-2.0*nu + 1.0)), hHat_4_2_2(-0.10647942749999*nu + 0.0354931424999967),
    hHat_4_3_3(0.268926437100239*I*delta*(2.0*nu - 1.0)), hHat_4_4_2(2.25374467927604*nu - 0.751248226425348),
    hHat_5_1_3(0.000176960830360287*I*delta*(-2.0*nu + 1.0)), hHat_5_3_3(0.0464469088412683*I*delta*(2.0*nu - 1.0)),
    hHat_5_5_3(-0.801376894396698*I*delta*(2.0*nu - 1.0)), hHat_spin_Symm_2_2_3(0.333333333333333*(-11.0*S_ell -
    5.0*Sigma_ell*delta)/pow(m, 2)), hHat_spin_Symm_2_1_2(0.5*I*Sigma_ell/pow(m, 2)),
    hHat_spin_Symm_3_2_3(0.563436169819011*(S_ell + Sigma_ell*delta)/pow(m, 2)),
    hHat_spin_Asymm_2_2_2(0.5*(-Sigma_lambda - I*Sigma_n)/pow(m, 2)),
    hHat_spin_Asymm_2_1_3(0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
    13.0*Sigma_n*delta)/pow(m, 2)), hHat_spin_Asymm_2_0_2(0.408248290463863*I*Sigma_n/pow(m, 2)),
    hHat_spin_Asymm_3_3_3(0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(m, 2)),
    hHat_spin_Asymm_3_1_3(0.17817416127495*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/pow(m, 2))
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_x_k, const double chi1_y_k, const double chi1_z_k,
    const double chi2_x_k, const double chi2_y_k, const double chi2_z_k,
    const double ellHat_x_k, const double ellHat_y_k, const double ellHat_z_k)
  {
    v = v_k;
    chi1_x = chi1_x_k;
    chi1_y = chi1_y_k;
    chi1_z = chi1_z_k;
    chi2_x = chi2_x_k;
    chi2_y = chi2_y_k;
    chi2_z = chi2_z_k;
    ellHat_x = ellHat_x_k;
    ellHat_y = ellHat_y_k;
    ellHat_z = ellHat_z_k;

    R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat);
    nHat = R*xHat*conjugate(R);
    lambdaHat = R*yHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    lambdaHat_x = lambdaHat[1];
    lambdaHat_y = lambdaHat[2];
    lambdaHat_z = lambdaHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1);
    chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1_ell = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi1_lambda = chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z;
    chi2_ell = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    chi2_lambda = chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z;
    S_ell = chi1_ell*pow(m1, 2) + chi2_ell*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    S_lambda = chi1_lambda*pow(m1, 2) + chi2_lambda*pow(m2, 2);
    Sigma_ell = m*(-chi1_ell*m1 + chi2_ell*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    Sigma_lambda = m*(-chi1_lambda*m1 + chi2_lambda*m2);
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    hHat_spin_Symm_2_2_3 = 0.333333333333333*(-11.0*S_ell - 5.0*Sigma_ell*delta)/pow(m, 2);
    hHat_spin_Symm_2_1_2 = 0.5*I*Sigma_ell/pow(m, 2);
    hHat_spin_Symm_3_2_3 = 0.563436169819011*(S_ell + Sigma_ell*delta)/pow(m, 2);
    hHat_spin_Asymm_2_2_2 = 0.5*(-Sigma_lambda - I*Sigma_n)/pow(m, 2);
    hHat_spin_Asymm_2_1_3 = 0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
      13.0*Sigma_n*delta)/pow(m, 2);
    hHat_spin_Asymm_2_0_2 = 0.408248290463863*I*Sigma_n/pow(m, 2);
    hHat_spin_Asymm_3_3_3 = 0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(m, 2);
    hHat_spin_Asymm_3_1_3 = 0.17817416127495*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/pow(m, 2);

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = (((conjugate(hHat_2_2_3) + conjugate(hHat_spin_Symm_2_2_3))*conjugate(v) +
      conjugate(hHat_2_2_2))*pow(conjugate(v), 2) + conjugate(hHat_2_2_0))*conjugate(rhOverM_coeff) -
      conjugate(hHat_spin_Asymm_2_2_2)*conjugate(rhOverM_coeff)*pow(conjugate(v), 2);
    Modes[i++] = ((conjugate(hHat_2_1_3)*conjugate(v) + conjugate(hHat_spin_Symm_2_1_2))*conjugate(v) +
      conjugate(hHat_2_1_1))*conjugate(rhOverM_coeff)*conjugate(v) -
      conjugate(hHat_spin_Asymm_2_1_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = hHat_2_0_0*rhOverM_coeff + hHat_spin_Asymm_2_0_2*rhOverM_coeff*pow(v, 2);
    Modes[i++] = hHat_spin_Asymm_2_1_3*rhOverM_coeff*pow(v, 3) + rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_2_1_3*v +
      hHat_spin_Symm_2_1_2));
    Modes[i++] = hHat_spin_Asymm_2_2_2*rhOverM_coeff*pow(v, 2) + rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 +
      v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3)));
    Modes[i++] = -(conjugate(hHat_3_3_1) + conjugate(hHat_3_3_3)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*conjugate(v) +
      conjugate(hHat_spin_Asymm_3_3_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -(conjugate(hHat_3_2_2) +
      conjugate(hHat_spin_Symm_3_2_3)*conjugate(v))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2);
    Modes[i++] = -(conjugate(hHat_3_1_1) + conjugate(hHat_3_1_3)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*conjugate(v) +
      conjugate(hHat_spin_Asymm_3_1_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = 0;
    Modes[i++] = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3) + rhOverM_coeff*v*(hHat_3_1_1 + hHat_3_1_3*pow(v, 2));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 + hHat_spin_Symm_3_2_3*v);
    Modes[i++] = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3) + rhOverM_coeff*v*(hHat_3_3_1 + hHat_3_3_3*pow(v, 2));
    Modes[i++] = conjugate(hHat_4_4_2)*conjugate(rhOverM_coeff)*pow(conjugate(v), 2);
    Modes[i++] = conjugate(hHat_4_3_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = conjugate(hHat_4_2_2)*conjugate(rhOverM_coeff)*pow(conjugate(v), 2);
    Modes[i++] = conjugate(hHat_4_1_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = hHat_4_0_0*rhOverM_coeff;
    Modes[i++] = hHat_4_1_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = hHat_4_2_2*rhOverM_coeff*pow(v, 2);
    Modes[i++] = hHat_4_3_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = hHat_4_4_2*rhOverM_coeff*pow(v, 2);
    Modes[i++] = -conjugate(hHat_5_5_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = 0;
    Modes[i++] = -conjugate(hHat_5_3_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = 0;
    Modes[i++] = -conjugate(hHat_5_1_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = 0;
    Modes[i++] = hHat_5_1_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = 0;
    Modes[i++] = hHat_5_3_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = 0;
    Modes[i++] = hHat_5_5_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;

    return Modes;
  }

}; // class WaveformModes_1p5PN : public WaveformModes_Base


class WaveformModes_2p0PN : public WaveformModes_Base {
private:
  const Quaternions::Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const Quaternions::Quaternion S_chi1, S_chi2;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z;
  const double m, delta, nu;
  Quaternions::Quaternion R, nHat, lambdaHat, ellHat;
  double nHat_x, nHat_y, nHat_z, lambdaHat_x, lambdaHat_y, lambdaHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternions::Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, chi1_ell, chi1_n, chi1_lambda, chi2_ell, chi2_n, chi2_lambda,
         S_ell, S_n, S_lambda, Sigma_ell, Sigma_n, Sigma_lambda, S1_ell, S1_n, S1_lambda, S2_ell, S2_n, S2_lambda;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> hHat_2_0_0, hHat_2_1_1, hHat_2_1_3, hHat_2_1_4, hHat_2_2_0, hHat_2_2_2, hHat_2_2_3,
                             hHat_2_2_4, hHat_3_1_1, hHat_3_1_3, hHat_3_1_4, hHat_3_2_2, hHat_3_2_4, hHat_3_3_1,
                             hHat_3_3_3, hHat_3_3_4, hHat_4_0_0, hHat_4_1_3, hHat_4_2_2, hHat_4_2_4, hHat_4_3_3,
                             hHat_4_4_2, hHat_4_4_4, hHat_5_1_3, hHat_5_2_4, hHat_5_3_3, hHat_5_4_4, hHat_5_5_3,
                             hHat_6_2_4, hHat_6_4_4, hHat_6_6_4;
  std::complex<double> hHat_spin_Symm_2_2_3, hHat_spin_Symm_2_2_4, hHat_spin_Symm_2_1_2, hHat_spin_Symm_2_1_4,
                       hHat_spin_Symm_2_0_4, hHat_spin_Symm_3_3_4, hHat_spin_Symm_3_2_3, hHat_spin_Symm_3_1_4,
                       hHat_spin_Symm_4_3_4, hHat_spin_Symm_4_1_4, hHat_spin_Asymm_2_2_2, hHat_spin_Asymm_2_2_4,
                       hHat_spin_Asymm_2_1_3, hHat_spin_Asymm_2_1_4, hHat_spin_Asymm_2_0_2, hHat_spin_Asymm_2_0_4,
                       hHat_spin_Asymm_3_3_3, hHat_spin_Asymm_3_2_4, hHat_spin_Asymm_3_1_3, hHat_spin_Asymm_3_0_4,
                       hHat_spin_Asymm_4_4_4, hHat_spin_Asymm_4_2_4, hHat_spin_Asymm_4_0_4;

public:
  WaveformModes_2p0PN(const Quaternions::Quaternion xHat_i, const Quaternions::Quaternion yHat_i, const
                      Quaternions::Quaternion zHat_i, const double m1_i, const double m2_i, const double v_i, const
                      Quaternions::Quaternion S_chi1_i, const Quaternions::Quaternion S_chi2_i, const double
                      rfrak_chi1_x_i, const double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double
                      rfrak_chi2_y_i, const double rfrak_frame_x_i, const double rfrak_frame_y_i, const double
                      rfrak_frame_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), S_chi1(S_chi1_i), S_chi2(S_chi2_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_frame_x(rfrak_frame_x_i), rfrak_frame_y(rfrak_frame_y_i),
    rfrak_frame_z(rfrak_frame_z_i), m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_frame_x*xHat +
    rfrak_frame_y*yHat + rfrak_frame_z*zHat)), nHat(R*xHat*conjugate(R)), lambdaHat(R*yHat*conjugate(R)),
    ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]), lambdaHat_x(lambdaHat[1]),
    lambdaHat_y(lambdaHat[2]), lambdaHat_z(lambdaHat[3]), ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]),
    R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)), R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)),
    chiVec1(S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)),
    chiVec2(S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]),
    chi1_z(chiVec1[3]), chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1_ell(chi1_x*ellHat_x +
    chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z),
    chi1_lambda(chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z), chi2_ell(chi2_x*ellHat_x +
    chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z),
    chi2_lambda(chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z), S_ell(chi1_ell*pow(m1, 2) +
    chi2_ell*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)), S_lambda(chi1_lambda*pow(m1, 2) +
    chi2_lambda*pow(m2, 2)), Sigma_ell(m*(-chi1_ell*m1 + chi2_ell*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)),
    Sigma_lambda(m*(-chi1_lambda*m1 + chi2_lambda*m2)), S1_ell(chi1_ell*pow(m1, 2)), S1_n(chi1_n*pow(m1, 2)),
    S1_lambda(chi1_lambda*pow(m1, 2)), S2_ell(chi2_ell*pow(m2, 2)), S2_n(chi2_n*pow(m2, 2)),
    S2_lambda(chi2_lambda*pow(m2, 2)), rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), hHat_2_0_0(-0.145802960879951),
    hHat_2_1_1(0.333333333333333*I*delta), hHat_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)),
    hHat_2_1_4(delta*(0.628764787039964 + 1.0471975511966*I)), hHat_2_2_0(1.00000000000000),
    hHat_2_2_2(1.30952380952381*nu - 2.54761904761905), hHat_2_2_3(6.28318530717959),
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
    + 0.903141370807658), hHat_spin_Symm_2_2_3(0.333333333333333*(-11.0*S_ell - 5.0*Sigma_ell*delta)/pow(m, 2)),
    hHat_spin_Symm_2_2_4(0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
    15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(m, 4)*nu)), hHat_spin_Symm_2_1_2(0.5*I*Sigma_ell/pow(m, 2)),
    hHat_spin_Symm_2_1_4(0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(m, 2)),
    hHat_spin_Symm_2_0_4(0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(m, 4)*nu)),
    hHat_spin_Symm_3_3_4(0.0143763658196324*I*delta*(193.0*S_ell + 78732.0*Sigma_ell*nu)/pow(m, 2)),
    hHat_spin_Symm_3_2_3(0.563436169819011*(S_ell + Sigma_ell*delta)/pow(m, 2)),
    hHat_spin_Symm_3_1_4(0.0556794253984217*I*delta*(S_ell + 60.0*Sigma_ell*nu)/pow(m, 2)),
    hHat_spin_Symm_4_3_4(0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(m, 2)),
    hHat_spin_Symm_4_1_4(0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(m, 2)),
    hHat_spin_Asymm_2_2_2(0.5*(-Sigma_lambda - I*Sigma_n)/pow(m, 2)),
    hHat_spin_Asymm_2_2_4(0.00396825396825397*(-54180.0*Sigma_lambda*delta*nu + 5880.0*I*Sigma_n*nu +
    delta*(29.0*S_lambda + 546.0*I*S_n))/pow(m, 2)), hHat_spin_Asymm_2_1_3(0.166666666666667*(4.0*I*S_lambda + 25.0*S_n
    + 4.0*I*Sigma_lambda*delta + 13.0*Sigma_n*delta)/pow(m, 2)), hHat_spin_Asymm_2_1_4(0.5*(-3.0*S1_ell*S2_n -
    3.0*S1_n*S2_ell)/(pow(m, 4)*nu)), hHat_spin_Asymm_2_0_2(0.408248290463863*I*Sigma_n/pow(m, 2)),
    hHat_spin_Asymm_2_0_4(0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(m, 2)),
    hHat_spin_Asymm_3_3_3(0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(m, 2)),
    hHat_spin_Asymm_3_2_4(0.0117382535378961*(-50796.0*Sigma_lambda*delta*nu - 8580.0*I*Sigma_n*nu +
    delta*(71.0*S_lambda - 300.0*I*S_n))/pow(m, 2)), hHat_spin_Asymm_3_1_3(0.17817416127495*(I*S_lambda + S_n +
    delta*(I*Sigma_lambda + Sigma_n))/pow(m, 2)), hHat_spin_Asymm_3_0_4(-0.064293062484205*delta*(11.0*S_lambda +
    2268.0*Sigma_lambda*nu)/pow(m, 2)), hHat_spin_Asymm_4_4_4(0.950798536569581*(-3.0*Sigma_lambda*nu + Sigma_lambda -
    I*Sigma_n*(3.0*nu - 1.0) + delta*(S_lambda + I*S_n))/pow(m, 2)),
    hHat_spin_Asymm_4_2_4(0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) - 42.0*I*Sigma_n*nu +
    delta*(13.0*S_lambda - 14.0*I*S_n))/pow(m, 2)), hHat_spin_Asymm_4_0_4(0.00841793787126842*I*(S_n*delta -
    3.0*Sigma_n*nu + Sigma_n)/pow(m, 2))
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_x_k, const double chi1_y_k, const double chi1_z_k,
    const double chi2_x_k, const double chi2_y_k, const double chi2_z_k,
    const double ellHat_x_k, const double ellHat_y_k, const double ellHat_z_k)
  {
    v = v_k;
    chi1_x = chi1_x_k;
    chi1_y = chi1_y_k;
    chi1_z = chi1_z_k;
    chi2_x = chi2_x_k;
    chi2_y = chi2_y_k;
    chi2_z = chi2_z_k;
    ellHat_x = ellHat_x_k;
    ellHat_y = ellHat_y_k;
    ellHat_z = ellHat_z_k;

    R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat);
    nHat = R*xHat*conjugate(R);
    lambdaHat = R*yHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    lambdaHat_x = lambdaHat[1];
    lambdaHat_y = lambdaHat[2];
    lambdaHat_z = lambdaHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1);
    chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1_ell = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi1_lambda = chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z;
    chi2_ell = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    chi2_lambda = chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z;
    S_ell = chi1_ell*pow(m1, 2) + chi2_ell*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    S_lambda = chi1_lambda*pow(m1, 2) + chi2_lambda*pow(m2, 2);
    Sigma_ell = m*(-chi1_ell*m1 + chi2_ell*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    Sigma_lambda = m*(-chi1_lambda*m1 + chi2_lambda*m2);
    S1_ell = chi1_ell*pow(m1, 2);
    S1_n = chi1_n*pow(m1, 2);
    S1_lambda = chi1_lambda*pow(m1, 2);
    S2_ell = chi2_ell*pow(m2, 2);
    S2_n = chi2_n*pow(m2, 2);
    S2_lambda = chi2_lambda*pow(m2, 2);
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    hHat_spin_Symm_2_2_3 = 0.333333333333333*(-11.0*S_ell - 5.0*Sigma_ell*delta)/pow(m, 2);
    hHat_spin_Symm_2_2_4 = 0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
      15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(m, 4)*nu);
    hHat_spin_Symm_2_1_2 = 0.5*I*Sigma_ell/pow(m, 2);
    hHat_spin_Symm_2_1_4 = 0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(m, 2);
    hHat_spin_Symm_2_0_4 = 0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(m, 4)*nu);
    hHat_spin_Symm_3_3_4 = 0.0143763658196324*I*delta*(193.0*S_ell + 78732.0*Sigma_ell*nu)/pow(m, 2);
    hHat_spin_Symm_3_2_3 = 0.563436169819011*(S_ell + Sigma_ell*delta)/pow(m, 2);
    hHat_spin_Symm_3_1_4 = 0.0556794253984217*I*delta*(S_ell + 60.0*Sigma_ell*nu)/pow(m, 2);
    hHat_spin_Symm_4_3_4 = 0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(m, 2);
    hHat_spin_Symm_4_1_4 = 0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(m, 2);
    hHat_spin_Asymm_2_2_2 = 0.5*(-Sigma_lambda - I*Sigma_n)/pow(m, 2);
    hHat_spin_Asymm_2_2_4 = 0.00396825396825397*(-54180.0*Sigma_lambda*delta*nu + 5880.0*I*Sigma_n*nu +
      delta*(29.0*S_lambda + 546.0*I*S_n))/pow(m, 2);
    hHat_spin_Asymm_2_1_3 = 0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
      13.0*Sigma_n*delta)/pow(m, 2);
    hHat_spin_Asymm_2_1_4 = 0.5*(-3.0*S1_ell*S2_n - 3.0*S1_n*S2_ell)/(pow(m, 4)*nu);
    hHat_spin_Asymm_2_0_2 = 0.408248290463863*I*Sigma_n/pow(m, 2);
    hHat_spin_Asymm_2_0_4 = 0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(m, 2);
    hHat_spin_Asymm_3_3_3 = 0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(m, 2);
    hHat_spin_Asymm_3_2_4 = 0.0117382535378961*(-50796.0*Sigma_lambda*delta*nu - 8580.0*I*Sigma_n*nu +
      delta*(71.0*S_lambda - 300.0*I*S_n))/pow(m, 2);
    hHat_spin_Asymm_3_1_3 = 0.17817416127495*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/pow(m, 2);
    hHat_spin_Asymm_3_0_4 = -0.064293062484205*delta*(11.0*S_lambda + 2268.0*Sigma_lambda*nu)/pow(m, 2);
    hHat_spin_Asymm_4_4_4 = 0.950798536569581*(-3.0*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3.0*nu - 1.0) +
      delta*(S_lambda + I*S_n))/pow(m, 2);
    hHat_spin_Asymm_4_2_4 = 0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) - 42.0*I*Sigma_n*nu +
      delta*(13.0*S_lambda - 14.0*I*S_n))/pow(m, 2);
    hHat_spin_Asymm_4_0_4 = 0.00841793787126842*I*(S_n*delta - 3.0*Sigma_n*nu + Sigma_n)/pow(m, 2);

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = ((((conjugate(hHat_2_2_4) + conjugate(hHat_spin_Symm_2_2_4))*conjugate(v) + conjugate(hHat_2_2_3) +
      conjugate(hHat_spin_Symm_2_2_3))*conjugate(v) + conjugate(hHat_2_2_2))*pow(conjugate(v), 2) +
      conjugate(hHat_2_2_0))*conjugate(rhOverM_coeff) - (conjugate(hHat_spin_Asymm_2_2_2) +
      conjugate(hHat_spin_Asymm_2_2_4)*pow(conjugate(v), 2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2);
    Modes[i++] = ((((conjugate(hHat_2_1_4) + conjugate(hHat_spin_Symm_2_1_4))*conjugate(v) +
      conjugate(hHat_2_1_3))*conjugate(v) + conjugate(hHat_spin_Symm_2_1_2))*conjugate(v) +
      conjugate(hHat_2_1_1))*conjugate(rhOverM_coeff)*conjugate(v) - (conjugate(hHat_spin_Asymm_2_1_3) +
      conjugate(hHat_spin_Asymm_2_1_4)*conjugate(v))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*pow(v, 2)) +
      rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_4*pow(v, 4));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v) + rhOverM_coeff*v*(hHat_2_1_1
      + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_spin_Symm_2_1_4))));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*pow(v, 2)) +
      rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 +
      hHat_spin_Symm_2_2_4))));
    Modes[i++] = -(((conjugate(hHat_3_3_4) + conjugate(hHat_spin_Symm_3_3_4))*conjugate(v) +
      conjugate(hHat_3_3_3))*pow(conjugate(v), 2) + conjugate(hHat_3_3_1))*conjugate(rhOverM_coeff)*conjugate(v) +
      conjugate(hHat_spin_Asymm_3_3_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -((conjugate(hHat_3_2_4)*conjugate(v) + conjugate(hHat_spin_Symm_3_2_3))*conjugate(v) +
      conjugate(hHat_3_2_2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2) +
      conjugate(hHat_spin_Asymm_3_2_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = -(((conjugate(hHat_3_1_4) + conjugate(hHat_spin_Symm_3_1_4))*conjugate(v) +
      conjugate(hHat_3_1_3))*pow(conjugate(v), 2) + conjugate(hHat_3_1_1))*conjugate(rhOverM_coeff)*conjugate(v) +
      conjugate(hHat_spin_Asymm_3_1_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = hHat_spin_Asymm_3_0_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3) + rhOverM_coeff*v*(hHat_3_1_1 + pow(v, 2)*(hHat_3_1_3 +
      v*(hHat_3_1_4 + hHat_spin_Symm_3_1_4)));
    Modes[i++] = hHat_spin_Asymm_3_2_4*rhOverM_coeff*pow(v, 4) + rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 + v*(hHat_3_2_4*v +
      hHat_spin_Symm_3_2_3));
    Modes[i++] = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3) + rhOverM_coeff*v*(hHat_3_3_1 + pow(v, 2)*(hHat_3_3_3 +
      v*(hHat_3_3_4 + hHat_spin_Symm_3_3_4)));
    Modes[i++] = (conjugate(hHat_4_4_2) + conjugate(hHat_4_4_4)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2) -
      conjugate(hHat_spin_Asymm_4_4_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = (conjugate(hHat_4_3_3) +
      conjugate(hHat_spin_Symm_4_3_4)*conjugate(v))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = (conjugate(hHat_4_2_2) + conjugate(hHat_4_2_4)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2) -
      conjugate(hHat_spin_Asymm_4_2_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = (conjugate(hHat_4_1_3) +
      conjugate(hHat_spin_Symm_4_1_4)*conjugate(v))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = hHat_4_0_0*rhOverM_coeff + hHat_spin_Asymm_4_0_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_4_1_3 + hHat_spin_Symm_4_1_4*v);
    Modes[i++] = hHat_spin_Asymm_4_2_4*rhOverM_coeff*pow(v, 4) + rhOverM_coeff*pow(v, 2)*(hHat_4_2_2 + hHat_4_2_4*pow(v,
      2));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_4_3_3 + hHat_spin_Symm_4_3_4*v);
    Modes[i++] = hHat_spin_Asymm_4_4_4*rhOverM_coeff*pow(v, 4) + rhOverM_coeff*pow(v, 2)*(hHat_4_4_2 + hHat_4_4_4*pow(v,
      2));
    Modes[i++] = -conjugate(hHat_5_5_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -conjugate(hHat_5_4_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = -conjugate(hHat_5_3_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -conjugate(hHat_5_2_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = -conjugate(hHat_5_1_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = 0;
    Modes[i++] = hHat_5_1_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = hHat_5_2_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = hHat_5_3_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = hHat_5_4_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = hHat_5_5_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = conjugate(hHat_6_6_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = 0;
    Modes[i++] = conjugate(hHat_6_4_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = 0;
    Modes[i++] = conjugate(hHat_6_2_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = hHat_6_2_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = 0;
    Modes[i++] = hHat_6_4_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = 0;
    Modes[i++] = hHat_6_6_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;

    return Modes;
  }

}; // class WaveformModes_2p0PN : public WaveformModes_Base


class WaveformModes_2p5PN : public WaveformModes_Base {
private:
  const Quaternions::Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const Quaternions::Quaternion S_chi1, S_chi2;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z;
  const double m, delta, nu;
  Quaternions::Quaternion R, nHat, lambdaHat, ellHat;
  double nHat_x, nHat_y, nHat_z, lambdaHat_x, lambdaHat_y, lambdaHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternions::Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, chi1_ell, chi1_n, chi1_lambda, chi2_ell, chi2_n, chi2_lambda,
         S_ell, S_n, S_lambda, Sigma_ell, Sigma_n, Sigma_lambda, S1_ell, S1_n, S1_lambda, S2_ell, S2_n, S2_lambda;
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
                       hHat_spin_Symm_2_0_4, hHat_spin_Symm_3_3_4, hHat_spin_Symm_3_2_3, hHat_spin_Symm_3_1_4,
                       hHat_spin_Symm_4_3_4, hHat_spin_Symm_4_1_4, hHat_spin_Asymm_2_2_2, hHat_spin_Asymm_2_2_4,
                       hHat_spin_Asymm_2_1_3, hHat_spin_Asymm_2_1_4, hHat_spin_Asymm_2_0_2, hHat_spin_Asymm_2_0_4,
                       hHat_spin_Asymm_3_3_3, hHat_spin_Asymm_3_2_4, hHat_spin_Asymm_3_1_3, hHat_spin_Asymm_3_0_4,
                       hHat_spin_Asymm_4_4_4, hHat_spin_Asymm_4_2_4, hHat_spin_Asymm_4_0_4;

public:
  WaveformModes_2p5PN(const Quaternions::Quaternion xHat_i, const Quaternions::Quaternion yHat_i, const
                      Quaternions::Quaternion zHat_i, const double m1_i, const double m2_i, const double v_i, const
                      Quaternions::Quaternion S_chi1_i, const Quaternions::Quaternion S_chi2_i, const double
                      rfrak_chi1_x_i, const double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double
                      rfrak_chi2_y_i, const double rfrak_frame_x_i, const double rfrak_frame_y_i, const double
                      rfrak_frame_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), S_chi1(S_chi1_i), S_chi2(S_chi2_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_frame_x(rfrak_frame_x_i), rfrak_frame_y(rfrak_frame_y_i),
    rfrak_frame_z(rfrak_frame_z_i), m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_frame_x*xHat +
    rfrak_frame_y*yHat + rfrak_frame_z*zHat)), nHat(R*xHat*conjugate(R)), lambdaHat(R*yHat*conjugate(R)),
    ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]), lambdaHat_x(lambdaHat[1]),
    lambdaHat_y(lambdaHat[2]), lambdaHat_z(lambdaHat[3]), ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]),
    R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)), R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)),
    chiVec1(S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)),
    chiVec2(S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]),
    chi1_z(chiVec1[3]), chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1_ell(chi1_x*ellHat_x +
    chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z),
    chi1_lambda(chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z), chi2_ell(chi2_x*ellHat_x +
    chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z),
    chi2_lambda(chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z), S_ell(chi1_ell*pow(m1, 2) +
    chi2_ell*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)), S_lambda(chi1_lambda*pow(m1, 2) +
    chi2_lambda*pow(m2, 2)), Sigma_ell(m*(-chi1_ell*m1 + chi2_ell*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)),
    Sigma_lambda(m*(-chi1_lambda*m1 + chi2_lambda*m2)), S1_ell(chi1_ell*pow(m1, 2)), S1_n(chi1_n*pow(m1, 2)),
    S1_lambda(chi1_lambda*pow(m1, 2)), S2_ell(chi2_ell*pow(m2, 2)), S2_n(chi2_n*pow(m2, 2)),
    S2_lambda(chi2_lambda*pow(m2, 2)), rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), hHat_2_0_0(-0.145802960879951),
    hHat_2_1_1(0.333333333333333*I*delta), hHat_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)),
    hHat_2_1_4(delta*(0.628764787039964 + 1.0471975511966*I)), hHat_2_1_5(0.000661375661375661*I*delta*(nu*(237.0*nu -
    2036.0) - 172.0)), hHat_2_2_0(1.00000000000000), hHat_2_2_2(1.30952380952381*nu - 2.54761904761905),
    hHat_2_2_3(6.28318530717959), hHat_2_2_4(0.000661375661375661*nu*(2047.0*nu - 7483.0) - 1.43716931216931),
    hHat_2_2_5(5.08638810581205*nu - 24.0*I*nu - 16.0071625682908), hHat_3_0_5(-0.370328039909021*I*nu),
    hHat_3_1_1(0.0222717701593687*I*delta), hHat_3_1_3(-0.0148478467729125*I*delta*(nu + 4.0)),
    hHat_3_1_4(delta*(0.0620557076072073 + 0.0699688295151131*I)),
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
    hHat_7_7_5(-1.05422444934392*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), hHat_spin_Symm_2_2_3(0.333333333333333*(-11.0*S_ell
    - 5.0*Sigma_ell*delta)/pow(m, 2)), hHat_spin_Symm_2_2_4(0.166666666666667*(12.0*S1_ell*S2_ell +
    10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n - 15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(m, 4)*nu)),
    hHat_spin_Symm_2_1_2(0.5*I*Sigma_ell/pow(m, 2)), hHat_spin_Symm_2_1_4(0.0238095238095238*I*(-86.0*S_ell*delta +
    Sigma_ell*(139.0*nu - 79.0))/pow(m, 2)), hHat_spin_Symm_2_0_4(0.816496580927726*(-S1_lambda*S2_lambda +
    S1_n*S2_n)/(pow(m, 4)*nu)), hHat_spin_Symm_3_3_4(0.0143763658196324*I*delta*(193.0*S_ell +
    78732.0*Sigma_ell*nu)/pow(m, 2)), hHat_spin_Symm_3_2_3(0.563436169819011*(S_ell + Sigma_ell*delta)/pow(m, 2)),
    hHat_spin_Symm_3_1_4(0.0556794253984217*I*delta*(S_ell + 60.0*Sigma_ell*nu)/pow(m, 2)),
    hHat_spin_Symm_4_3_4(0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(m, 2)),
    hHat_spin_Symm_4_1_4(0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(m, 2)),
    hHat_spin_Asymm_2_2_2(0.5*(-Sigma_lambda - I*Sigma_n)/pow(m, 2)),
    hHat_spin_Asymm_2_2_4(0.00396825396825397*(-54180.0*Sigma_lambda*delta*nu + 5880.0*I*Sigma_n*nu +
    delta*(29.0*S_lambda + 546.0*I*S_n))/pow(m, 2)), hHat_spin_Asymm_2_1_3(0.166666666666667*(4.0*I*S_lambda + 25.0*S_n
    + 4.0*I*Sigma_lambda*delta + 13.0*Sigma_n*delta)/pow(m, 2)), hHat_spin_Asymm_2_1_4(0.5*(-3.0*S1_ell*S2_n -
    3.0*S1_n*S2_ell)/(pow(m, 4)*nu)), hHat_spin_Asymm_2_0_2(0.408248290463863*I*Sigma_n/pow(m, 2)),
    hHat_spin_Asymm_2_0_4(0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(m, 2)),
    hHat_spin_Asymm_3_3_3(0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(m, 2)),
    hHat_spin_Asymm_3_2_4(0.0117382535378961*(-50796.0*Sigma_lambda*delta*nu - 8580.0*I*Sigma_n*nu +
    delta*(71.0*S_lambda - 300.0*I*S_n))/pow(m, 2)), hHat_spin_Asymm_3_1_3(0.17817416127495*(I*S_lambda + S_n +
    delta*(I*Sigma_lambda + Sigma_n))/pow(m, 2)), hHat_spin_Asymm_3_0_4(-0.064293062484205*delta*(11.0*S_lambda +
    2268.0*Sigma_lambda*nu)/pow(m, 2)), hHat_spin_Asymm_4_4_4(0.950798536569581*(-3.0*Sigma_lambda*nu + Sigma_lambda -
    I*Sigma_n*(3.0*nu - 1.0) + delta*(S_lambda + I*S_n))/pow(m, 2)),
    hHat_spin_Asymm_4_2_4(0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) - 42.0*I*Sigma_n*nu +
    delta*(13.0*S_lambda - 14.0*I*S_n))/pow(m, 2)), hHat_spin_Asymm_4_0_4(0.00841793787126842*I*(S_n*delta -
    3.0*Sigma_n*nu + Sigma_n)/pow(m, 2))
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_x_k, const double chi1_y_k, const double chi1_z_k,
    const double chi2_x_k, const double chi2_y_k, const double chi2_z_k,
    const double ellHat_x_k, const double ellHat_y_k, const double ellHat_z_k)
  {
    v = v_k;
    chi1_x = chi1_x_k;
    chi1_y = chi1_y_k;
    chi1_z = chi1_z_k;
    chi2_x = chi2_x_k;
    chi2_y = chi2_y_k;
    chi2_z = chi2_z_k;
    ellHat_x = ellHat_x_k;
    ellHat_y = ellHat_y_k;
    ellHat_z = ellHat_z_k;

    R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat);
    nHat = R*xHat*conjugate(R);
    lambdaHat = R*yHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    lambdaHat_x = lambdaHat[1];
    lambdaHat_y = lambdaHat[2];
    lambdaHat_z = lambdaHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1);
    chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1_ell = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi1_lambda = chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z;
    chi2_ell = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    chi2_lambda = chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z;
    S_ell = chi1_ell*pow(m1, 2) + chi2_ell*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    S_lambda = chi1_lambda*pow(m1, 2) + chi2_lambda*pow(m2, 2);
    Sigma_ell = m*(-chi1_ell*m1 + chi2_ell*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    Sigma_lambda = m*(-chi1_lambda*m1 + chi2_lambda*m2);
    S1_ell = chi1_ell*pow(m1, 2);
    S1_n = chi1_n*pow(m1, 2);
    S1_lambda = chi1_lambda*pow(m1, 2);
    S2_ell = chi2_ell*pow(m2, 2);
    S2_n = chi2_n*pow(m2, 2);
    S2_lambda = chi2_lambda*pow(m2, 2);
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    hHat_spin_Symm_2_2_3 = 0.333333333333333*(-11.0*S_ell - 5.0*Sigma_ell*delta)/pow(m, 2);
    hHat_spin_Symm_2_2_4 = 0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
      15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(m, 4)*nu);
    hHat_spin_Symm_2_1_2 = 0.5*I*Sigma_ell/pow(m, 2);
    hHat_spin_Symm_2_1_4 = 0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(m, 2);
    hHat_spin_Symm_2_0_4 = 0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(m, 4)*nu);
    hHat_spin_Symm_3_3_4 = 0.0143763658196324*I*delta*(193.0*S_ell + 78732.0*Sigma_ell*nu)/pow(m, 2);
    hHat_spin_Symm_3_2_3 = 0.563436169819011*(S_ell + Sigma_ell*delta)/pow(m, 2);
    hHat_spin_Symm_3_1_4 = 0.0556794253984217*I*delta*(S_ell + 60.0*Sigma_ell*nu)/pow(m, 2);
    hHat_spin_Symm_4_3_4 = 0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(m, 2);
    hHat_spin_Symm_4_1_4 = 0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(m, 2);
    hHat_spin_Asymm_2_2_2 = 0.5*(-Sigma_lambda - I*Sigma_n)/pow(m, 2);
    hHat_spin_Asymm_2_2_4 = 0.00396825396825397*(-54180.0*Sigma_lambda*delta*nu + 5880.0*I*Sigma_n*nu +
      delta*(29.0*S_lambda + 546.0*I*S_n))/pow(m, 2);
    hHat_spin_Asymm_2_1_3 = 0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
      13.0*Sigma_n*delta)/pow(m, 2);
    hHat_spin_Asymm_2_1_4 = 0.5*(-3.0*S1_ell*S2_n - 3.0*S1_n*S2_ell)/(pow(m, 4)*nu);
    hHat_spin_Asymm_2_0_2 = 0.408248290463863*I*Sigma_n/pow(m, 2);
    hHat_spin_Asymm_2_0_4 = 0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(m, 2);
    hHat_spin_Asymm_3_3_3 = 0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(m, 2);
    hHat_spin_Asymm_3_2_4 = 0.0117382535378961*(-50796.0*Sigma_lambda*delta*nu - 8580.0*I*Sigma_n*nu +
      delta*(71.0*S_lambda - 300.0*I*S_n))/pow(m, 2);
    hHat_spin_Asymm_3_1_3 = 0.17817416127495*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/pow(m, 2);
    hHat_spin_Asymm_3_0_4 = -0.064293062484205*delta*(11.0*S_lambda + 2268.0*Sigma_lambda*nu)/pow(m, 2);
    hHat_spin_Asymm_4_4_4 = 0.950798536569581*(-3.0*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3.0*nu - 1.0) +
      delta*(S_lambda + I*S_n))/pow(m, 2);
    hHat_spin_Asymm_4_2_4 = 0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) - 42.0*I*Sigma_n*nu +
      delta*(13.0*S_lambda - 14.0*I*S_n))/pow(m, 2);
    hHat_spin_Asymm_4_0_4 = 0.00841793787126842*I*(S_n*delta - 3.0*Sigma_n*nu + Sigma_n)/pow(m, 2);

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = ((((conjugate(hHat_2_2_4) + conjugate(hHat_2_2_5)*conjugate(v) +
      conjugate(hHat_spin_Symm_2_2_4))*conjugate(v) + conjugate(hHat_2_2_3) +
      conjugate(hHat_spin_Symm_2_2_3))*conjugate(v) + conjugate(hHat_2_2_2))*pow(conjugate(v), 2) +
      conjugate(hHat_2_2_0))*conjugate(rhOverM_coeff) - (conjugate(hHat_spin_Asymm_2_2_2) +
      conjugate(hHat_spin_Asymm_2_2_4)*pow(conjugate(v), 2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2);
    Modes[i++] = ((((conjugate(hHat_2_1_4) + conjugate(hHat_2_1_5)*conjugate(v) +
      conjugate(hHat_spin_Symm_2_1_4))*conjugate(v) + conjugate(hHat_2_1_3))*conjugate(v) +
      conjugate(hHat_spin_Symm_2_1_2))*conjugate(v) + conjugate(hHat_2_1_1))*conjugate(rhOverM_coeff)*conjugate(v) -
      (conjugate(hHat_spin_Asymm_2_1_3) +
      conjugate(hHat_spin_Asymm_2_1_4)*conjugate(v))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*pow(v, 2)) +
      rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_4*pow(v, 4));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v) + rhOverM_coeff*v*(hHat_2_1_1
      + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_2_1_5*v + hHat_spin_Symm_2_1_4))));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*pow(v, 2)) +
      rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 +
      hHat_2_2_5*v + hHat_spin_Symm_2_2_4))));
    Modes[i++] = -(((conjugate(hHat_3_3_4) + conjugate(hHat_3_3_5)*conjugate(v) +
      conjugate(hHat_spin_Symm_3_3_4))*conjugate(v) + conjugate(hHat_3_3_3))*pow(conjugate(v), 2) +
      conjugate(hHat_3_3_1))*conjugate(rhOverM_coeff)*conjugate(v) +
      conjugate(hHat_spin_Asymm_3_3_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -(((conjugate(hHat_3_2_4) + conjugate(hHat_3_2_5)*conjugate(v))*conjugate(v) +
      conjugate(hHat_spin_Symm_3_2_3))*conjugate(v) + conjugate(hHat_3_2_2))*conjugate(rhOverM_coeff)*pow(conjugate(v),
      2) + conjugate(hHat_spin_Asymm_3_2_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = -(((conjugate(hHat_3_1_4) + conjugate(hHat_3_1_5)*conjugate(v) +
      conjugate(hHat_spin_Symm_3_1_4))*conjugate(v) + conjugate(hHat_3_1_3))*pow(conjugate(v), 2) +
      conjugate(hHat_3_1_1))*conjugate(rhOverM_coeff)*conjugate(v) +
      conjugate(hHat_spin_Asymm_3_1_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = hHat_3_0_5*rhOverM_coeff*pow(v, 5) + hHat_spin_Asymm_3_0_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3) + rhOverM_coeff*v*(hHat_3_1_1 + pow(v, 2)*(hHat_3_1_3 +
      v*(hHat_3_1_4 + hHat_3_1_5*v + hHat_spin_Symm_3_1_4)));
    Modes[i++] = hHat_spin_Asymm_3_2_4*rhOverM_coeff*pow(v, 4) + rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 +
      v*(hHat_spin_Symm_3_2_3 + v*(hHat_3_2_4 + hHat_3_2_5*v)));
    Modes[i++] = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3) + rhOverM_coeff*v*(hHat_3_3_1 + pow(v, 2)*(hHat_3_3_3 +
      v*(hHat_3_3_4 + hHat_3_3_5*v + hHat_spin_Symm_3_3_4)));
    Modes[i++] = ((conjugate(hHat_4_4_4) + conjugate(hHat_4_4_5)*conjugate(v))*pow(conjugate(v), 2) +
      conjugate(hHat_4_4_2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2) -
      conjugate(hHat_spin_Asymm_4_4_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = ((conjugate(hHat_4_3_5)*conjugate(v) + conjugate(hHat_spin_Symm_4_3_4))*conjugate(v) +
      conjugate(hHat_4_3_3))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = ((conjugate(hHat_4_2_4) + conjugate(hHat_4_2_5)*conjugate(v))*pow(conjugate(v), 2) +
      conjugate(hHat_4_2_2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2) -
      conjugate(hHat_spin_Asymm_4_2_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = ((conjugate(hHat_4_1_5)*conjugate(v) + conjugate(hHat_spin_Symm_4_1_4))*conjugate(v) +
      conjugate(hHat_4_1_3))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = hHat_4_0_0*rhOverM_coeff + hHat_spin_Asymm_4_0_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_4_1_3 + v*(hHat_4_1_5*v + hHat_spin_Symm_4_1_4));
    Modes[i++] = hHat_spin_Asymm_4_2_4*rhOverM_coeff*pow(v, 4) + rhOverM_coeff*pow(v, 2)*(hHat_4_2_2 + pow(v,
      2)*(hHat_4_2_4 + hHat_4_2_5*v));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_4_3_3 + v*(hHat_4_3_5*v + hHat_spin_Symm_4_3_4));
    Modes[i++] = hHat_spin_Asymm_4_4_4*rhOverM_coeff*pow(v, 4) + rhOverM_coeff*pow(v, 2)*(hHat_4_4_2 + pow(v,
      2)*(hHat_4_4_4 + hHat_4_4_5*v));
    Modes[i++] = -(conjugate(hHat_5_5_3) + conjugate(hHat_5_5_5)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -conjugate(hHat_5_4_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = -(conjugate(hHat_5_3_3) + conjugate(hHat_5_3_5)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -conjugate(hHat_5_2_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = -(conjugate(hHat_5_1_3) + conjugate(hHat_5_1_5)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_5_1_3 + hHat_5_1_5*pow(v, 2));
    Modes[i++] = hHat_5_2_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_5_3_3 + hHat_5_3_5*pow(v, 2));
    Modes[i++] = hHat_5_4_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_5_5_3 + hHat_5_5_5*pow(v, 2));
    Modes[i++] = conjugate(hHat_6_6_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = conjugate(hHat_6_5_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = conjugate(hHat_6_4_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = conjugate(hHat_6_3_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = conjugate(hHat_6_2_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = conjugate(hHat_6_1_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = 0;
    Modes[i++] = hHat_6_1_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = hHat_6_2_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = hHat_6_3_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = hHat_6_4_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = hHat_6_5_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = hHat_6_6_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = -conjugate(hHat_7_7_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = 0;
    Modes[i++] = -conjugate(hHat_7_5_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = 0;
    Modes[i++] = -conjugate(hHat_7_3_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = 0;
    Modes[i++] = -conjugate(hHat_7_1_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = 0;
    Modes[i++] = hHat_7_1_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = 0;
    Modes[i++] = hHat_7_3_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = 0;
    Modes[i++] = hHat_7_5_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = 0;
    Modes[i++] = hHat_7_7_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;

    return Modes;
  }

}; // class WaveformModes_2p5PN : public WaveformModes_Base


class WaveformModes_3p0PN : public WaveformModes_Base {
private:
  const Quaternions::Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const Quaternions::Quaternion S_chi1, S_chi2;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z;
  const double m, delta, nu;
  Quaternions::Quaternion R, nHat, lambdaHat, ellHat;
  double nHat_x, nHat_y, nHat_z, lambdaHat_x, lambdaHat_y, lambdaHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternions::Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, chi1_ell, chi1_n, chi1_lambda, chi2_ell, chi2_n, chi2_lambda,
         S_ell, S_n, S_lambda, Sigma_ell, Sigma_n, Sigma_lambda, S1_ell, S1_n, S1_lambda, S2_ell, S2_n, S2_lambda, logv;
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
                       hHat_spin_Symm_2_0_4, hHat_spin_Symm_3_3_4, hHat_spin_Symm_3_2_3, hHat_spin_Symm_3_1_4,
                       hHat_spin_Symm_4_3_4, hHat_spin_Symm_4_1_4, hHat_spin_Asymm_2_2_2, hHat_spin_Asymm_2_2_4,
                       hHat_spin_Asymm_2_1_3, hHat_spin_Asymm_2_1_4, hHat_spin_Asymm_2_0_2, hHat_spin_Asymm_2_0_4,
                       hHat_spin_Asymm_3_3_3, hHat_spin_Asymm_3_2_4, hHat_spin_Asymm_3_1_3, hHat_spin_Asymm_3_0_4,
                       hHat_spin_Asymm_4_4_4, hHat_spin_Asymm_4_2_4, hHat_spin_Asymm_4_0_4;

public:
  WaveformModes_3p0PN(const Quaternions::Quaternion xHat_i, const Quaternions::Quaternion yHat_i, const
                      Quaternions::Quaternion zHat_i, const double m1_i, const double m2_i, const double v_i, const
                      Quaternions::Quaternion S_chi1_i, const Quaternions::Quaternion S_chi2_i, const double
                      rfrak_chi1_x_i, const double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double
                      rfrak_chi2_y_i, const double rfrak_frame_x_i, const double rfrak_frame_y_i, const double
                      rfrak_frame_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), S_chi1(S_chi1_i), S_chi2(S_chi2_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_frame_x(rfrak_frame_x_i), rfrak_frame_y(rfrak_frame_y_i),
    rfrak_frame_z(rfrak_frame_z_i), m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_frame_x*xHat +
    rfrak_frame_y*yHat + rfrak_frame_z*zHat)), nHat(R*xHat*conjugate(R)), lambdaHat(R*yHat*conjugate(R)),
    ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]), lambdaHat_x(lambdaHat[1]),
    lambdaHat_y(lambdaHat[2]), lambdaHat_z(lambdaHat[3]), ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]),
    R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)), R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)),
    chiVec1(S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)),
    chiVec2(S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]),
    chi1_z(chiVec1[3]), chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1_ell(chi1_x*ellHat_x +
    chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z),
    chi1_lambda(chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z), chi2_ell(chi2_x*ellHat_x +
    chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z),
    chi2_lambda(chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z), S_ell(chi1_ell*pow(m1, 2) +
    chi2_ell*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)), S_lambda(chi1_lambda*pow(m1, 2) +
    chi2_lambda*pow(m2, 2)), Sigma_ell(m*(-chi1_ell*m1 + chi2_ell*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)),
    Sigma_lambda(m*(-chi1_lambda*m1 + chi2_lambda*m2)), S1_ell(chi1_ell*pow(m1, 2)), S1_n(chi1_n*pow(m1, 2)),
    S1_lambda(chi1_lambda*pow(m1, 2)), S2_ell(chi2_ell*pow(m2, 2)), S2_n(chi2_n*pow(m2, 2)),
    S2_lambda(chi2_lambda*pow(m2, 2)), logv(log(v)), rhOverM_coeff(6.34132367616962*nu*pow(v, 2)),
    hHat_2_0_0(-0.145802960879951), hHat_2_1_1(0.333333333333333*I*delta),
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
    1.0, 2) - 1.26086295798513), hHat_spin_Symm_2_2_3(0.333333333333333*(-11.0*S_ell - 5.0*Sigma_ell*delta)/pow(m, 2)),
    hHat_spin_Symm_2_2_4(0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
    15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(m, 4)*nu)), hHat_spin_Symm_2_1_2(0.5*I*Sigma_ell/pow(m, 2)),
    hHat_spin_Symm_2_1_4(0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(m, 2)),
    hHat_spin_Symm_2_0_4(0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(m, 4)*nu)),
    hHat_spin_Symm_3_3_4(0.0143763658196324*I*delta*(193.0*S_ell + 78732.0*Sigma_ell*nu)/pow(m, 2)),
    hHat_spin_Symm_3_2_3(0.563436169819011*(S_ell + Sigma_ell*delta)/pow(m, 2)),
    hHat_spin_Symm_3_1_4(0.0556794253984217*I*delta*(S_ell + 60.0*Sigma_ell*nu)/pow(m, 2)),
    hHat_spin_Symm_4_3_4(0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(m, 2)),
    hHat_spin_Symm_4_1_4(0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(m, 2)),
    hHat_spin_Asymm_2_2_2(0.5*(-Sigma_lambda - I*Sigma_n)/pow(m, 2)),
    hHat_spin_Asymm_2_2_4(0.00396825396825397*(-54180.0*Sigma_lambda*delta*nu + 5880.0*I*Sigma_n*nu +
    delta*(29.0*S_lambda + 546.0*I*S_n))/pow(m, 2)), hHat_spin_Asymm_2_1_3(0.166666666666667*(4.0*I*S_lambda + 25.0*S_n
    + 4.0*I*Sigma_lambda*delta + 13.0*Sigma_n*delta)/pow(m, 2)), hHat_spin_Asymm_2_1_4(0.5*(-3.0*S1_ell*S2_n -
    3.0*S1_n*S2_ell)/(pow(m, 4)*nu)), hHat_spin_Asymm_2_0_2(0.408248290463863*I*Sigma_n/pow(m, 2)),
    hHat_spin_Asymm_2_0_4(0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(m, 2)),
    hHat_spin_Asymm_3_3_3(0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(m, 2)),
    hHat_spin_Asymm_3_2_4(0.0117382535378961*(-50796.0*Sigma_lambda*delta*nu - 8580.0*I*Sigma_n*nu +
    delta*(71.0*S_lambda - 300.0*I*S_n))/pow(m, 2)), hHat_spin_Asymm_3_1_3(0.17817416127495*(I*S_lambda + S_n +
    delta*(I*Sigma_lambda + Sigma_n))/pow(m, 2)), hHat_spin_Asymm_3_0_4(-0.064293062484205*delta*(11.0*S_lambda +
    2268.0*Sigma_lambda*nu)/pow(m, 2)), hHat_spin_Asymm_4_4_4(0.950798536569581*(-3.0*Sigma_lambda*nu + Sigma_lambda -
    I*Sigma_n*(3.0*nu - 1.0) + delta*(S_lambda + I*S_n))/pow(m, 2)),
    hHat_spin_Asymm_4_2_4(0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) - 42.0*I*Sigma_n*nu +
    delta*(13.0*S_lambda - 14.0*I*S_n))/pow(m, 2)), hHat_spin_Asymm_4_0_4(0.00841793787126842*I*(S_n*delta -
    3.0*Sigma_n*nu + Sigma_n)/pow(m, 2))
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_x_k, const double chi1_y_k, const double chi1_z_k,
    const double chi2_x_k, const double chi2_y_k, const double chi2_z_k,
    const double ellHat_x_k, const double ellHat_y_k, const double ellHat_z_k)
  {
    v = v_k;
    chi1_x = chi1_x_k;
    chi1_y = chi1_y_k;
    chi1_z = chi1_z_k;
    chi2_x = chi2_x_k;
    chi2_y = chi2_y_k;
    chi2_z = chi2_z_k;
    ellHat_x = ellHat_x_k;
    ellHat_y = ellHat_y_k;
    ellHat_z = ellHat_z_k;

    R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat);
    nHat = R*xHat*conjugate(R);
    lambdaHat = R*yHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    lambdaHat_x = lambdaHat[1];
    lambdaHat_y = lambdaHat[2];
    lambdaHat_z = lambdaHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1);
    chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1_ell = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi1_lambda = chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z;
    chi2_ell = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    chi2_lambda = chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z;
    S_ell = chi1_ell*pow(m1, 2) + chi2_ell*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    S_lambda = chi1_lambda*pow(m1, 2) + chi2_lambda*pow(m2, 2);
    Sigma_ell = m*(-chi1_ell*m1 + chi2_ell*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    Sigma_lambda = m*(-chi1_lambda*m1 + chi2_lambda*m2);
    S1_ell = chi1_ell*pow(m1, 2);
    S1_n = chi1_n*pow(m1, 2);
    S1_lambda = chi1_lambda*pow(m1, 2);
    S2_ell = chi2_ell*pow(m2, 2);
    S2_n = chi2_n*pow(m2, 2);
    S2_lambda = chi2_lambda*pow(m2, 2);
    logv = log(v);
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    hHat_spin_Symm_2_2_3 = 0.333333333333333*(-11.0*S_ell - 5.0*Sigma_ell*delta)/pow(m, 2);
    hHat_spin_Symm_2_2_4 = 0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
      15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(m, 4)*nu);
    hHat_spin_Symm_2_1_2 = 0.5*I*Sigma_ell/pow(m, 2);
    hHat_spin_Symm_2_1_4 = 0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(m, 2);
    hHat_spin_Symm_2_0_4 = 0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(m, 4)*nu);
    hHat_spin_Symm_3_3_4 = 0.0143763658196324*I*delta*(193.0*S_ell + 78732.0*Sigma_ell*nu)/pow(m, 2);
    hHat_spin_Symm_3_2_3 = 0.563436169819011*(S_ell + Sigma_ell*delta)/pow(m, 2);
    hHat_spin_Symm_3_1_4 = 0.0556794253984217*I*delta*(S_ell + 60.0*Sigma_ell*nu)/pow(m, 2);
    hHat_spin_Symm_4_3_4 = 0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(m, 2);
    hHat_spin_Symm_4_1_4 = 0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(m, 2);
    hHat_spin_Asymm_2_2_2 = 0.5*(-Sigma_lambda - I*Sigma_n)/pow(m, 2);
    hHat_spin_Asymm_2_2_4 = 0.00396825396825397*(-54180.0*Sigma_lambda*delta*nu + 5880.0*I*Sigma_n*nu +
      delta*(29.0*S_lambda + 546.0*I*S_n))/pow(m, 2);
    hHat_spin_Asymm_2_1_3 = 0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
      13.0*Sigma_n*delta)/pow(m, 2);
    hHat_spin_Asymm_2_1_4 = 0.5*(-3.0*S1_ell*S2_n - 3.0*S1_n*S2_ell)/(pow(m, 4)*nu);
    hHat_spin_Asymm_2_0_2 = 0.408248290463863*I*Sigma_n/pow(m, 2);
    hHat_spin_Asymm_2_0_4 = 0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(m, 2);
    hHat_spin_Asymm_3_3_3 = 0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(m, 2);
    hHat_spin_Asymm_3_2_4 = 0.0117382535378961*(-50796.0*Sigma_lambda*delta*nu - 8580.0*I*Sigma_n*nu +
      delta*(71.0*S_lambda - 300.0*I*S_n))/pow(m, 2);
    hHat_spin_Asymm_3_1_3 = 0.17817416127495*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/pow(m, 2);
    hHat_spin_Asymm_3_0_4 = -0.064293062484205*delta*(11.0*S_lambda + 2268.0*Sigma_lambda*nu)/pow(m, 2);
    hHat_spin_Asymm_4_4_4 = 0.950798536569581*(-3.0*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3.0*nu - 1.0) +
      delta*(S_lambda + I*S_n))/pow(m, 2);
    hHat_spin_Asymm_4_2_4 = 0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) - 42.0*I*Sigma_n*nu +
      delta*(13.0*S_lambda - 14.0*I*S_n))/pow(m, 2);
    hHat_spin_Asymm_4_0_4 = 0.00841793787126842*I*(S_n*delta - 3.0*Sigma_n*nu + Sigma_n)/pow(m, 2);

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = ((((((conjugate(hHat_2_2_6) + conjugate(hHat_2_2_lnv_6)*conjugate(logv))*conjugate(v) +
      conjugate(hHat_2_2_5))*conjugate(v) + conjugate(hHat_2_2_4) + conjugate(hHat_spin_Symm_2_2_4))*conjugate(v) +
      conjugate(hHat_2_2_3) + conjugate(hHat_spin_Symm_2_2_3))*conjugate(v) + conjugate(hHat_2_2_2))*pow(conjugate(v),
      2) + conjugate(hHat_2_2_0))*conjugate(rhOverM_coeff) - (conjugate(hHat_spin_Asymm_2_2_2) +
      conjugate(hHat_spin_Asymm_2_2_4)*pow(conjugate(v), 2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2);
    Modes[i++] = (((((conjugate(hHat_2_1_5) + conjugate(hHat_2_1_6)*conjugate(v))*conjugate(v) + conjugate(hHat_2_1_4) +
      conjugate(hHat_spin_Symm_2_1_4))*conjugate(v) + conjugate(hHat_2_1_3))*conjugate(v) +
      conjugate(hHat_spin_Symm_2_1_2))*conjugate(v) + conjugate(hHat_2_1_1))*conjugate(rhOverM_coeff)*conjugate(v) -
      (conjugate(hHat_spin_Asymm_2_1_3) +
      conjugate(hHat_spin_Asymm_2_1_4)*conjugate(v))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*pow(v, 2)) +
      rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_4*pow(v, 4));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v) + rhOverM_coeff*v*(hHat_2_1_1
      + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_spin_Symm_2_1_4 + v*(hHat_2_1_5 +
      hHat_2_1_6*v)))));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*pow(v, 2)) +
      rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 +
      hHat_spin_Symm_2_2_4 + v*(hHat_2_2_5 + v*(hHat_2_2_6 + hHat_2_2_lnv_6*logv))))));
    Modes[i++] = -((((conjugate(hHat_3_3_5) + conjugate(hHat_3_3_6)*conjugate(v))*conjugate(v) + conjugate(hHat_3_3_4) +
      conjugate(hHat_spin_Symm_3_3_4))*conjugate(v) + conjugate(hHat_3_3_3))*pow(conjugate(v), 2) +
      conjugate(hHat_3_3_1))*conjugate(rhOverM_coeff)*conjugate(v) +
      conjugate(hHat_spin_Asymm_3_3_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -((((conjugate(hHat_3_2_5) + conjugate(hHat_3_2_6)*conjugate(v))*conjugate(v) +
      conjugate(hHat_3_2_4))*conjugate(v) + conjugate(hHat_spin_Symm_3_2_3))*conjugate(v) +
      conjugate(hHat_3_2_2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2) +
      conjugate(hHat_spin_Asymm_3_2_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = -((((conjugate(hHat_3_1_5) + conjugate(hHat_3_1_6)*conjugate(v))*conjugate(v) + conjugate(hHat_3_1_4) +
      conjugate(hHat_spin_Symm_3_1_4))*conjugate(v) + conjugate(hHat_3_1_3))*pow(conjugate(v), 2) +
      conjugate(hHat_3_1_1))*conjugate(rhOverM_coeff)*conjugate(v) +
      conjugate(hHat_spin_Asymm_3_1_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = hHat_3_0_5*rhOverM_coeff*pow(v, 5) + hHat_spin_Asymm_3_0_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3) + rhOverM_coeff*v*(hHat_3_1_1 + pow(v, 2)*(hHat_3_1_3 +
      v*(hHat_3_1_4 + hHat_spin_Symm_3_1_4 + v*(hHat_3_1_5 + hHat_3_1_6*v))));
    Modes[i++] = hHat_spin_Asymm_3_2_4*rhOverM_coeff*pow(v, 4) + rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 +
      v*(hHat_spin_Symm_3_2_3 + v*(hHat_3_2_4 + v*(hHat_3_2_5 + hHat_3_2_6*v))));
    Modes[i++] = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3) + rhOverM_coeff*v*(hHat_3_3_1 + pow(v, 2)*(hHat_3_3_3 +
      v*(hHat_3_3_4 + hHat_spin_Symm_3_3_4 + v*(hHat_3_3_5 + hHat_3_3_6*v))));
    Modes[i++] = (((conjugate(hHat_4_4_5) + conjugate(hHat_4_4_6)*conjugate(v))*conjugate(v) +
      conjugate(hHat_4_4_4))*pow(conjugate(v), 2) + conjugate(hHat_4_4_2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2)
      - conjugate(hHat_spin_Asymm_4_4_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = (((conjugate(hHat_4_3_5) + conjugate(hHat_4_3_6)*conjugate(v))*conjugate(v) +
      conjugate(hHat_spin_Symm_4_3_4))*conjugate(v) + conjugate(hHat_4_3_3))*conjugate(rhOverM_coeff)*pow(conjugate(v),
      3);
    Modes[i++] = (((conjugate(hHat_4_2_5) + conjugate(hHat_4_2_6)*conjugate(v))*conjugate(v) +
      conjugate(hHat_4_2_4))*pow(conjugate(v), 2) + conjugate(hHat_4_2_2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2)
      - conjugate(hHat_spin_Asymm_4_2_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = (((conjugate(hHat_4_1_5) + conjugate(hHat_4_1_6)*conjugate(v))*conjugate(v) +
      conjugate(hHat_spin_Symm_4_1_4))*conjugate(v) + conjugate(hHat_4_1_3))*conjugate(rhOverM_coeff)*pow(conjugate(v),
      3);
    Modes[i++] = hHat_4_0_0*rhOverM_coeff + hHat_spin_Asymm_4_0_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_4_1_3 + v*(hHat_spin_Symm_4_1_4 + v*(hHat_4_1_5 + hHat_4_1_6*v)));
    Modes[i++] = hHat_spin_Asymm_4_2_4*rhOverM_coeff*pow(v, 4) + rhOverM_coeff*pow(v, 2)*(hHat_4_2_2 + pow(v,
      2)*(hHat_4_2_4 + v*(hHat_4_2_5 + hHat_4_2_6*v)));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_4_3_3 + v*(hHat_spin_Symm_4_3_4 + v*(hHat_4_3_5 + hHat_4_3_6*v)));
    Modes[i++] = hHat_spin_Asymm_4_4_4*rhOverM_coeff*pow(v, 4) + rhOverM_coeff*pow(v, 2)*(hHat_4_4_2 + pow(v,
      2)*(hHat_4_4_4 + v*(hHat_4_4_5 + hHat_4_4_6*v)));
    Modes[i++] = -((conjugate(hHat_5_5_5) + conjugate(hHat_5_5_6)*conjugate(v))*pow(conjugate(v), 2) +
      conjugate(hHat_5_5_3))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -(conjugate(hHat_5_4_4) + conjugate(hHat_5_4_6)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = -((conjugate(hHat_5_3_5) + conjugate(hHat_5_3_6)*conjugate(v))*pow(conjugate(v), 2) +
      conjugate(hHat_5_3_3))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -(conjugate(hHat_5_2_4) + conjugate(hHat_5_2_6)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = -((conjugate(hHat_5_1_5) + conjugate(hHat_5_1_6)*conjugate(v))*pow(conjugate(v), 2) +
      conjugate(hHat_5_1_3))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_5_1_3 + pow(v, 2)*(hHat_5_1_5 + hHat_5_1_6*v));
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(hHat_5_2_4 + hHat_5_2_6*pow(v, 2));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_5_3_3 + pow(v, 2)*(hHat_5_3_5 + hHat_5_3_6*v));
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(hHat_5_4_4 + hHat_5_4_6*pow(v, 2));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_5_5_3 + pow(v, 2)*(hHat_5_5_5 + hHat_5_5_6*v));
    Modes[i++] = (conjugate(hHat_6_6_4) + conjugate(hHat_6_6_6)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = conjugate(hHat_6_5_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = (conjugate(hHat_6_4_4) + conjugate(hHat_6_4_6)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = conjugate(hHat_6_3_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = (conjugate(hHat_6_2_4) + conjugate(hHat_6_2_6)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = conjugate(hHat_6_1_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = 0;
    Modes[i++] = hHat_6_1_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(hHat_6_2_4 + hHat_6_2_6*pow(v, 2));
    Modes[i++] = hHat_6_3_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(hHat_6_4_4 + hHat_6_4_6*pow(v, 2));
    Modes[i++] = hHat_6_5_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(hHat_6_6_4 + hHat_6_6_6*pow(v, 2));
    Modes[i++] = -conjugate(hHat_7_7_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = -conjugate(hHat_7_6_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = -conjugate(hHat_7_5_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = -conjugate(hHat_7_4_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = -conjugate(hHat_7_3_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = -conjugate(hHat_7_2_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = -conjugate(hHat_7_1_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = 0;
    Modes[i++] = hHat_7_1_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = hHat_7_2_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = hHat_7_3_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = hHat_7_4_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = hHat_7_5_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = hHat_7_6_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = hHat_7_7_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = conjugate(hHat_8_8_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = 0;
    Modes[i++] = conjugate(hHat_8_6_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = 0;
    Modes[i++] = conjugate(hHat_8_4_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = 0;
    Modes[i++] = conjugate(hHat_8_2_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = hHat_8_2_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = 0;
    Modes[i++] = hHat_8_4_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = 0;
    Modes[i++] = hHat_8_6_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = 0;
    Modes[i++] = hHat_8_8_6*rhOverM_coeff*pow(v, 6);

    return Modes;
  }

}; // class WaveformModes_3p0PN : public WaveformModes_Base


class WaveformModes_3p5PN : public WaveformModes_Base {
private:
  const Quaternions::Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const Quaternions::Quaternion S_chi1, S_chi2;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z;
  const double m, delta, nu;
  Quaternions::Quaternion R, nHat, lambdaHat, ellHat;
  double nHat_x, nHat_y, nHat_z, lambdaHat_x, lambdaHat_y, lambdaHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternions::Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, chi1_ell, chi1_n, chi1_lambda, chi2_ell, chi2_n, chi2_lambda,
         S_ell, S_n, S_lambda, Sigma_ell, Sigma_n, Sigma_lambda, S1_ell, S1_n, S1_lambda, S2_ell, S2_n, S2_lambda, logv;
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
                       hHat_spin_Symm_2_0_4, hHat_spin_Symm_3_3_4, hHat_spin_Symm_3_2_3, hHat_spin_Symm_3_1_4,
                       hHat_spin_Symm_4_3_4, hHat_spin_Symm_4_1_4, hHat_spin_Asymm_2_2_2, hHat_spin_Asymm_2_2_4,
                       hHat_spin_Asymm_2_1_3, hHat_spin_Asymm_2_1_4, hHat_spin_Asymm_2_0_2, hHat_spin_Asymm_2_0_4,
                       hHat_spin_Asymm_3_3_3, hHat_spin_Asymm_3_2_4, hHat_spin_Asymm_3_1_3, hHat_spin_Asymm_3_0_4,
                       hHat_spin_Asymm_4_4_4, hHat_spin_Asymm_4_2_4, hHat_spin_Asymm_4_0_4;

public:
  WaveformModes_3p5PN(const Quaternions::Quaternion xHat_i, const Quaternions::Quaternion yHat_i, const
                      Quaternions::Quaternion zHat_i, const double m1_i, const double m2_i, const double v_i, const
                      Quaternions::Quaternion S_chi1_i, const Quaternions::Quaternion S_chi2_i, const double
                      rfrak_chi1_x_i, const double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double
                      rfrak_chi2_y_i, const double rfrak_frame_x_i, const double rfrak_frame_y_i, const double
                      rfrak_frame_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), S_chi1(S_chi1_i), S_chi2(S_chi2_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_frame_x(rfrak_frame_x_i), rfrak_frame_y(rfrak_frame_y_i),
    rfrak_frame_z(rfrak_frame_z_i), m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_frame_x*xHat +
    rfrak_frame_y*yHat + rfrak_frame_z*zHat)), nHat(R*xHat*conjugate(R)), lambdaHat(R*yHat*conjugate(R)),
    ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]), lambdaHat_x(lambdaHat[1]),
    lambdaHat_y(lambdaHat[2]), lambdaHat_z(lambdaHat[3]), ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]),
    R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)), R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)),
    chiVec1(S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)),
    chiVec2(S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]),
    chi1_z(chiVec1[3]), chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1_ell(chi1_x*ellHat_x +
    chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z),
    chi1_lambda(chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z), chi2_ell(chi2_x*ellHat_x +
    chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z),
    chi2_lambda(chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z), S_ell(chi1_ell*pow(m1, 2) +
    chi2_ell*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)), S_lambda(chi1_lambda*pow(m1, 2) +
    chi2_lambda*pow(m2, 2)), Sigma_ell(m*(-chi1_ell*m1 + chi2_ell*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)),
    Sigma_lambda(m*(-chi1_lambda*m1 + chi2_lambda*m2)), S1_ell(chi1_ell*pow(m1, 2)), S1_n(chi1_n*pow(m1, 2)),
    S1_lambda(chi1_lambda*pow(m1, 2)), S2_ell(chi2_ell*pow(m2, 2)), S2_n(chi2_n*pow(m2, 2)),
    S2_lambda(chi2_lambda*pow(m2, 2)), logv(log(v)), rhOverM_coeff(6.34132367616962*nu*pow(v, 2)),
    hHat_2_0_0(-0.145802960879951), hHat_2_1_1(0.333333333333333*I*delta),
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
    1.0, 2) - 1.26086295798513), hHat_spin_Symm_2_2_3(0.333333333333333*(-11.0*S_ell - 5.0*Sigma_ell*delta)/pow(m, 2)),
    hHat_spin_Symm_2_2_4(0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
    15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(m, 4)*nu)), hHat_spin_Symm_2_1_2(0.5*I*Sigma_ell/pow(m, 2)),
    hHat_spin_Symm_2_1_4(0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(m, 2)),
    hHat_spin_Symm_2_0_4(0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(m, 4)*nu)),
    hHat_spin_Symm_3_3_4(0.0143763658196324*I*delta*(193.0*S_ell + 78732.0*Sigma_ell*nu)/pow(m, 2)),
    hHat_spin_Symm_3_2_3(0.563436169819011*(S_ell + Sigma_ell*delta)/pow(m, 2)),
    hHat_spin_Symm_3_1_4(0.0556794253984217*I*delta*(S_ell + 60.0*Sigma_ell*nu)/pow(m, 2)),
    hHat_spin_Symm_4_3_4(0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(m, 2)),
    hHat_spin_Symm_4_1_4(0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(m, 2)),
    hHat_spin_Asymm_2_2_2(0.5*(-Sigma_lambda - I*Sigma_n)/pow(m, 2)),
    hHat_spin_Asymm_2_2_4(0.00396825396825397*(-54180.0*Sigma_lambda*delta*nu + 5880.0*I*Sigma_n*nu +
    delta*(29.0*S_lambda + 546.0*I*S_n))/pow(m, 2)), hHat_spin_Asymm_2_1_3(0.166666666666667*(4.0*I*S_lambda + 25.0*S_n
    + 4.0*I*Sigma_lambda*delta + 13.0*Sigma_n*delta)/pow(m, 2)), hHat_spin_Asymm_2_1_4(0.5*(-3.0*S1_ell*S2_n -
    3.0*S1_n*S2_ell)/(pow(m, 4)*nu)), hHat_spin_Asymm_2_0_2(0.408248290463863*I*Sigma_n/pow(m, 2)),
    hHat_spin_Asymm_2_0_4(0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(m, 2)),
    hHat_spin_Asymm_3_3_3(0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(m, 2)),
    hHat_spin_Asymm_3_2_4(0.0117382535378961*(-50796.0*Sigma_lambda*delta*nu - 8580.0*I*Sigma_n*nu +
    delta*(71.0*S_lambda - 300.0*I*S_n))/pow(m, 2)), hHat_spin_Asymm_3_1_3(0.17817416127495*(I*S_lambda + S_n +
    delta*(I*Sigma_lambda + Sigma_n))/pow(m, 2)), hHat_spin_Asymm_3_0_4(-0.064293062484205*delta*(11.0*S_lambda +
    2268.0*Sigma_lambda*nu)/pow(m, 2)), hHat_spin_Asymm_4_4_4(0.950798536569581*(-3.0*Sigma_lambda*nu + Sigma_lambda -
    I*Sigma_n*(3.0*nu - 1.0) + delta*(S_lambda + I*S_n))/pow(m, 2)),
    hHat_spin_Asymm_4_2_4(0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) - 42.0*I*Sigma_n*nu +
    delta*(13.0*S_lambda - 14.0*I*S_n))/pow(m, 2)), hHat_spin_Asymm_4_0_4(0.00841793787126842*I*(S_n*delta -
    3.0*Sigma_n*nu + Sigma_n)/pow(m, 2))
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_x_k, const double chi1_y_k, const double chi1_z_k,
    const double chi2_x_k, const double chi2_y_k, const double chi2_z_k,
    const double ellHat_x_k, const double ellHat_y_k, const double ellHat_z_k)
  {
    v = v_k;
    chi1_x = chi1_x_k;
    chi1_y = chi1_y_k;
    chi1_z = chi1_z_k;
    chi2_x = chi2_x_k;
    chi2_y = chi2_y_k;
    chi2_z = chi2_z_k;
    ellHat_x = ellHat_x_k;
    ellHat_y = ellHat_y_k;
    ellHat_z = ellHat_z_k;

    R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat);
    nHat = R*xHat*conjugate(R);
    lambdaHat = R*yHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    lambdaHat_x = lambdaHat[1];
    lambdaHat_y = lambdaHat[2];
    lambdaHat_z = lambdaHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1);
    chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1_ell = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi1_lambda = chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z;
    chi2_ell = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    chi2_lambda = chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z;
    S_ell = chi1_ell*pow(m1, 2) + chi2_ell*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    S_lambda = chi1_lambda*pow(m1, 2) + chi2_lambda*pow(m2, 2);
    Sigma_ell = m*(-chi1_ell*m1 + chi2_ell*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    Sigma_lambda = m*(-chi1_lambda*m1 + chi2_lambda*m2);
    S1_ell = chi1_ell*pow(m1, 2);
    S1_n = chi1_n*pow(m1, 2);
    S1_lambda = chi1_lambda*pow(m1, 2);
    S2_ell = chi2_ell*pow(m2, 2);
    S2_n = chi2_n*pow(m2, 2);
    S2_lambda = chi2_lambda*pow(m2, 2);
    logv = log(v);
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    hHat_spin_Symm_2_2_3 = 0.333333333333333*(-11.0*S_ell - 5.0*Sigma_ell*delta)/pow(m, 2);
    hHat_spin_Symm_2_2_4 = 0.166666666666667*(12.0*S1_ell*S2_ell + 10.0*S1_lambda*S2_lambda - 15.0*I*S1_lambda*S2_n -
      15.0*I*S1_n*S2_lambda - 22.0*S1_n*S2_n)/(pow(m, 4)*nu);
    hHat_spin_Symm_2_1_2 = 0.5*I*Sigma_ell/pow(m, 2);
    hHat_spin_Symm_2_1_4 = 0.0238095238095238*I*(-86.0*S_ell*delta + Sigma_ell*(139.0*nu - 79.0))/pow(m, 2);
    hHat_spin_Symm_2_0_4 = 0.816496580927726*(-S1_lambda*S2_lambda + S1_n*S2_n)/(pow(m, 4)*nu);
    hHat_spin_Symm_3_3_4 = 0.0143763658196324*I*delta*(193.0*S_ell + 78732.0*Sigma_ell*nu)/pow(m, 2);
    hHat_spin_Symm_3_2_3 = 0.563436169819011*(S_ell + Sigma_ell*delta)/pow(m, 2);
    hHat_spin_Symm_3_1_4 = 0.0556794253984217*I*delta*(S_ell + 60.0*Sigma_ell*nu)/pow(m, 2);
    hHat_spin_Symm_4_3_4 = 0.672316092750596*I*(-S_ell*delta + 3.0*Sigma_ell*nu - Sigma_ell)/pow(m, 2);
    hHat_spin_Symm_4_1_4 = 0.00941154065526303*I*(S_ell*delta - 3.0*Sigma_ell*nu + Sigma_ell)/pow(m, 2);
    hHat_spin_Asymm_2_2_2 = 0.5*(-Sigma_lambda - I*Sigma_n)/pow(m, 2);
    hHat_spin_Asymm_2_2_4 = 0.00396825396825397*(-54180.0*Sigma_lambda*delta*nu + 5880.0*I*Sigma_n*nu +
      delta*(29.0*S_lambda + 546.0*I*S_n))/pow(m, 2);
    hHat_spin_Asymm_2_1_3 = 0.166666666666667*(4.0*I*S_lambda + 25.0*S_n + 4.0*I*Sigma_lambda*delta +
      13.0*Sigma_n*delta)/pow(m, 2);
    hHat_spin_Asymm_2_1_4 = 0.5*(-3.0*S1_ell*S2_n - 3.0*S1_n*S2_ell)/(pow(m, 4)*nu);
    hHat_spin_Asymm_2_0_2 = 0.408248290463863*I*Sigma_n/pow(m, 2);
    hHat_spin_Asymm_2_0_4 = 0.0194403947839935*I*(255.0*S_n*delta - Sigma_n*(506.0*nu - 45.0))/pow(m, 2);
    hHat_spin_Asymm_3_3_3 = 0.690065559342354*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/pow(m, 2);
    hHat_spin_Asymm_3_2_4 = 0.0117382535378961*(-50796.0*Sigma_lambda*delta*nu - 8580.0*I*Sigma_n*nu +
      delta*(71.0*S_lambda - 300.0*I*S_n))/pow(m, 2);
    hHat_spin_Asymm_3_1_3 = 0.17817416127495*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/pow(m, 2);
    hHat_spin_Asymm_3_0_4 = -0.064293062484205*delta*(11.0*S_lambda + 2268.0*Sigma_lambda*nu)/pow(m, 2);
    hHat_spin_Asymm_4_4_4 = 0.950798536569581*(-3.0*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3.0*nu - 1.0) +
      delta*(S_lambda + I*S_n))/pow(m, 2);
    hHat_spin_Asymm_4_2_4 = 0.0133099284374987*(-13.0*Sigma_lambda*(3.0*nu - 1.0) - 42.0*I*Sigma_n*nu +
      delta*(13.0*S_lambda - 14.0*I*S_n))/pow(m, 2);
    hHat_spin_Asymm_4_0_4 = 0.00841793787126842*I*(S_n*delta - 3.0*Sigma_n*nu + Sigma_n)/pow(m, 2);

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = ((((((conjugate(hHat_2_2_6) + conjugate(hHat_2_2_7)*conjugate(v) +
      conjugate(hHat_2_2_lnv_6)*conjugate(logv))*conjugate(v) + conjugate(hHat_2_2_5))*conjugate(v) +
      conjugate(hHat_2_2_4) + conjugate(hHat_spin_Symm_2_2_4))*conjugate(v) + conjugate(hHat_2_2_3) +
      conjugate(hHat_spin_Symm_2_2_3))*conjugate(v) + conjugate(hHat_2_2_2))*pow(conjugate(v), 2) +
      conjugate(hHat_2_2_0))*conjugate(rhOverM_coeff) - (conjugate(hHat_spin_Asymm_2_2_2) +
      conjugate(hHat_spin_Asymm_2_2_4)*pow(conjugate(v), 2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2);
    Modes[i++] = (((((conjugate(hHat_2_1_5) + conjugate(hHat_2_1_6)*conjugate(v))*conjugate(v) + conjugate(hHat_2_1_4) +
      conjugate(hHat_spin_Symm_2_1_4))*conjugate(v) + conjugate(hHat_2_1_3))*conjugate(v) +
      conjugate(hHat_spin_Symm_2_1_2))*conjugate(v) + conjugate(hHat_2_1_1))*conjugate(rhOverM_coeff)*conjugate(v) -
      (conjugate(hHat_spin_Asymm_2_1_3) +
      conjugate(hHat_spin_Asymm_2_1_4)*conjugate(v))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*pow(v, 2)) +
      rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_4*pow(v, 4));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v) + rhOverM_coeff*v*(hHat_2_1_1
      + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_spin_Symm_2_1_4 + v*(hHat_2_1_5 +
      hHat_2_1_6*v)))));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*pow(v, 2)) +
      rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 +
      hHat_spin_Symm_2_2_4 + v*(hHat_2_2_5 + v*(hHat_2_2_6 + hHat_2_2_7*v + hHat_2_2_lnv_6*logv))))));
    Modes[i++] = -((((conjugate(hHat_3_3_5) + conjugate(hHat_3_3_6)*conjugate(v))*conjugate(v) + conjugate(hHat_3_3_4) +
      conjugate(hHat_spin_Symm_3_3_4))*conjugate(v) + conjugate(hHat_3_3_3))*pow(conjugate(v), 2) +
      conjugate(hHat_3_3_1))*conjugate(rhOverM_coeff)*conjugate(v) +
      conjugate(hHat_spin_Asymm_3_3_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -((((conjugate(hHat_3_2_5) + conjugate(hHat_3_2_6)*conjugate(v))*conjugate(v) +
      conjugate(hHat_3_2_4))*conjugate(v) + conjugate(hHat_spin_Symm_3_2_3))*conjugate(v) +
      conjugate(hHat_3_2_2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2) +
      conjugate(hHat_spin_Asymm_3_2_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = -((((conjugate(hHat_3_1_5) + conjugate(hHat_3_1_6)*conjugate(v))*conjugate(v) + conjugate(hHat_3_1_4) +
      conjugate(hHat_spin_Symm_3_1_4))*conjugate(v) + conjugate(hHat_3_1_3))*pow(conjugate(v), 2) +
      conjugate(hHat_3_1_1))*conjugate(rhOverM_coeff)*conjugate(v) +
      conjugate(hHat_spin_Asymm_3_1_3)*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = hHat_3_0_5*rhOverM_coeff*pow(v, 5) + hHat_spin_Asymm_3_0_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3) + rhOverM_coeff*v*(hHat_3_1_1 + pow(v, 2)*(hHat_3_1_3 +
      v*(hHat_3_1_4 + hHat_spin_Symm_3_1_4 + v*(hHat_3_1_5 + hHat_3_1_6*v))));
    Modes[i++] = hHat_spin_Asymm_3_2_4*rhOverM_coeff*pow(v, 4) + rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 +
      v*(hHat_spin_Symm_3_2_3 + v*(hHat_3_2_4 + v*(hHat_3_2_5 + hHat_3_2_6*v))));
    Modes[i++] = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3) + rhOverM_coeff*v*(hHat_3_3_1 + pow(v, 2)*(hHat_3_3_3 +
      v*(hHat_3_3_4 + hHat_spin_Symm_3_3_4 + v*(hHat_3_3_5 + hHat_3_3_6*v))));
    Modes[i++] = (((conjugate(hHat_4_4_5) + conjugate(hHat_4_4_6)*conjugate(v))*conjugate(v) +
      conjugate(hHat_4_4_4))*pow(conjugate(v), 2) + conjugate(hHat_4_4_2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2)
      - conjugate(hHat_spin_Asymm_4_4_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = (((conjugate(hHat_4_3_5) + conjugate(hHat_4_3_6)*conjugate(v))*conjugate(v) +
      conjugate(hHat_spin_Symm_4_3_4))*conjugate(v) + conjugate(hHat_4_3_3))*conjugate(rhOverM_coeff)*pow(conjugate(v),
      3);
    Modes[i++] = (((conjugate(hHat_4_2_5) + conjugate(hHat_4_2_6)*conjugate(v))*conjugate(v) +
      conjugate(hHat_4_2_4))*pow(conjugate(v), 2) + conjugate(hHat_4_2_2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 2)
      - conjugate(hHat_spin_Asymm_4_2_4)*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = (((conjugate(hHat_4_1_5) + conjugate(hHat_4_1_6)*conjugate(v))*conjugate(v) +
      conjugate(hHat_spin_Symm_4_1_4))*conjugate(v) + conjugate(hHat_4_1_3))*conjugate(rhOverM_coeff)*pow(conjugate(v),
      3);
    Modes[i++] = hHat_4_0_0*rhOverM_coeff + hHat_spin_Asymm_4_0_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_4_1_3 + v*(hHat_spin_Symm_4_1_4 + v*(hHat_4_1_5 + hHat_4_1_6*v)));
    Modes[i++] = hHat_spin_Asymm_4_2_4*rhOverM_coeff*pow(v, 4) + rhOverM_coeff*pow(v, 2)*(hHat_4_2_2 + pow(v,
      2)*(hHat_4_2_4 + v*(hHat_4_2_5 + hHat_4_2_6*v)));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_4_3_3 + v*(hHat_spin_Symm_4_3_4 + v*(hHat_4_3_5 + hHat_4_3_6*v)));
    Modes[i++] = hHat_spin_Asymm_4_4_4*rhOverM_coeff*pow(v, 4) + rhOverM_coeff*pow(v, 2)*(hHat_4_4_2 + pow(v,
      2)*(hHat_4_4_4 + v*(hHat_4_4_5 + hHat_4_4_6*v)));
    Modes[i++] = -((conjugate(hHat_5_5_5) + conjugate(hHat_5_5_6)*conjugate(v))*pow(conjugate(v), 2) +
      conjugate(hHat_5_5_3))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -(conjugate(hHat_5_4_4) + conjugate(hHat_5_4_6)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = -((conjugate(hHat_5_3_5) + conjugate(hHat_5_3_6)*conjugate(v))*pow(conjugate(v), 2) +
      conjugate(hHat_5_3_3))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = -(conjugate(hHat_5_2_4) + conjugate(hHat_5_2_6)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = -((conjugate(hHat_5_1_5) + conjugate(hHat_5_1_6)*conjugate(v))*pow(conjugate(v), 2) +
      conjugate(hHat_5_1_3))*conjugate(rhOverM_coeff)*pow(conjugate(v), 3);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_5_1_3 + pow(v, 2)*(hHat_5_1_5 + hHat_5_1_6*v));
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(hHat_5_2_4 + hHat_5_2_6*pow(v, 2));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_5_3_3 + pow(v, 2)*(hHat_5_3_5 + hHat_5_3_6*v));
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(hHat_5_4_4 + hHat_5_4_6*pow(v, 2));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(hHat_5_5_3 + pow(v, 2)*(hHat_5_5_5 + hHat_5_5_6*v));
    Modes[i++] = (conjugate(hHat_6_6_4) + conjugate(hHat_6_6_6)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = conjugate(hHat_6_5_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = (conjugate(hHat_6_4_4) + conjugate(hHat_6_4_6)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = conjugate(hHat_6_3_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = (conjugate(hHat_6_2_4) + conjugate(hHat_6_2_6)*pow(conjugate(v),
      2))*conjugate(rhOverM_coeff)*pow(conjugate(v), 4);
    Modes[i++] = conjugate(hHat_6_1_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = 0;
    Modes[i++] = hHat_6_1_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(hHat_6_2_4 + hHat_6_2_6*pow(v, 2));
    Modes[i++] = hHat_6_3_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(hHat_6_4_4 + hHat_6_4_6*pow(v, 2));
    Modes[i++] = hHat_6_5_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(hHat_6_6_4 + hHat_6_6_6*pow(v, 2));
    Modes[i++] = -conjugate(hHat_7_7_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = -conjugate(hHat_7_6_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = -conjugate(hHat_7_5_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = -conjugate(hHat_7_4_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = -conjugate(hHat_7_3_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = -conjugate(hHat_7_2_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = -conjugate(hHat_7_1_5)*conjugate(rhOverM_coeff)*pow(conjugate(v), 5);
    Modes[i++] = 0;
    Modes[i++] = hHat_7_1_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = hHat_7_2_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = hHat_7_3_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = hHat_7_4_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = hHat_7_5_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = hHat_7_6_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = hHat_7_7_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = conjugate(hHat_8_8_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = 0;
    Modes[i++] = conjugate(hHat_8_6_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = 0;
    Modes[i++] = conjugate(hHat_8_4_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = 0;
    Modes[i++] = conjugate(hHat_8_2_6)*conjugate(rhOverM_coeff)*pow(conjugate(v), 6);
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = hHat_8_2_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = 0;
    Modes[i++] = hHat_8_4_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = 0;
    Modes[i++] = hHat_8_6_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = 0;
    Modes[i++] = hHat_8_8_6*rhOverM_coeff*pow(v, 6);

    return Modes;
  }

}; // class WaveformModes_3p5PN : public WaveformModes_Base
