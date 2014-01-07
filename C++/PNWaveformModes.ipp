// File produced automatically by WaveformModeCodeGen.ipynb

class WaveformModes {
private:
  virtual std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_x_k, const double chi1_y_k, const double chi1_z_k,
    const double chi2_x_k, const double chi2_y_k, const double chi2_z_k,
    const double ellHat_x_k, const double ellHat_y_k, const double ellHat_z_k);
public:
  virtual std::vector<std::complex<double> > operator()(
    const double v_k, const std::vector<double> chi1, const std::vector<double> chi2)
  {
    return this->operator()(v_k, chi1[0], chi1[1], chi1[2], chi2[0], chi2[1], chi2[2], 0.0, 0.0, 1.0);
  }
};

const std::complex<double> I(0,1.0);

class WaveformModes_0PN : public WaveformModes {
private:
  const double m1, m2;
  double v;
  const double m, nu;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> rhOverM_2_0_0, rhOverM_2_2_0, rhOverM_4_0_0;

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

    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = conjugate(rhOverM_2_2_0)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_2_0_0*rhOverM_coeff;
    Modes[i++] = 0;
    Modes[i++] = rhOverM_2_2_0*rhOverM_coeff;
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
    Modes[i++] = rhOverM_4_0_0*rhOverM_coeff;
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

public:
  WaveformModes_0PN(const double m1_i, const double m2_i, const double v_i) :
    m1(m1_i), m2(m2_i), v(v_i), m(m1 + m2), nu(m1*m2/pow(m, 2)), rhOverM_coeff(6.34132367616962*nu*pow(v, 2)),
    rhOverM_2_0_0(-0.145802960879951), rhOverM_2_2_0(1.00000000000000), rhOverM_4_0_0(-0.00140298964521140)
  { }

}; // class WaveformModes_0PN : public WaveformModes


class WaveformModes_0p50PN : public WaveformModes {
private:
  const double m1, m2;
  double v;
  const double m, delta, nu;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> rhOverM_2_0_0, rhOverM_2_1_1, rhOverM_2_2_0, rhOverM_3_1_1, rhOverM_3_3_1, rhOverM_4_0_0;

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

    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = conjugate(rhOverM_2_2_0)*conjugate(rhOverM_coeff);
    Modes[i++] = v*conjugate(rhOverM_2_1_1)*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_2_0_0*rhOverM_coeff;
    Modes[i++] = rhOverM_2_1_1*rhOverM_coeff*v;
    Modes[i++] = rhOverM_2_2_0*rhOverM_coeff;
    Modes[i++] = -v*conjugate(rhOverM_3_3_1)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = -v*conjugate(rhOverM_3_1_1)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_3_1_1*rhOverM_coeff*v;
    Modes[i++] = 0;
    Modes[i++] = rhOverM_3_3_1*rhOverM_coeff*v;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = rhOverM_4_0_0*rhOverM_coeff;
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

public:
  WaveformModes_0p50PN(const double m1_i, const double m2_i, const double v_i) :
    m1(m1_i), m2(m2_i), v(v_i), m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)),
    rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), rhOverM_2_0_0(-0.145802960879951),
    rhOverM_2_1_1(0.333333333333333*I*delta), rhOverM_2_2_0(1.00000000000000),
    rhOverM_3_1_1(0.0222717701593687*I*delta), rhOverM_3_3_1(-0.776323754260148*I*delta),
    rhOverM_4_0_0(-0.00140298964521140)
  { }

}; // class WaveformModes_0p50PN : public WaveformModes


class WaveformModes_1p0PN : public WaveformModes {
private:
  const double m1, m2;
  double v, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, ellHat_x, ellHat_y, ellHat_z;
  const double m, delta, nu;
  double chi1_l, chi2_l, chi_s_l, chi_a_l;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> rhOverM_2_0_0, rhOverM_2_1_1, rhOverM_2_2_0, rhOverM_2_2_2, rhOverM_3_1_1, rhOverM_3_2_2,
                             rhOverM_3_3_1, rhOverM_4_0_0, rhOverM_4_2_2, rhOverM_4_4_2;
  std::complex<double> rhOverM_2_1_SO_2;

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

    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    rhOverM_2_1_SO_2 = 0.5*I*(chi_a_l + chi_s_l*delta);

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = (pow(v, 2)*conjugate(rhOverM_2_2_2) + conjugate(rhOverM_2_2_0))*conjugate(rhOverM_coeff);
    Modes[i++] = v*(v*conjugate(rhOverM_2_1_SO_2) + conjugate(rhOverM_2_1_1))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_2_0_0*rhOverM_coeff;
    Modes[i++] = rhOverM_coeff*v*(rhOverM_2_1_1 + rhOverM_2_1_SO_2*v);
    Modes[i++] = rhOverM_coeff*(rhOverM_2_2_0 + rhOverM_2_2_2*pow(v, 2));
    Modes[i++] = -v*conjugate(rhOverM_3_3_1)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 2)*conjugate(rhOverM_3_2_2)*conjugate(rhOverM_coeff);
    Modes[i++] = -v*conjugate(rhOverM_3_1_1)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_3_1_1*rhOverM_coeff*v;
    Modes[i++] = rhOverM_3_2_2*rhOverM_coeff*pow(v, 2);
    Modes[i++] = rhOverM_3_3_1*rhOverM_coeff*v;
    Modes[i++] = pow(v, 2)*conjugate(rhOverM_4_4_2)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = pow(v, 2)*conjugate(rhOverM_4_2_2)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_4_0_0*rhOverM_coeff;
    Modes[i++] = 0;
    Modes[i++] = rhOverM_4_2_2*rhOverM_coeff*pow(v, 2);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_4_4_2*rhOverM_coeff*pow(v, 2);
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

public:
  WaveformModes_1p0PN(const double m1_i, const double m2_i, const double v_i, const double chi1_x_i, const double
                      chi1_y_i, const double chi1_z_i, const double chi2_x_i, const double chi2_y_i, const double
                      chi2_z_i, const double ellHat_x_i, const double ellHat_y_i, const double ellHat_z_i) :
    m1(m1_i), m2(m2_i), v(v_i), chi1_x(chi1_x_i), chi1_y(chi1_y_i), chi1_z(chi1_z_i), chi2_x(chi2_x_i),
    chi2_y(chi2_y_i), chi2_z(chi2_z_i), ellHat_x(ellHat_x_i), ellHat_y(ellHat_y_i), ellHat_z(ellHat_z_i), m(m1 + m2),
    delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z),
    chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi_s_l(0.5*chi1_l + 0.5*chi2_l), chi_a_l(0.5*chi1_l -
    0.5*chi2_l), rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), rhOverM_2_0_0(-0.145802960879951),
    rhOverM_2_1_1(0.333333333333333*I*delta), rhOverM_2_2_0(1.00000000000000), rhOverM_2_2_2(1.30952380952381*nu -
    2.54761904761905), rhOverM_3_1_1(0.0222717701593687*I*delta), rhOverM_3_2_2(-0.845154254728516*nu +
    0.281718084909506), rhOverM_3_3_1(-0.776323754260148*I*delta), rhOverM_4_0_0(-0.00140298964521140),
    rhOverM_4_2_2(-0.10647942749999*nu + 0.0354931424999967), rhOverM_4_4_2(2.25374467927604*nu - 0.751248226425348),
    rhOverM_2_1_SO_2(0.5*I*(chi_a_l + chi_s_l*delta))
  { }

}; // class WaveformModes_1p0PN : public WaveformModes


class WaveformModes_1p5PN : public WaveformModes {
private:
  const double m1, m2;
  double v, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, ellHat_x, ellHat_y, ellHat_z;
  const double m, delta, nu;
  double chi1_l, chi2_l, chi_s_l, chi_a_l;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> rhOverM_2_0_0, rhOverM_2_1_1, rhOverM_2_1_3, rhOverM_2_2_0, rhOverM_2_2_2, rhOverM_2_2_3,
                             rhOverM_3_1_1, rhOverM_3_1_3, rhOverM_3_2_2, rhOverM_3_3_1, rhOverM_3_3_3, rhOverM_4_0_0,
                             rhOverM_4_1_3, rhOverM_4_2_2, rhOverM_4_3_3, rhOverM_4_4_2, rhOverM_5_1_3, rhOverM_5_3_3,
                             rhOverM_5_5_3;
  std::complex<double> rhOverM_2_1_SO_2, rhOverM_2_2_SO_3, rhOverM_3_2_SO_3;

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

    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    rhOverM_2_1_SO_2 = 0.5*I*(chi_a_l + chi_s_l*delta);
    rhOverM_2_2_SO_3 = -1.33333333333333*chi_a_l*delta + 1.33333333333333*chi_s_l*nu - 1.33333333333333*chi_s_l;
    rhOverM_3_2_SO_3 = 1.12687233963802*chi_s_l*nu;

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = (pow(v, 2)*(v*(conjugate(rhOverM_2_2_3) + conjugate(rhOverM_2_2_SO_3)) + conjugate(rhOverM_2_2_2)) +
      conjugate(rhOverM_2_2_0))*conjugate(rhOverM_coeff);
    Modes[i++] = v*(v*(v*conjugate(rhOverM_2_1_3) + conjugate(rhOverM_2_1_SO_2)) +
      conjugate(rhOverM_2_1_1))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_2_0_0*rhOverM_coeff;
    Modes[i++] = rhOverM_coeff*v*(rhOverM_2_1_1 + v*(rhOverM_2_1_3*v + rhOverM_2_1_SO_2));
    Modes[i++] = rhOverM_coeff*(rhOverM_2_2_0 + pow(v, 2)*(rhOverM_2_2_2 + v*(rhOverM_2_2_3 + rhOverM_2_2_SO_3)));
    Modes[i++] = -v*(pow(v, 2)*conjugate(rhOverM_3_3_3) + conjugate(rhOverM_3_3_1))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 2)*(v*conjugate(rhOverM_3_2_SO_3) + conjugate(rhOverM_3_2_2))*conjugate(rhOverM_coeff);
    Modes[i++] = -v*(pow(v, 2)*conjugate(rhOverM_3_1_3) + conjugate(rhOverM_3_1_1))*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_coeff*v*(rhOverM_3_1_1 + rhOverM_3_1_3*pow(v, 2));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_3_2_2 + rhOverM_3_2_SO_3*v);
    Modes[i++] = rhOverM_coeff*v*(rhOverM_3_3_1 + rhOverM_3_3_3*pow(v, 2));
    Modes[i++] = pow(v, 2)*conjugate(rhOverM_4_4_2)*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 3)*conjugate(rhOverM_4_3_3)*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 2)*conjugate(rhOverM_4_2_2)*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 3)*conjugate(rhOverM_4_1_3)*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_4_0_0*rhOverM_coeff;
    Modes[i++] = rhOverM_4_1_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = rhOverM_4_2_2*rhOverM_coeff*pow(v, 2);
    Modes[i++] = rhOverM_4_3_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = rhOverM_4_4_2*rhOverM_coeff*pow(v, 2);
    Modes[i++] = -pow(v, 3)*conjugate(rhOverM_5_5_3)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = -pow(v, 3)*conjugate(rhOverM_5_3_3)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = -pow(v, 3)*conjugate(rhOverM_5_1_3)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_5_1_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_5_3_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_5_5_3*rhOverM_coeff*pow(v, 3);
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

public:
  WaveformModes_1p5PN(const double m1_i, const double m2_i, const double v_i, const double chi1_x_i, const double
                      chi1_y_i, const double chi1_z_i, const double chi2_x_i, const double chi2_y_i, const double
                      chi2_z_i, const double ellHat_x_i, const double ellHat_y_i, const double ellHat_z_i) :
    m1(m1_i), m2(m2_i), v(v_i), chi1_x(chi1_x_i), chi1_y(chi1_y_i), chi1_z(chi1_z_i), chi2_x(chi2_x_i),
    chi2_y(chi2_y_i), chi2_z(chi2_z_i), ellHat_x(ellHat_x_i), ellHat_y(ellHat_y_i), ellHat_z(ellHat_z_i), m(m1 + m2),
    delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z),
    chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi_s_l(0.5*chi1_l + 0.5*chi2_l), chi_a_l(0.5*chi1_l -
    0.5*chi2_l), rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), rhOverM_2_0_0(-0.145802960879951),
    rhOverM_2_1_1(0.333333333333333*I*delta), rhOverM_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)),
    rhOverM_2_2_0(1.00000000000000), rhOverM_2_2_2(1.30952380952381*nu - 2.54761904761905),
    rhOverM_2_2_3(6.28318530717959), rhOverM_3_1_1(0.0222717701593687*I*delta),
    rhOverM_3_1_3(-0.0148478467729125*I*delta*(nu + 4.0)), rhOverM_3_2_2(-0.845154254728516*nu + 0.281718084909506),
    rhOverM_3_3_1(-0.776323754260148*I*delta), rhOverM_3_3_3(-1.5526475085203*I*delta*(nu - 2.0)),
    rhOverM_4_0_0(-0.00140298964521140), rhOverM_4_1_3(0.00376461626210521*I*delta*(-2.0*nu + 1.0)),
    rhOverM_4_2_2(-0.10647942749999*nu + 0.0354931424999967), rhOverM_4_3_3(0.268926437100239*I*delta*(2.0*nu - 1.0)),
    rhOverM_4_4_2(2.25374467927604*nu - 0.751248226425348), rhOverM_5_1_3(0.000176960830360287*I*delta*(-2.0*nu + 1.0)),
    rhOverM_5_3_3(0.0464469088412683*I*delta*(2.0*nu - 1.0)), rhOverM_5_5_3(-0.801376894396698*I*delta*(2.0*nu - 1.0)),
    rhOverM_2_1_SO_2(0.5*I*(chi_a_l + chi_s_l*delta)), rhOverM_2_2_SO_3(-1.33333333333333*chi_a_l*delta +
    1.33333333333333*chi_s_l*nu - 1.33333333333333*chi_s_l), rhOverM_3_2_SO_3(1.12687233963802*chi_s_l*nu)
  { }

}; // class WaveformModes_1p5PN : public WaveformModes


class WaveformModes_2p0PN : public WaveformModes {
private:
  const double m1, m2;
  double v, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, ellHat_x, ellHat_y, ellHat_z;
  const double m, delta, nu;
  double chi1chi2, chi1_l, chi2_l, chi_s_l, chi_a_l;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> rhOverM_2_0_0, rhOverM_2_1_1, rhOverM_2_1_3, rhOverM_2_1_4, rhOverM_2_2_0, rhOverM_2_2_2,
                             rhOverM_2_2_3, rhOverM_2_2_4, rhOverM_3_1_1, rhOverM_3_1_3, rhOverM_3_1_4, rhOverM_3_2_2,
                             rhOverM_3_2_4, rhOverM_3_3_1, rhOverM_3_3_3, rhOverM_3_3_4, rhOverM_4_0_0, rhOverM_4_1_3,
                             rhOverM_4_2_2, rhOverM_4_2_4, rhOverM_4_3_3, rhOverM_4_4_2, rhOverM_4_4_4, rhOverM_5_1_3,
                             rhOverM_5_2_4, rhOverM_5_3_3, rhOverM_5_4_4, rhOverM_5_5_3, rhOverM_6_2_4, rhOverM_6_4_4,
                             rhOverM_6_6_4;
  std::complex<double> rhOverM_2_1_SO_2, rhOverM_2_1_SO_4, rhOverM_2_2_SO_3, rhOverM_3_1_SO_4, rhOverM_3_2_SO_3,
                       rhOverM_3_3_SO_4, rhOverM_4_1_SO_4, rhOverM_4_3_SO_4, rhOverM_2_2_SQ_4;

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

    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    rhOverM_2_1_SO_2 = 0.5*I*(chi_a_l + chi_s_l*delta);
    rhOverM_2_1_SO_4 = 0.0238095238095238*chi_a_l*(-205.0*nu + 7.0) + 0.0238095238095238*chi_s_l*delta*(-33.0*nu + 7.0);
    rhOverM_2_2_SO_3 = -1.33333333333333*chi_a_l*delta + 1.33333333333333*chi_s_l*nu - 1.33333333333333*chi_s_l;
    rhOverM_3_1_SO_4 = 0.0111358850796843*chi_a_l*(-11.0*nu + 4.0) + 0.0111358850796843*chi_s_l*delta*(-13.0*nu + 4.0);
    rhOverM_3_2_SO_3 = 1.12687233963802*chi_s_l*nu;
    rhOverM_3_3_SO_4 = -0.388161877130074*chi_a_l*(19.0*nu - 4.0) - 0.388161877130074*chi_s_l*delta*(5.0*nu - 4.0);
    rhOverM_4_1_SO_4 = 0.00941154065526303*nu*(chi_a_l - chi_s_l*delta);
    rhOverM_4_3_SO_4 = 0.672316092750596*nu*(chi_a_l - chi_s_l*delta);
    rhOverM_2_2_SQ_4 = 2.0*chi1chi2*nu;

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = (pow(v, 2)*(v*(v*(conjugate(rhOverM_2_2_4) + conjugate(rhOverM_2_2_SQ_4)) + conjugate(rhOverM_2_2_3) +
      conjugate(rhOverM_2_2_SO_3)) + conjugate(rhOverM_2_2_2)) + conjugate(rhOverM_2_2_0))*conjugate(rhOverM_coeff);
    Modes[i++] = v*(v*(v*(v*(conjugate(rhOverM_2_1_4) + conjugate(rhOverM_2_1_SO_4)) + conjugate(rhOverM_2_1_3)) +
      conjugate(rhOverM_2_1_SO_2)) + conjugate(rhOverM_2_1_1))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_2_0_0*rhOverM_coeff;
    Modes[i++] = rhOverM_coeff*v*(rhOverM_2_1_1 + v*(rhOverM_2_1_SO_2 + v*(rhOverM_2_1_3 + v*(rhOverM_2_1_4 +
      rhOverM_2_1_SO_4))));
    Modes[i++] = rhOverM_coeff*(rhOverM_2_2_0 + pow(v, 2)*(rhOverM_2_2_2 + v*(rhOverM_2_2_3 + rhOverM_2_2_SO_3 +
      v*(rhOverM_2_2_4 + rhOverM_2_2_SQ_4))));
    Modes[i++] = -v*(pow(v, 2)*(v*(conjugate(rhOverM_3_3_4) + conjugate(rhOverM_3_3_SO_4)) + conjugate(rhOverM_3_3_3)) +
      conjugate(rhOverM_3_3_1))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 2)*(v*(v*conjugate(rhOverM_3_2_4) + conjugate(rhOverM_3_2_SO_3)) +
      conjugate(rhOverM_3_2_2))*conjugate(rhOverM_coeff);
    Modes[i++] = -v*(pow(v, 2)*(v*(conjugate(rhOverM_3_1_4) + conjugate(rhOverM_3_1_SO_4)) + conjugate(rhOverM_3_1_3)) +
      conjugate(rhOverM_3_1_1))*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_coeff*v*(rhOverM_3_1_1 + pow(v, 2)*(rhOverM_3_1_3 + v*(rhOverM_3_1_4 + rhOverM_3_1_SO_4)));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_3_2_2 + v*(rhOverM_3_2_4*v + rhOverM_3_2_SO_3));
    Modes[i++] = rhOverM_coeff*v*(rhOverM_3_3_1 + pow(v, 2)*(rhOverM_3_3_3 + v*(rhOverM_3_3_4 + rhOverM_3_3_SO_4)));
    Modes[i++] = pow(v, 2)*(pow(v, 2)*conjugate(rhOverM_4_4_4) + conjugate(rhOverM_4_4_2))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 3)*(v*conjugate(rhOverM_4_3_SO_4) + conjugate(rhOverM_4_3_3))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 2)*(pow(v, 2)*conjugate(rhOverM_4_2_4) + conjugate(rhOverM_4_2_2))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 3)*(v*conjugate(rhOverM_4_1_SO_4) + conjugate(rhOverM_4_1_3))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_4_0_0*rhOverM_coeff;
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_4_1_3 + rhOverM_4_1_SO_4*v);
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_4_2_2 + rhOverM_4_2_4*pow(v, 2));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_4_3_3 + rhOverM_4_3_SO_4*v);
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_4_4_2 + rhOverM_4_4_4*pow(v, 2));
    Modes[i++] = -pow(v, 3)*conjugate(rhOverM_5_5_3)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 4)*conjugate(rhOverM_5_4_4)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 3)*conjugate(rhOverM_5_3_3)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 4)*conjugate(rhOverM_5_2_4)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 3)*conjugate(rhOverM_5_1_3)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_5_1_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = rhOverM_5_2_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = rhOverM_5_3_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = rhOverM_5_4_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = rhOverM_5_5_3*rhOverM_coeff*pow(v, 3);
    Modes[i++] = pow(v, 4)*conjugate(rhOverM_6_6_4)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = pow(v, 4)*conjugate(rhOverM_6_4_4)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = pow(v, 4)*conjugate(rhOverM_6_2_4)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = rhOverM_6_2_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_6_4_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_6_6_4*rhOverM_coeff*pow(v, 4);
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

public:
  WaveformModes_2p0PN(const double m1_i, const double m2_i, const double v_i, const double chi1_x_i, const double
                      chi1_y_i, const double chi1_z_i, const double chi2_x_i, const double chi2_y_i, const double
                      chi2_z_i, const double ellHat_x_i, const double ellHat_y_i, const double ellHat_z_i) :
    m1(m1_i), m2(m2_i), v(v_i), chi1_x(chi1_x_i), chi1_y(chi1_y_i), chi1_z(chi1_z_i), chi2_x(chi2_x_i),
    chi2_y(chi2_y_i), chi2_z(chi2_z_i), ellHat_x(ellHat_x_i), ellHat_y(ellHat_y_i), ellHat_z(ellHat_z_i), m(m1 + m2),
    delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z),
    chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y +
    chi2_z*ellHat_z), chi_s_l(0.5*chi1_l + 0.5*chi2_l), chi_a_l(0.5*chi1_l - 0.5*chi2_l),
    rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), rhOverM_2_0_0(-0.145802960879951),
    rhOverM_2_1_1(0.333333333333333*I*delta), rhOverM_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)),
    rhOverM_2_1_4(delta*(0.628764787039964 + 1.0471975511966*I)), rhOverM_2_2_0(1.00000000000000),
    rhOverM_2_2_2(1.30952380952381*nu - 2.54761904761905), rhOverM_2_2_3(6.28318530717959),
    rhOverM_2_2_4(0.000661375661375661*nu*(2047.0*nu - 7483.0) - 1.43716931216931),
    rhOverM_3_1_1(0.0222717701593687*I*delta), rhOverM_3_1_3(-0.0148478467729125*I*delta*(nu + 4.0)),
    rhOverM_3_1_4(delta*(0.0620557076072073 + 0.0699688295151131*I)), rhOverM_3_2_2(-0.845154254728516*nu +
    0.281718084909506), rhOverM_3_2_4(0.00313020094343895*nu*(-365.0*nu + 725.0) - 0.604128782083717),
    rhOverM_3_3_1(-0.776323754260148*I*delta), rhOverM_3_3_3(-1.5526475085203*I*delta*(nu - 2.0)),
    rhOverM_3_3_4(delta*(-1.37192659820446 - 7.31667900957279*I)), rhOverM_4_0_0(-0.00140298964521140),
    rhOverM_4_1_3(0.00376461626210521*I*delta*(-2.0*nu + 1.0)), rhOverM_4_2_2(-0.10647942749999*nu +
    0.0354931424999967), rhOverM_4_2_4(0.000107554977272717*nu*(-285.0*nu + 4025.0) - 0.141004575204532),
    rhOverM_4_3_3(0.268926437100239*I*delta*(2.0*nu - 1.0)), rhOverM_4_4_2(2.25374467927604*nu - 0.751248226425348),
    rhOverM_4_4_4(0.0113825488852325*nu*(525.0*nu - 1273.0) + 4.04991089336574),
    rhOverM_5_1_3(0.000176960830360287*I*delta*(-2.0*nu + 1.0)), rhOverM_5_2_4(0.00998814611056655*nu*(5.0*nu - 5.0) +
    0.00998814611056655), rhOverM_5_3_3(0.0464469088412683*I*delta*(2.0*nu - 1.0)),
    rhOverM_5_4_4(-0.276799624590764*nu*(5.0*nu - 5.0) - 0.276799624590764),
    rhOverM_5_5_3(-0.801376894396698*I*delta*(2.0*nu - 1.0)), rhOverM_6_2_4(0.000835250737974468*nu*(5.0*nu - 5.0) +
    0.000835250737974468), rhOverM_6_4_4(-0.058558165806266*nu*(5.0*nu - 5.0) - 0.058558165806266),
    rhOverM_6_6_4(0.903141370807658*nu*(5.0*nu - 5.0) + 0.903141370807658), rhOverM_2_1_SO_2(0.5*I*(chi_a_l +
    chi_s_l*delta)), rhOverM_2_1_SO_4(0.0238095238095238*chi_a_l*(-205.0*nu + 7.0) +
    0.0238095238095238*chi_s_l*delta*(-33.0*nu + 7.0)), rhOverM_2_2_SO_3(-1.33333333333333*chi_a_l*delta +
    1.33333333333333*chi_s_l*nu - 1.33333333333333*chi_s_l), rhOverM_3_1_SO_4(0.0111358850796843*chi_a_l*(-11.0*nu +
    4.0) + 0.0111358850796843*chi_s_l*delta*(-13.0*nu + 4.0)), rhOverM_3_2_SO_3(1.12687233963802*chi_s_l*nu),
    rhOverM_3_3_SO_4(-0.388161877130074*chi_a_l*(19.0*nu - 4.0) - 0.388161877130074*chi_s_l*delta*(5.0*nu - 4.0)),
    rhOverM_4_1_SO_4(0.00941154065526303*nu*(chi_a_l - chi_s_l*delta)), rhOverM_4_3_SO_4(0.672316092750596*nu*(chi_a_l -
    chi_s_l*delta)), rhOverM_2_2_SQ_4(2.0*chi1chi2*nu)
  { }

}; // class WaveformModes_2p0PN : public WaveformModes


class WaveformModes_2p5PN : public WaveformModes {
private:
  const double m1, m2;
  double v, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, ellHat_x, ellHat_y, ellHat_z;
  const double m, delta, nu;
  double chi1chi2, chi1_l, chi2_l, chi_s_l, chi_a_l;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> rhOverM_2_0_0, rhOverM_2_1_1, rhOverM_2_1_3, rhOverM_2_1_4, rhOverM_2_1_5, rhOverM_2_2_0,
                             rhOverM_2_2_2, rhOverM_2_2_3, rhOverM_2_2_4, rhOverM_2_2_5, rhOverM_3_0_5, rhOverM_3_1_1,
                             rhOverM_3_1_3, rhOverM_3_1_4, rhOverM_3_1_5, rhOverM_3_2_2, rhOverM_3_2_4, rhOverM_3_2_5,
                             rhOverM_3_3_1, rhOverM_3_3_3, rhOverM_3_3_4, rhOverM_3_3_5, rhOverM_4_0_0, rhOverM_4_1_3,
                             rhOverM_4_1_5, rhOverM_4_2_2, rhOverM_4_2_4, rhOverM_4_2_5, rhOverM_4_3_3, rhOverM_4_3_5,
                             rhOverM_4_4_2, rhOverM_4_4_4, rhOverM_4_4_5, rhOverM_5_1_3, rhOverM_5_1_5, rhOverM_5_2_4,
                             rhOverM_5_3_3, rhOverM_5_3_5, rhOverM_5_4_4, rhOverM_5_5_3, rhOverM_5_5_5, rhOverM_6_1_5,
                             rhOverM_6_2_4, rhOverM_6_3_5, rhOverM_6_4_4, rhOverM_6_5_5, rhOverM_6_6_4, rhOverM_7_1_5,
                             rhOverM_7_3_5, rhOverM_7_5_5, rhOverM_7_7_5;
  std::complex<double> rhOverM_2_1_SO_2, rhOverM_2_1_SO_4, rhOverM_2_2_SO_3, rhOverM_3_1_SO_4, rhOverM_3_2_SO_3,
                       rhOverM_3_3_SO_4, rhOverM_4_1_SO_4, rhOverM_4_3_SO_4, rhOverM_2_2_SQ_4;

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

    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    rhOverM_2_1_SO_2 = 0.5*I*(chi_a_l + chi_s_l*delta);
    rhOverM_2_1_SO_4 = 0.0238095238095238*chi_a_l*(-205.0*nu + 7.0) + 0.0238095238095238*chi_s_l*delta*(-33.0*nu + 7.0);
    rhOverM_2_2_SO_3 = -1.33333333333333*chi_a_l*delta + 1.33333333333333*chi_s_l*nu - 1.33333333333333*chi_s_l;
    rhOverM_3_1_SO_4 = 0.0111358850796843*chi_a_l*(-11.0*nu + 4.0) + 0.0111358850796843*chi_s_l*delta*(-13.0*nu + 4.0);
    rhOverM_3_2_SO_3 = 1.12687233963802*chi_s_l*nu;
    rhOverM_3_3_SO_4 = -0.388161877130074*chi_a_l*(19.0*nu - 4.0) - 0.388161877130074*chi_s_l*delta*(5.0*nu - 4.0);
    rhOverM_4_1_SO_4 = 0.00941154065526303*nu*(chi_a_l - chi_s_l*delta);
    rhOverM_4_3_SO_4 = 0.672316092750596*nu*(chi_a_l - chi_s_l*delta);
    rhOverM_2_2_SQ_4 = 2.0*chi1chi2*nu;

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = (pow(v, 2)*(v*(v*(v*conjugate(rhOverM_2_2_5) + conjugate(rhOverM_2_2_4) + conjugate(rhOverM_2_2_SQ_4))
      + conjugate(rhOverM_2_2_3) + conjugate(rhOverM_2_2_SO_3)) + conjugate(rhOverM_2_2_2)) +
      conjugate(rhOverM_2_2_0))*conjugate(rhOverM_coeff);
    Modes[i++] = v*(v*(v*(v*(v*conjugate(rhOverM_2_1_5) + conjugate(rhOverM_2_1_4) + conjugate(rhOverM_2_1_SO_4)) +
      conjugate(rhOverM_2_1_3)) + conjugate(rhOverM_2_1_SO_2)) + conjugate(rhOverM_2_1_1))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_2_0_0*rhOverM_coeff;
    Modes[i++] = rhOverM_coeff*v*(rhOverM_2_1_1 + v*(rhOverM_2_1_SO_2 + v*(rhOverM_2_1_3 + v*(rhOverM_2_1_4 +
      rhOverM_2_1_5*v + rhOverM_2_1_SO_4))));
    Modes[i++] = rhOverM_coeff*(rhOverM_2_2_0 + pow(v, 2)*(rhOverM_2_2_2 + v*(rhOverM_2_2_3 + rhOverM_2_2_SO_3 +
      v*(rhOverM_2_2_4 + rhOverM_2_2_5*v + rhOverM_2_2_SQ_4))));
    Modes[i++] = -v*(pow(v, 2)*(v*(v*conjugate(rhOverM_3_3_5) + conjugate(rhOverM_3_3_4) + conjugate(rhOverM_3_3_SO_4))
      + conjugate(rhOverM_3_3_3)) + conjugate(rhOverM_3_3_1))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 2)*(v*(v*(v*conjugate(rhOverM_3_2_5) + conjugate(rhOverM_3_2_4)) + conjugate(rhOverM_3_2_SO_3))
      + conjugate(rhOverM_3_2_2))*conjugate(rhOverM_coeff);
    Modes[i++] = -v*(pow(v, 2)*(v*(v*conjugate(rhOverM_3_1_5) + conjugate(rhOverM_3_1_4) + conjugate(rhOverM_3_1_SO_4))
      + conjugate(rhOverM_3_1_3)) + conjugate(rhOverM_3_1_1))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_3_0_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*v*(rhOverM_3_1_1 + pow(v, 2)*(rhOverM_3_1_3 + v*(rhOverM_3_1_4 + rhOverM_3_1_5*v +
      rhOverM_3_1_SO_4)));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_3_2_2 + v*(rhOverM_3_2_SO_3 + v*(rhOverM_3_2_4 + rhOverM_3_2_5*v)));
    Modes[i++] = rhOverM_coeff*v*(rhOverM_3_3_1 + pow(v, 2)*(rhOverM_3_3_3 + v*(rhOverM_3_3_4 + rhOverM_3_3_5*v +
      rhOverM_3_3_SO_4)));
    Modes[i++] = pow(v, 2)*(pow(v, 2)*(v*conjugate(rhOverM_4_4_5) + conjugate(rhOverM_4_4_4)) +
      conjugate(rhOverM_4_4_2))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 3)*(v*(v*conjugate(rhOverM_4_3_5) + conjugate(rhOverM_4_3_SO_4)) +
      conjugate(rhOverM_4_3_3))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 2)*(pow(v, 2)*(v*conjugate(rhOverM_4_2_5) + conjugate(rhOverM_4_2_4)) +
      conjugate(rhOverM_4_2_2))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 3)*(v*(v*conjugate(rhOverM_4_1_5) + conjugate(rhOverM_4_1_SO_4)) +
      conjugate(rhOverM_4_1_3))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_4_0_0*rhOverM_coeff;
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_4_1_3 + v*(rhOverM_4_1_5*v + rhOverM_4_1_SO_4));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_4_2_2 + pow(v, 2)*(rhOverM_4_2_4 + rhOverM_4_2_5*v));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_4_3_3 + v*(rhOverM_4_3_5*v + rhOverM_4_3_SO_4));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_4_4_2 + pow(v, 2)*(rhOverM_4_4_4 + rhOverM_4_4_5*v));
    Modes[i++] = -pow(v, 3)*(pow(v, 2)*conjugate(rhOverM_5_5_5) + conjugate(rhOverM_5_5_3))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 4)*conjugate(rhOverM_5_4_4)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 3)*(pow(v, 2)*conjugate(rhOverM_5_3_5) + conjugate(rhOverM_5_3_3))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 4)*conjugate(rhOverM_5_2_4)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 3)*(pow(v, 2)*conjugate(rhOverM_5_1_5) + conjugate(rhOverM_5_1_3))*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_5_1_3 + rhOverM_5_1_5*pow(v, 2));
    Modes[i++] = rhOverM_5_2_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_5_3_3 + rhOverM_5_3_5*pow(v, 2));
    Modes[i++] = rhOverM_5_4_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_5_5_3 + rhOverM_5_5_5*pow(v, 2));
    Modes[i++] = pow(v, 4)*conjugate(rhOverM_6_6_4)*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 5)*conjugate(rhOverM_6_5_5)*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 4)*conjugate(rhOverM_6_4_4)*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 5)*conjugate(rhOverM_6_3_5)*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 4)*conjugate(rhOverM_6_2_4)*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 5)*conjugate(rhOverM_6_1_5)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_6_1_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_6_2_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = rhOverM_6_3_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_6_4_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = rhOverM_6_5_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_6_6_4*rhOverM_coeff*pow(v, 4);
    Modes[i++] = -pow(v, 5)*conjugate(rhOverM_7_7_5)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = -pow(v, 5)*conjugate(rhOverM_7_5_5)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = -pow(v, 5)*conjugate(rhOverM_7_3_5)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = -pow(v, 5)*conjugate(rhOverM_7_1_5)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_7_1_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_7_3_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_7_5_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_7_7_5*rhOverM_coeff*pow(v, 5);
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

public:
  WaveformModes_2p5PN(const double m1_i, const double m2_i, const double v_i, const double chi1_x_i, const double
                      chi1_y_i, const double chi1_z_i, const double chi2_x_i, const double chi2_y_i, const double
                      chi2_z_i, const double ellHat_x_i, const double ellHat_y_i, const double ellHat_z_i) :
    m1(m1_i), m2(m2_i), v(v_i), chi1_x(chi1_x_i), chi1_y(chi1_y_i), chi1_z(chi1_z_i), chi2_x(chi2_x_i),
    chi2_y(chi2_y_i), chi2_z(chi2_z_i), ellHat_x(ellHat_x_i), ellHat_y(ellHat_y_i), ellHat_z(ellHat_z_i), m(m1 + m2),
    delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z),
    chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y +
    chi2_z*ellHat_z), chi_s_l(0.5*chi1_l + 0.5*chi2_l), chi_a_l(0.5*chi1_l - 0.5*chi2_l),
    rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), rhOverM_2_0_0(-0.145802960879951),
    rhOverM_2_1_1(0.333333333333333*I*delta), rhOverM_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)),
    rhOverM_2_1_4(delta*(0.628764787039964 + 1.0471975511966*I)),
    rhOverM_2_1_5(0.000661375661375661*I*delta*(nu*(237.0*nu - 2036.0) - 172.0)), rhOverM_2_2_0(1.00000000000000),
    rhOverM_2_2_2(1.30952380952381*nu - 2.54761904761905), rhOverM_2_2_3(6.28318530717959),
    rhOverM_2_2_4(0.000661375661375661*nu*(2047.0*nu - 7483.0) - 1.43716931216931), rhOverM_2_2_5(5.08638810581205*nu -
    24.0*I*nu - 16.0071625682908), rhOverM_3_0_5(-0.370328039909021*I*nu), rhOverM_3_1_1(0.0222717701593687*I*delta),
    rhOverM_3_1_3(-0.0148478467729125*I*delta*(nu + 4.0)), rhOverM_3_1_4(delta*(0.0620557076072073 +
    0.0699688295151131*I)), rhOverM_3_1_5(-0.000112483687673579*I*delta*(nu*(247.0*nu + 272.0) - 607.0)),
    rhOverM_3_2_2(-0.845154254728516*nu + 0.281718084909506), rhOverM_3_2_4(0.00313020094343895*nu*(-365.0*nu + 725.0) -
    0.604128782083717), rhOverM_3_2_5(-5.31026079561053*nu + 3.71867872080547*I*nu + 1.77008693187018 -
    0.845154254728517*I), rhOverM_3_3_1(-0.776323754260148*I*delta), rhOverM_3_3_3(-1.5526475085203*I*delta*(nu - 2.0)),
    rhOverM_3_3_4(delta*(-1.37192659820446 - 7.31667900957279*I)),
    rhOverM_3_3_5(-0.00235249622503075*I*delta*(nu*(887.0*nu - 3676.0) + 369.0)), rhOverM_4_0_0(-0.00140298964521140),
    rhOverM_4_1_3(0.00376461626210521*I*delta*(-2.0*nu + 1.0)), rhOverM_4_1_5(-2.85198201674637e-5*I*delta*(nu*(332.0*nu
    - 1011.0) + 404.0)), rhOverM_4_2_2(-0.10647942749999*nu + 0.0354931424999967),
    rhOverM_4_2_4(0.000107554977272717*nu*(-285.0*nu + 4025.0) - 0.141004575204532),
    rhOverM_4_2_5(0.00709862849999933*nu*(-94.2477796076938 + 84.0*I) + 0.22300999146161 - 0.149071198499986*I),
    rhOverM_4_3_3(0.268926437100239*I*delta*(2.0*nu - 1.0)), rhOverM_4_3_5(0.00203732149318363*I*delta*(nu*(524.0*nu -
    1267.0) + 468.0)), rhOverM_4_4_2(2.25374467927604*nu - 0.751248226425348),
    rhOverM_4_4_4(0.0113825488852325*nu*(525.0*nu - 1273.0) + 4.04991089336574), rhOverM_4_4_5(28.3213909099228*nu +
    0.0187812056606337*I*(-527.578706662453*nu + 114.192902220818) - 9.44046363664094),
    rhOverM_5_1_3(0.000176960830360287*I*delta*(-2.0*nu + 1.0)), rhOverM_5_1_5(-4.5374571887253e-6*I*delta*(nu*(4.0*nu -
    352.0) + 179.0)), rhOverM_5_2_4(0.00998814611056655*nu*(5.0*nu - 5.0) + 0.00998814611056655),
    rhOverM_5_3_3(0.0464469088412683*I*delta*(2.0*nu - 1.0)), rhOverM_5_3_5(0.00119094638054534*I*delta*(8.0*nu*(11.0*nu
    - 58.0) + 207.0)), rhOverM_5_4_4(-0.276799624590764*nu*(5.0*nu - 5.0) - 0.276799624590764),
    rhOverM_5_5_3(-0.801376894396698*I*delta*(2.0*nu - 1.0)),
    rhOverM_5_5_5(-0.0205481254973512*I*delta*(16.0*nu*(16.0*nu - 43.0) + 263.0)),
    rhOverM_6_1_5(2.35829888333555e-5*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), rhOverM_6_2_4(0.000835250737974468*nu*(5.0*nu
    - 5.0) + 0.000835250737974468), rhOverM_6_3_5(-0.0163097621781264*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    rhOverM_6_4_4(-0.058558165806266*nu*(5.0*nu - 5.0) - 0.058558165806266), rhOverM_6_5_5(0.299357979653564*I*delta*(nu
    - 1.0)*(3.0*nu - 1.0)), rhOverM_6_6_4(0.903141370807658*nu*(5.0*nu - 5.0) + 0.903141370807658),
    rhOverM_7_1_5(8.17593033339979e-7*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), rhOverM_7_3_5(-0.0018582230503756*I*delta*(nu
    - 1.0)*(3.0*nu - 1.0)), rhOverM_7_5_5(0.0733861624905401*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    rhOverM_7_7_5(-1.05422444934392*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), rhOverM_2_1_SO_2(0.5*I*(chi_a_l +
    chi_s_l*delta)), rhOverM_2_1_SO_4(0.0238095238095238*chi_a_l*(-205.0*nu + 7.0) +
    0.0238095238095238*chi_s_l*delta*(-33.0*nu + 7.0)), rhOverM_2_2_SO_3(-1.33333333333333*chi_a_l*delta +
    1.33333333333333*chi_s_l*nu - 1.33333333333333*chi_s_l), rhOverM_3_1_SO_4(0.0111358850796843*chi_a_l*(-11.0*nu +
    4.0) + 0.0111358850796843*chi_s_l*delta*(-13.0*nu + 4.0)), rhOverM_3_2_SO_3(1.12687233963802*chi_s_l*nu),
    rhOverM_3_3_SO_4(-0.388161877130074*chi_a_l*(19.0*nu - 4.0) - 0.388161877130074*chi_s_l*delta*(5.0*nu - 4.0)),
    rhOverM_4_1_SO_4(0.00941154065526303*nu*(chi_a_l - chi_s_l*delta)), rhOverM_4_3_SO_4(0.672316092750596*nu*(chi_a_l -
    chi_s_l*delta)), rhOverM_2_2_SQ_4(2.0*chi1chi2*nu)
  { }

}; // class WaveformModes_2p5PN : public WaveformModes


class WaveformModes_3p0PN : public WaveformModes {
private:
  const double m1, m2;
  double v, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, ellHat_x, ellHat_y, ellHat_z;
  const double m, delta, nu;
  double chi1chi2, chi1_l, chi2_l, chi_s_l, chi_a_l, logv;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> rhOverM_2_0_0, rhOverM_2_1_1, rhOverM_2_1_3, rhOverM_2_1_4, rhOverM_2_1_5, rhOverM_2_1_6,
                             rhOverM_2_2_0, rhOverM_2_2_2, rhOverM_2_2_3, rhOverM_2_2_4, rhOverM_2_2_5, rhOverM_2_2_6,
                             rhOverM_2_2_lnv_6, rhOverM_3_0_5, rhOverM_3_1_1, rhOverM_3_1_3, rhOverM_3_1_4,
                             rhOverM_3_1_5, rhOverM_3_1_6, rhOverM_3_2_2, rhOverM_3_2_4, rhOverM_3_2_5, rhOverM_3_2_6,
                             rhOverM_3_3_1, rhOverM_3_3_3, rhOverM_3_3_4, rhOverM_3_3_5, rhOverM_3_3_6, rhOverM_4_0_0,
                             rhOverM_4_1_3, rhOverM_4_1_5, rhOverM_4_1_6, rhOverM_4_2_2, rhOverM_4_2_4, rhOverM_4_2_5,
                             rhOverM_4_2_6, rhOverM_4_3_3, rhOverM_4_3_5, rhOverM_4_3_6, rhOverM_4_4_2, rhOverM_4_4_4,
                             rhOverM_4_4_5, rhOverM_4_4_6, rhOverM_5_1_3, rhOverM_5_1_5, rhOverM_5_1_6, rhOverM_5_2_4,
                             rhOverM_5_2_6, rhOverM_5_3_3, rhOverM_5_3_5, rhOverM_5_3_6, rhOverM_5_4_4, rhOverM_5_4_6,
                             rhOverM_5_5_3, rhOverM_5_5_5, rhOverM_5_5_6, rhOverM_6_1_5, rhOverM_6_2_4, rhOverM_6_2_6,
                             rhOverM_6_3_5, rhOverM_6_4_4, rhOverM_6_4_6, rhOverM_6_5_5, rhOverM_6_6_4, rhOverM_6_6_6,
                             rhOverM_7_1_5, rhOverM_7_2_6, rhOverM_7_3_5, rhOverM_7_4_6, rhOverM_7_5_5, rhOverM_7_6_6,
                             rhOverM_7_7_5, rhOverM_8_2_6, rhOverM_8_4_6, rhOverM_8_6_6, rhOverM_8_8_6;
  std::complex<double> rhOverM_2_1_SO_2, rhOverM_2_1_SO_4, rhOverM_2_2_SO_3, rhOverM_3_1_SO_4, rhOverM_3_2_SO_3,
                       rhOverM_3_3_SO_4, rhOverM_4_1_SO_4, rhOverM_4_3_SO_4, rhOverM_2_2_SQ_4;

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

    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    logv = log(v);
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    rhOverM_2_1_SO_2 = 0.5*I*(chi_a_l + chi_s_l*delta);
    rhOverM_2_1_SO_4 = 0.0238095238095238*chi_a_l*(-205.0*nu + 7.0) + 0.0238095238095238*chi_s_l*delta*(-33.0*nu + 7.0);
    rhOverM_2_2_SO_3 = -1.33333333333333*chi_a_l*delta + 1.33333333333333*chi_s_l*nu - 1.33333333333333*chi_s_l;
    rhOverM_3_1_SO_4 = 0.0111358850796843*chi_a_l*(-11.0*nu + 4.0) + 0.0111358850796843*chi_s_l*delta*(-13.0*nu + 4.0);
    rhOverM_3_2_SO_3 = 1.12687233963802*chi_s_l*nu;
    rhOverM_3_3_SO_4 = -0.388161877130074*chi_a_l*(19.0*nu - 4.0) - 0.388161877130074*chi_s_l*delta*(5.0*nu - 4.0);
    rhOverM_4_1_SO_4 = 0.00941154065526303*nu*(chi_a_l - chi_s_l*delta);
    rhOverM_4_3_SO_4 = 0.672316092750596*nu*(chi_a_l - chi_s_l*delta);
    rhOverM_2_2_SQ_4 = 2.0*chi1chi2*nu;

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = (pow(v, 2)*(v*(v*(v*(v*(logv*conjugate(rhOverM_2_2_lnv_6) + conjugate(rhOverM_2_2_6)) +
      conjugate(rhOverM_2_2_5)) + conjugate(rhOverM_2_2_4) + conjugate(rhOverM_2_2_SQ_4)) + conjugate(rhOverM_2_2_3) +
      conjugate(rhOverM_2_2_SO_3)) + conjugate(rhOverM_2_2_2)) + conjugate(rhOverM_2_2_0))*conjugate(rhOverM_coeff);
    Modes[i++] = v*(v*(v*(v*(v*(v*conjugate(rhOverM_2_1_6) + conjugate(rhOverM_2_1_5)) + conjugate(rhOverM_2_1_4) +
      conjugate(rhOverM_2_1_SO_4)) + conjugate(rhOverM_2_1_3)) + conjugate(rhOverM_2_1_SO_2)) +
      conjugate(rhOverM_2_1_1))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_2_0_0*rhOverM_coeff;
    Modes[i++] = rhOverM_coeff*v*(rhOverM_2_1_1 + v*(rhOverM_2_1_SO_2 + v*(rhOverM_2_1_3 + v*(rhOverM_2_1_4 +
      rhOverM_2_1_SO_4 + v*(rhOverM_2_1_5 + rhOverM_2_1_6*v)))));
    Modes[i++] = rhOverM_coeff*(rhOverM_2_2_0 + pow(v, 2)*(rhOverM_2_2_2 + v*(rhOverM_2_2_3 + rhOverM_2_2_SO_3 +
      v*(rhOverM_2_2_4 + rhOverM_2_2_SQ_4 + v*(rhOverM_2_2_5 + v*(logv*rhOverM_2_2_lnv_6 + rhOverM_2_2_6))))));
    Modes[i++] = -v*(pow(v, 2)*(v*(v*(v*conjugate(rhOverM_3_3_6) + conjugate(rhOverM_3_3_5)) + conjugate(rhOverM_3_3_4)
      + conjugate(rhOverM_3_3_SO_4)) + conjugate(rhOverM_3_3_3)) + conjugate(rhOverM_3_3_1))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 2)*(v*(v*(v*(v*conjugate(rhOverM_3_2_6) + conjugate(rhOverM_3_2_5)) + conjugate(rhOverM_3_2_4))
      + conjugate(rhOverM_3_2_SO_3)) + conjugate(rhOverM_3_2_2))*conjugate(rhOverM_coeff);
    Modes[i++] = -v*(pow(v, 2)*(v*(v*(v*conjugate(rhOverM_3_1_6) + conjugate(rhOverM_3_1_5)) + conjugate(rhOverM_3_1_4)
      + conjugate(rhOverM_3_1_SO_4)) + conjugate(rhOverM_3_1_3)) + conjugate(rhOverM_3_1_1))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_3_0_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*v*(rhOverM_3_1_1 + pow(v, 2)*(rhOverM_3_1_3 + v*(rhOverM_3_1_4 + rhOverM_3_1_SO_4 +
      v*(rhOverM_3_1_5 + rhOverM_3_1_6*v))));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_3_2_2 + v*(rhOverM_3_2_SO_3 + v*(rhOverM_3_2_4 + v*(rhOverM_3_2_5 +
      rhOverM_3_2_6*v))));
    Modes[i++] = rhOverM_coeff*v*(rhOverM_3_3_1 + pow(v, 2)*(rhOverM_3_3_3 + v*(rhOverM_3_3_4 + rhOverM_3_3_SO_4 +
      v*(rhOverM_3_3_5 + rhOverM_3_3_6*v))));
    Modes[i++] = pow(v, 2)*(pow(v, 2)*(v*(v*conjugate(rhOverM_4_4_6) + conjugate(rhOverM_4_4_5)) +
      conjugate(rhOverM_4_4_4)) + conjugate(rhOverM_4_4_2))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 3)*(v*(v*(v*conjugate(rhOverM_4_3_6) + conjugate(rhOverM_4_3_5)) + conjugate(rhOverM_4_3_SO_4))
      + conjugate(rhOverM_4_3_3))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 2)*(pow(v, 2)*(v*(v*conjugate(rhOverM_4_2_6) + conjugate(rhOverM_4_2_5)) +
      conjugate(rhOverM_4_2_4)) + conjugate(rhOverM_4_2_2))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 3)*(v*(v*(v*conjugate(rhOverM_4_1_6) + conjugate(rhOverM_4_1_5)) + conjugate(rhOverM_4_1_SO_4))
      + conjugate(rhOverM_4_1_3))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_4_0_0*rhOverM_coeff;
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_4_1_3 + v*(rhOverM_4_1_SO_4 + v*(rhOverM_4_1_5 + rhOverM_4_1_6*v)));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_4_2_2 + pow(v, 2)*(rhOverM_4_2_4 + v*(rhOverM_4_2_5 +
      rhOverM_4_2_6*v)));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_4_3_3 + v*(rhOverM_4_3_SO_4 + v*(rhOverM_4_3_5 + rhOverM_4_3_6*v)));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_4_4_2 + pow(v, 2)*(rhOverM_4_4_4 + v*(rhOverM_4_4_5 +
      rhOverM_4_4_6*v)));
    Modes[i++] = -pow(v, 3)*(pow(v, 2)*(v*conjugate(rhOverM_5_5_6) + conjugate(rhOverM_5_5_5)) +
      conjugate(rhOverM_5_5_3))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 4)*(pow(v, 2)*conjugate(rhOverM_5_4_6) + conjugate(rhOverM_5_4_4))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 3)*(pow(v, 2)*(v*conjugate(rhOverM_5_3_6) + conjugate(rhOverM_5_3_5)) +
      conjugate(rhOverM_5_3_3))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 4)*(pow(v, 2)*conjugate(rhOverM_5_2_6) + conjugate(rhOverM_5_2_4))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 3)*(pow(v, 2)*(v*conjugate(rhOverM_5_1_6) + conjugate(rhOverM_5_1_5)) +
      conjugate(rhOverM_5_1_3))*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_5_1_3 + pow(v, 2)*(rhOverM_5_1_5 + rhOverM_5_1_6*v));
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(rhOverM_5_2_4 + rhOverM_5_2_6*pow(v, 2));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_5_3_3 + pow(v, 2)*(rhOverM_5_3_5 + rhOverM_5_3_6*v));
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(rhOverM_5_4_4 + rhOverM_5_4_6*pow(v, 2));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_5_5_3 + pow(v, 2)*(rhOverM_5_5_5 + rhOverM_5_5_6*v));
    Modes[i++] = pow(v, 4)*(pow(v, 2)*conjugate(rhOverM_6_6_6) + conjugate(rhOverM_6_6_4))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 5)*conjugate(rhOverM_6_5_5)*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 4)*(pow(v, 2)*conjugate(rhOverM_6_4_6) + conjugate(rhOverM_6_4_4))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 5)*conjugate(rhOverM_6_3_5)*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 4)*(pow(v, 2)*conjugate(rhOverM_6_2_6) + conjugate(rhOverM_6_2_4))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 5)*conjugate(rhOverM_6_1_5)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_6_1_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(rhOverM_6_2_4 + rhOverM_6_2_6*pow(v, 2));
    Modes[i++] = rhOverM_6_3_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(rhOverM_6_4_4 + rhOverM_6_4_6*pow(v, 2));
    Modes[i++] = rhOverM_6_5_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(rhOverM_6_6_4 + rhOverM_6_6_6*pow(v, 2));
    Modes[i++] = -pow(v, 5)*conjugate(rhOverM_7_7_5)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 6)*conjugate(rhOverM_7_6_6)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 5)*conjugate(rhOverM_7_5_5)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 6)*conjugate(rhOverM_7_4_6)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 5)*conjugate(rhOverM_7_3_5)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 6)*conjugate(rhOverM_7_2_6)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 5)*conjugate(rhOverM_7_1_5)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_7_1_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_7_2_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = rhOverM_7_3_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_7_4_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = rhOverM_7_5_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_7_6_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = rhOverM_7_7_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = pow(v, 6)*conjugate(rhOverM_8_8_6)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = pow(v, 6)*conjugate(rhOverM_8_6_6)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = pow(v, 6)*conjugate(rhOverM_8_4_6)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = pow(v, 6)*conjugate(rhOverM_8_2_6)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = rhOverM_8_2_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_8_4_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_8_6_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_8_8_6*rhOverM_coeff*pow(v, 6);

    return Modes;
  }

public:
  WaveformModes_3p0PN(const double m1_i, const double m2_i, const double v_i, const double chi1_x_i, const double
                      chi1_y_i, const double chi1_z_i, const double chi2_x_i, const double chi2_y_i, const double
                      chi2_z_i, const double ellHat_x_i, const double ellHat_y_i, const double ellHat_z_i) :
    m1(m1_i), m2(m2_i), v(v_i), chi1_x(chi1_x_i), chi1_y(chi1_y_i), chi1_z(chi1_z_i), chi2_x(chi2_x_i),
    chi2_y(chi2_y_i), chi2_z(chi2_z_i), ellHat_x(ellHat_x_i), ellHat_y(ellHat_y_i), ellHat_z(ellHat_z_i), m(m1 + m2),
    delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z),
    chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y +
    chi2_z*ellHat_z), chi_s_l(0.5*chi1_l + 0.5*chi2_l), chi_a_l(0.5*chi1_l - 0.5*chi2_l), logv(log(v)),
    rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), rhOverM_2_0_0(-0.145802960879951),
    rhOverM_2_1_1(0.333333333333333*I*delta), rhOverM_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)),
    rhOverM_2_1_4(delta*(0.628764787039964 + 1.0471975511966*I)),
    rhOverM_2_1_5(0.000661375661375661*I*delta*(nu*(237.0*nu - 2036.0) - 172.0)),
    rhOverM_2_1_6(0.00595238095238095*delta*(nu*(722.635532333439 + 37.6991118430775*I) - 64.1340082780763 -
    106.814150222053*I)), rhOverM_2_2_0(1.00000000000000), rhOverM_2_2_2(1.30952380952381*nu - 2.54761904761905),
    rhOverM_2_2_3(6.28318530717959), rhOverM_2_2_4(0.000661375661375661*nu*(2047.0*nu - 7483.0) - 1.43716931216931),
    rhOverM_2_2_5(5.08638810581205*nu - 24.0*I*nu - 16.0071625682908),
    rhOverM_2_2_6(1.00208433541767e-5*nu*(nu*(114635.0*nu - 729396.0) - 834555.0) + 4.21514354629858*nu +
    32.3588011610077 + 12.8057300546327*I), rhOverM_2_2_lnv_6(-8.15238095238095),
    rhOverM_3_0_5(-0.370328039909021*I*nu), rhOverM_3_1_1(0.0222717701593687*I*delta),
    rhOverM_3_1_3(-0.0148478467729125*I*delta*(nu + 4.0)), rhOverM_3_1_4(delta*(0.0620557076072073 +
    0.0699688295151131*I)), rhOverM_3_1_5(-0.000112483687673579*I*delta*(nu*(247.0*nu + 272.0) - 607.0)),
    rhOverM_3_1_6(0.000742392338645623*delta*(-46.5203026391962*nu - 15.707963267949*I*(7.0*nu + 16.0) -
    222.903548889591)), rhOverM_3_2_2(-0.845154254728516*nu + 0.281718084909506),
    rhOverM_3_2_4(0.00313020094343895*nu*(-365.0*nu + 725.0) - 0.604128782083717), rhOverM_3_2_5(-5.31026079561053*nu +
    3.71867872080547*I*nu + 1.77008693187018 - 0.845154254728517*I),
    rhOverM_3_2_6(7.11409305327034e-5*nu*(nu*(-16023.0*nu + 100026.0) - 17387.0) - 0.103225490202953),
    rhOverM_3_3_1(-0.776323754260148*I*delta), rhOverM_3_3_3(-1.5526475085203*I*delta*(nu - 2.0)),
    rhOverM_3_3_4(delta*(-1.37192659820446 - 7.31667900957279*I)),
    rhOverM_3_3_5(-0.00235249622503075*I*delta*(nu*(887.0*nu - 3676.0) + 369.0)),
    rhOverM_3_3_6(0.000319474795991831*delta*(-87338.4780856744*nu - 11451.1052223348*I*(3.0*nu - 8.0) +
    17177.2748951318)), rhOverM_4_0_0(-0.00140298964521140), rhOverM_4_1_3(0.00376461626210521*I*delta*(-2.0*nu + 1.0)),
    rhOverM_4_1_5(-2.85198201674637e-5*I*delta*(nu*(332.0*nu - 1011.0) + 404.0)),
    rhOverM_4_1_6(0.00012548720873684*delta*(-1744.17766166719*nu - 94.2477796076938*I*(2.0*nu - 1.0) +
    105.588830833597)), rhOverM_4_2_2(-0.10647942749999*nu + 0.0354931424999967),
    rhOverM_4_2_4(0.000107554977272717*nu*(-285.0*nu + 4025.0) - 0.141004575204532),
    rhOverM_4_2_5(0.00709862849999933*nu*(-94.2477796076938 + 84.0*I) + 0.22300999146161 - 0.149071198499986*I),
    rhOverM_4_2_6(1.37890996503484e-7*nu*(115.0*nu*(3363.0*nu + 34822.0) - 5460759.0) + 0.184032298439331),
    rhOverM_4_3_3(0.268926437100239*I*delta*(2.0*nu - 1.0)), rhOverM_4_3_5(0.00203732149318363*I*delta*(nu*(524.0*nu -
    1267.0) + 468.0)), rhOverM_4_3_6(0.000332007947037332*delta*(12359.8791491886*nu + 7634.0701482232*I*(2.0*nu - 1.0)
    - 3213.43957459432)), rhOverM_4_4_2(2.25374467927604*nu - 0.751248226425348),
    rhOverM_4_4_4(0.0113825488852325*nu*(525.0*nu - 1273.0) + 4.04991089336574), rhOverM_4_4_5(28.3213909099228*nu +
    0.0187812056606337*I*(-527.578706662453*nu + 114.192902220818) - 9.44046363664094),
    rhOverM_4_4_6(2.91860227826476e-6*nu*(5.0*nu*(678291.0*nu - 3231338.0) + 9793071.0) - 4.0101757911199),
    rhOverM_5_1_3(0.000176960830360287*I*delta*(-2.0*nu + 1.0)), rhOverM_5_1_5(-4.5374571887253e-6*I*delta*(nu*(4.0*nu -
    352.0) + 179.0)), rhOverM_5_1_6(2.52801186228981e-6*delta*(-8958.08121055678*nu - 219.911485751286*I*(2.0*nu - 1.0)
    + 278.040605278392)), rhOverM_5_2_4(0.00998814611056655*nu*(5.0*nu - 5.0) + 0.00998814611056655),
    rhOverM_5_2_6(7.68318931582042e-5*nu*(35.0*nu*(33.0*nu - 118.0) + 3079.0) - 0.0429270763059624),
    rhOverM_5_3_3(0.0464469088412683*I*delta*(2.0*nu - 1.0)), rhOverM_5_3_5(0.00119094638054534*I*delta*(8.0*nu*(11.0*nu
    - 58.0) + 207.0)), rhOverM_5_3_6(9.10188297888856e-7*delta*(923537.386398884*nu + 480946.419338061*I*(2.0*nu - 1.0)
    - 271701.693199442)), rhOverM_5_4_4(-0.276799624590764*nu*(5.0*nu - 5.0) - 0.276799624590764),
    rhOverM_5_4_6(-0.00212922788146741*nu*(5.0*nu*(339.0*nu - 1042.0) + 3619.0) + 1.35388475720164),
    rhOverM_5_5_3(-0.801376894396698*I*delta*(2.0*nu - 1.0)),
    rhOverM_5_5_5(-0.0205481254973512*I*delta*(16.0*nu*(16.0*nu - 43.0) + 263.0)),
    rhOverM_5_5_6(1.83171861576388e-5*delta*(-679921.609610114*nu - 687223.392972767*I*(2.0*nu - 1.0) +
    164747.804805057)), rhOverM_6_1_5(2.35829888333555e-5*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    rhOverM_6_2_4(0.000835250737974468*nu*(5.0*nu - 5.0) + 0.000835250737974468),
    rhOverM_6_2_6(0.000417625368987234*nu*(nu*(7.0*nu - 64.0) + 59.0) - 0.00483252212685228),
    rhOverM_6_3_5(-0.0163097621781264*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), rhOverM_6_4_4(-0.058558165806266*nu*(5.0*nu -
    5.0) - 0.058558165806266), rhOverM_6_4_6(-0.029279082903133*nu*(nu*(19.0*nu - 88.0) + 71.0) + 0.388993529998767),
    rhOverM_6_5_5(0.299357979653564*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), rhOverM_6_6_4(0.903141370807658*nu*(5.0*nu -
    5.0) + 0.903141370807658), rhOverM_6_6_6(0.451570685403829*nu*(nu*(39.0*nu - 128.0) + 91.0) - 7.2896410643761),
    rhOverM_7_1_5(8.17593033339979e-7*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), rhOverM_7_2_6(-0.00134580482328584*nu*pow(nu -
    1.0, 2) + 0.000192257831897977), rhOverM_7_3_5(-0.0018582230503756*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    rhOverM_7_4_6(0.161597034311329*nu*pow(nu - 1.0, 2) - 0.0230852906159042),
    rhOverM_7_5_5(0.0733861624905401*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), rhOverM_7_6_6(-2.3464301109844*nu*pow(nu - 1.0,
    2) + 0.3352043015692), rhOverM_7_7_5(-1.05422444934392*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    rhOverM_8_2_6(-8.42775671401151e-5*nu*pow(nu - 1.0, 2) + 1.20396524485879e-5),
    rhOverM_8_4_6(0.0226281108784145*nu*pow(nu - 1.0, 2) - 0.00323258726834493),
    rhOverM_8_6_6(-0.645290686836342*nu*pow(nu - 1.0, 2) + 0.0921843838337631), rhOverM_8_8_6(8.82604070589592*nu*pow(nu
    - 1.0, 2) - 1.26086295798513), rhOverM_2_1_SO_2(0.5*I*(chi_a_l + chi_s_l*delta)),
    rhOverM_2_1_SO_4(0.0238095238095238*chi_a_l*(-205.0*nu + 7.0) + 0.0238095238095238*chi_s_l*delta*(-33.0*nu + 7.0)),
    rhOverM_2_2_SO_3(-1.33333333333333*chi_a_l*delta + 1.33333333333333*chi_s_l*nu - 1.33333333333333*chi_s_l),
    rhOverM_3_1_SO_4(0.0111358850796843*chi_a_l*(-11.0*nu + 4.0) + 0.0111358850796843*chi_s_l*delta*(-13.0*nu + 4.0)),
    rhOverM_3_2_SO_3(1.12687233963802*chi_s_l*nu), rhOverM_3_3_SO_4(-0.388161877130074*chi_a_l*(19.0*nu - 4.0) -
    0.388161877130074*chi_s_l*delta*(5.0*nu - 4.0)), rhOverM_4_1_SO_4(0.00941154065526303*nu*(chi_a_l - chi_s_l*delta)),
    rhOverM_4_3_SO_4(0.672316092750596*nu*(chi_a_l - chi_s_l*delta)), rhOverM_2_2_SQ_4(2.0*chi1chi2*nu)
  { }

}; // class WaveformModes_3p0PN : public WaveformModes


class WaveformModes_3p5PN : public WaveformModes {
private:
  const double m1, m2;
  double v, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, ellHat_x, ellHat_y, ellHat_z;
  const double m, delta, nu;
  double chi1chi2, chi1_l, chi2_l, chi_s_l, chi_a_l, logv;
  std::complex<double> rhOverM_coeff;
  const std::complex<double> rhOverM_2_0_0, rhOverM_2_1_1, rhOverM_2_1_3, rhOverM_2_1_4, rhOverM_2_1_5, rhOverM_2_1_6,
                             rhOverM_2_2_0, rhOverM_2_2_2, rhOverM_2_2_3, rhOverM_2_2_4, rhOverM_2_2_5, rhOverM_2_2_6,
                             rhOverM_2_2_lnv_6, rhOverM_2_2_7, rhOverM_3_0_5, rhOverM_3_1_1, rhOverM_3_1_3,
                             rhOverM_3_1_4, rhOverM_3_1_5, rhOverM_3_1_6, rhOverM_3_2_2, rhOverM_3_2_4, rhOverM_3_2_5,
                             rhOverM_3_2_6, rhOverM_3_3_1, rhOverM_3_3_3, rhOverM_3_3_4, rhOverM_3_3_5, rhOverM_3_3_6,
                             rhOverM_4_0_0, rhOverM_4_1_3, rhOverM_4_1_5, rhOverM_4_1_6, rhOverM_4_2_2, rhOverM_4_2_4,
                             rhOverM_4_2_5, rhOverM_4_2_6, rhOverM_4_3_3, rhOverM_4_3_5, rhOverM_4_3_6, rhOverM_4_4_2,
                             rhOverM_4_4_4, rhOverM_4_4_5, rhOverM_4_4_6, rhOverM_5_1_3, rhOverM_5_1_5, rhOverM_5_1_6,
                             rhOverM_5_2_4, rhOverM_5_2_6, rhOverM_5_3_3, rhOverM_5_3_5, rhOverM_5_3_6, rhOverM_5_4_4,
                             rhOverM_5_4_6, rhOverM_5_5_3, rhOverM_5_5_5, rhOverM_5_5_6, rhOverM_6_1_5, rhOverM_6_2_4,
                             rhOverM_6_2_6, rhOverM_6_3_5, rhOverM_6_4_4, rhOverM_6_4_6, rhOverM_6_5_5, rhOverM_6_6_4,
                             rhOverM_6_6_6, rhOverM_7_1_5, rhOverM_7_2_6, rhOverM_7_3_5, rhOverM_7_4_6, rhOverM_7_5_5,
                             rhOverM_7_6_6, rhOverM_7_7_5, rhOverM_8_2_6, rhOverM_8_4_6, rhOverM_8_6_6, rhOverM_8_8_6;
  std::complex<double> rhOverM_2_1_SO_2, rhOverM_2_1_SO_4, rhOverM_2_2_SO_3, rhOverM_3_1_SO_4, rhOverM_3_2_SO_3,
                       rhOverM_3_3_SO_4, rhOverM_4_1_SO_4, rhOverM_4_3_SO_4, rhOverM_2_2_SQ_4;

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

    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    logv = log(v);
    rhOverM_coeff = 6.34132367616962*nu*pow(v, 2);
    rhOverM_2_1_SO_2 = 0.5*I*(chi_a_l + chi_s_l*delta);
    rhOverM_2_1_SO_4 = 0.0238095238095238*chi_a_l*(-205.0*nu + 7.0) + 0.0238095238095238*chi_s_l*delta*(-33.0*nu + 7.0);
    rhOverM_2_2_SO_3 = -1.33333333333333*chi_a_l*delta + 1.33333333333333*chi_s_l*nu - 1.33333333333333*chi_s_l;
    rhOverM_3_1_SO_4 = 0.0111358850796843*chi_a_l*(-11.0*nu + 4.0) + 0.0111358850796843*chi_s_l*delta*(-13.0*nu + 4.0);
    rhOverM_3_2_SO_3 = 1.12687233963802*chi_s_l*nu;
    rhOverM_3_3_SO_4 = -0.388161877130074*chi_a_l*(19.0*nu - 4.0) - 0.388161877130074*chi_s_l*delta*(5.0*nu - 4.0);
    rhOverM_4_1_SO_4 = 0.00941154065526303*nu*(chi_a_l - chi_s_l*delta);
    rhOverM_4_3_SO_4 = 0.672316092750596*nu*(chi_a_l - chi_s_l*delta);
    rhOverM_2_2_SQ_4 = 2.0*chi1chi2*nu;

    unsigned int i=0;
    std::vector<std::complex<double> > Modes(77);
    Modes[i++] = (pow(v, 2)*(v*(v*(v*(v*(logv*conjugate(rhOverM_2_2_lnv_6) + v*conjugate(rhOverM_2_2_7) +
      conjugate(rhOverM_2_2_6)) + conjugate(rhOverM_2_2_5)) + conjugate(rhOverM_2_2_4) + conjugate(rhOverM_2_2_SQ_4)) +
      conjugate(rhOverM_2_2_3) + conjugate(rhOverM_2_2_SO_3)) + conjugate(rhOverM_2_2_2)) +
      conjugate(rhOverM_2_2_0))*conjugate(rhOverM_coeff);
    Modes[i++] = v*(v*(v*(v*(v*(v*conjugate(rhOverM_2_1_6) + conjugate(rhOverM_2_1_5)) + conjugate(rhOverM_2_1_4) +
      conjugate(rhOverM_2_1_SO_4)) + conjugate(rhOverM_2_1_3)) + conjugate(rhOverM_2_1_SO_2)) +
      conjugate(rhOverM_2_1_1))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_2_0_0*rhOverM_coeff;
    Modes[i++] = rhOverM_coeff*v*(rhOverM_2_1_1 + v*(rhOverM_2_1_SO_2 + v*(rhOverM_2_1_3 + v*(rhOverM_2_1_4 +
      rhOverM_2_1_SO_4 + v*(rhOverM_2_1_5 + rhOverM_2_1_6*v)))));
    Modes[i++] = rhOverM_coeff*(rhOverM_2_2_0 + pow(v, 2)*(rhOverM_2_2_2 + v*(rhOverM_2_2_3 + rhOverM_2_2_SO_3 +
      v*(rhOverM_2_2_4 + rhOverM_2_2_SQ_4 + v*(rhOverM_2_2_5 + v*(logv*rhOverM_2_2_lnv_6 + rhOverM_2_2_6 +
      rhOverM_2_2_7*v))))));
    Modes[i++] = -v*(pow(v, 2)*(v*(v*(v*conjugate(rhOverM_3_3_6) + conjugate(rhOverM_3_3_5)) + conjugate(rhOverM_3_3_4)
      + conjugate(rhOverM_3_3_SO_4)) + conjugate(rhOverM_3_3_3)) + conjugate(rhOverM_3_3_1))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 2)*(v*(v*(v*(v*conjugate(rhOverM_3_2_6) + conjugate(rhOverM_3_2_5)) + conjugate(rhOverM_3_2_4))
      + conjugate(rhOverM_3_2_SO_3)) + conjugate(rhOverM_3_2_2))*conjugate(rhOverM_coeff);
    Modes[i++] = -v*(pow(v, 2)*(v*(v*(v*conjugate(rhOverM_3_1_6) + conjugate(rhOverM_3_1_5)) + conjugate(rhOverM_3_1_4)
      + conjugate(rhOverM_3_1_SO_4)) + conjugate(rhOverM_3_1_3)) + conjugate(rhOverM_3_1_1))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_3_0_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*v*(rhOverM_3_1_1 + pow(v, 2)*(rhOverM_3_1_3 + v*(rhOverM_3_1_4 + rhOverM_3_1_SO_4 +
      v*(rhOverM_3_1_5 + rhOverM_3_1_6*v))));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_3_2_2 + v*(rhOverM_3_2_SO_3 + v*(rhOverM_3_2_4 + v*(rhOverM_3_2_5 +
      rhOverM_3_2_6*v))));
    Modes[i++] = rhOverM_coeff*v*(rhOverM_3_3_1 + pow(v, 2)*(rhOverM_3_3_3 + v*(rhOverM_3_3_4 + rhOverM_3_3_SO_4 +
      v*(rhOverM_3_3_5 + rhOverM_3_3_6*v))));
    Modes[i++] = pow(v, 2)*(pow(v, 2)*(v*(v*conjugate(rhOverM_4_4_6) + conjugate(rhOverM_4_4_5)) +
      conjugate(rhOverM_4_4_4)) + conjugate(rhOverM_4_4_2))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 3)*(v*(v*(v*conjugate(rhOverM_4_3_6) + conjugate(rhOverM_4_3_5)) + conjugate(rhOverM_4_3_SO_4))
      + conjugate(rhOverM_4_3_3))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 2)*(pow(v, 2)*(v*(v*conjugate(rhOverM_4_2_6) + conjugate(rhOverM_4_2_5)) +
      conjugate(rhOverM_4_2_4)) + conjugate(rhOverM_4_2_2))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 3)*(v*(v*(v*conjugate(rhOverM_4_1_6) + conjugate(rhOverM_4_1_5)) + conjugate(rhOverM_4_1_SO_4))
      + conjugate(rhOverM_4_1_3))*conjugate(rhOverM_coeff);
    Modes[i++] = rhOverM_4_0_0*rhOverM_coeff;
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_4_1_3 + v*(rhOverM_4_1_SO_4 + v*(rhOverM_4_1_5 + rhOverM_4_1_6*v)));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_4_2_2 + pow(v, 2)*(rhOverM_4_2_4 + v*(rhOverM_4_2_5 +
      rhOverM_4_2_6*v)));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_4_3_3 + v*(rhOverM_4_3_SO_4 + v*(rhOverM_4_3_5 + rhOverM_4_3_6*v)));
    Modes[i++] = rhOverM_coeff*pow(v, 2)*(rhOverM_4_4_2 + pow(v, 2)*(rhOverM_4_4_4 + v*(rhOverM_4_4_5 +
      rhOverM_4_4_6*v)));
    Modes[i++] = -pow(v, 3)*(pow(v, 2)*(v*conjugate(rhOverM_5_5_6) + conjugate(rhOverM_5_5_5)) +
      conjugate(rhOverM_5_5_3))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 4)*(pow(v, 2)*conjugate(rhOverM_5_4_6) + conjugate(rhOverM_5_4_4))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 3)*(pow(v, 2)*(v*conjugate(rhOverM_5_3_6) + conjugate(rhOverM_5_3_5)) +
      conjugate(rhOverM_5_3_3))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 4)*(pow(v, 2)*conjugate(rhOverM_5_2_6) + conjugate(rhOverM_5_2_4))*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 3)*(pow(v, 2)*(v*conjugate(rhOverM_5_1_6) + conjugate(rhOverM_5_1_5)) +
      conjugate(rhOverM_5_1_3))*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_5_1_3 + pow(v, 2)*(rhOverM_5_1_5 + rhOverM_5_1_6*v));
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(rhOverM_5_2_4 + rhOverM_5_2_6*pow(v, 2));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_5_3_3 + pow(v, 2)*(rhOverM_5_3_5 + rhOverM_5_3_6*v));
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(rhOverM_5_4_4 + rhOverM_5_4_6*pow(v, 2));
    Modes[i++] = rhOverM_coeff*pow(v, 3)*(rhOverM_5_5_3 + pow(v, 2)*(rhOverM_5_5_5 + rhOverM_5_5_6*v));
    Modes[i++] = pow(v, 4)*(pow(v, 2)*conjugate(rhOverM_6_6_6) + conjugate(rhOverM_6_6_4))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 5)*conjugate(rhOverM_6_5_5)*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 4)*(pow(v, 2)*conjugate(rhOverM_6_4_6) + conjugate(rhOverM_6_4_4))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 5)*conjugate(rhOverM_6_3_5)*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 4)*(pow(v, 2)*conjugate(rhOverM_6_2_6) + conjugate(rhOverM_6_2_4))*conjugate(rhOverM_coeff);
    Modes[i++] = pow(v, 5)*conjugate(rhOverM_6_1_5)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_6_1_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(rhOverM_6_2_4 + rhOverM_6_2_6*pow(v, 2));
    Modes[i++] = rhOverM_6_3_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(rhOverM_6_4_4 + rhOverM_6_4_6*pow(v, 2));
    Modes[i++] = rhOverM_6_5_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_coeff*pow(v, 4)*(rhOverM_6_6_4 + rhOverM_6_6_6*pow(v, 2));
    Modes[i++] = -pow(v, 5)*conjugate(rhOverM_7_7_5)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 6)*conjugate(rhOverM_7_6_6)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 5)*conjugate(rhOverM_7_5_5)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 6)*conjugate(rhOverM_7_4_6)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 5)*conjugate(rhOverM_7_3_5)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 6)*conjugate(rhOverM_7_2_6)*conjugate(rhOverM_coeff);
    Modes[i++] = -pow(v, 5)*conjugate(rhOverM_7_1_5)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_7_1_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_7_2_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = rhOverM_7_3_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_7_4_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = rhOverM_7_5_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = rhOverM_7_6_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = rhOverM_7_7_5*rhOverM_coeff*pow(v, 5);
    Modes[i++] = pow(v, 6)*conjugate(rhOverM_8_8_6)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = pow(v, 6)*conjugate(rhOverM_8_6_6)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = pow(v, 6)*conjugate(rhOverM_8_4_6)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = pow(v, 6)*conjugate(rhOverM_8_2_6)*conjugate(rhOverM_coeff);
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = 0;
    Modes[i++] = rhOverM_8_2_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_8_4_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_8_6_6*rhOverM_coeff*pow(v, 6);
    Modes[i++] = 0;
    Modes[i++] = rhOverM_8_8_6*rhOverM_coeff*pow(v, 6);

    return Modes;
  }

public:
  WaveformModes_3p5PN(const double m1_i, const double m2_i, const double v_i, const double chi1_x_i, const double
                      chi1_y_i, const double chi1_z_i, const double chi2_x_i, const double chi2_y_i, const double
                      chi2_z_i, const double ellHat_x_i, const double ellHat_y_i, const double ellHat_z_i) :
    m1(m1_i), m2(m2_i), v(v_i), chi1_x(chi1_x_i), chi1_y(chi1_y_i), chi1_z(chi1_z_i), chi2_x(chi2_x_i),
    chi2_y(chi2_y_i), chi2_z(chi2_z_i), ellHat_x(ellHat_x_i), ellHat_y(ellHat_y_i), ellHat_z(ellHat_z_i), m(m1 + m2),
    delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z),
    chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y +
    chi2_z*ellHat_z), chi_s_l(0.5*chi1_l + 0.5*chi2_l), chi_a_l(0.5*chi1_l - 0.5*chi2_l), logv(log(v)),
    rhOverM_coeff(6.34132367616962*nu*pow(v, 2)), rhOverM_2_0_0(-0.145802960879951),
    rhOverM_2_1_1(0.333333333333333*I*delta), rhOverM_2_1_3(0.0119047619047619*I*delta*(20.0*nu - 17.0)),
    rhOverM_2_1_4(delta*(0.628764787039964 + 1.0471975511966*I)),
    rhOverM_2_1_5(0.000661375661375661*I*delta*(nu*(237.0*nu - 2036.0) - 172.0)),
    rhOverM_2_1_6(0.00595238095238095*delta*(nu*(722.635532333439 + 37.6991118430775*I) - 64.1340082780763 -
    106.814150222053*I)), rhOverM_2_2_0(1.00000000000000), rhOverM_2_2_2(1.30952380952381*nu - 2.54761904761905),
    rhOverM_2_2_3(6.28318530717959), rhOverM_2_2_4(0.000661375661375661*nu*(2047.0*nu - 7483.0) - 1.43716931216931),
    rhOverM_2_2_5(5.08638810581205*nu - 24.0*I*nu - 16.0071625682908),
    rhOverM_2_2_6(1.00208433541767e-5*nu*(nu*(114635.0*nu - 729396.0) - 834555.0) + 4.21514354629858*nu +
    32.3588011610077 + 12.8057300546327*I), rhOverM_2_2_lnv_6(-8.15238095238095),
    rhOverM_2_2_7(0.00831109167616348*nu*(560.0*nu - 2459.0) - 0.00017636684303351*I*nu*(24396.0*nu - 501655.0) -
    9.03000110615162), rhOverM_3_0_5(-0.370328039909021*I*nu), rhOverM_3_1_1(0.0222717701593687*I*delta),
    rhOverM_3_1_3(-0.0148478467729125*I*delta*(nu + 4.0)), rhOverM_3_1_4(delta*(0.0620557076072073 +
    0.0699688295151131*I)), rhOverM_3_1_5(-0.000112483687673579*I*delta*(nu*(247.0*nu + 272.0) - 607.0)),
    rhOverM_3_1_6(0.000742392338645623*delta*(-46.5203026391962*nu - 15.707963267949*I*(7.0*nu + 16.0) -
    222.903548889591)), rhOverM_3_2_2(-0.845154254728516*nu + 0.281718084909506),
    rhOverM_3_2_4(0.00313020094343895*nu*(-365.0*nu + 725.0) - 0.604128782083717), rhOverM_3_2_5(-5.31026079561053*nu +
    3.71867872080547*I*nu + 1.77008693187018 - 0.845154254728517*I),
    rhOverM_3_2_6(7.11409305327034e-5*nu*(nu*(-16023.0*nu + 100026.0) - 17387.0) - 0.103225490202953),
    rhOverM_3_3_1(-0.776323754260148*I*delta), rhOverM_3_3_3(-1.5526475085203*I*delta*(nu - 2.0)),
    rhOverM_3_3_4(delta*(-1.37192659820446 - 7.31667900957279*I)),
    rhOverM_3_3_5(-0.00235249622503075*I*delta*(nu*(887.0*nu - 3676.0) + 369.0)),
    rhOverM_3_3_6(0.000319474795991831*delta*(-87338.4780856744*nu - 11451.1052223348*I*(3.0*nu - 8.0) +
    17177.2748951318)), rhOverM_4_0_0(-0.00140298964521140), rhOverM_4_1_3(0.00376461626210521*I*delta*(-2.0*nu + 1.0)),
    rhOverM_4_1_5(-2.85198201674637e-5*I*delta*(nu*(332.0*nu - 1011.0) + 404.0)),
    rhOverM_4_1_6(0.00012548720873684*delta*(-1744.17766166719*nu - 94.2477796076938*I*(2.0*nu - 1.0) +
    105.588830833597)), rhOverM_4_2_2(-0.10647942749999*nu + 0.0354931424999967),
    rhOverM_4_2_4(0.000107554977272717*nu*(-285.0*nu + 4025.0) - 0.141004575204532),
    rhOverM_4_2_5(0.00709862849999933*nu*(-94.2477796076938 + 84.0*I) + 0.22300999146161 - 0.149071198499986*I),
    rhOverM_4_2_6(1.37890996503484e-7*nu*(115.0*nu*(3363.0*nu + 34822.0) - 5460759.0) + 0.184032298439331),
    rhOverM_4_3_3(0.268926437100239*I*delta*(2.0*nu - 1.0)), rhOverM_4_3_5(0.00203732149318363*I*delta*(nu*(524.0*nu -
    1267.0) + 468.0)), rhOverM_4_3_6(0.000332007947037332*delta*(12359.8791491886*nu + 7634.0701482232*I*(2.0*nu - 1.0)
    - 3213.43957459432)), rhOverM_4_4_2(2.25374467927604*nu - 0.751248226425348),
    rhOverM_4_4_4(0.0113825488852325*nu*(525.0*nu - 1273.0) + 4.04991089336574), rhOverM_4_4_5(28.3213909099228*nu +
    0.0187812056606337*I*(-527.578706662453*nu + 114.192902220818) - 9.44046363664094),
    rhOverM_4_4_6(2.91860227826476e-6*nu*(5.0*nu*(678291.0*nu - 3231338.0) + 9793071.0) - 4.0101757911199),
    rhOverM_5_1_3(0.000176960830360287*I*delta*(-2.0*nu + 1.0)), rhOverM_5_1_5(-4.5374571887253e-6*I*delta*(nu*(4.0*nu -
    352.0) + 179.0)), rhOverM_5_1_6(2.52801186228981e-6*delta*(-8958.08121055678*nu - 219.911485751286*I*(2.0*nu - 1.0)
    + 278.040605278392)), rhOverM_5_2_4(0.00998814611056655*nu*(5.0*nu - 5.0) + 0.00998814611056655),
    rhOverM_5_2_6(7.68318931582042e-5*nu*(35.0*nu*(33.0*nu - 118.0) + 3079.0) - 0.0429270763059624),
    rhOverM_5_3_3(0.0464469088412683*I*delta*(2.0*nu - 1.0)), rhOverM_5_3_5(0.00119094638054534*I*delta*(8.0*nu*(11.0*nu
    - 58.0) + 207.0)), rhOverM_5_3_6(9.10188297888856e-7*delta*(923537.386398884*nu + 480946.419338061*I*(2.0*nu - 1.0)
    - 271701.693199442)), rhOverM_5_4_4(-0.276799624590764*nu*(5.0*nu - 5.0) - 0.276799624590764),
    rhOverM_5_4_6(-0.00212922788146741*nu*(5.0*nu*(339.0*nu - 1042.0) + 3619.0) + 1.35388475720164),
    rhOverM_5_5_3(-0.801376894396698*I*delta*(2.0*nu - 1.0)),
    rhOverM_5_5_5(-0.0205481254973512*I*delta*(16.0*nu*(16.0*nu - 43.0) + 263.0)),
    rhOverM_5_5_6(1.83171861576388e-5*delta*(-679921.609610114*nu - 687223.392972767*I*(2.0*nu - 1.0) +
    164747.804805057)), rhOverM_6_1_5(2.35829888333555e-5*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    rhOverM_6_2_4(0.000835250737974468*nu*(5.0*nu - 5.0) + 0.000835250737974468),
    rhOverM_6_2_6(0.000417625368987234*nu*(nu*(7.0*nu - 64.0) + 59.0) - 0.00483252212685228),
    rhOverM_6_3_5(-0.0163097621781264*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), rhOverM_6_4_4(-0.058558165806266*nu*(5.0*nu -
    5.0) - 0.058558165806266), rhOverM_6_4_6(-0.029279082903133*nu*(nu*(19.0*nu - 88.0) + 71.0) + 0.388993529998767),
    rhOverM_6_5_5(0.299357979653564*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), rhOverM_6_6_4(0.903141370807658*nu*(5.0*nu -
    5.0) + 0.903141370807658), rhOverM_6_6_6(0.451570685403829*nu*(nu*(39.0*nu - 128.0) + 91.0) - 7.2896410643761),
    rhOverM_7_1_5(8.17593033339979e-7*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), rhOverM_7_2_6(-0.00134580482328584*nu*pow(nu -
    1.0, 2) + 0.000192257831897977), rhOverM_7_3_5(-0.0018582230503756*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    rhOverM_7_4_6(0.161597034311329*nu*pow(nu - 1.0, 2) - 0.0230852906159042),
    rhOverM_7_5_5(0.0733861624905401*I*delta*(nu - 1.0)*(3.0*nu - 1.0)), rhOverM_7_6_6(-2.3464301109844*nu*pow(nu - 1.0,
    2) + 0.3352043015692), rhOverM_7_7_5(-1.05422444934392*I*delta*(nu - 1.0)*(3.0*nu - 1.0)),
    rhOverM_8_2_6(-8.42775671401151e-5*nu*pow(nu - 1.0, 2) + 1.20396524485879e-5),
    rhOverM_8_4_6(0.0226281108784145*nu*pow(nu - 1.0, 2) - 0.00323258726834493),
    rhOverM_8_6_6(-0.645290686836342*nu*pow(nu - 1.0, 2) + 0.0921843838337631), rhOverM_8_8_6(8.82604070589592*nu*pow(nu
    - 1.0, 2) - 1.26086295798513), rhOverM_2_1_SO_2(0.5*I*(chi_a_l + chi_s_l*delta)),
    rhOverM_2_1_SO_4(0.0238095238095238*chi_a_l*(-205.0*nu + 7.0) + 0.0238095238095238*chi_s_l*delta*(-33.0*nu + 7.0)),
    rhOverM_2_2_SO_3(-1.33333333333333*chi_a_l*delta + 1.33333333333333*chi_s_l*nu - 1.33333333333333*chi_s_l),
    rhOverM_3_1_SO_4(0.0111358850796843*chi_a_l*(-11.0*nu + 4.0) + 0.0111358850796843*chi_s_l*delta*(-13.0*nu + 4.0)),
    rhOverM_3_2_SO_3(1.12687233963802*chi_s_l*nu), rhOverM_3_3_SO_4(-0.388161877130074*chi_a_l*(19.0*nu - 4.0) -
    0.388161877130074*chi_s_l*delta*(5.0*nu - 4.0)), rhOverM_4_1_SO_4(0.00941154065526303*nu*(chi_a_l - chi_s_l*delta)),
    rhOverM_4_3_SO_4(0.672316092750596*nu*(chi_a_l - chi_s_l*delta)), rhOverM_2_2_SQ_4(2.0*chi1chi2*nu)
  { }

}; // class WaveformModes_3p5PN : public WaveformModes
