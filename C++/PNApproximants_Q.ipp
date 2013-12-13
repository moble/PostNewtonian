// File produced automatically by CodeOutput_Q.ipynb

class TaylorTn_0PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, Fcal_coeff;
  const double Fcal_0, E_0;
  double Phi;

public:
  TaylorTn_0PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    E_0(1.00000000000000), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*ellHat*(-0.75*delta + 0.5*nu + 0.75);
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*ellHat*(0.75*delta + 0.5*nu + 0.75);
  }
  Quaternion OmegaVec() {
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_0 = 1.00000000000000;
    return a_ell_0*gamma_PN_0*nHat*pow(v, 6)/pow(m, 3) + ellHat*pow(v, 3)/m;
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_0*Fcal_coeff;
    const double dEdv = -E_0*m*nu*v;
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = 1.0*Fcal_0*Fcal_coeff/(E_0*m*nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = 1.0*E_0*m*nu*v/(Fcal_0*Fcal_coeff);
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_0PN_Q : public TaylorTn_Q


class TaylorTn_0p50PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, Fcal_coeff;
  const double Fcal_0, E_0;
  double Phi;

public:
  TaylorTn_0p50PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    E_0(1.00000000000000), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*(-chiVec2*pow(m2, 2)*v/pow(m, 2) + ellHat*(-0.75*delta + 0.5*nu + 0.75) + nHat*v*(3.0*chi1_n*nu
      + 3.0*chi2_n*pow(m2, 2)/pow(m, 2)));
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*(-chiVec1*pow(m1, 2)*v/pow(m, 2) + ellHat*(0.75*delta + 0.5*nu + 0.75) +
      nHat*v*(3.0*chi1_n*pow(m1, 2)/pow(m, 2) + 3.0*chi2_n*nu));
  }
  Quaternion OmegaVec() {
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_0 = 1.00000000000000;
    return a_ell_0*gamma_PN_0*nHat*pow(v, 6)/pow(m, 3) + ellHat*pow(v, 3)/m;
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_0*Fcal_coeff;
    const double dEdv = -E_0*m*nu*v;
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = 1.0*Fcal_0*Fcal_coeff/(E_0*m*nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = 1.0*E_0*m*nu*v/(Fcal_0*Fcal_coeff);
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_0p50PN_Q : public TaylorTn_Q


class TaylorTn_1p0PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, Fcal_coeff;
  const double Fcal_0, Fcal_2, E_0, E_2;
  double Phi;

public:
  TaylorTn_1p0PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    Fcal_2(-2.91666666666667*nu - 3.71130952380952), E_0(1.00000000000000), E_2(-0.0833333333333333*nu - 0.75), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*(-chiVec2*pow(m2, 2)*v/pow(m, 2) + ellHat*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu -
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*chi2_n*pow(m2,
      2)/pow(m, 2)));
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*(-chiVec1*pow(m1, 2)*v/pow(m, 2) + ellHat*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu +
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*pow(m1, 2)/pow(m, 2) +
      3.0*chi2_n*nu));
  }
  Quaternion OmegaVec() {
    double gamma_PN_2 = -0.333333333333333*nu + 1.0;
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_0 = 1.00000000000000;
    double a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0);
    return ellHat*pow(v, 3)/m + nHat*pow(v, 6)*(a_ell_0 + a_ell_2*pow(v, 2))*(gamma_PN_0 + gamma_PN_2*pow(v, 2))/pow(m,
      3);
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_coeff*(Fcal_0 + Fcal_2*pow(v, 2));
    const double dEdv = -m*nu*v*(E_0 + 2.0*E_2*pow(v, 2));
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = -2.0*Fcal_coeff*(-0.5*Fcal_0/(E_0*m) + pow(v, 2)*(-0.5*Fcal_2/m +
      1.0*E_2*Fcal_0/(E_0*m))/E_0)/(nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = -0.5*nu*v*(-2.0*E_0*m/Fcal_0 + m*pow(v, 2)*(2.0*E_0*Fcal_2/Fcal_0 - 4.0*E_2)/Fcal_0)/Fcal_coeff;
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_1p0PN_Q : public TaylorTn_Q


class TaylorTn_1p5PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, Fcal_coeff;
  const double Fcal_0, Fcal_2, Fcal_3;
  double Fcal_SO_3;
  const double E_0, E_2;
  double E_SO_3;
  double Phi;

public:
  TaylorTn_1p5PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    Fcal_2(-2.91666666666667*nu - 3.71130952380952), Fcal_3(12.5663706143592), Fcal_SO_3((-4.0*S_l -
    1.25*Sigma_l*delta)/pow(m, 2)), E_0(1.00000000000000), E_2(-0.0833333333333333*nu - 0.75),
    E_SO_3((4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2)), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
    Fcal_SO_3 = (-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2);
    E_SO_3 = (4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*(-chiVec2*pow(m2, 2)*v/pow(m, 2) + ellHat*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu -
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*chi2_n*pow(m2,
      2)/pow(m, 2)));
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*(-chiVec1*pow(m1, 2)*v/pow(m, 2) + ellHat*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu +
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*pow(m1, 2)/pow(m, 2) +
      3.0*chi2_n*nu));
  }
  Quaternion OmegaVec() {
    double gamma_PN_2 = -0.333333333333333*nu + 1.0;
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_3 = (1.66666666666667*S_l + Sigma_l*delta)/pow(m, 2);
    double gamma_PN_0 = 1.00000000000000;
    double a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0);
    return ellHat*pow(v, 3)/m + nHat*pow(v, 6)*(a_ell_0 + a_ell_2*pow(v, 2))*(gamma_PN_0 + pow(v, 2)*(gamma_PN_2 +
      gamma_PN_3*v))/pow(m, 3);
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_coeff*(Fcal_0 + pow(v, 2)*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3)));
    const double dEdv = -0.5*m*nu*v*(2.0*E_0 + pow(v, 2)*(4.0*E_2 + 5.0*E_SO_3*v));
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = -2.0*Fcal_coeff*(pow(v, 2)*(v*((-0.5*Fcal_3 - 0.5*Fcal_SO_3)/m +
      1.25*E_SO_3*Fcal_0/(E_0*m))/E_0 + (-0.5*Fcal_2/m + 1.0*E_2*Fcal_0/(E_0*m))/E_0) - 0.5*Fcal_0/(E_0*m))/(nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = -0.5*nu*v*(-2.0*E_0*m/Fcal_0 + pow(v, 2)*(m*v*(E_0*(2.0*Fcal_3 + 2.0*Fcal_SO_3)/Fcal_0 -
      5.0*E_SO_3)/Fcal_0 + m*(2.0*E_0*Fcal_2/Fcal_0 - 4.0*E_2)/Fcal_0))/Fcal_coeff;
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_1p5PN_Q : public TaylorTn_Q


class TaylorTn_2p0PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, Fcal_coeff;
  const double Fcal_0, Fcal_2, Fcal_3, Fcal_4;
  double Fcal_SQ_4, Fcal_SO_3;
  const double E_0, E_2, E_4;
  double E_SQ_4, E_SO_3;
  double Phi;

public:
  TaylorTn_2p0PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    Fcal_2(-2.91666666666667*nu - 3.71130952380952), Fcal_3(12.5663706143592), Fcal_4(3.61111111111111*pow(nu, 2) +
    18.3948412698413*nu - 4.92846119929453), Fcal_SQ_4(chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu -
    0.463541666666667) - 2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu -
    0.463541666666667) + chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
    2)*(0.0416666666666667*nu + 2.98958333333333)), Fcal_SO_3((-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2)),
    E_0(1.00000000000000), E_2(-0.0833333333333333*nu - 0.75), E_4(-0.0416666666666667*pow(nu, 2) + 2.375*nu - 3.375),
    E_SQ_4(-1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 +
    6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0)), E_SO_3((4.66666666666667*S_l +
    2.0*Sigma_l*delta)/pow(m, 2)), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
    Fcal_SQ_4 = chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
      2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
      chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
      2)*(0.0416666666666667*nu + 2.98958333333333);
    Fcal_SO_3 = (-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2);
    E_SQ_4 = -1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 +
      6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0);
    E_SO_3 = (4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*(-chiVec2*pow(m2, 2)*v/pow(m, 2) + ellHat*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu -
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu +
      3.0*chi2_n*pow(m2, 2)/pow(m, 2)));
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*(-chiVec1*pow(m1, 2)*v/pow(m, 2) + ellHat*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu +
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*pow(m1,
      2)/pow(m, 2) + 3.0*chi2_n*nu));
  }
  Quaternion OmegaVec() {
    double gamma_PN_2 = -0.333333333333333*nu + 1.0;
    double a_ell_4 = S_n*(5.77777777777778*pow(nu, 2) + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*pow(nu, 2) +
      9.125*nu + 1.5);
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_3 = (1.66666666666667*S_l + Sigma_l*delta)/pow(m, 2);
    double gamma_PN_0 = 1.00000000000000;
    double a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0);
    double gamma_PN_4 = -5.41666666666667*nu + 1.0;
    return ellHat*pow(v, 3)/m + nHat*pow(v, 6)*(a_ell_0 + pow(v, 2)*(a_ell_2 + a_ell_4*pow(v, 2)))*(gamma_PN_0 + pow(v,
      2)*(gamma_PN_2 + v*(gamma_PN_3 + gamma_PN_4*v)))/pow(m, 3);
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_coeff*(Fcal_0 + pow(v, 2)*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4))));
    const double dEdv = -0.5*m*nu*v*(2.0*E_0 + pow(v, 2)*(4.0*E_2 + v*(5.0*E_SO_3 + 6.0*v*(E_4 + E_SQ_4))));
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = -2.0*Fcal_coeff*(pow(v, 2)*(v*(v*((-0.5*Fcal_4 - 0.5*Fcal_SQ_4)/m + ((1.0*E_2*Fcal_2 +
      1.5*E_4*Fcal_0 + 1.5*E_SQ_4*Fcal_0)/m - 2.0*pow(E_2, 2)*Fcal_0/(E_0*m))/E_0)/E_0 + ((-0.5*Fcal_3 -
      0.5*Fcal_SO_3)/m + 1.25*E_SO_3*Fcal_0/(E_0*m))/E_0) + (-0.5*Fcal_2/m + 1.0*E_2*Fcal_0/(E_0*m))/E_0) -
      0.5*Fcal_0/(E_0*m))/(nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = -0.5*nu*v*(-2.0*E_0*m/Fcal_0 + pow(v, 2)*(v*(m*v*(-6.0*E_4 - 6.0*E_SQ_4 + (E_0*(2.0*Fcal_4 +
      2.0*Fcal_SQ_4) - 2.0*E_0*pow(Fcal_2, 2)/Fcal_0 + 4.0*E_2*Fcal_2)/Fcal_0)/Fcal_0 + m*(E_0*(2.0*Fcal_3 +
      2.0*Fcal_SO_3)/Fcal_0 - 5.0*E_SO_3)/Fcal_0) + m*(2.0*E_0*Fcal_2/Fcal_0 - 4.0*E_2)/Fcal_0))/Fcal_coeff;
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_2p0PN_Q : public TaylorTn_Q


class TaylorTn_2p5PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, Fcal_coeff;
  const double Fcal_0, Fcal_2, Fcal_3, Fcal_4, Fcal_5;
  double Fcal_SQ_4, Fcal_SO_3, Fcal_SO_5;
  const double E_0, E_2, E_4;
  double E_SQ_4, E_SO_3, E_SO_5;
  double Phi;

public:
  TaylorTn_2p5PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    Fcal_2(-2.91666666666667*nu - 3.71130952380952), Fcal_3(12.5663706143592), Fcal_4(3.61111111111111*pow(nu, 2) +
    18.3948412698413*nu - 4.92846119929453), Fcal_5(-76.3145215434521*nu - 38.2928354546934),
    Fcal_SQ_4(chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
    2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
    chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
    2)*(0.0416666666666667*nu + 2.98958333333333)), Fcal_SO_3((-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2)),
    Fcal_SO_5((S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2)), E_0(1.00000000000000),
    E_2(-0.0833333333333333*nu - 0.75), E_4(-0.0416666666666667*pow(nu, 2) + 2.375*nu - 3.375), E_SQ_4(-1.5*pow(chi_a_l,
    2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 + 6.0*pow(chi_a_l, 2)) +
    0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0)), E_SO_3((4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2)),
    E_SO_5((S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2)), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
    Fcal_SQ_4 = chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
      2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
      chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
      2)*(0.0416666666666667*nu + 2.98958333333333);
    Fcal_SO_3 = (-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_5 = (S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2);
    E_SQ_4 = -1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 +
      6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0);
    E_SO_3 = (4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2);
    E_SO_5 = (S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*(-chiVec2*pow(m2, 2)*v/pow(m, 2) + ellHat*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu -
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu +
      3.0*chi2_n*pow(m2, 2)/pow(m, 2)));
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*(-chiVec1*pow(m1, 2)*v/pow(m, 2) + ellHat*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu +
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*pow(m1,
      2)/pow(m, 2) + 3.0*chi2_n*nu));
  }
  Quaternion OmegaVec() {
    double gamma_PN_2 = -0.333333333333333*nu + 1.0;
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_3 = (1.66666666666667*S_l + Sigma_l*delta)/pow(m, 2);
    double gamma_PN_5 = (S_l*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_l*delta)/pow(m, 2);
    double gamma_PN_0 = 1.00000000000000;
    double a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0);
    double gamma_PN_4 = -5.41666666666667*nu + 1.0;
    double a_ell_4 = S_n*(5.77777777777778*pow(nu, 2) + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*pow(nu, 2) +
      9.125*nu + 1.5);
    return ellHat*pow(v, 3)/m + nHat*pow(v, 6)*(a_ell_0 + pow(v, 2)*(a_ell_2 + a_ell_4*pow(v, 2)))*(gamma_PN_0 + pow(v,
      2)*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + gamma_PN_5*v))))/pow(m, 3);
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_coeff*(Fcal_0 + pow(v, 2)*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 +
      v*(Fcal_5 + Fcal_SO_5)))));
    const double dEdv = -0.5*m*nu*v*(2.0*E_0 + pow(v, 2)*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 7.0*E_SO_5*v +
      6.0*E_SQ_4))));
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = -2.0*Fcal_coeff*(pow(v, 2)*(v*(v*(v*((-0.5*Fcal_5 - 0.5*Fcal_SO_5)/m + ((E_2*(1.0*Fcal_3 +
      1.0*Fcal_SO_3) + 1.25*E_SO_3*Fcal_2 + 1.75*E_SO_5*Fcal_0)/m - 5.0*E_2*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0 +
      ((-0.5*Fcal_4 - 0.5*Fcal_SQ_4)/m + ((1.0*E_2*Fcal_2 + 1.5*E_4*Fcal_0 + 1.5*E_SQ_4*Fcal_0)/m - 2.0*pow(E_2,
      2)*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_3 - 0.5*Fcal_SO_3)/m + 1.25*E_SO_3*Fcal_0/(E_0*m))/E_0) +
      (-0.5*Fcal_2/m + 1.0*E_2*Fcal_0/(E_0*m))/E_0) - 0.5*Fcal_0/(E_0*m))/(nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = -0.5*nu*v*(-2.0*E_0*m/Fcal_0 + pow(v, 2)*(v*(v*(m*v*(-7.0*E_SO_5 + (E_0*(2.0*Fcal_5 +
      2.0*Fcal_SO_5) + E_0*Fcal_2*(-4.0*Fcal_3 - 4.0*Fcal_SO_3)/Fcal_0 + E_2*(4.0*Fcal_3 + 4.0*Fcal_SO_3) +
      5.0*E_SO_3*Fcal_2)/Fcal_0)/Fcal_0 + m*(-6.0*E_4 - 6.0*E_SQ_4 + (E_0*(2.0*Fcal_4 + 2.0*Fcal_SQ_4) -
      2.0*E_0*pow(Fcal_2, 2)/Fcal_0 + 4.0*E_2*Fcal_2)/Fcal_0)/Fcal_0) + m*(E_0*(2.0*Fcal_3 + 2.0*Fcal_SO_3)/Fcal_0 -
      5.0*E_SO_3)/Fcal_0) + m*(2.0*E_0*Fcal_2/Fcal_0 - 4.0*E_2)/Fcal_0))/Fcal_coeff;
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_2p5PN_Q : public TaylorTn_Q


class TaylorTn_3p0PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, logv, Fcal_coeff;
  const double Fcal_0, Fcal_2, Fcal_3, Fcal_4, Fcal_5, Fcal_6, Fcal_lnv_6;
  double Fcal_SQ_4, Fcal_SO_3, Fcal_SO_5, Fcal_SO_6;
  const double E_0, E_2, E_4, E_6;
  double E_SQ_4, E_SO_3, E_SO_5;
  double Phi;

public:
  TaylorTn_3p0PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), logv(log(v)), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    Fcal_2(-2.91666666666667*nu - 3.71130952380952), Fcal_3(12.5663706143592), Fcal_4(3.61111111111111*pow(nu, 2) +
    18.3948412698413*nu - 4.92846119929453), Fcal_5(-76.3145215434521*nu - 38.2928354546934),
    Fcal_6(-2.39197530864198*pow(nu, 3) - 31.2179232804233*pow(nu, 2) - 8.87205344238227*nu + 115.731716675611),
    Fcal_lnv_6(-16.3047619047619), Fcal_SQ_4(chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu -
    0.463541666666667) - 2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu -
    0.463541666666667) + chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
    2)*(0.0416666666666667*nu + 2.98958333333333)), Fcal_SO_3((-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2)),
    Fcal_SO_5((S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2)),
    Fcal_SO_6((-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2)), E_0(1.00000000000000),
    E_2(-0.0833333333333333*nu - 0.75), E_4(-0.0416666666666667*pow(nu, 2) + 2.375*nu - 3.375),
    E_6(-0.00675154320987654*pow(nu, 3) - 1.61458333333333*pow(nu, 2) + 38.7246294907293*nu - 10.546875),
    E_SQ_4(-1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 +
    6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0)), E_SO_3((4.66666666666667*S_l +
    2.0*Sigma_l*delta)/pow(m, 2)), E_SO_5((S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu +
    3.0))/pow(m, 2)), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    logv = log(v);
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
    Fcal_SQ_4 = chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
      2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
      chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
      2)*(0.0416666666666667*nu + 2.98958333333333);
    Fcal_SO_3 = (-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_5 = (S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2);
    Fcal_SO_6 = (-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2);
    E_SQ_4 = -1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 +
      6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0);
    E_SO_3 = (4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2);
    E_SO_5 = (S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*(-chiVec2*pow(m2, 2)*v/pow(m, 2) + ellHat*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu -
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu +
      3.0*chi2_n*pow(m2, 2)/pow(m, 2)));
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*(-chiVec1*pow(m1, 2)*v/pow(m, 2) + ellHat*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu +
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*pow(m1,
      2)/pow(m, 2) + 3.0*chi2_n*nu));
  }
  Quaternion OmegaVec() {
    double gamma_PN_2 = -0.333333333333333*nu + 1.0;
    double gamma_PN_6 = 0.0123456790123457*pow(nu, 3) + 6.36111111111111*pow(nu, 2) - 2.98177812235564*nu + 1.0;
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_3 = (1.66666666666667*S_l + Sigma_l*delta)/pow(m, 2);
    double gamma_PN_5 = (S_l*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_l*delta)/pow(m, 2);
    double gamma_PN_0 = 1.00000000000000;
    double a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0);
    double gamma_PN_4 = -5.41666666666667*nu + 1.0;
    double a_ell_4 = S_n*(5.77777777777778*pow(nu, 2) + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*pow(nu, 2) +
      9.125*nu + 1.5);
    return ellHat*pow(v, 3)/m + nHat*pow(v, 6)*(a_ell_0 + pow(v, 2)*(a_ell_2 + a_ell_4*pow(v, 2)))*(gamma_PN_0 + pow(v,
      2)*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + gamma_PN_6*v)))))/pow(m, 3);
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_coeff*(Fcal_0 + pow(v, 2)*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 +
      v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv))))));
    const double dEdv = -0.5*m*nu*v*(2.0*E_0 + pow(v, 2)*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 +
      v*(8.0*E_6*v + 7.0*E_SO_5)))));
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = -2.0*Fcal_coeff*(pow(v, 2)*(v*(v*(v*(v*((-0.5*Fcal_6 - 0.5*Fcal_SO_6 - 0.5*Fcal_lnv_6*logv)/m
      + ((E_2*(1.0*Fcal_4 + 1.0*Fcal_SQ_4) + 1.5*E_4*Fcal_2 + 2.0*E_6*Fcal_0 + E_SO_3*(1.25*Fcal_3 + 1.25*Fcal_SO_3) +
      1.5*E_SQ_4*Fcal_2)/m + ((E_2*(-2.0*E_2*Fcal_2 - 6.0*E_4*Fcal_0 - 6.0*E_SQ_4*Fcal_0) - 3.125*pow(E_SO_3,
      2)*Fcal_0)/m + 4.0*pow(E_2, 3)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0 + ((-0.5*Fcal_5 - 0.5*Fcal_SO_5)/m +
      ((E_2*(1.0*Fcal_3 + 1.0*Fcal_SO_3) + 1.25*E_SO_3*Fcal_2 + 1.75*E_SO_5*Fcal_0)/m -
      5.0*E_2*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_4 - 0.5*Fcal_SQ_4)/m + ((1.0*E_2*Fcal_2 + 1.5*E_4*Fcal_0 +
      1.5*E_SQ_4*Fcal_0)/m - 2.0*pow(E_2, 2)*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_3 - 0.5*Fcal_SO_3)/m +
      1.25*E_SO_3*Fcal_0/(E_0*m))/E_0) + (-0.5*Fcal_2/m + 1.0*E_2*Fcal_0/(E_0*m))/E_0) - 0.5*Fcal_0/(E_0*m))/(nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = -0.5*nu*v*(-2.0*E_0*m/Fcal_0 + pow(v, 2)*(v*(v*(v*(m*v*(-8.0*E_6 + (E_0*(2.0*Fcal_6 +
      2.0*Fcal_SO_6 + 2.0*Fcal_lnv_6*logv) + E_2*(4.0*Fcal_4 + 4.0*Fcal_SQ_4) + 6.0*E_4*Fcal_2 + E_SO_3*(5.0*Fcal_3 +
      5.0*Fcal_SO_3) + 6.0*E_SQ_4*Fcal_2 + (E_0*(Fcal_2*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) + Fcal_3*(-2.0*Fcal_3 -
      4.0*Fcal_SO_3) - 2.0*pow(Fcal_SO_3, 2)) + 2.0*E_0*pow(Fcal_2, 3)/Fcal_0 - 4.0*E_2*pow(Fcal_2,
      2))/Fcal_0)/Fcal_0)/Fcal_0 + m*(-7.0*E_SO_5 + (E_0*(2.0*Fcal_5 + 2.0*Fcal_SO_5) + E_0*Fcal_2*(-4.0*Fcal_3 -
      4.0*Fcal_SO_3)/Fcal_0 + E_2*(4.0*Fcal_3 + 4.0*Fcal_SO_3) + 5.0*E_SO_3*Fcal_2)/Fcal_0)/Fcal_0) + m*(-6.0*E_4 -
      6.0*E_SQ_4 + (E_0*(2.0*Fcal_4 + 2.0*Fcal_SQ_4) - 2.0*E_0*pow(Fcal_2, 2)/Fcal_0 + 4.0*E_2*Fcal_2)/Fcal_0)/Fcal_0) +
      m*(E_0*(2.0*Fcal_3 + 2.0*Fcal_SO_3)/Fcal_0 - 5.0*E_SO_3)/Fcal_0) + m*(2.0*E_0*Fcal_2/Fcal_0 -
      4.0*E_2)/Fcal_0))/Fcal_coeff;
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_3p0PN_Q : public TaylorTn_Q


class TaylorTn_3p5PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, logv, Fcal_coeff;
  const double Fcal_0, Fcal_2, Fcal_3, Fcal_4, Fcal_5, Fcal_6, Fcal_lnv_6, Fcal_7;
  double Fcal_SQ_4, Fcal_SO_3, Fcal_SO_5, Fcal_SO_6, Fcal_SO_7;
  const double E_0, E_2, E_4, E_6;
  double E_SQ_4, E_SO_3, E_SO_5, E_SO_7;
  double Phi;

public:
  TaylorTn_3p5PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), logv(log(v)), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    Fcal_2(-2.91666666666667*nu - 3.71130952380952), Fcal_3(12.5663706143592), Fcal_4(3.61111111111111*pow(nu, 2) +
    18.3948412698413*nu - 4.92846119929453), Fcal_5(-76.3145215434521*nu - 38.2928354546934),
    Fcal_6(-2.39197530864198*pow(nu, 3) - 31.2179232804233*pow(nu, 2) - 8.87205344238227*nu + 115.731716675611),
    Fcal_lnv_6(-16.3047619047619), Fcal_7(200.905057974359*pow(nu, 2) + 390.417427312002*nu - 101.509595959742),
    Fcal_SQ_4(chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
    2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
    chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
    2)*(0.0416666666666667*nu + 2.98958333333333)), Fcal_SO_3((-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2)),
    Fcal_SO_5((S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2)),
    Fcal_SO_6((-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2)),
    Fcal_SO_7((S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
    Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762))/pow(m, 2)),
    E_0(1.00000000000000), E_2(-0.0833333333333333*nu - 0.75), E_4(-0.0416666666666667*pow(nu, 2) + 2.375*nu - 3.375),
    E_6(-0.00675154320987654*pow(nu, 3) - 1.61458333333333*pow(nu, 2) + 38.7246294907293*nu - 10.546875),
    E_SQ_4(-1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 +
    6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0)), E_SO_3((4.66666666666667*S_l +
    2.0*Sigma_l*delta)/pow(m, 2)), E_SO_5((S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu +
    3.0))/pow(m, 2)), E_SO_7((S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75) + Sigma_l*delta*(1.25*pow(nu, 2) -
    39.0*nu + 6.75))/pow(m, 2)), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    logv = log(v);
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
    Fcal_SQ_4 = chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
      2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
      chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
      2)*(0.0416666666666667*nu + 2.98958333333333);
    Fcal_SO_3 = (-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_5 = (S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2);
    Fcal_SO_6 = (-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_7 = (S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
      Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762))/pow(m, 2);
    E_SQ_4 = -1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 +
      6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0);
    E_SO_3 = (4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2);
    E_SO_5 = (S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2);
    E_SO_7 = (S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75) + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu +
      6.75))/pow(m, 2);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*(-chiVec2*pow(m2, 2)*v/pow(m, 2) + ellHat*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu -
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu +
      3.0*chi2_n*pow(m2, 2)/pow(m, 2)));
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*(-chiVec1*pow(m1, 2)*v/pow(m, 2) + ellHat*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu +
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*pow(m1,
      2)/pow(m, 2) + 3.0*chi2_n*nu));
  }
  Quaternion OmegaVec() {
    double gamma_PN_2 = -0.333333333333333*nu + 1.0;
    double gamma_PN_6 = 0.0123456790123457*pow(nu, 3) + 6.36111111111111*pow(nu, 2) - 2.98177812235564*nu + 1.0;
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_3 = (1.66666666666667*S_l + Sigma_l*delta)/pow(m, 2);
    double gamma_PN_7 = (S_l*(-6.0*pow(nu, 2) - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_l*delta*pow(nu, 2) +
      Sigma_l*delta*(-10.1666666666667*nu + 3.0))/pow(m, 2);
    double gamma_PN_5 = (S_l*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_l*delta)/pow(m, 2);
    double gamma_PN_0 = 1.00000000000000;
    double a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0);
    double gamma_PN_4 = -5.41666666666667*nu + 1.0;
    double a_ell_4 = S_n*(5.77777777777778*pow(nu, 2) + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*pow(nu, 2) +
      9.125*nu + 1.5);
    return ellHat*pow(v, 3)/m + nHat*pow(v, 6)*(a_ell_0 + pow(v, 2)*(a_ell_2 + a_ell_4*pow(v, 2)))*(gamma_PN_0 + pow(v,
      2)*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/pow(m, 3);
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_coeff*(Fcal_0 + pow(v, 2)*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 +
      v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7)))))));
    const double dEdv = -0.5*m*nu*v*(2.0*E_0 + pow(v, 2)*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 +
      v*(7.0*E_SO_5 + v*(8.0*E_6 + 9.0*E_SO_7*v))))));
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = -2.0*Fcal_coeff*(pow(v, 2)*(v*(v*(v*(v*(v*((-0.5*Fcal_7 - 0.5*Fcal_SO_7)/m +
      ((E_2*(1.0*Fcal_5 + 1.0*Fcal_SO_5) + E_4*(1.5*Fcal_3 + 1.5*Fcal_SO_3) + E_SO_3*(1.25*Fcal_4 + 1.25*Fcal_SQ_4) +
      1.75*E_SO_5*Fcal_2 + 2.25*E_SO_7*Fcal_0 + E_SQ_4*(1.5*Fcal_3 + 1.5*Fcal_SO_3))/m + ((E_2*(E_2*(-2.0*Fcal_3 -
      2.0*Fcal_SO_3) - 5.0*E_SO_3*Fcal_2 - 7.0*E_SO_5*Fcal_0) - 7.5*E_4*E_SO_3*Fcal_0 - 7.5*E_SO_3*E_SQ_4*Fcal_0)/m +
      15.0*pow(E_2, 2)*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0)/E_0 + ((-0.5*Fcal_6 - 0.5*Fcal_SO_6 - 0.5*Fcal_lnv_6*logv)/m +
      ((E_2*(1.0*Fcal_4 + 1.0*Fcal_SQ_4) + 1.5*E_4*Fcal_2 + 2.0*E_6*Fcal_0 + E_SO_3*(1.25*Fcal_3 + 1.25*Fcal_SO_3) +
      1.5*E_SQ_4*Fcal_2)/m + ((E_2*(-2.0*E_2*Fcal_2 - 6.0*E_4*Fcal_0 - 6.0*E_SQ_4*Fcal_0) - 3.125*pow(E_SO_3,
      2)*Fcal_0)/m + 4.0*pow(E_2, 3)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0) + ((-0.5*Fcal_5 - 0.5*Fcal_SO_5)/m +
      ((E_2*(1.0*Fcal_3 + 1.0*Fcal_SO_3) + 1.25*E_SO_3*Fcal_2 + 1.75*E_SO_5*Fcal_0)/m -
      5.0*E_2*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_4 - 0.5*Fcal_SQ_4)/m + ((1.0*E_2*Fcal_2 + 1.5*E_4*Fcal_0 +
      1.5*E_SQ_4*Fcal_0)/m - 2.0*pow(E_2, 2)*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_3 - 0.5*Fcal_SO_3)/m +
      1.25*E_SO_3*Fcal_0/(E_0*m))/E_0) + (-0.5*Fcal_2/m + 1.0*E_2*Fcal_0/(E_0*m))/E_0) - 0.5*Fcal_0/(E_0*m))/(nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = -0.5*nu*v*(-2.0*E_0*m/Fcal_0 + pow(v, 2)*(v*(v*(v*(v*(m*v*(-9.0*E_SO_7 + (E_0*(2.0*Fcal_7 +
      2.0*Fcal_SO_7) + E_2*(4.0*Fcal_5 + 4.0*Fcal_SO_5) + E_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) + E_SO_3*(5.0*Fcal_4 +
      5.0*Fcal_SQ_4) + 7.0*E_SO_5*Fcal_2 + E_SQ_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) + (E_0*(Fcal_2*(-4.0*Fcal_5 -
      4.0*Fcal_SO_5) + Fcal_3*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) - 4.0*Fcal_4*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SQ_4) +
      E_0*pow(Fcal_2, 2)*(6.0*Fcal_3 + 6.0*Fcal_SO_3)/Fcal_0 + E_2*Fcal_2*(-8.0*Fcal_3 - 8.0*Fcal_SO_3) -
      5.0*E_SO_3*pow(Fcal_2, 2))/Fcal_0)/Fcal_0)/Fcal_0 + m*(-8.0*E_6 + (E_0*(2.0*Fcal_6 + 2.0*Fcal_SO_6 +
      2.0*Fcal_lnv_6*logv) + E_2*(4.0*Fcal_4 + 4.0*Fcal_SQ_4) + 6.0*E_4*Fcal_2 + E_SO_3*(5.0*Fcal_3 + 5.0*Fcal_SO_3) +
      6.0*E_SQ_4*Fcal_2 + (E_0*(Fcal_2*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) + Fcal_3*(-2.0*Fcal_3 - 4.0*Fcal_SO_3) -
      2.0*pow(Fcal_SO_3, 2)) + 2.0*E_0*pow(Fcal_2, 3)/Fcal_0 - 4.0*E_2*pow(Fcal_2, 2))/Fcal_0)/Fcal_0)/Fcal_0) +
      m*(-7.0*E_SO_5 + (E_0*(2.0*Fcal_5 + 2.0*Fcal_SO_5) + E_0*Fcal_2*(-4.0*Fcal_3 - 4.0*Fcal_SO_3)/Fcal_0 +
      E_2*(4.0*Fcal_3 + 4.0*Fcal_SO_3) + 5.0*E_SO_3*Fcal_2)/Fcal_0)/Fcal_0) + m*(-6.0*E_4 - 6.0*E_SQ_4 +
      (E_0*(2.0*Fcal_4 + 2.0*Fcal_SQ_4) - 2.0*E_0*pow(Fcal_2, 2)/Fcal_0 + 4.0*E_2*Fcal_2)/Fcal_0)/Fcal_0) +
      m*(E_0*(2.0*Fcal_3 + 2.0*Fcal_SO_3)/Fcal_0 - 5.0*E_SO_3)/Fcal_0) + m*(2.0*E_0*Fcal_2/Fcal_0 -
      4.0*E_2)/Fcal_0))/Fcal_coeff;
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_3p5PN_Q : public TaylorTn_Q


class TaylorTn_4p0PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, logv, Fcal_coeff;
  const double Fcal_0, Fcal_2, Fcal_3, Fcal_4, Fcal_5, Fcal_6, Fcal_lnv_6, Fcal_7, Fcal_8, Fcal_lnv_8;
  double Fcal_SQ_4, Fcal_SO_3, Fcal_SO_5, Fcal_SO_6, Fcal_SO_7, Fcal_SO_8;
  const double E_0, E_2, E_4, E_6, E_8, E_lnv_8;
  double E_SQ_4, E_SO_3, E_SO_5, E_SO_7;
  double Phi;

public:
  TaylorTn_4p0PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), logv(log(v)), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    Fcal_2(-2.91666666666667*nu - 3.71130952380952), Fcal_3(12.5663706143592), Fcal_4(3.61111111111111*pow(nu, 2) +
    18.3948412698413*nu - 4.92846119929453), Fcal_5(-76.3145215434521*nu - 38.2928354546934),
    Fcal_6(-2.39197530864198*pow(nu, 3) - 31.2179232804233*pow(nu, 2) - 8.87205344238227*nu + 115.731716675611),
    Fcal_lnv_6(-16.3047619047619), Fcal_7(200.905057974359*pow(nu, 2) + 390.417427312002*nu - 101.509595959742),
    Fcal_8(-117.504390722677), Fcal_lnv_8(52.7430839002268), Fcal_SQ_4(chi1chi1*(-0.463541666666667*delta +
    0.927083333333333*nu - 0.463541666666667) - 2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta +
    0.927083333333333*nu - 0.463541666666667) + chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) +
    5.97916666666667*chi_s_l*delta) + pow(chi_s_l, 2)*(0.0416666666666667*nu + 2.98958333333333)), Fcal_SO_3((-4.0*S_l -
    1.25*Sigma_l*delta)/pow(m, 2)), Fcal_SO_5((S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu -
    0.8125))/pow(m, 2)), Fcal_SO_6((-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2)),
    Fcal_SO_7((S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
    Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762))/pow(m, 2)),
    Fcal_SO_8((3.14159265358979*S_l*(192.763888888889*nu - 36.3020833333333) +
    3.14159265358979*Sigma_l*delta*(64.7733134920635*nu - 10.6592261904762))/pow(m, 2)), E_0(1.00000000000000),
    E_2(-0.0833333333333333*nu - 0.75), E_4(-0.0416666666666667*pow(nu, 2) + 2.375*nu - 3.375),
    E_6(-0.00675154320987654*pow(nu, 3) - 1.61458333333333*pow(nu, 2) + 38.7246294907293*nu - 10.546875),
    E_8(0.0024755658436214*pow(nu, 4) + 0.174189814814815*pow(nu, 3) - 90.1327990262052*pow(nu, 2) + 153.88379682994*nu
    - 31.0078125), E_lnv_8(59.7333333333333*nu), E_SQ_4(-1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2
    + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 + 6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0)),
    E_SO_3((4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2)), E_SO_5((S_l*(-6.77777777777778*nu + 11.0) +
    Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2)), E_SO_7((S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75)
    + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu + 6.75))/pow(m, 2)), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    logv = log(v);
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
    Fcal_SQ_4 = chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
      2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
      chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
      2)*(0.0416666666666667*nu + 2.98958333333333);
    Fcal_SO_3 = (-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_5 = (S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2);
    Fcal_SO_6 = (-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_7 = (S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
      Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762))/pow(m, 2);
    Fcal_SO_8 = (3.14159265358979*S_l*(192.763888888889*nu - 36.3020833333333) +
      3.14159265358979*Sigma_l*delta*(64.7733134920635*nu - 10.6592261904762))/pow(m, 2);
    E_SQ_4 = -1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 +
      6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0);
    E_SO_3 = (4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2);
    E_SO_5 = (S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2);
    E_SO_7 = (S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75) + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu +
      6.75))/pow(m, 2);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*(-chiVec2*pow(m2, 2)*v/pow(m, 2) + ellHat*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu -
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu +
      3.0*chi2_n*pow(m2, 2)/pow(m, 2)));
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*(-chiVec1*pow(m1, 2)*v/pow(m, 2) + ellHat*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu +
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*pow(m1,
      2)/pow(m, 2) + 3.0*chi2_n*nu));
  }
  Quaternion OmegaVec() {
    double gamma_PN_2 = -0.333333333333333*nu + 1.0;
    double gamma_PN_6 = 0.0123456790123457*pow(nu, 3) + 6.36111111111111*pow(nu, 2) - 2.98177812235564*nu + 1.0;
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_3 = (1.66666666666667*S_l + Sigma_l*delta)/pow(m, 2);
    double gamma_PN_7 = (S_l*(-6.0*pow(nu, 2) - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_l*delta*pow(nu, 2) +
      Sigma_l*delta*(-10.1666666666667*nu + 3.0))/pow(m, 2);
    double gamma_PN_5 = (S_l*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_l*delta)/pow(m, 2);
    double gamma_PN_0 = 1.00000000000000;
    double a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0);
    double gamma_PN_4 = -5.41666666666667*nu + 1.0;
    double a_ell_4 = S_n*(5.77777777777778*pow(nu, 2) + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*pow(nu, 2) +
      9.125*nu + 1.5);
    return ellHat*pow(v, 3)/m + nHat*pow(v, 6)*(a_ell_0 + pow(v, 2)*(a_ell_2 + a_ell_4*pow(v, 2)))*(gamma_PN_0 + pow(v,
      2)*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/pow(m, 3);
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_coeff*(Fcal_0 + pow(v, 2)*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 +
      v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 +
      Fcal_lnv_8*logv))))))));
    const double dEdv = -0.5*m*nu*v*(2.0*E_0 + pow(v, 2)*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 +
      v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0)))))))));
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = -2.0*Fcal_coeff*(pow(v, 2)*(v*(v*(v*(v*(v*(v*((-0.5*Fcal_8 - 0.5*Fcal_SO_8 -
      0.5*Fcal_lnv_8*logv)/m + ((E_2*(1.0*Fcal_6 + 1.0*Fcal_SO_6 + 1.0*Fcal_lnv_6*logv) + E_4*(1.5*Fcal_4 +
      1.5*Fcal_SQ_4) + 2.0*E_6*Fcal_2 + 2.5*E_8*Fcal_0 + E_SO_3*(1.25*Fcal_5 + 1.25*Fcal_SO_5) + E_SO_5*(1.75*Fcal_3 +
      1.75*Fcal_SO_3) + E_SQ_4*(1.5*Fcal_4 + 1.5*Fcal_SQ_4) + E_lnv_8*Fcal_0*(2.5*logv + 0.25))/m +
      ((E_2*(E_2*(-2.0*Fcal_4 - 2.0*Fcal_SQ_4) - 6.0*E_4*Fcal_2 - 8.0*E_6*Fcal_0 + E_SO_3*(-5.0*Fcal_3 - 5.0*Fcal_SO_3)
      - 6.0*E_SQ_4*Fcal_2) + E_4*(-4.5*E_4*Fcal_0 - 9.0*E_SQ_4*Fcal_0) + E_SO_3*(-3.125*E_SO_3*Fcal_2 -
      8.75*E_SO_5*Fcal_0) - 4.5*pow(E_SQ_4, 2)*Fcal_0)/m + (E_2*(E_2*(4.0*E_2*Fcal_2 + 18.0*E_4*Fcal_0 +
      18.0*E_SQ_4*Fcal_0) + 18.75*pow(E_SO_3, 2)*Fcal_0)/m - 8.0*pow(E_2, 4)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0 +
      ((-0.5*Fcal_7 - 0.5*Fcal_SO_7)/m + ((E_2*(1.0*Fcal_5 + 1.0*Fcal_SO_5) + E_4*(1.5*Fcal_3 + 1.5*Fcal_SO_3) +
      E_SO_3*(1.25*Fcal_4 + 1.25*Fcal_SQ_4) + 1.75*E_SO_5*Fcal_2 + 2.25*E_SO_7*Fcal_0 + E_SQ_4*(1.5*Fcal_3 +
      1.5*Fcal_SO_3))/m + ((E_2*(E_2*(-2.0*Fcal_3 - 2.0*Fcal_SO_3) - 5.0*E_SO_3*Fcal_2 - 7.0*E_SO_5*Fcal_0) -
      7.5*E_4*E_SO_3*Fcal_0 - 7.5*E_SO_3*E_SQ_4*Fcal_0)/m + 15.0*pow(E_2, 2)*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0)/E_0) +
      ((-0.5*Fcal_6 - 0.5*Fcal_SO_6 - 0.5*Fcal_lnv_6*logv)/m + ((E_2*(1.0*Fcal_4 + 1.0*Fcal_SQ_4) + 1.5*E_4*Fcal_2 +
      2.0*E_6*Fcal_0 + E_SO_3*(1.25*Fcal_3 + 1.25*Fcal_SO_3) + 1.5*E_SQ_4*Fcal_2)/m + ((E_2*(-2.0*E_2*Fcal_2 -
      6.0*E_4*Fcal_0 - 6.0*E_SQ_4*Fcal_0) - 3.125*pow(E_SO_3, 2)*Fcal_0)/m + 4.0*pow(E_2,
      3)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0) + ((-0.5*Fcal_5 - 0.5*Fcal_SO_5)/m + ((E_2*(1.0*Fcal_3 + 1.0*Fcal_SO_3) +
      1.25*E_SO_3*Fcal_2 + 1.75*E_SO_5*Fcal_0)/m - 5.0*E_2*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_4 -
      0.5*Fcal_SQ_4)/m + ((1.0*E_2*Fcal_2 + 1.5*E_4*Fcal_0 + 1.5*E_SQ_4*Fcal_0)/m - 2.0*pow(E_2,
      2)*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_3 - 0.5*Fcal_SO_3)/m + 1.25*E_SO_3*Fcal_0/(E_0*m))/E_0) +
      (-0.5*Fcal_2/m + 1.0*E_2*Fcal_0/(E_0*m))/E_0) - 0.5*Fcal_0/(E_0*m))/(nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = -0.5*nu*v*(-2.0*E_0*m/Fcal_0 + pow(v, 2)*(v*(v*(v*(v*(v*(m*v*(-10.0*E_8 + E_lnv_8*(-10.0*logv -
      1.0) + (E_0*(2.0*Fcal_8 + 2.0*Fcal_SO_8 + 2.0*Fcal_lnv_8*logv) + E_2*(4.0*Fcal_6 + 4.0*Fcal_SO_6 +
      4.0*Fcal_lnv_6*logv) + E_4*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + 8.0*E_6*Fcal_2 + E_SO_3*(5.0*Fcal_5 + 5.0*Fcal_SO_5) +
      E_SO_5*(7.0*Fcal_3 + 7.0*Fcal_SO_3) + E_SQ_4*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + (E_0*(Fcal_2*(-4.0*Fcal_6 -
      4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) + Fcal_3*(-4.0*Fcal_5 - 4.0*Fcal_SO_5) + Fcal_4*(-2.0*Fcal_4 - 4.0*Fcal_SQ_4)
      - 4.0*Fcal_5*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SO_5 - 2.0*pow(Fcal_SQ_4, 2)) + E_2*(Fcal_2*(-8.0*Fcal_4 -
      8.0*Fcal_SQ_4) + Fcal_3*(-4.0*Fcal_3 - 8.0*Fcal_SO_3) - 4.0*pow(Fcal_SO_3, 2)) - 6.0*E_4*pow(Fcal_2, 2) +
      E_SO_3*Fcal_2*(-10.0*Fcal_3 - 10.0*Fcal_SO_3) - 6.0*E_SQ_4*pow(Fcal_2, 2) + (E_0*Fcal_2*(Fcal_2*(6.0*Fcal_4 +
      6.0*Fcal_SQ_4) + Fcal_3*(6.0*Fcal_3 + 12.0*Fcal_SO_3) + 6.0*pow(Fcal_SO_3, 2)) - 2.0*E_0*pow(Fcal_2, 4)/Fcal_0 +
      4.0*E_2*pow(Fcal_2, 3))/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0 + m*(-9.0*E_SO_7 + (E_0*(2.0*Fcal_7 + 2.0*Fcal_SO_7) +
      E_2*(4.0*Fcal_5 + 4.0*Fcal_SO_5) + E_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) + E_SO_3*(5.0*Fcal_4 + 5.0*Fcal_SQ_4) +
      7.0*E_SO_5*Fcal_2 + E_SQ_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) + (E_0*(Fcal_2*(-4.0*Fcal_5 - 4.0*Fcal_SO_5) +
      Fcal_3*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) - 4.0*Fcal_4*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SQ_4) + E_0*pow(Fcal_2,
      2)*(6.0*Fcal_3 + 6.0*Fcal_SO_3)/Fcal_0 + E_2*Fcal_2*(-8.0*Fcal_3 - 8.0*Fcal_SO_3) - 5.0*E_SO_3*pow(Fcal_2,
      2))/Fcal_0)/Fcal_0)/Fcal_0) + m*(-8.0*E_6 + (E_0*(2.0*Fcal_6 + 2.0*Fcal_SO_6 + 2.0*Fcal_lnv_6*logv) +
      E_2*(4.0*Fcal_4 + 4.0*Fcal_SQ_4) + 6.0*E_4*Fcal_2 + E_SO_3*(5.0*Fcal_3 + 5.0*Fcal_SO_3) + 6.0*E_SQ_4*Fcal_2 +
      (E_0*(Fcal_2*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) + Fcal_3*(-2.0*Fcal_3 - 4.0*Fcal_SO_3) - 2.0*pow(Fcal_SO_3, 2)) +
      2.0*E_0*pow(Fcal_2, 3)/Fcal_0 - 4.0*E_2*pow(Fcal_2, 2))/Fcal_0)/Fcal_0)/Fcal_0) + m*(-7.0*E_SO_5 +
      (E_0*(2.0*Fcal_5 + 2.0*Fcal_SO_5) + E_0*Fcal_2*(-4.0*Fcal_3 - 4.0*Fcal_SO_3)/Fcal_0 + E_2*(4.0*Fcal_3 +
      4.0*Fcal_SO_3) + 5.0*E_SO_3*Fcal_2)/Fcal_0)/Fcal_0) + m*(-6.0*E_4 - 6.0*E_SQ_4 + (E_0*(2.0*Fcal_4 + 2.0*Fcal_SQ_4)
      - 2.0*E_0*pow(Fcal_2, 2)/Fcal_0 + 4.0*E_2*Fcal_2)/Fcal_0)/Fcal_0) + m*(E_0*(2.0*Fcal_3 + 2.0*Fcal_SO_3)/Fcal_0 -
      5.0*E_SO_3)/Fcal_0) + m*(2.0*E_0*Fcal_2/Fcal_0 - 4.0*E_2)/Fcal_0))/Fcal_coeff;
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_4p0PN_Q : public TaylorTn_Q


class TaylorTn_4p5PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, logv, Fcal_coeff;
  const double Fcal_0, Fcal_2, Fcal_3, Fcal_4, Fcal_5, Fcal_6, Fcal_lnv_6, Fcal_7, Fcal_8, Fcal_lnv_8, Fcal_9,
               Fcal_lnv_9;
  double Fcal_SQ_4, Fcal_SO_3, Fcal_SO_5, Fcal_SO_6, Fcal_SO_7, Fcal_SO_8;
  const double E_0, E_2, E_4, E_6, E_8, E_lnv_8;
  double E_SQ_4, E_SO_3, E_SO_5, E_SO_7;
  double Phi;

public:
  TaylorTn_4p5PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), logv(log(v)), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    Fcal_2(-2.91666666666667*nu - 3.71130952380952), Fcal_3(12.5663706143592), Fcal_4(3.61111111111111*pow(nu, 2) +
    18.3948412698413*nu - 4.92846119929453), Fcal_5(-76.3145215434521*nu - 38.2928354546934),
    Fcal_6(-2.39197530864198*pow(nu, 3) - 31.2179232804233*pow(nu, 2) - 8.87205344238227*nu + 115.731716675611),
    Fcal_lnv_6(-16.3047619047619), Fcal_7(200.905057974359*pow(nu, 2) + 390.417427312002*nu - 101.509595959742),
    Fcal_8(-117.504390722677), Fcal_lnv_8(52.7430839002268), Fcal_9(719.128342233430), Fcal_lnv_9(-204.891680874123),
    Fcal_SQ_4(chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
    2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
    chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
    2)*(0.0416666666666667*nu + 2.98958333333333)), Fcal_SO_3((-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2)),
    Fcal_SO_5((S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2)),
    Fcal_SO_6((-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2)),
    Fcal_SO_7((S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
    Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762))/pow(m, 2)),
    Fcal_SO_8((3.14159265358979*S_l*(192.763888888889*nu - 36.3020833333333) +
    3.14159265358979*Sigma_l*delta*(64.7733134920635*nu - 10.6592261904762))/pow(m, 2)), E_0(1.00000000000000),
    E_2(-0.0833333333333333*nu - 0.75), E_4(-0.0416666666666667*pow(nu, 2) + 2.375*nu - 3.375),
    E_6(-0.00675154320987654*pow(nu, 3) - 1.61458333333333*pow(nu, 2) + 38.7246294907293*nu - 10.546875),
    E_8(0.0024755658436214*pow(nu, 4) + 0.174189814814815*pow(nu, 3) - 90.1327990262052*pow(nu, 2) + 153.88379682994*nu
    - 31.0078125), E_lnv_8(59.7333333333333*nu), E_SQ_4(-1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2
    + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 + 6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0)),
    E_SO_3((4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2)), E_SO_5((S_l*(-6.77777777777778*nu + 11.0) +
    Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2)), E_SO_7((S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75)
    + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu + 6.75))/pow(m, 2)), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    logv = log(v);
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
    Fcal_SQ_4 = chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
      2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
      chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
      2)*(0.0416666666666667*nu + 2.98958333333333);
    Fcal_SO_3 = (-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_5 = (S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2);
    Fcal_SO_6 = (-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_7 = (S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
      Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762))/pow(m, 2);
    Fcal_SO_8 = (3.14159265358979*S_l*(192.763888888889*nu - 36.3020833333333) +
      3.14159265358979*Sigma_l*delta*(64.7733134920635*nu - 10.6592261904762))/pow(m, 2);
    E_SQ_4 = -1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 +
      6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0);
    E_SO_3 = (4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2);
    E_SO_5 = (S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2);
    E_SO_7 = (S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75) + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu +
      6.75))/pow(m, 2);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*(-chiVec2*pow(m2, 2)*v/pow(m, 2) + ellHat*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu -
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu +
      3.0*chi2_n*pow(m2, 2)/pow(m, 2)));
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*(-chiVec1*pow(m1, 2)*v/pow(m, 2) + ellHat*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu +
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*pow(m1,
      2)/pow(m, 2) + 3.0*chi2_n*nu));
  }
  Quaternion OmegaVec() {
    double gamma_PN_2 = -0.333333333333333*nu + 1.0;
    double gamma_PN_6 = 0.0123456790123457*pow(nu, 3) + 6.36111111111111*pow(nu, 2) - 2.98177812235564*nu + 1.0;
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_3 = (1.66666666666667*S_l + Sigma_l*delta)/pow(m, 2);
    double gamma_PN_7 = (S_l*(-6.0*pow(nu, 2) - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_l*delta*pow(nu, 2) +
      Sigma_l*delta*(-10.1666666666667*nu + 3.0))/pow(m, 2);
    double gamma_PN_5 = (S_l*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_l*delta)/pow(m, 2);
    double gamma_PN_0 = 1.00000000000000;
    double a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0);
    double gamma_PN_4 = -5.41666666666667*nu + 1.0;
    double a_ell_4 = S_n*(5.77777777777778*pow(nu, 2) + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*pow(nu, 2) +
      9.125*nu + 1.5);
    return ellHat*pow(v, 3)/m + nHat*pow(v, 6)*(a_ell_0 + pow(v, 2)*(a_ell_2 + a_ell_4*pow(v, 2)))*(gamma_PN_0 + pow(v,
      2)*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/pow(m, 3);
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_coeff*(Fcal_0 + pow(v, 2)*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 +
      v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 +
      Fcal_lnv_8*logv + v*(Fcal_9 + Fcal_lnv_9*logv)))))))));
    const double dEdv = -0.5*m*nu*v*(2.0*E_0 + pow(v, 2)*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 +
      v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0)))))))));
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = -2.0*Fcal_coeff*(pow(v, 2)*(v*(v*(v*(v*(v*(v*(v*((-0.5*Fcal_9 - 0.5*Fcal_lnv_9*logv)/m +
      ((E_2*(1.0*Fcal_7 + 1.0*Fcal_SO_7) + E_4*(1.5*Fcal_5 + 1.5*Fcal_SO_5) + E_6*(2.0*Fcal_3 + 2.0*Fcal_SO_3) +
      E_SO_3*(1.25*Fcal_6 + 1.25*Fcal_SO_6 + 1.25*Fcal_lnv_6*logv) + E_SO_5*(1.75*Fcal_4 + 1.75*Fcal_SQ_4) +
      2.25*E_SO_7*Fcal_2 + E_SQ_4*(1.5*Fcal_5 + 1.5*Fcal_SO_5))/m + ((E_2*(E_2*(-2.0*Fcal_5 - 2.0*Fcal_SO_5) +
      E_4*(-6.0*Fcal_3 - 6.0*Fcal_SO_3) + E_SO_3*(-5.0*Fcal_4 - 5.0*Fcal_SQ_4) - 7.0*E_SO_5*Fcal_2 - 9.0*E_SO_7*Fcal_0 +
      E_SQ_4*(-6.0*Fcal_3 - 6.0*Fcal_SO_3)) + E_4*(-7.5*E_SO_3*Fcal_2 - 10.5*E_SO_5*Fcal_0) - 10.0*E_6*E_SO_3*Fcal_0 +
      E_SO_3*(E_SO_3*(-3.125*Fcal_3 - 3.125*Fcal_SO_3) - 7.5*E_SQ_4*Fcal_2) - 10.5*E_SO_5*E_SQ_4*Fcal_0)/m +
      ((E_2*(E_2*(E_2*(4.0*Fcal_3 + 4.0*Fcal_SO_3) + 15.0*E_SO_3*Fcal_2 + 21.0*E_SO_5*Fcal_0) + 45.0*E_4*E_SO_3*Fcal_0 +
      45.0*E_SO_3*E_SQ_4*Fcal_0) + 7.8125*pow(E_SO_3, 3)*Fcal_0)/m - 40.0*pow(E_2,
      3)*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0 + ((-0.5*Fcal_8 - 0.5*Fcal_SO_8 - 0.5*Fcal_lnv_8*logv)/m +
      ((E_2*(1.0*Fcal_6 + 1.0*Fcal_SO_6 + 1.0*Fcal_lnv_6*logv) + E_4*(1.5*Fcal_4 + 1.5*Fcal_SQ_4) + 2.0*E_6*Fcal_2 +
      2.5*E_8*Fcal_0 + E_SO_3*(1.25*Fcal_5 + 1.25*Fcal_SO_5) + E_SO_5*(1.75*Fcal_3 + 1.75*Fcal_SO_3) +
      E_SQ_4*(1.5*Fcal_4 + 1.5*Fcal_SQ_4) + E_lnv_8*Fcal_0*(2.5*logv + 0.25))/m + ((E_2*(E_2*(-2.0*Fcal_4 -
      2.0*Fcal_SQ_4) - 6.0*E_4*Fcal_2 - 8.0*E_6*Fcal_0 + E_SO_3*(-5.0*Fcal_3 - 5.0*Fcal_SO_3) - 6.0*E_SQ_4*Fcal_2) +
      E_4*(-4.5*E_4*Fcal_0 - 9.0*E_SQ_4*Fcal_0) + E_SO_3*(-3.125*E_SO_3*Fcal_2 - 8.75*E_SO_5*Fcal_0) - 4.5*pow(E_SQ_4,
      2)*Fcal_0)/m + (E_2*(E_2*(4.0*E_2*Fcal_2 + 18.0*E_4*Fcal_0 + 18.0*E_SQ_4*Fcal_0) + 18.75*pow(E_SO_3, 2)*Fcal_0)/m
      - 8.0*pow(E_2, 4)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0) + ((-0.5*Fcal_7 - 0.5*Fcal_SO_7)/m + ((E_2*(1.0*Fcal_5 +
      1.0*Fcal_SO_5) + E_4*(1.5*Fcal_3 + 1.5*Fcal_SO_3) + E_SO_3*(1.25*Fcal_4 + 1.25*Fcal_SQ_4) + 1.75*E_SO_5*Fcal_2 +
      2.25*E_SO_7*Fcal_0 + E_SQ_4*(1.5*Fcal_3 + 1.5*Fcal_SO_3))/m + ((E_2*(E_2*(-2.0*Fcal_3 - 2.0*Fcal_SO_3) -
      5.0*E_SO_3*Fcal_2 - 7.0*E_SO_5*Fcal_0) - 7.5*E_4*E_SO_3*Fcal_0 - 7.5*E_SO_3*E_SQ_4*Fcal_0)/m + 15.0*pow(E_2,
      2)*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0)/E_0) + ((-0.5*Fcal_6 - 0.5*Fcal_SO_6 - 0.5*Fcal_lnv_6*logv)/m +
      ((E_2*(1.0*Fcal_4 + 1.0*Fcal_SQ_4) + 1.5*E_4*Fcal_2 + 2.0*E_6*Fcal_0 + E_SO_3*(1.25*Fcal_3 + 1.25*Fcal_SO_3) +
      1.5*E_SQ_4*Fcal_2)/m + ((E_2*(-2.0*E_2*Fcal_2 - 6.0*E_4*Fcal_0 - 6.0*E_SQ_4*Fcal_0) - 3.125*pow(E_SO_3,
      2)*Fcal_0)/m + 4.0*pow(E_2, 3)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0) + ((-0.5*Fcal_5 - 0.5*Fcal_SO_5)/m +
      ((E_2*(1.0*Fcal_3 + 1.0*Fcal_SO_3) + 1.25*E_SO_3*Fcal_2 + 1.75*E_SO_5*Fcal_0)/m -
      5.0*E_2*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_4 - 0.5*Fcal_SQ_4)/m + ((1.0*E_2*Fcal_2 + 1.5*E_4*Fcal_0 +
      1.5*E_SQ_4*Fcal_0)/m - 2.0*pow(E_2, 2)*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_3 - 0.5*Fcal_SO_3)/m +
      1.25*E_SO_3*Fcal_0/(E_0*m))/E_0) + (-0.5*Fcal_2/m + 1.0*E_2*Fcal_0/(E_0*m))/E_0) - 0.5*Fcal_0/(E_0*m))/(nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = -0.5*nu*v*(-2.0*E_0*m/Fcal_0 + pow(v, 2)*(v*(v*(v*(v*(v*(v*(m*(-10.0*E_8 + E_lnv_8*(-10.0*logv -
      1.0) + (E_0*(2.0*Fcal_8 + 2.0*Fcal_SO_8 + 2.0*Fcal_lnv_8*logv) + E_2*(4.0*Fcal_6 + 4.0*Fcal_SO_6 +
      4.0*Fcal_lnv_6*logv) + E_4*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + 8.0*E_6*Fcal_2 + E_SO_3*(5.0*Fcal_5 + 5.0*Fcal_SO_5) +
      E_SO_5*(7.0*Fcal_3 + 7.0*Fcal_SO_3) + E_SQ_4*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + (E_0*(Fcal_2*(-4.0*Fcal_6 -
      4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) + Fcal_3*(-4.0*Fcal_5 - 4.0*Fcal_SO_5) + Fcal_4*(-2.0*Fcal_4 - 4.0*Fcal_SQ_4)
      - 4.0*Fcal_5*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SO_5 - 2.0*pow(Fcal_SQ_4, 2)) + E_2*(Fcal_2*(-8.0*Fcal_4 -
      8.0*Fcal_SQ_4) + Fcal_3*(-4.0*Fcal_3 - 8.0*Fcal_SO_3) - 4.0*pow(Fcal_SO_3, 2)) - 6.0*E_4*pow(Fcal_2, 2) +
      E_SO_3*Fcal_2*(-10.0*Fcal_3 - 10.0*Fcal_SO_3) - 6.0*E_SQ_4*pow(Fcal_2, 2) + (E_0*Fcal_2*(Fcal_2*(6.0*Fcal_4 +
      6.0*Fcal_SQ_4) + Fcal_3*(6.0*Fcal_3 + 12.0*Fcal_SO_3) + 6.0*pow(Fcal_SO_3, 2)) - 2.0*E_0*pow(Fcal_2, 4)/Fcal_0 +
      4.0*E_2*pow(Fcal_2, 3))/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0 + m*v*(E_0*(2.0*Fcal_9 + 2.0*Fcal_lnv_9*logv) +
      E_2*(4.0*Fcal_7 + 4.0*Fcal_SO_7) + E_4*(6.0*Fcal_5 + 6.0*Fcal_SO_5) + E_6*(8.0*Fcal_3 + 8.0*Fcal_SO_3) +
      E_SO_3*(5.0*Fcal_6 + 5.0*Fcal_SO_6 + 5.0*Fcal_lnv_6*logv) + E_SO_5*(7.0*Fcal_4 + 7.0*Fcal_SQ_4) +
      9.0*E_SO_7*Fcal_2 + E_SQ_4*(6.0*Fcal_5 + 6.0*Fcal_SO_5) + (E_0*(Fcal_2*(-4.0*Fcal_7 - 4.0*Fcal_SO_7) +
      Fcal_3*(-4.0*Fcal_6 - 4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) + Fcal_4*(-4.0*Fcal_5 - 4.0*Fcal_SO_5) -
      4.0*Fcal_5*Fcal_SQ_4 - 4.0*Fcal_6*Fcal_SO_3 + Fcal_SO_3*(-4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) -
      4.0*Fcal_SO_5*Fcal_SQ_4) + E_2*(Fcal_2*(-8.0*Fcal_5 - 8.0*Fcal_SO_5) + Fcal_3*(-8.0*Fcal_4 - 8.0*Fcal_SQ_4) -
      8.0*Fcal_4*Fcal_SO_3 - 8.0*Fcal_SO_3*Fcal_SQ_4) + E_4*Fcal_2*(-12.0*Fcal_3 - 12.0*Fcal_SO_3) +
      E_SO_3*(Fcal_2*(-10.0*Fcal_4 - 10.0*Fcal_SQ_4) + Fcal_3*(-5.0*Fcal_3 - 10.0*Fcal_SO_3) - 5.0*pow(Fcal_SO_3, 2)) -
      7.0*E_SO_5*pow(Fcal_2, 2) + E_SQ_4*Fcal_2*(-12.0*Fcal_3 - 12.0*Fcal_SO_3) + (E_0*(Fcal_2*(Fcal_2*(6.0*Fcal_5 +
      6.0*Fcal_SO_5) + Fcal_3*(12.0*Fcal_4 + 12.0*Fcal_SQ_4) + 12.0*Fcal_4*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SQ_4) +
      Fcal_3*(Fcal_3*(2.0*Fcal_3 + 6.0*Fcal_SO_3) + 6.0*pow(Fcal_SO_3, 2)) + 2.0*pow(Fcal_SO_3, 3)) + E_0*pow(Fcal_2,
      3)*(-8.0*Fcal_3 - 8.0*Fcal_SO_3)/Fcal_0 + E_2*pow(Fcal_2, 2)*(12.0*Fcal_3 + 12.0*Fcal_SO_3) +
      5.0*E_SO_3*pow(Fcal_2, 3))/Fcal_0)/Fcal_0)/pow(Fcal_0, 2)) + m*(-9.0*E_SO_7 + (E_0*(2.0*Fcal_7 + 2.0*Fcal_SO_7) +
      E_2*(4.0*Fcal_5 + 4.0*Fcal_SO_5) + E_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) + E_SO_3*(5.0*Fcal_4 + 5.0*Fcal_SQ_4) +
      7.0*E_SO_5*Fcal_2 + E_SQ_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) + (E_0*(Fcal_2*(-4.0*Fcal_5 - 4.0*Fcal_SO_5) +
      Fcal_3*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) - 4.0*Fcal_4*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SQ_4) + E_0*pow(Fcal_2,
      2)*(6.0*Fcal_3 + 6.0*Fcal_SO_3)/Fcal_0 + E_2*Fcal_2*(-8.0*Fcal_3 - 8.0*Fcal_SO_3) - 5.0*E_SO_3*pow(Fcal_2,
      2))/Fcal_0)/Fcal_0)/Fcal_0) + m*(-8.0*E_6 + (E_0*(2.0*Fcal_6 + 2.0*Fcal_SO_6 + 2.0*Fcal_lnv_6*logv) +
      E_2*(4.0*Fcal_4 + 4.0*Fcal_SQ_4) + 6.0*E_4*Fcal_2 + E_SO_3*(5.0*Fcal_3 + 5.0*Fcal_SO_3) + 6.0*E_SQ_4*Fcal_2 +
      (E_0*(Fcal_2*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) + Fcal_3*(-2.0*Fcal_3 - 4.0*Fcal_SO_3) - 2.0*pow(Fcal_SO_3, 2)) +
      2.0*E_0*pow(Fcal_2, 3)/Fcal_0 - 4.0*E_2*pow(Fcal_2, 2))/Fcal_0)/Fcal_0)/Fcal_0) + m*(-7.0*E_SO_5 +
      (E_0*(2.0*Fcal_5 + 2.0*Fcal_SO_5) + E_0*Fcal_2*(-4.0*Fcal_3 - 4.0*Fcal_SO_3)/Fcal_0 + E_2*(4.0*Fcal_3 +
      4.0*Fcal_SO_3) + 5.0*E_SO_3*Fcal_2)/Fcal_0)/Fcal_0) + m*(-6.0*E_4 - 6.0*E_SQ_4 + (E_0*(2.0*Fcal_4 + 2.0*Fcal_SQ_4)
      - 2.0*E_0*pow(Fcal_2, 2)/Fcal_0 + 4.0*E_2*Fcal_2)/Fcal_0)/Fcal_0) + m*(E_0*(2.0*Fcal_3 + 2.0*Fcal_SO_3)/Fcal_0 -
      5.0*E_SO_3)/Fcal_0) + m*(2.0*E_0*Fcal_2/Fcal_0 - 4.0*E_2)/Fcal_0))/Fcal_coeff;
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_4p5PN_Q : public TaylorTn_Q


class TaylorTn_5p0PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, logv, Fcal_coeff;
  const double Fcal_0, Fcal_2, Fcal_3, Fcal_4, Fcal_5, Fcal_6, Fcal_lnv_6, Fcal_7, Fcal_8, Fcal_lnv_8, Fcal_9,
               Fcal_lnv_9, Fcal_10, Fcal_lnv_10;
  double Fcal_SQ_4, Fcal_SO_3, Fcal_SO_5, Fcal_SO_6, Fcal_SO_7, Fcal_SO_8;
  const double E_0, E_2, E_4, E_6, E_8, E_lnv_8, E_10, E_lnv_10;
  double E_SQ_4, E_SO_3, E_SO_5, E_SO_7;
  double Phi;

public:
  TaylorTn_5p0PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), logv(log(v)), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    Fcal_2(-2.91666666666667*nu - 3.71130952380952), Fcal_3(12.5663706143592), Fcal_4(3.61111111111111*pow(nu, 2) +
    18.3948412698413*nu - 4.92846119929453), Fcal_5(-76.3145215434521*nu - 38.2928354546934),
    Fcal_6(-2.39197530864198*pow(nu, 3) - 31.2179232804233*pow(nu, 2) - 8.87205344238227*nu + 115.731716675611),
    Fcal_lnv_6(-16.3047619047619), Fcal_7(200.905057974359*pow(nu, 2) + 390.417427312002*nu - 101.509595959742),
    Fcal_8(-117.504390722677), Fcal_lnv_8(52.7430839002268), Fcal_9(719.128342233430), Fcal_lnv_9(-204.891680874123),
    Fcal_10(-1216.90699131704), Fcal_lnv_10(116.639876594109), Fcal_SQ_4(chi1chi1*(-0.463541666666667*delta +
    0.927083333333333*nu - 0.463541666666667) - 2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta +
    0.927083333333333*nu - 0.463541666666667) + chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) +
    5.97916666666667*chi_s_l*delta) + pow(chi_s_l, 2)*(0.0416666666666667*nu + 2.98958333333333)), Fcal_SO_3((-4.0*S_l -
    1.25*Sigma_l*delta)/pow(m, 2)), Fcal_SO_5((S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu -
    0.8125))/pow(m, 2)), Fcal_SO_6((-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2)),
    Fcal_SO_7((S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
    Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762))/pow(m, 2)),
    Fcal_SO_8((3.14159265358979*S_l*(192.763888888889*nu - 36.3020833333333) +
    3.14159265358979*Sigma_l*delta*(64.7733134920635*nu - 10.6592261904762))/pow(m, 2)), E_0(1.00000000000000),
    E_2(-0.0833333333333333*nu - 0.75), E_4(-0.0416666666666667*pow(nu, 2) + 2.375*nu - 3.375),
    E_6(-0.00675154320987654*pow(nu, 3) - 1.61458333333333*pow(nu, 2) + 38.7246294907293*nu - 10.546875),
    E_8(0.0024755658436214*pow(nu, 4) + 0.174189814814815*pow(nu, 3) - 90.1327990262052*pow(nu, 2) + 153.88379682994*nu
    - 31.0078125), E_lnv_8(59.7333333333333*nu), E_10(0.001953125*pow(nu, 5) + 0.107421875*pow(nu, 4) +
    83.4293954895551*pow(nu, 3) - 71.3641942901125*pow(nu, 2) - 62.1764445589718*nu - 89.701171875),
    E_lnv_10(-262.4*pow(nu, 2) - 285.028571428571*nu), E_SQ_4(-1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) -
    delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 + 6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta
    - 2.0*nu + 1.0)), E_SO_3((4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2)), E_SO_5((S_l*(-6.77777777777778*nu +
    11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2)), E_SO_7((S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu
    + 33.75) + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu + 6.75))/pow(m, 2)), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    logv = log(v);
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
    Fcal_SQ_4 = chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
      2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
      chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
      2)*(0.0416666666666667*nu + 2.98958333333333);
    Fcal_SO_3 = (-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_5 = (S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2);
    Fcal_SO_6 = (-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_7 = (S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
      Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762))/pow(m, 2);
    Fcal_SO_8 = (3.14159265358979*S_l*(192.763888888889*nu - 36.3020833333333) +
      3.14159265358979*Sigma_l*delta*(64.7733134920635*nu - 10.6592261904762))/pow(m, 2);
    E_SQ_4 = -1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 +
      6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0);
    E_SO_3 = (4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2);
    E_SO_5 = (S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2);
    E_SO_7 = (S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75) + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu +
      6.75))/pow(m, 2);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*(-chiVec2*pow(m2, 2)*v/pow(m, 2) + ellHat*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu -
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu +
      3.0*chi2_n*pow(m2, 2)/pow(m, 2)));
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*(-chiVec1*pow(m1, 2)*v/pow(m, 2) + ellHat*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu +
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*pow(m1,
      2)/pow(m, 2) + 3.0*chi2_n*nu));
  }
  Quaternion OmegaVec() {
    double gamma_PN_2 = -0.333333333333333*nu + 1.0;
    double gamma_PN_6 = 0.0123456790123457*pow(nu, 3) + 6.36111111111111*pow(nu, 2) - 2.98177812235564*nu + 1.0;
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_3 = (1.66666666666667*S_l + Sigma_l*delta)/pow(m, 2);
    double gamma_PN_7 = (S_l*(-6.0*pow(nu, 2) - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_l*delta*pow(nu, 2) +
      Sigma_l*delta*(-10.1666666666667*nu + 3.0))/pow(m, 2);
    double gamma_PN_5 = (S_l*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_l*delta)/pow(m, 2);
    double gamma_PN_0 = 1.00000000000000;
    double a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0);
    double gamma_PN_4 = -5.41666666666667*nu + 1.0;
    double a_ell_4 = S_n*(5.77777777777778*pow(nu, 2) + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*pow(nu, 2) +
      9.125*nu + 1.5);
    return ellHat*pow(v, 3)/m + nHat*pow(v, 6)*(a_ell_0 + pow(v, 2)*(a_ell_2 + a_ell_4*pow(v, 2)))*(gamma_PN_0 + pow(v,
      2)*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/pow(m, 3);
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_coeff*(Fcal_0 + pow(v, 2)*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 +
      v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 +
      Fcal_lnv_8*logv + v*(Fcal_9 + Fcal_lnv_9*logv + v*(Fcal_10 + Fcal_lnv_10*logv))))))))));
    const double dEdv = -0.5*m*nu*v*(2.0*E_0 + pow(v, 2)*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 +
      v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0) + pow(v, 2)*(12.0*E_10 +
      E_lnv_10*(12.0*logv + 1.0))))))))));
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = -2.0*Fcal_coeff*(pow(v, 2)*(v*(v*(v*(v*(v*(v*(v*(v*((-0.5*Fcal_10 - 0.5*Fcal_lnv_10*logv)/m +
      ((3.0*E_10*Fcal_0 + E_2*(1.0*Fcal_8 + 1.0*Fcal_SO_8 + 1.0*Fcal_lnv_8*logv) + E_4*(1.5*Fcal_6 + 1.5*Fcal_SO_6 +
      1.5*Fcal_lnv_6*logv) + E_6*(2.0*Fcal_4 + 2.0*Fcal_SQ_4) + 2.5*E_8*Fcal_2 + E_SO_3*(1.25*Fcal_7 + 1.25*Fcal_SO_7) +
      E_SO_5*(1.75*Fcal_5 + 1.75*Fcal_SO_5) + E_SO_7*(2.25*Fcal_3 + 2.25*Fcal_SO_3) + E_SQ_4*(1.5*Fcal_6 + 1.5*Fcal_SO_6
      + 1.5*Fcal_lnv_6*logv) + E_lnv_10*Fcal_0*(3.0*logv + 0.25) + E_lnv_8*Fcal_2*(2.5*logv + 0.25))/m +
      ((E_2*(E_2*(-2.0*Fcal_6 - 2.0*Fcal_SO_6 - 2.0*Fcal_lnv_6*logv) + E_4*(-6.0*Fcal_4 - 6.0*Fcal_SQ_4) -
      8.0*E_6*Fcal_2 - 10.0*E_8*Fcal_0 + E_SO_3*(-5.0*Fcal_5 - 5.0*Fcal_SO_5) + E_SO_5*(-7.0*Fcal_3 - 7.0*Fcal_SO_3) +
      E_SQ_4*(-6.0*Fcal_4 - 6.0*Fcal_SQ_4) + E_lnv_8*Fcal_0*(-10.0*logv - 1.0)) + E_4*(-4.5*E_4*Fcal_2 - 12.0*E_6*Fcal_0
      + E_SO_3*(-7.5*Fcal_3 - 7.5*Fcal_SO_3) - 9.0*E_SQ_4*Fcal_2) - 12.0*E_6*E_SQ_4*Fcal_0 +
      E_SO_3*(E_SO_3*(-3.125*Fcal_4 - 3.125*Fcal_SQ_4) - 8.75*E_SO_5*Fcal_2 - 11.25*E_SO_7*Fcal_0 + E_SQ_4*(-7.5*Fcal_3
      - 7.5*Fcal_SO_3)) - 6.125*pow(E_SO_5, 2)*Fcal_0 - 4.5*pow(E_SQ_4, 2)*Fcal_2)/m + ((E_2*(E_2*(E_2*(4.0*Fcal_4 +
      4.0*Fcal_SQ_4) + 18.0*E_4*Fcal_2 + 24.0*E_6*Fcal_0 + E_SO_3*(15.0*Fcal_3 + 15.0*Fcal_SO_3) + 18.0*E_SQ_4*Fcal_2) +
      E_4*(27.0*E_4*Fcal_0 + 54.0*E_SQ_4*Fcal_0) + E_SO_3*(18.75*E_SO_3*Fcal_2 + 52.5*E_SO_5*Fcal_0) + 27.0*pow(E_SQ_4,
      2)*Fcal_0) + 28.125*E_4*pow(E_SO_3, 2)*Fcal_0 + 28.125*pow(E_SO_3, 2)*E_SQ_4*Fcal_0)/m + (pow(E_2,
      2)*(E_2*(-8.0*E_2*Fcal_2 - 48.0*E_4*Fcal_0 - 48.0*E_SQ_4*Fcal_0) - 75.0*pow(E_SO_3, 2)*Fcal_0)/m + 16.0*pow(E_2,
      5)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0)/E_0 + ((-0.5*Fcal_9 - 0.5*Fcal_lnv_9*logv)/m + ((E_2*(1.0*Fcal_7 +
      1.0*Fcal_SO_7) + E_4*(1.5*Fcal_5 + 1.5*Fcal_SO_5) + E_6*(2.0*Fcal_3 + 2.0*Fcal_SO_3) + E_SO_3*(1.25*Fcal_6 +
      1.25*Fcal_SO_6 + 1.25*Fcal_lnv_6*logv) + E_SO_5*(1.75*Fcal_4 + 1.75*Fcal_SQ_4) + 2.25*E_SO_7*Fcal_2 +
      E_SQ_4*(1.5*Fcal_5 + 1.5*Fcal_SO_5))/m + ((E_2*(E_2*(-2.0*Fcal_5 - 2.0*Fcal_SO_5) + E_4*(-6.0*Fcal_3 -
      6.0*Fcal_SO_3) + E_SO_3*(-5.0*Fcal_4 - 5.0*Fcal_SQ_4) - 7.0*E_SO_5*Fcal_2 - 9.0*E_SO_7*Fcal_0 +
      E_SQ_4*(-6.0*Fcal_3 - 6.0*Fcal_SO_3)) + E_4*(-7.5*E_SO_3*Fcal_2 - 10.5*E_SO_5*Fcal_0) - 10.0*E_6*E_SO_3*Fcal_0 +
      E_SO_3*(E_SO_3*(-3.125*Fcal_3 - 3.125*Fcal_SO_3) - 7.5*E_SQ_4*Fcal_2) - 10.5*E_SO_5*E_SQ_4*Fcal_0)/m +
      ((E_2*(E_2*(E_2*(4.0*Fcal_3 + 4.0*Fcal_SO_3) + 15.0*E_SO_3*Fcal_2 + 21.0*E_SO_5*Fcal_0) + 45.0*E_4*E_SO_3*Fcal_0 +
      45.0*E_SO_3*E_SQ_4*Fcal_0) + 7.8125*pow(E_SO_3, 3)*Fcal_0)/m - 40.0*pow(E_2,
      3)*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0) + ((-0.5*Fcal_8 - 0.5*Fcal_SO_8 - 0.5*Fcal_lnv_8*logv)/m +
      ((E_2*(1.0*Fcal_6 + 1.0*Fcal_SO_6 + 1.0*Fcal_lnv_6*logv) + E_4*(1.5*Fcal_4 + 1.5*Fcal_SQ_4) + 2.0*E_6*Fcal_2 +
      2.5*E_8*Fcal_0 + E_SO_3*(1.25*Fcal_5 + 1.25*Fcal_SO_5) + E_SO_5*(1.75*Fcal_3 + 1.75*Fcal_SO_3) +
      E_SQ_4*(1.5*Fcal_4 + 1.5*Fcal_SQ_4) + E_lnv_8*Fcal_0*(2.5*logv + 0.25))/m + ((E_2*(E_2*(-2.0*Fcal_4 -
      2.0*Fcal_SQ_4) - 6.0*E_4*Fcal_2 - 8.0*E_6*Fcal_0 + E_SO_3*(-5.0*Fcal_3 - 5.0*Fcal_SO_3) - 6.0*E_SQ_4*Fcal_2) +
      E_4*(-4.5*E_4*Fcal_0 - 9.0*E_SQ_4*Fcal_0) + E_SO_3*(-3.125*E_SO_3*Fcal_2 - 8.75*E_SO_5*Fcal_0) - 4.5*pow(E_SQ_4,
      2)*Fcal_0)/m + (E_2*(E_2*(4.0*E_2*Fcal_2 + 18.0*E_4*Fcal_0 + 18.0*E_SQ_4*Fcal_0) + 18.75*pow(E_SO_3, 2)*Fcal_0)/m
      - 8.0*pow(E_2, 4)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0) + ((-0.5*Fcal_7 - 0.5*Fcal_SO_7)/m + ((E_2*(1.0*Fcal_5 +
      1.0*Fcal_SO_5) + E_4*(1.5*Fcal_3 + 1.5*Fcal_SO_3) + E_SO_3*(1.25*Fcal_4 + 1.25*Fcal_SQ_4) + 1.75*E_SO_5*Fcal_2 +
      2.25*E_SO_7*Fcal_0 + E_SQ_4*(1.5*Fcal_3 + 1.5*Fcal_SO_3))/m + ((E_2*(E_2*(-2.0*Fcal_3 - 2.0*Fcal_SO_3) -
      5.0*E_SO_3*Fcal_2 - 7.0*E_SO_5*Fcal_0) - 7.5*E_4*E_SO_3*Fcal_0 - 7.5*E_SO_3*E_SQ_4*Fcal_0)/m + 15.0*pow(E_2,
      2)*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0)/E_0) + ((-0.5*Fcal_6 - 0.5*Fcal_SO_6 - 0.5*Fcal_lnv_6*logv)/m +
      ((E_2*(1.0*Fcal_4 + 1.0*Fcal_SQ_4) + 1.5*E_4*Fcal_2 + 2.0*E_6*Fcal_0 + E_SO_3*(1.25*Fcal_3 + 1.25*Fcal_SO_3) +
      1.5*E_SQ_4*Fcal_2)/m + ((E_2*(-2.0*E_2*Fcal_2 - 6.0*E_4*Fcal_0 - 6.0*E_SQ_4*Fcal_0) - 3.125*pow(E_SO_3,
      2)*Fcal_0)/m + 4.0*pow(E_2, 3)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0) + ((-0.5*Fcal_5 - 0.5*Fcal_SO_5)/m +
      ((E_2*(1.0*Fcal_3 + 1.0*Fcal_SO_3) + 1.25*E_SO_3*Fcal_2 + 1.75*E_SO_5*Fcal_0)/m -
      5.0*E_2*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_4 - 0.5*Fcal_SQ_4)/m + ((1.0*E_2*Fcal_2 + 1.5*E_4*Fcal_0 +
      1.5*E_SQ_4*Fcal_0)/m - 2.0*pow(E_2, 2)*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_3 - 0.5*Fcal_SO_3)/m +
      1.25*E_SO_3*Fcal_0/(E_0*m))/E_0) + (-0.5*Fcal_2/m + 1.0*E_2*Fcal_0/(E_0*m))/E_0) - 0.5*Fcal_0/(E_0*m))/(nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = -0.5*nu*v*(-2.0*E_0*m/Fcal_0 + pow(v, 2)*(v*(v*(v*(v*(v*(v*(v*(m*v*(-12.0*E_10 +
      E_lnv_10*(-12.0*logv - 1.0) + (E_0*(2.0*Fcal_10 + 2.0*Fcal_lnv_10*logv) + E_2*(4.0*Fcal_8 + 4.0*Fcal_SO_8 +
      4.0*Fcal_lnv_8*logv) + E_4*(6.0*Fcal_6 + 6.0*Fcal_SO_6 + 6.0*Fcal_lnv_6*logv) + E_6*(8.0*Fcal_4 + 8.0*Fcal_SQ_4) +
      10.0*E_8*Fcal_2 + E_SO_3*(5.0*Fcal_7 + 5.0*Fcal_SO_7) + E_SO_5*(7.0*Fcal_5 + 7.0*Fcal_SO_5) + E_SO_7*(9.0*Fcal_3 +
      9.0*Fcal_SO_3) + E_SQ_4*(6.0*Fcal_6 + 6.0*Fcal_SO_6 + 6.0*Fcal_lnv_6*logv) + E_lnv_8*Fcal_2*(10.0*logv + 1.0) +
      (E_0*(Fcal_2*(-4.0*Fcal_8 - 4.0*Fcal_SO_8 - 4.0*Fcal_lnv_8*logv) + Fcal_3*(-4.0*Fcal_7 - 4.0*Fcal_SO_7) +
      Fcal_4*(-4.0*Fcal_6 - 4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) + Fcal_5*(-2.0*Fcal_5 - 4.0*Fcal_SO_5) -
      4.0*Fcal_6*Fcal_SQ_4 - 4.0*Fcal_7*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SO_7 - 2.0*pow(Fcal_SO_5, 2) -
      4.0*Fcal_SO_6*Fcal_SQ_4 - 4.0*Fcal_SQ_4*Fcal_lnv_6*logv) + E_2*(Fcal_2*(-8.0*Fcal_6 - 8.0*Fcal_SO_6 -
      8.0*Fcal_lnv_6*logv) + Fcal_3*(-8.0*Fcal_5 - 8.0*Fcal_SO_5) + Fcal_4*(-4.0*Fcal_4 - 8.0*Fcal_SQ_4) -
      8.0*Fcal_5*Fcal_SO_3 - 8.0*Fcal_SO_3*Fcal_SO_5 - 4.0*pow(Fcal_SQ_4, 2)) + E_4*(Fcal_2*(-12.0*Fcal_4 -
      12.0*Fcal_SQ_4) + Fcal_3*(-6.0*Fcal_3 - 12.0*Fcal_SO_3) - 6.0*pow(Fcal_SO_3, 2)) - 8.0*E_6*pow(Fcal_2, 2) +
      E_SO_3*(Fcal_2*(-10.0*Fcal_5 - 10.0*Fcal_SO_5) + Fcal_3*(-10.0*Fcal_4 - 10.0*Fcal_SQ_4) - 10.0*Fcal_4*Fcal_SO_3 -
      10.0*Fcal_SO_3*Fcal_SQ_4) + E_SO_5*Fcal_2*(-14.0*Fcal_3 - 14.0*Fcal_SO_3) + E_SQ_4*(Fcal_2*(-12.0*Fcal_4 -
      12.0*Fcal_SQ_4) + Fcal_3*(-6.0*Fcal_3 - 12.0*Fcal_SO_3) - 6.0*pow(Fcal_SO_3, 2)) +
      (E_0*(Fcal_2*(Fcal_2*(6.0*Fcal_6 + 6.0*Fcal_SO_6 + 6.0*Fcal_lnv_6*logv) + Fcal_3*(12.0*Fcal_5 + 12.0*Fcal_SO_5) +
      Fcal_4*(6.0*Fcal_4 + 12.0*Fcal_SQ_4) + 12.0*Fcal_5*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SO_5 + 6.0*pow(Fcal_SQ_4, 2)) +
      Fcal_3*(Fcal_3*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + 12.0*Fcal_4*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SQ_4) +
      6.0*Fcal_4*pow(Fcal_SO_3, 2) + 6.0*pow(Fcal_SO_3, 2)*Fcal_SQ_4) + E_2*Fcal_2*(Fcal_2*(12.0*Fcal_4 +
      12.0*Fcal_SQ_4) + Fcal_3*(12.0*Fcal_3 + 24.0*Fcal_SO_3) + 12.0*pow(Fcal_SO_3, 2)) + 6.0*E_4*pow(Fcal_2, 3) +
      E_SO_3*pow(Fcal_2, 2)*(15.0*Fcal_3 + 15.0*Fcal_SO_3) + 6.0*E_SQ_4*pow(Fcal_2, 3) + (E_0*pow(Fcal_2,
      2)*(Fcal_2*(-8.0*Fcal_4 - 8.0*Fcal_SQ_4) + Fcal_3*(-12.0*Fcal_3 - 24.0*Fcal_SO_3) - 12.0*pow(Fcal_SO_3, 2)) +
      2.0*E_0*pow(Fcal_2, 5)/Fcal_0 - 4.0*E_2*pow(Fcal_2, 4))/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0 +
      m*(E_0*(2.0*Fcal_9 + 2.0*Fcal_lnv_9*logv) + E_2*(4.0*Fcal_7 + 4.0*Fcal_SO_7) + E_4*(6.0*Fcal_5 + 6.0*Fcal_SO_5) +
      E_6*(8.0*Fcal_3 + 8.0*Fcal_SO_3) + E_SO_3*(5.0*Fcal_6 + 5.0*Fcal_SO_6 + 5.0*Fcal_lnv_6*logv) + E_SO_5*(7.0*Fcal_4
      + 7.0*Fcal_SQ_4) + 9.0*E_SO_7*Fcal_2 + E_SQ_4*(6.0*Fcal_5 + 6.0*Fcal_SO_5) + (E_0*(Fcal_2*(-4.0*Fcal_7 -
      4.0*Fcal_SO_7) + Fcal_3*(-4.0*Fcal_6 - 4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) + Fcal_4*(-4.0*Fcal_5 - 4.0*Fcal_SO_5)
      - 4.0*Fcal_5*Fcal_SQ_4 - 4.0*Fcal_6*Fcal_SO_3 + Fcal_SO_3*(-4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) -
      4.0*Fcal_SO_5*Fcal_SQ_4) + E_2*(Fcal_2*(-8.0*Fcal_5 - 8.0*Fcal_SO_5) + Fcal_3*(-8.0*Fcal_4 - 8.0*Fcal_SQ_4) -
      8.0*Fcal_4*Fcal_SO_3 - 8.0*Fcal_SO_3*Fcal_SQ_4) + E_4*Fcal_2*(-12.0*Fcal_3 - 12.0*Fcal_SO_3) +
      E_SO_3*(Fcal_2*(-10.0*Fcal_4 - 10.0*Fcal_SQ_4) + Fcal_3*(-5.0*Fcal_3 - 10.0*Fcal_SO_3) - 5.0*pow(Fcal_SO_3, 2)) -
      7.0*E_SO_5*pow(Fcal_2, 2) + E_SQ_4*Fcal_2*(-12.0*Fcal_3 - 12.0*Fcal_SO_3) + (E_0*(Fcal_2*(Fcal_2*(6.0*Fcal_5 +
      6.0*Fcal_SO_5) + Fcal_3*(12.0*Fcal_4 + 12.0*Fcal_SQ_4) + 12.0*Fcal_4*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SQ_4) +
      Fcal_3*(Fcal_3*(2.0*Fcal_3 + 6.0*Fcal_SO_3) + 6.0*pow(Fcal_SO_3, 2)) + 2.0*pow(Fcal_SO_3, 3)) + E_0*pow(Fcal_2,
      3)*(-8.0*Fcal_3 - 8.0*Fcal_SO_3)/Fcal_0 + E_2*pow(Fcal_2, 2)*(12.0*Fcal_3 + 12.0*Fcal_SO_3) +
      5.0*E_SO_3*pow(Fcal_2, 3))/Fcal_0)/Fcal_0)/pow(Fcal_0, 2)) + m*(-10.0*E_8 + E_lnv_8*(-10.0*logv - 1.0) +
      (E_0*(2.0*Fcal_8 + 2.0*Fcal_SO_8 + 2.0*Fcal_lnv_8*logv) + E_2*(4.0*Fcal_6 + 4.0*Fcal_SO_6 + 4.0*Fcal_lnv_6*logv) +
      E_4*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + 8.0*E_6*Fcal_2 + E_SO_3*(5.0*Fcal_5 + 5.0*Fcal_SO_5) + E_SO_5*(7.0*Fcal_3 +
      7.0*Fcal_SO_3) + E_SQ_4*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + (E_0*(Fcal_2*(-4.0*Fcal_6 - 4.0*Fcal_SO_6 -
      4.0*Fcal_lnv_6*logv) + Fcal_3*(-4.0*Fcal_5 - 4.0*Fcal_SO_5) + Fcal_4*(-2.0*Fcal_4 - 4.0*Fcal_SQ_4) -
      4.0*Fcal_5*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SO_5 - 2.0*pow(Fcal_SQ_4, 2)) + E_2*(Fcal_2*(-8.0*Fcal_4 -
      8.0*Fcal_SQ_4) + Fcal_3*(-4.0*Fcal_3 - 8.0*Fcal_SO_3) - 4.0*pow(Fcal_SO_3, 2)) - 6.0*E_4*pow(Fcal_2, 2) +
      E_SO_3*Fcal_2*(-10.0*Fcal_3 - 10.0*Fcal_SO_3) - 6.0*E_SQ_4*pow(Fcal_2, 2) + (E_0*Fcal_2*(Fcal_2*(6.0*Fcal_4 +
      6.0*Fcal_SQ_4) + Fcal_3*(6.0*Fcal_3 + 12.0*Fcal_SO_3) + 6.0*pow(Fcal_SO_3, 2)) - 2.0*E_0*pow(Fcal_2, 4)/Fcal_0 +
      4.0*E_2*pow(Fcal_2, 3))/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0) + m*(-9.0*E_SO_7 + (E_0*(2.0*Fcal_7 + 2.0*Fcal_SO_7) +
      E_2*(4.0*Fcal_5 + 4.0*Fcal_SO_5) + E_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) + E_SO_3*(5.0*Fcal_4 + 5.0*Fcal_SQ_4) +
      7.0*E_SO_5*Fcal_2 + E_SQ_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) + (E_0*(Fcal_2*(-4.0*Fcal_5 - 4.0*Fcal_SO_5) +
      Fcal_3*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) - 4.0*Fcal_4*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SQ_4) + E_0*pow(Fcal_2,
      2)*(6.0*Fcal_3 + 6.0*Fcal_SO_3)/Fcal_0 + E_2*Fcal_2*(-8.0*Fcal_3 - 8.0*Fcal_SO_3) - 5.0*E_SO_3*pow(Fcal_2,
      2))/Fcal_0)/Fcal_0)/Fcal_0) + m*(-8.0*E_6 + (E_0*(2.0*Fcal_6 + 2.0*Fcal_SO_6 + 2.0*Fcal_lnv_6*logv) +
      E_2*(4.0*Fcal_4 + 4.0*Fcal_SQ_4) + 6.0*E_4*Fcal_2 + E_SO_3*(5.0*Fcal_3 + 5.0*Fcal_SO_3) + 6.0*E_SQ_4*Fcal_2 +
      (E_0*(Fcal_2*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) + Fcal_3*(-2.0*Fcal_3 - 4.0*Fcal_SO_3) - 2.0*pow(Fcal_SO_3, 2)) +
      2.0*E_0*pow(Fcal_2, 3)/Fcal_0 - 4.0*E_2*pow(Fcal_2, 2))/Fcal_0)/Fcal_0)/Fcal_0) + m*(-7.0*E_SO_5 +
      (E_0*(2.0*Fcal_5 + 2.0*Fcal_SO_5) + E_0*Fcal_2*(-4.0*Fcal_3 - 4.0*Fcal_SO_3)/Fcal_0 + E_2*(4.0*Fcal_3 +
      4.0*Fcal_SO_3) + 5.0*E_SO_3*Fcal_2)/Fcal_0)/Fcal_0) + m*(-6.0*E_4 - 6.0*E_SQ_4 + (E_0*(2.0*Fcal_4 + 2.0*Fcal_SQ_4)
      - 2.0*E_0*pow(Fcal_2, 2)/Fcal_0 + 4.0*E_2*Fcal_2)/Fcal_0)/Fcal_0) + m*(E_0*(2.0*Fcal_3 + 2.0*Fcal_SO_3)/Fcal_0 -
      5.0*E_SO_3)/Fcal_0) + m*(2.0*E_0*Fcal_2/Fcal_0 - 4.0*E_2)/Fcal_0))/Fcal_coeff;
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_5p0PN_Q : public TaylorTn_Q


class TaylorTn_5p5PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, logv, Fcal_coeff;
  const double Fcal_0, Fcal_2, Fcal_3, Fcal_4, Fcal_5, Fcal_6, Fcal_lnv_6, Fcal_7, Fcal_8, Fcal_lnv_8, Fcal_9,
               Fcal_lnv_9, Fcal_10, Fcal_lnv_10, Fcal_11, Fcal_lnv_11;
  double Fcal_SQ_4, Fcal_SO_3, Fcal_SO_5, Fcal_SO_6, Fcal_SO_7, Fcal_SO_8;
  const double E_0, E_2, E_4, E_6, E_8, E_lnv_8, E_10, E_lnv_10, E_11;
  double E_SQ_4, E_SO_3, E_SO_5, E_SO_7;
  double Phi;

public:
  TaylorTn_5p5PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), logv(log(v)), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    Fcal_2(-2.91666666666667*nu - 3.71130952380952), Fcal_3(12.5663706143592), Fcal_4(3.61111111111111*pow(nu, 2) +
    18.3948412698413*nu - 4.92846119929453), Fcal_5(-76.3145215434521*nu - 38.2928354546934),
    Fcal_6(-2.39197530864198*pow(nu, 3) - 31.2179232804233*pow(nu, 2) - 8.87205344238227*nu + 115.731716675611),
    Fcal_lnv_6(-16.3047619047619), Fcal_7(200.905057974359*pow(nu, 2) + 390.417427312002*nu - 101.509595959742),
    Fcal_8(-117.504390722677), Fcal_lnv_8(52.7430839002268), Fcal_9(719.128342233430), Fcal_lnv_9(-204.891680874123),
    Fcal_10(-1216.90699131704), Fcal_lnv_10(116.639876594109), Fcal_11(958.934970119567), Fcal_lnv_11(473.624478174231),
    Fcal_SQ_4(chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
    2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
    chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
    2)*(0.0416666666666667*nu + 2.98958333333333)), Fcal_SO_3((-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2)),
    Fcal_SO_5((S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2)),
    Fcal_SO_6((-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2)),
    Fcal_SO_7((S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
    Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762))/pow(m, 2)),
    Fcal_SO_8((3.14159265358979*S_l*(192.763888888889*nu - 36.3020833333333) +
    3.14159265358979*Sigma_l*delta*(64.7733134920635*nu - 10.6592261904762))/pow(m, 2)), E_0(1.00000000000000),
    E_2(-0.0833333333333333*nu - 0.75), E_4(-0.0416666666666667*pow(nu, 2) + 2.375*nu - 3.375),
    E_6(-0.00675154320987654*pow(nu, 3) - 1.61458333333333*pow(nu, 2) + 38.7246294907293*nu - 10.546875),
    E_8(0.0024755658436214*pow(nu, 4) + 0.174189814814815*pow(nu, 3) - 90.1327990262052*pow(nu, 2) + 153.88379682994*nu
    - 31.0078125), E_lnv_8(59.7333333333333*nu), E_10(0.001953125*pow(nu, 5) + 0.107421875*pow(nu, 4) +
    83.4293954895551*pow(nu, 3) - 71.3641942901125*pow(nu, 2) - 62.1764445589718*nu - 89.701171875),
    E_lnv_10(-262.4*pow(nu, 2) - 285.028571428571*nu), E_11(273.188907832164*nu), E_SQ_4(-1.5*pow(chi_a_l, 2) -
    1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 + 6.0*pow(chi_a_l, 2)) +
    0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0)), E_SO_3((4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2)),
    E_SO_5((S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2)),
    E_SO_7((S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75) + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu +
    6.75))/pow(m, 2)), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    logv = log(v);
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
    Fcal_SQ_4 = chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
      2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
      chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
      2)*(0.0416666666666667*nu + 2.98958333333333);
    Fcal_SO_3 = (-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_5 = (S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2);
    Fcal_SO_6 = (-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_7 = (S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
      Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762))/pow(m, 2);
    Fcal_SO_8 = (3.14159265358979*S_l*(192.763888888889*nu - 36.3020833333333) +
      3.14159265358979*Sigma_l*delta*(64.7733134920635*nu - 10.6592261904762))/pow(m, 2);
    E_SQ_4 = -1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 +
      6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0);
    E_SO_3 = (4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2);
    E_SO_5 = (S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2);
    E_SO_7 = (S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75) + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu +
      6.75))/pow(m, 2);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*(-chiVec2*pow(m2, 2)*v/pow(m, 2) + ellHat*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu -
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu +
      3.0*chi2_n*pow(m2, 2)/pow(m, 2)));
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*(-chiVec1*pow(m1, 2)*v/pow(m, 2) + ellHat*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu +
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*pow(m1,
      2)/pow(m, 2) + 3.0*chi2_n*nu));
  }
  Quaternion OmegaVec() {
    double gamma_PN_2 = -0.333333333333333*nu + 1.0;
    double gamma_PN_6 = 0.0123456790123457*pow(nu, 3) + 6.36111111111111*pow(nu, 2) - 2.98177812235564*nu + 1.0;
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_3 = (1.66666666666667*S_l + Sigma_l*delta)/pow(m, 2);
    double gamma_PN_7 = (S_l*(-6.0*pow(nu, 2) - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_l*delta*pow(nu, 2) +
      Sigma_l*delta*(-10.1666666666667*nu + 3.0))/pow(m, 2);
    double gamma_PN_5 = (S_l*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_l*delta)/pow(m, 2);
    double gamma_PN_0 = 1.00000000000000;
    double a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0);
    double gamma_PN_4 = -5.41666666666667*nu + 1.0;
    double a_ell_4 = S_n*(5.77777777777778*pow(nu, 2) + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*pow(nu, 2) +
      9.125*nu + 1.5);
    return ellHat*pow(v, 3)/m + nHat*pow(v, 6)*(a_ell_0 + pow(v, 2)*(a_ell_2 + a_ell_4*pow(v, 2)))*(gamma_PN_0 + pow(v,
      2)*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/pow(m, 3);
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_coeff*(Fcal_0 + pow(v, 2)*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 +
      v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 +
      Fcal_lnv_8*logv + v*(Fcal_9 + Fcal_lnv_9*logv + v*(Fcal_10 + Fcal_lnv_10*logv + v*(Fcal_11 +
      Fcal_lnv_11*logv)))))))))));
    const double dEdv = -0.5*m*nu*v*(2.0*E_0 + pow(v, 2)*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 +
      v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0) + pow(v, 2)*(12.0*E_10 +
      13.0*E_11*v + E_lnv_10*(12.0*logv + 1.0))))))))));
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = -2.0*Fcal_coeff*(pow(v, 2)*(v*(v*(v*(v*(v*(v*(v*(v*(v*((-0.5*Fcal_11 -
      0.5*Fcal_lnv_11*logv)/m + ((3.25*E_11*Fcal_0 + E_2*(1.0*Fcal_9 + 1.0*Fcal_lnv_9*logv) + E_4*(1.5*Fcal_7 +
      1.5*Fcal_SO_7) + E_6*(2.0*Fcal_5 + 2.0*Fcal_SO_5) + E_8*(2.5*Fcal_3 + 2.5*Fcal_SO_3) + E_SO_3*(1.25*Fcal_8 +
      1.25*Fcal_SO_8 + 1.25*Fcal_lnv_8*logv) + E_SO_5*(1.75*Fcal_6 + 1.75*Fcal_SO_6 + 1.75*Fcal_lnv_6*logv) +
      E_SO_7*(2.25*Fcal_4 + 2.25*Fcal_SQ_4) + E_SQ_4*(1.5*Fcal_7 + 1.5*Fcal_SO_7) + E_lnv_8*(Fcal_3*(2.5*logv + 0.25) +
      Fcal_SO_3*(2.5*logv + 0.25)))/m + ((E_2*(E_2*(-2.0*Fcal_7 - 2.0*Fcal_SO_7) + E_4*(-6.0*Fcal_5 - 6.0*Fcal_SO_5) +
      E_6*(-8.0*Fcal_3 - 8.0*Fcal_SO_3) + E_SO_3*(-5.0*Fcal_6 - 5.0*Fcal_SO_6 - 5.0*Fcal_lnv_6*logv) +
      E_SO_5*(-7.0*Fcal_4 - 7.0*Fcal_SQ_4) - 9.0*E_SO_7*Fcal_2 + E_SQ_4*(-6.0*Fcal_5 - 6.0*Fcal_SO_5)) +
      E_4*(E_4*(-4.5*Fcal_3 - 4.5*Fcal_SO_3) + E_SO_3*(-7.5*Fcal_4 - 7.5*Fcal_SQ_4) - 10.5*E_SO_5*Fcal_2 -
      13.5*E_SO_7*Fcal_0 + E_SQ_4*(-9.0*Fcal_3 - 9.0*Fcal_SO_3)) + E_6*(-10.0*E_SO_3*Fcal_2 - 14.0*E_SO_5*Fcal_0) -
      12.5*E_8*E_SO_3*Fcal_0 + E_SO_3*(E_SO_3*(-3.125*Fcal_5 - 3.125*Fcal_SO_5) + E_SO_5*(-8.75*Fcal_3 - 8.75*Fcal_SO_3)
      + E_SQ_4*(-7.5*Fcal_4 - 7.5*Fcal_SQ_4) + E_lnv_8*Fcal_0*(-12.5*logv - 1.25)) - 10.5*E_SO_5*E_SQ_4*Fcal_2 -
      13.5*E_SO_7*E_SQ_4*Fcal_0 + pow(E_SQ_4, 2)*(-4.5*Fcal_3 - 4.5*Fcal_SO_3))/m + ((E_2*(E_2*(E_2*(4.0*Fcal_5 +
      4.0*Fcal_SO_5) + E_4*(18.0*Fcal_3 + 18.0*Fcal_SO_3) + E_SO_3*(15.0*Fcal_4 + 15.0*Fcal_SQ_4) + 21.0*E_SO_5*Fcal_2 +
      27.0*E_SO_7*Fcal_0 + E_SQ_4*(18.0*Fcal_3 + 18.0*Fcal_SO_3)) + E_4*(45.0*E_SO_3*Fcal_2 + 63.0*E_SO_5*Fcal_0) +
      60.0*E_6*E_SO_3*Fcal_0 + E_SO_3*(E_SO_3*(18.75*Fcal_3 + 18.75*Fcal_SO_3) + 45.0*E_SQ_4*Fcal_2) +
      63.0*E_SO_5*E_SQ_4*Fcal_0) + E_4*(33.75*E_4*E_SO_3*Fcal_0 + 67.5*E_SO_3*E_SQ_4*Fcal_0) +
      E_SO_3*(E_SO_3*(7.8125*E_SO_3*Fcal_2 + 32.8125*E_SO_5*Fcal_0) + 33.75*pow(E_SQ_4, 2)*Fcal_0))/m +
      (E_2*(E_2*(E_2*(E_2*(-8.0*Fcal_3 - 8.0*Fcal_SO_3) - 40.0*E_SO_3*Fcal_2 - 56.0*E_SO_5*Fcal_0) -
      180.0*E_4*E_SO_3*Fcal_0 - 180.0*E_SO_3*E_SQ_4*Fcal_0) - 62.5*pow(E_SO_3, 3)*Fcal_0)/m + 100.0*pow(E_2,
      4)*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0)/E_0 + ((-0.5*Fcal_10 - 0.5*Fcal_lnv_10*logv)/m + ((3.0*E_10*Fcal_0 +
      E_2*(1.0*Fcal_8 + 1.0*Fcal_SO_8 + 1.0*Fcal_lnv_8*logv) + E_4*(1.5*Fcal_6 + 1.5*Fcal_SO_6 + 1.5*Fcal_lnv_6*logv) +
      E_6*(2.0*Fcal_4 + 2.0*Fcal_SQ_4) + 2.5*E_8*Fcal_2 + E_SO_3*(1.25*Fcal_7 + 1.25*Fcal_SO_7) + E_SO_5*(1.75*Fcal_5 +
      1.75*Fcal_SO_5) + E_SO_7*(2.25*Fcal_3 + 2.25*Fcal_SO_3) + E_SQ_4*(1.5*Fcal_6 + 1.5*Fcal_SO_6 +
      1.5*Fcal_lnv_6*logv) + E_lnv_10*Fcal_0*(3.0*logv + 0.25) + E_lnv_8*Fcal_2*(2.5*logv + 0.25))/m +
      ((E_2*(E_2*(-2.0*Fcal_6 - 2.0*Fcal_SO_6 - 2.0*Fcal_lnv_6*logv) + E_4*(-6.0*Fcal_4 - 6.0*Fcal_SQ_4) -
      8.0*E_6*Fcal_2 - 10.0*E_8*Fcal_0 + E_SO_3*(-5.0*Fcal_5 - 5.0*Fcal_SO_5) + E_SO_5*(-7.0*Fcal_3 - 7.0*Fcal_SO_3) +
      E_SQ_4*(-6.0*Fcal_4 - 6.0*Fcal_SQ_4) + E_lnv_8*Fcal_0*(-10.0*logv - 1.0)) + E_4*(-4.5*E_4*Fcal_2 - 12.0*E_6*Fcal_0
      + E_SO_3*(-7.5*Fcal_3 - 7.5*Fcal_SO_3) - 9.0*E_SQ_4*Fcal_2) - 12.0*E_6*E_SQ_4*Fcal_0 +
      E_SO_3*(E_SO_3*(-3.125*Fcal_4 - 3.125*Fcal_SQ_4) - 8.75*E_SO_5*Fcal_2 - 11.25*E_SO_7*Fcal_0 + E_SQ_4*(-7.5*Fcal_3
      - 7.5*Fcal_SO_3)) - 6.125*pow(E_SO_5, 2)*Fcal_0 - 4.5*pow(E_SQ_4, 2)*Fcal_2)/m + ((E_2*(E_2*(E_2*(4.0*Fcal_4 +
      4.0*Fcal_SQ_4) + 18.0*E_4*Fcal_2 + 24.0*E_6*Fcal_0 + E_SO_3*(15.0*Fcal_3 + 15.0*Fcal_SO_3) + 18.0*E_SQ_4*Fcal_2) +
      E_4*(27.0*E_4*Fcal_0 + 54.0*E_SQ_4*Fcal_0) + E_SO_3*(18.75*E_SO_3*Fcal_2 + 52.5*E_SO_5*Fcal_0) + 27.0*pow(E_SQ_4,
      2)*Fcal_0) + 28.125*E_4*pow(E_SO_3, 2)*Fcal_0 + 28.125*pow(E_SO_3, 2)*E_SQ_4*Fcal_0)/m + (pow(E_2,
      2)*(E_2*(-8.0*E_2*Fcal_2 - 48.0*E_4*Fcal_0 - 48.0*E_SQ_4*Fcal_0) - 75.0*pow(E_SO_3, 2)*Fcal_0)/m + 16.0*pow(E_2,
      5)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0)/E_0) + ((-0.5*Fcal_9 - 0.5*Fcal_lnv_9*logv)/m + ((E_2*(1.0*Fcal_7 +
      1.0*Fcal_SO_7) + E_4*(1.5*Fcal_5 + 1.5*Fcal_SO_5) + E_6*(2.0*Fcal_3 + 2.0*Fcal_SO_3) + E_SO_3*(1.25*Fcal_6 +
      1.25*Fcal_SO_6 + 1.25*Fcal_lnv_6*logv) + E_SO_5*(1.75*Fcal_4 + 1.75*Fcal_SQ_4) + 2.25*E_SO_7*Fcal_2 +
      E_SQ_4*(1.5*Fcal_5 + 1.5*Fcal_SO_5))/m + ((E_2*(E_2*(-2.0*Fcal_5 - 2.0*Fcal_SO_5) + E_4*(-6.0*Fcal_3 -
      6.0*Fcal_SO_3) + E_SO_3*(-5.0*Fcal_4 - 5.0*Fcal_SQ_4) - 7.0*E_SO_5*Fcal_2 - 9.0*E_SO_7*Fcal_0 +
      E_SQ_4*(-6.0*Fcal_3 - 6.0*Fcal_SO_3)) + E_4*(-7.5*E_SO_3*Fcal_2 - 10.5*E_SO_5*Fcal_0) - 10.0*E_6*E_SO_3*Fcal_0 +
      E_SO_3*(E_SO_3*(-3.125*Fcal_3 - 3.125*Fcal_SO_3) - 7.5*E_SQ_4*Fcal_2) - 10.5*E_SO_5*E_SQ_4*Fcal_0)/m +
      ((E_2*(E_2*(E_2*(4.0*Fcal_3 + 4.0*Fcal_SO_3) + 15.0*E_SO_3*Fcal_2 + 21.0*E_SO_5*Fcal_0) + 45.0*E_4*E_SO_3*Fcal_0 +
      45.0*E_SO_3*E_SQ_4*Fcal_0) + 7.8125*pow(E_SO_3, 3)*Fcal_0)/m - 40.0*pow(E_2,
      3)*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0) + ((-0.5*Fcal_8 - 0.5*Fcal_SO_8 - 0.5*Fcal_lnv_8*logv)/m +
      ((E_2*(1.0*Fcal_6 + 1.0*Fcal_SO_6 + 1.0*Fcal_lnv_6*logv) + E_4*(1.5*Fcal_4 + 1.5*Fcal_SQ_4) + 2.0*E_6*Fcal_2 +
      2.5*E_8*Fcal_0 + E_SO_3*(1.25*Fcal_5 + 1.25*Fcal_SO_5) + E_SO_5*(1.75*Fcal_3 + 1.75*Fcal_SO_3) +
      E_SQ_4*(1.5*Fcal_4 + 1.5*Fcal_SQ_4) + E_lnv_8*Fcal_0*(2.5*logv + 0.25))/m + ((E_2*(E_2*(-2.0*Fcal_4 -
      2.0*Fcal_SQ_4) - 6.0*E_4*Fcal_2 - 8.0*E_6*Fcal_0 + E_SO_3*(-5.0*Fcal_3 - 5.0*Fcal_SO_3) - 6.0*E_SQ_4*Fcal_2) +
      E_4*(-4.5*E_4*Fcal_0 - 9.0*E_SQ_4*Fcal_0) + E_SO_3*(-3.125*E_SO_3*Fcal_2 - 8.75*E_SO_5*Fcal_0) - 4.5*pow(E_SQ_4,
      2)*Fcal_0)/m + (E_2*(E_2*(4.0*E_2*Fcal_2 + 18.0*E_4*Fcal_0 + 18.0*E_SQ_4*Fcal_0) + 18.75*pow(E_SO_3, 2)*Fcal_0)/m
      - 8.0*pow(E_2, 4)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0) + ((-0.5*Fcal_7 - 0.5*Fcal_SO_7)/m + ((E_2*(1.0*Fcal_5 +
      1.0*Fcal_SO_5) + E_4*(1.5*Fcal_3 + 1.5*Fcal_SO_3) + E_SO_3*(1.25*Fcal_4 + 1.25*Fcal_SQ_4) + 1.75*E_SO_5*Fcal_2 +
      2.25*E_SO_7*Fcal_0 + E_SQ_4*(1.5*Fcal_3 + 1.5*Fcal_SO_3))/m + ((E_2*(E_2*(-2.0*Fcal_3 - 2.0*Fcal_SO_3) -
      5.0*E_SO_3*Fcal_2 - 7.0*E_SO_5*Fcal_0) - 7.5*E_4*E_SO_3*Fcal_0 - 7.5*E_SO_3*E_SQ_4*Fcal_0)/m + 15.0*pow(E_2,
      2)*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0)/E_0) + ((-0.5*Fcal_6 - 0.5*Fcal_SO_6 - 0.5*Fcal_lnv_6*logv)/m +
      ((E_2*(1.0*Fcal_4 + 1.0*Fcal_SQ_4) + 1.5*E_4*Fcal_2 + 2.0*E_6*Fcal_0 + E_SO_3*(1.25*Fcal_3 + 1.25*Fcal_SO_3) +
      1.5*E_SQ_4*Fcal_2)/m + ((E_2*(-2.0*E_2*Fcal_2 - 6.0*E_4*Fcal_0 - 6.0*E_SQ_4*Fcal_0) - 3.125*pow(E_SO_3,
      2)*Fcal_0)/m + 4.0*pow(E_2, 3)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0) + ((-0.5*Fcal_5 - 0.5*Fcal_SO_5)/m +
      ((E_2*(1.0*Fcal_3 + 1.0*Fcal_SO_3) + 1.25*E_SO_3*Fcal_2 + 1.75*E_SO_5*Fcal_0)/m -
      5.0*E_2*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_4 - 0.5*Fcal_SQ_4)/m + ((1.0*E_2*Fcal_2 + 1.5*E_4*Fcal_0 +
      1.5*E_SQ_4*Fcal_0)/m - 2.0*pow(E_2, 2)*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_3 - 0.5*Fcal_SO_3)/m +
      1.25*E_SO_3*Fcal_0/(E_0*m))/E_0) + (-0.5*Fcal_2/m + 1.0*E_2*Fcal_0/(E_0*m))/E_0) - 0.5*Fcal_0/(E_0*m))/(nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = -0.5*nu*v*(-2.0*E_0*m/Fcal_0 + pow(v, 2)*(v*(v*(v*(v*(v*(v*(v*(v*(m*v*(-13.0*E_11 +
      (E_0*(2.0*Fcal_11 + 2.0*Fcal_lnv_11*logv) + E_2*(4.0*Fcal_9 + 4.0*Fcal_lnv_9*logv) + E_4*(6.0*Fcal_7 +
      6.0*Fcal_SO_7) + E_6*(8.0*Fcal_5 + 8.0*Fcal_SO_5) + E_8*(10.0*Fcal_3 + 10.0*Fcal_SO_3) + E_SO_3*(5.0*Fcal_8 +
      5.0*Fcal_SO_8 + 5.0*Fcal_lnv_8*logv) + E_SO_5*(7.0*Fcal_6 + 7.0*Fcal_SO_6 + 7.0*Fcal_lnv_6*logv) +
      E_SO_7*(9.0*Fcal_4 + 9.0*Fcal_SQ_4) + E_SQ_4*(6.0*Fcal_7 + 6.0*Fcal_SO_7) + E_lnv_8*(Fcal_3*(10.0*logv + 1.0) +
      Fcal_SO_3*(10.0*logv + 1.0)) + (E_0*(Fcal_2*(-4.0*Fcal_9 - 4.0*Fcal_lnv_9*logv) + Fcal_3*(-4.0*Fcal_8 -
      4.0*Fcal_SO_8 - 4.0*Fcal_lnv_8*logv) + Fcal_4*(-4.0*Fcal_7 - 4.0*Fcal_SO_7) + Fcal_5*(-4.0*Fcal_6 - 4.0*Fcal_SO_6
      - 4.0*Fcal_lnv_6*logv) - 4.0*Fcal_6*Fcal_SO_5 - 4.0*Fcal_7*Fcal_SQ_4 - 4.0*Fcal_8*Fcal_SO_3 +
      Fcal_SO_3*(-4.0*Fcal_SO_8 - 4.0*Fcal_lnv_8*logv) + Fcal_SO_5*(-4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) -
      4.0*Fcal_SO_7*Fcal_SQ_4) + E_2*(Fcal_2*(-8.0*Fcal_7 - 8.0*Fcal_SO_7) + Fcal_3*(-8.0*Fcal_6 - 8.0*Fcal_SO_6 -
      8.0*Fcal_lnv_6*logv) + Fcal_4*(-8.0*Fcal_5 - 8.0*Fcal_SO_5) - 8.0*Fcal_5*Fcal_SQ_4 - 8.0*Fcal_6*Fcal_SO_3 +
      Fcal_SO_3*(-8.0*Fcal_SO_6 - 8.0*Fcal_lnv_6*logv) - 8.0*Fcal_SO_5*Fcal_SQ_4) + E_4*(Fcal_2*(-12.0*Fcal_5 -
      12.0*Fcal_SO_5) + Fcal_3*(-12.0*Fcal_4 - 12.0*Fcal_SQ_4) - 12.0*Fcal_4*Fcal_SO_3 - 12.0*Fcal_SO_3*Fcal_SQ_4) +
      E_6*Fcal_2*(-16.0*Fcal_3 - 16.0*Fcal_SO_3) + E_SO_3*(Fcal_2*(-10.0*Fcal_6 - 10.0*Fcal_SO_6 - 10.0*Fcal_lnv_6*logv)
      + Fcal_3*(-10.0*Fcal_5 - 10.0*Fcal_SO_5) + Fcal_4*(-5.0*Fcal_4 - 10.0*Fcal_SQ_4) - 10.0*Fcal_5*Fcal_SO_3 -
      10.0*Fcal_SO_3*Fcal_SO_5 - 5.0*pow(Fcal_SQ_4, 2)) + E_SO_5*(Fcal_2*(-14.0*Fcal_4 - 14.0*Fcal_SQ_4) +
      Fcal_3*(-7.0*Fcal_3 - 14.0*Fcal_SO_3) - 7.0*pow(Fcal_SO_3, 2)) - 9.0*E_SO_7*pow(Fcal_2, 2) +
      E_SQ_4*(Fcal_2*(-12.0*Fcal_5 - 12.0*Fcal_SO_5) + Fcal_3*(-12.0*Fcal_4 - 12.0*Fcal_SQ_4) - 12.0*Fcal_4*Fcal_SO_3 -
      12.0*Fcal_SO_3*Fcal_SQ_4) + (E_0*(Fcal_2*(Fcal_2*(6.0*Fcal_7 + 6.0*Fcal_SO_7) + Fcal_3*(12.0*Fcal_6 +
      12.0*Fcal_SO_6 + 12.0*Fcal_lnv_6*logv) + Fcal_4*(12.0*Fcal_5 + 12.0*Fcal_SO_5) + 12.0*Fcal_5*Fcal_SQ_4 +
      12.0*Fcal_6*Fcal_SO_3 + Fcal_SO_3*(12.0*Fcal_SO_6 + 12.0*Fcal_lnv_6*logv) + 12.0*Fcal_SO_5*Fcal_SQ_4) +
      Fcal_3*(Fcal_3*(6.0*Fcal_5 + 6.0*Fcal_SO_5) + Fcal_4*(6.0*Fcal_4 + 12.0*Fcal_SQ_4) + 12.0*Fcal_5*Fcal_SO_3 +
      12.0*Fcal_SO_3*Fcal_SO_5 + 6.0*pow(Fcal_SQ_4, 2)) + Fcal_4*(6.0*Fcal_4*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SQ_4) +
      6.0*Fcal_5*pow(Fcal_SO_3, 2) + Fcal_SO_3*(6.0*Fcal_SO_3*Fcal_SO_5 + 6.0*pow(Fcal_SQ_4, 2))) +
      E_2*(Fcal_2*(Fcal_2*(12.0*Fcal_5 + 12.0*Fcal_SO_5) + Fcal_3*(24.0*Fcal_4 + 24.0*Fcal_SQ_4) + 24.0*Fcal_4*Fcal_SO_3
      + 24.0*Fcal_SO_3*Fcal_SQ_4) + Fcal_3*(Fcal_3*(4.0*Fcal_3 + 12.0*Fcal_SO_3) + 12.0*pow(Fcal_SO_3, 2)) +
      4.0*pow(Fcal_SO_3, 3)) + E_4*pow(Fcal_2, 2)*(18.0*Fcal_3 + 18.0*Fcal_SO_3) + E_SO_3*Fcal_2*(Fcal_2*(15.0*Fcal_4 +
      15.0*Fcal_SQ_4) + Fcal_3*(15.0*Fcal_3 + 30.0*Fcal_SO_3) + 15.0*pow(Fcal_SO_3, 2)) + 7.0*E_SO_5*pow(Fcal_2, 3) +
      E_SQ_4*pow(Fcal_2, 2)*(18.0*Fcal_3 + 18.0*Fcal_SO_3) + (E_0*Fcal_2*(Fcal_2*(Fcal_2*(-8.0*Fcal_5 - 8.0*Fcal_SO_5) +
      Fcal_3*(-24.0*Fcal_4 - 24.0*Fcal_SQ_4) - 24.0*Fcal_4*Fcal_SO_3 - 24.0*Fcal_SO_3*Fcal_SQ_4) +
      Fcal_3*(Fcal_3*(-8.0*Fcal_3 - 24.0*Fcal_SO_3) - 24.0*pow(Fcal_SO_3, 2)) - 8.0*pow(Fcal_SO_3, 3)) + E_0*pow(Fcal_2,
      4)*(10.0*Fcal_3 + 10.0*Fcal_SO_3)/Fcal_0 + E_2*pow(Fcal_2, 3)*(-16.0*Fcal_3 - 16.0*Fcal_SO_3) -
      5.0*E_SO_3*pow(Fcal_2, 4))/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0 + m*(-12.0*E_10 + E_lnv_10*(-12.0*logv - 1.0) +
      (E_0*(2.0*Fcal_10 + 2.0*Fcal_lnv_10*logv) + E_2*(4.0*Fcal_8 + 4.0*Fcal_SO_8 + 4.0*Fcal_lnv_8*logv) +
      E_4*(6.0*Fcal_6 + 6.0*Fcal_SO_6 + 6.0*Fcal_lnv_6*logv) + E_6*(8.0*Fcal_4 + 8.0*Fcal_SQ_4) + 10.0*E_8*Fcal_2 +
      E_SO_3*(5.0*Fcal_7 + 5.0*Fcal_SO_7) + E_SO_5*(7.0*Fcal_5 + 7.0*Fcal_SO_5) + E_SO_7*(9.0*Fcal_3 + 9.0*Fcal_SO_3) +
      E_SQ_4*(6.0*Fcal_6 + 6.0*Fcal_SO_6 + 6.0*Fcal_lnv_6*logv) + E_lnv_8*Fcal_2*(10.0*logv + 1.0) +
      (E_0*(Fcal_2*(-4.0*Fcal_8 - 4.0*Fcal_SO_8 - 4.0*Fcal_lnv_8*logv) + Fcal_3*(-4.0*Fcal_7 - 4.0*Fcal_SO_7) +
      Fcal_4*(-4.0*Fcal_6 - 4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) + Fcal_5*(-2.0*Fcal_5 - 4.0*Fcal_SO_5) -
      4.0*Fcal_6*Fcal_SQ_4 - 4.0*Fcal_7*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SO_7 - 2.0*pow(Fcal_SO_5, 2) -
      4.0*Fcal_SO_6*Fcal_SQ_4 - 4.0*Fcal_SQ_4*Fcal_lnv_6*logv) + E_2*(Fcal_2*(-8.0*Fcal_6 - 8.0*Fcal_SO_6 -
      8.0*Fcal_lnv_6*logv) + Fcal_3*(-8.0*Fcal_5 - 8.0*Fcal_SO_5) + Fcal_4*(-4.0*Fcal_4 - 8.0*Fcal_SQ_4) -
      8.0*Fcal_5*Fcal_SO_3 - 8.0*Fcal_SO_3*Fcal_SO_5 - 4.0*pow(Fcal_SQ_4, 2)) + E_4*(Fcal_2*(-12.0*Fcal_4 -
      12.0*Fcal_SQ_4) + Fcal_3*(-6.0*Fcal_3 - 12.0*Fcal_SO_3) - 6.0*pow(Fcal_SO_3, 2)) - 8.0*E_6*pow(Fcal_2, 2) +
      E_SO_3*(Fcal_2*(-10.0*Fcal_5 - 10.0*Fcal_SO_5) + Fcal_3*(-10.0*Fcal_4 - 10.0*Fcal_SQ_4) - 10.0*Fcal_4*Fcal_SO_3 -
      10.0*Fcal_SO_3*Fcal_SQ_4) + E_SO_5*Fcal_2*(-14.0*Fcal_3 - 14.0*Fcal_SO_3) + E_SQ_4*(Fcal_2*(-12.0*Fcal_4 -
      12.0*Fcal_SQ_4) + Fcal_3*(-6.0*Fcal_3 - 12.0*Fcal_SO_3) - 6.0*pow(Fcal_SO_3, 2)) +
      (E_0*(Fcal_2*(Fcal_2*(6.0*Fcal_6 + 6.0*Fcal_SO_6 + 6.0*Fcal_lnv_6*logv) + Fcal_3*(12.0*Fcal_5 + 12.0*Fcal_SO_5) +
      Fcal_4*(6.0*Fcal_4 + 12.0*Fcal_SQ_4) + 12.0*Fcal_5*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SO_5 + 6.0*pow(Fcal_SQ_4, 2)) +
      Fcal_3*(Fcal_3*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + 12.0*Fcal_4*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SQ_4) +
      6.0*Fcal_4*pow(Fcal_SO_3, 2) + 6.0*pow(Fcal_SO_3, 2)*Fcal_SQ_4) + E_2*Fcal_2*(Fcal_2*(12.0*Fcal_4 +
      12.0*Fcal_SQ_4) + Fcal_3*(12.0*Fcal_3 + 24.0*Fcal_SO_3) + 12.0*pow(Fcal_SO_3, 2)) + 6.0*E_4*pow(Fcal_2, 3) +
      E_SO_3*pow(Fcal_2, 2)*(15.0*Fcal_3 + 15.0*Fcal_SO_3) + 6.0*E_SQ_4*pow(Fcal_2, 3) + (E_0*pow(Fcal_2,
      2)*(Fcal_2*(-8.0*Fcal_4 - 8.0*Fcal_SQ_4) + Fcal_3*(-12.0*Fcal_3 - 24.0*Fcal_SO_3) - 12.0*pow(Fcal_SO_3, 2)) +
      2.0*E_0*pow(Fcal_2, 5)/Fcal_0 - 4.0*E_2*pow(Fcal_2, 4))/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0) +
      m*(E_0*(2.0*Fcal_9 + 2.0*Fcal_lnv_9*logv) + E_2*(4.0*Fcal_7 + 4.0*Fcal_SO_7) + E_4*(6.0*Fcal_5 + 6.0*Fcal_SO_5) +
      E_6*(8.0*Fcal_3 + 8.0*Fcal_SO_3) + E_SO_3*(5.0*Fcal_6 + 5.0*Fcal_SO_6 + 5.0*Fcal_lnv_6*logv) + E_SO_5*(7.0*Fcal_4
      + 7.0*Fcal_SQ_4) + 9.0*E_SO_7*Fcal_2 + E_SQ_4*(6.0*Fcal_5 + 6.0*Fcal_SO_5) + (E_0*(Fcal_2*(-4.0*Fcal_7 -
      4.0*Fcal_SO_7) + Fcal_3*(-4.0*Fcal_6 - 4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) + Fcal_4*(-4.0*Fcal_5 - 4.0*Fcal_SO_5)
      - 4.0*Fcal_5*Fcal_SQ_4 - 4.0*Fcal_6*Fcal_SO_3 + Fcal_SO_3*(-4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) -
      4.0*Fcal_SO_5*Fcal_SQ_4) + E_2*(Fcal_2*(-8.0*Fcal_5 - 8.0*Fcal_SO_5) + Fcal_3*(-8.0*Fcal_4 - 8.0*Fcal_SQ_4) -
      8.0*Fcal_4*Fcal_SO_3 - 8.0*Fcal_SO_3*Fcal_SQ_4) + E_4*Fcal_2*(-12.0*Fcal_3 - 12.0*Fcal_SO_3) +
      E_SO_3*(Fcal_2*(-10.0*Fcal_4 - 10.0*Fcal_SQ_4) + Fcal_3*(-5.0*Fcal_3 - 10.0*Fcal_SO_3) - 5.0*pow(Fcal_SO_3, 2)) -
      7.0*E_SO_5*pow(Fcal_2, 2) + E_SQ_4*Fcal_2*(-12.0*Fcal_3 - 12.0*Fcal_SO_3) + (E_0*(Fcal_2*(Fcal_2*(6.0*Fcal_5 +
      6.0*Fcal_SO_5) + Fcal_3*(12.0*Fcal_4 + 12.0*Fcal_SQ_4) + 12.0*Fcal_4*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SQ_4) +
      Fcal_3*(Fcal_3*(2.0*Fcal_3 + 6.0*Fcal_SO_3) + 6.0*pow(Fcal_SO_3, 2)) + 2.0*pow(Fcal_SO_3, 3)) + E_0*pow(Fcal_2,
      3)*(-8.0*Fcal_3 - 8.0*Fcal_SO_3)/Fcal_0 + E_2*pow(Fcal_2, 2)*(12.0*Fcal_3 + 12.0*Fcal_SO_3) +
      5.0*E_SO_3*pow(Fcal_2, 3))/Fcal_0)/Fcal_0)/pow(Fcal_0, 2)) + m*(-10.0*E_8 + E_lnv_8*(-10.0*logv - 1.0) +
      (E_0*(2.0*Fcal_8 + 2.0*Fcal_SO_8 + 2.0*Fcal_lnv_8*logv) + E_2*(4.0*Fcal_6 + 4.0*Fcal_SO_6 + 4.0*Fcal_lnv_6*logv) +
      E_4*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + 8.0*E_6*Fcal_2 + E_SO_3*(5.0*Fcal_5 + 5.0*Fcal_SO_5) + E_SO_5*(7.0*Fcal_3 +
      7.0*Fcal_SO_3) + E_SQ_4*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + (E_0*(Fcal_2*(-4.0*Fcal_6 - 4.0*Fcal_SO_6 -
      4.0*Fcal_lnv_6*logv) + Fcal_3*(-4.0*Fcal_5 - 4.0*Fcal_SO_5) + Fcal_4*(-2.0*Fcal_4 - 4.0*Fcal_SQ_4) -
      4.0*Fcal_5*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SO_5 - 2.0*pow(Fcal_SQ_4, 2)) + E_2*(Fcal_2*(-8.0*Fcal_4 -
      8.0*Fcal_SQ_4) + Fcal_3*(-4.0*Fcal_3 - 8.0*Fcal_SO_3) - 4.0*pow(Fcal_SO_3, 2)) - 6.0*E_4*pow(Fcal_2, 2) +
      E_SO_3*Fcal_2*(-10.0*Fcal_3 - 10.0*Fcal_SO_3) - 6.0*E_SQ_4*pow(Fcal_2, 2) + (E_0*Fcal_2*(Fcal_2*(6.0*Fcal_4 +
      6.0*Fcal_SQ_4) + Fcal_3*(6.0*Fcal_3 + 12.0*Fcal_SO_3) + 6.0*pow(Fcal_SO_3, 2)) - 2.0*E_0*pow(Fcal_2, 4)/Fcal_0 +
      4.0*E_2*pow(Fcal_2, 3))/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0) + m*(-9.0*E_SO_7 + (E_0*(2.0*Fcal_7 + 2.0*Fcal_SO_7) +
      E_2*(4.0*Fcal_5 + 4.0*Fcal_SO_5) + E_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) + E_SO_3*(5.0*Fcal_4 + 5.0*Fcal_SQ_4) +
      7.0*E_SO_5*Fcal_2 + E_SQ_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) + (E_0*(Fcal_2*(-4.0*Fcal_5 - 4.0*Fcal_SO_5) +
      Fcal_3*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) - 4.0*Fcal_4*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SQ_4) + E_0*pow(Fcal_2,
      2)*(6.0*Fcal_3 + 6.0*Fcal_SO_3)/Fcal_0 + E_2*Fcal_2*(-8.0*Fcal_3 - 8.0*Fcal_SO_3) - 5.0*E_SO_3*pow(Fcal_2,
      2))/Fcal_0)/Fcal_0)/Fcal_0) + m*(-8.0*E_6 + (E_0*(2.0*Fcal_6 + 2.0*Fcal_SO_6 + 2.0*Fcal_lnv_6*logv) +
      E_2*(4.0*Fcal_4 + 4.0*Fcal_SQ_4) + 6.0*E_4*Fcal_2 + E_SO_3*(5.0*Fcal_3 + 5.0*Fcal_SO_3) + 6.0*E_SQ_4*Fcal_2 +
      (E_0*(Fcal_2*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) + Fcal_3*(-2.0*Fcal_3 - 4.0*Fcal_SO_3) - 2.0*pow(Fcal_SO_3, 2)) +
      2.0*E_0*pow(Fcal_2, 3)/Fcal_0 - 4.0*E_2*pow(Fcal_2, 2))/Fcal_0)/Fcal_0)/Fcal_0) + m*(-7.0*E_SO_5 +
      (E_0*(2.0*Fcal_5 + 2.0*Fcal_SO_5) + E_0*Fcal_2*(-4.0*Fcal_3 - 4.0*Fcal_SO_3)/Fcal_0 + E_2*(4.0*Fcal_3 +
      4.0*Fcal_SO_3) + 5.0*E_SO_3*Fcal_2)/Fcal_0)/Fcal_0) + m*(-6.0*E_4 - 6.0*E_SQ_4 + (E_0*(2.0*Fcal_4 + 2.0*Fcal_SQ_4)
      - 2.0*E_0*pow(Fcal_2, 2)/Fcal_0 + 4.0*E_2*Fcal_2)/Fcal_0)/Fcal_0) + m*(E_0*(2.0*Fcal_3 + 2.0*Fcal_SO_3)/Fcal_0 -
      5.0*E_SO_3)/Fcal_0) + m*(2.0*E_0*Fcal_2/Fcal_0 - 4.0*E_2)/Fcal_0))/Fcal_coeff;
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_5p5PN_Q : public TaylorTn_Q


class TaylorTn_6p0PN_Q : public TaylorTn_Q {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1, m2;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m, delta, nu;
  Quaternion R, nHat, ellHat;
  double nHat_x, nHat_y, nHat_z, ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chiVec1, chiVec2;
  double chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z;
  const double chi1chi1;
  double chi1chi2;
  const double chi2chi2;
  double chi1_l, chi1_n, chi2_l, chi2_n, S_l, S_n, Sigma_l, Sigma_n, chi_s_l, chi_a_l, logv, Fcal_coeff;
  const double Fcal_0, Fcal_2, Fcal_3, Fcal_4, Fcal_5, Fcal_6, Fcal_lnv_6, Fcal_7, Fcal_8, Fcal_lnv_8, Fcal_9,
               Fcal_lnv_9, Fcal_10, Fcal_lnv_10, Fcal_11, Fcal_lnv_11, Fcal_12, Fcal_lnv_12, Fcal_lnv2_12;
  double Fcal_SQ_4, Fcal_SO_3, Fcal_SO_5, Fcal_SO_6, Fcal_SO_7, Fcal_SO_8;
  const double E_0, E_2, E_4, E_6, E_8, E_lnv_8, E_10, E_lnv_10, E_11, E_12, E_lnv_12;
  double E_SQ_4, E_SO_3, E_SO_5, E_SO_7;
  double Phi;

public:
  TaylorTn_6p0PN_Q(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           m2_i, const double v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const
           double rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i,
           const double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), m2(m2_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m(m1 + m2), delta((m1 - m2)/m), nu(m1*m2/pow(m, 2)), R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)),
    nHat(R*xHat*conjugate(R)), ellHat(R*zHat*conjugate(R)), nHat_x(nHat[1]), nHat_y(nHat[2]), nHat_z(nHat[3]),
    ellHat_x(ellHat[1]), ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chiVec1(chi1Mag*R_S1*zHat*conjugate(R_S1)),
    chiVec2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi1_x(chiVec1[1]), chi1_y(chiVec1[2]), chi1_z(chiVec1[3]),
    chi2_x(chiVec2[1]), chi2_y(chiVec2[2]), chi2_z(chiVec2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z,
    2)), chi1chi2(chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y +
    chi2_z*nHat_z), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)),
    Sigma_l(m*(-chi1_l*m1 + chi2_l*m2)), Sigma_n(m*(-chi1_n*m1 + chi2_n*m2)), chi_s_l(0.5*chi1_l + 0.5*chi2_l),
    chi_a_l(0.5*chi1_l - 0.5*chi2_l), logv(log(v)), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    Fcal_2(-2.91666666666667*nu - 3.71130952380952), Fcal_3(12.5663706143592), Fcal_4(3.61111111111111*pow(nu, 2) +
    18.3948412698413*nu - 4.92846119929453), Fcal_5(-76.3145215434521*nu - 38.2928354546934),
    Fcal_6(-2.39197530864198*pow(nu, 3) - 31.2179232804233*pow(nu, 2) - 8.87205344238227*nu + 115.731716675611),
    Fcal_lnv_6(-16.3047619047619), Fcal_7(200.905057974359*pow(nu, 2) + 390.417427312002*nu - 101.509595959742),
    Fcal_8(-117.504390722677), Fcal_lnv_8(52.7430839002268), Fcal_9(719.128342233430), Fcal_lnv_9(-204.891680874123),
    Fcal_10(-1216.90699131704), Fcal_lnv_10(116.639876594109), Fcal_11(958.934970119567), Fcal_lnv_11(473.624478174231),
    Fcal_12(2034.78998213038), Fcal_lnv_12(-1900.72932518810), Fcal_lnv2_12(132.922630385488),
    Fcal_SQ_4(chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
    2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
    chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
    2)*(0.0416666666666667*nu + 2.98958333333333)), Fcal_SO_3((-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2)),
    Fcal_SO_5((S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2)),
    Fcal_SO_6((-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2)),
    Fcal_SO_7((S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
    Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762))/pow(m, 2)),
    Fcal_SO_8((3.14159265358979*S_l*(192.763888888889*nu - 36.3020833333333) +
    3.14159265358979*Sigma_l*delta*(64.7733134920635*nu - 10.6592261904762))/pow(m, 2)), E_0(1.00000000000000),
    E_2(-0.0833333333333333*nu - 0.75), E_4(-0.0416666666666667*pow(nu, 2) + 2.375*nu - 3.375),
    E_6(-0.00675154320987654*pow(nu, 3) - 1.61458333333333*pow(nu, 2) + 38.7246294907293*nu - 10.546875),
    E_8(0.0024755658436214*pow(nu, 4) + 0.174189814814815*pow(nu, 3) - 90.1327990262052*pow(nu, 2) + 153.88379682994*nu
    - 31.0078125), E_lnv_8(59.7333333333333*nu), E_10(0.001953125*pow(nu, 5) + 0.107421875*pow(nu, 4) +
    83.4293954895551*pow(nu, 3) - 71.3641942901125*pow(nu, 2) - 62.1764445589718*nu - 89.701171875),
    E_lnv_10(-262.4*pow(nu, 2) - 285.028571428571*nu), E_11(273.188907832164*nu), E_12(0.000404407912284713*pow(nu, 6) +
    0.0207328639403292*pow(nu, 5) - 33.3947463172536*pow(nu, 4) - 701.460497654495*pow(nu, 3) + 2936.51869895088*pow(nu,
    2) + 3015.56151979737*nu - 258.4248046875), E_lnv_12(362.42962962963*pow(nu, 3) + 757.219047619048*pow(nu, 2) -
    462.61822457378*nu), E_SQ_4(-1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l)
    + nu*(chi1chi2 + 6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0)),
    E_SO_3((4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2)), E_SO_5((S_l*(-6.77777777777778*nu + 11.0) +
    Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2)), E_SO_7((S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75)
    + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu + 6.75))/pow(m, 2)), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_z = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    ellHat = R*zHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chiVec1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chiVec2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi1_x = chiVec1[1];
    chi1_y = chiVec1[2];
    chi1_z = chiVec1[3];
    chi2_x = chiVec2[1];
    chi2_y = chiVec2[2];
    chi2_z = chiVec2[3];
    chi1chi2 = chi1_x*chi2_x + chi1_y*chi2_y + chi1_z*chi2_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = m*(-chi1_l*m1 + chi2_l*m2);
    Sigma_n = m*(-chi1_n*m1 + chi2_n*m2);
    chi_s_l = 0.5*chi1_l + 0.5*chi2_l;
    chi_a_l = 0.5*chi1_l - 0.5*chi2_l;
    logv = log(v);
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
    Fcal_SQ_4 = chi1chi1*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
      2.14583333333333*chi1chi2*nu + chi2chi2*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) +
      chi_a_l*(chi_a_l*(-12.0*nu + 2.98958333333333) + 5.97916666666667*chi_s_l*delta) + pow(chi_s_l,
      2)*(0.0416666666666667*nu + 2.98958333333333);
    Fcal_SO_3 = (-4.0*S_l - 1.25*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_5 = (S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125))/pow(m, 2);
    Fcal_SO_6 = (-50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta)/pow(m, 2);
    Fcal_SO_7 = (S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
      Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762))/pow(m, 2);
    Fcal_SO_8 = (3.14159265358979*S_l*(192.763888888889*nu - 36.3020833333333) +
      3.14159265358979*Sigma_l*delta*(64.7733134920635*nu - 10.6592261904762))/pow(m, 2);
    E_SQ_4 = -1.5*pow(chi_a_l, 2) - 1.5*pow(chi_s_l, 2) - delta*(0.5*chi2chi2 + 3.0*chi_a_l*chi_s_l) + nu*(chi1chi2 +
      6.0*pow(chi_a_l, 2)) + 0.25*(chi1chi1 + chi2chi2)*(delta - 2.0*nu + 1.0);
    E_SO_3 = (4.66666666666667*S_l + 2.0*Sigma_l*delta)/pow(m, 2);
    E_SO_5 = (S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0))/pow(m, 2);
    E_SO_7 = (S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75) + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu +
      6.75))/pow(m, 2);
  }

  Quaternion OmegaVec_chiVec_1() {
    double Omega1_coeff = pow(v, 5)/m;
    return Omega1_coeff*(-chiVec2*pow(m2, 2)*v/pow(m, 2) + ellHat*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu -
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu +
      3.0*chi2_n*pow(m2, 2)/pow(m, 2)));
  }
  Quaternion OmegaVec_chiVec_2() {
    double Omega2_coeff = pow(v, 5)/m;
    return Omega2_coeff*(-chiVec1*pow(m1, 2)*v/pow(m, 2) + ellHat*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu +
      0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*pow(m1,
      2)/pow(m, 2) + 3.0*chi2_n*nu));
  }
  Quaternion OmegaVec() {
    double gamma_PN_2 = -0.333333333333333*nu + 1.0;
    double gamma_PN_6 = 0.0123456790123457*pow(nu, 3) + 6.36111111111111*pow(nu, 2) - 2.98177812235564*nu + 1.0;
    double a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    double gamma_PN_3 = (1.66666666666667*S_l + Sigma_l*delta)/pow(m, 2);
    double gamma_PN_7 = (S_l*(-6.0*pow(nu, 2) - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_l*delta*pow(nu, 2) +
      Sigma_l*delta*(-10.1666666666667*nu + 3.0))/pow(m, 2);
    double gamma_PN_5 = (S_l*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_l*delta)/pow(m, 2);
    double gamma_PN_0 = 1.00000000000000;
    double a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0);
    double gamma_PN_4 = -5.41666666666667*nu + 1.0;
    double a_ell_4 = S_n*(5.77777777777778*pow(nu, 2) + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*pow(nu, 2) +
      9.125*nu + 1.5);
    return ellHat*pow(v, 3)/m + nHat*pow(v, 6)*(a_ell_0 + pow(v, 2)*(a_ell_2 + a_ell_4*pow(v, 2)))*(gamma_PN_0 + pow(v,
      2)*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/pow(m, 3);
  }

  int TaylorT1(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_coeff*(Fcal_0 + pow(v, 2)*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 +
      v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 +
      Fcal_lnv_8*logv + v*(Fcal_9 + Fcal_lnv_9*logv + v*(Fcal_10 + Fcal_lnv_10*logv + v*(Fcal_11 + Fcal_lnv_11*logv +
      v*(Fcal_12 + Fcal_lnv2_12*pow(logv, 2) + Fcal_lnv_12*logv))))))))))));
    const double dEdv = -0.5*m*nu*v*(2.0*E_0 + pow(v, 2)*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 +
      v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0) + pow(v, 2)*(12.0*E_10 +
      E_lnv_10*(12.0*logv + 1.0) + v*(13.0*E_11 + v*(14.0*E_12 + E_lnv_12*(14.0*logv + 1.0))))))))))));
    const double Absorption = 0;
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = -2.0*Fcal_coeff*(pow(v, 2)*(v*(v*(v*(v*(v*(v*(v*(v*(v*(v*((-0.5*Fcal_12 -
      0.5*Fcal_lnv2_12*pow(logv, 2) - 0.5*Fcal_lnv_12*logv)/m + ((3.0*E_10*Fcal_2 + 3.5*E_12*Fcal_0 + E_2*(1.0*Fcal_10 +
      1.0*Fcal_lnv_10*logv) + E_4*(1.5*Fcal_8 + 1.5*Fcal_SO_8 + 1.5*Fcal_lnv_8*logv) + E_6*(2.0*Fcal_6 + 2.0*Fcal_SO_6 +
      2.0*Fcal_lnv_6*logv) + E_8*(2.5*Fcal_4 + 2.5*Fcal_SQ_4) + E_SO_3*(1.25*Fcal_9 + 1.25*Fcal_lnv_9*logv) +
      E_SO_5*(1.75*Fcal_7 + 1.75*Fcal_SO_7) + E_SO_7*(2.25*Fcal_5 + 2.25*Fcal_SO_5) + E_SQ_4*(1.5*Fcal_8 + 1.5*Fcal_SO_8
      + 1.5*Fcal_lnv_8*logv) + E_lnv_10*Fcal_2*(3.0*logv + 0.25) + E_lnv_12*Fcal_0*(3.5*logv + 0.25) +
      E_lnv_8*(Fcal_4*(2.5*logv + 0.25) + Fcal_SQ_4*(2.5*logv + 0.25)))/m + ((E_2*(-12.0*E_10*Fcal_0 + E_2*(-2.0*Fcal_8
      - 2.0*Fcal_SO_8 - 2.0*Fcal_lnv_8*logv) + E_4*(-6.0*Fcal_6 - 6.0*Fcal_SO_6 - 6.0*Fcal_lnv_6*logv) +
      E_6*(-8.0*Fcal_4 - 8.0*Fcal_SQ_4) - 10.0*E_8*Fcal_2 + E_SO_3*(-5.0*Fcal_7 - 5.0*Fcal_SO_7) + E_SO_5*(-7.0*Fcal_5 -
      7.0*Fcal_SO_5) + E_SO_7*(-9.0*Fcal_3 - 9.0*Fcal_SO_3) + E_SQ_4*(-6.0*Fcal_6 - 6.0*Fcal_SO_6 - 6.0*Fcal_lnv_6*logv)
      + E_lnv_10*Fcal_0*(-12.0*logv - 1.0) + E_lnv_8*Fcal_2*(-10.0*logv - 1.0)) + E_4*(E_4*(-4.5*Fcal_4 - 4.5*Fcal_SQ_4)
      - 12.0*E_6*Fcal_2 - 15.0*E_8*Fcal_0 + E_SO_3*(-7.5*Fcal_5 - 7.5*Fcal_SO_5) + E_SO_5*(-10.5*Fcal_3 -
      10.5*Fcal_SO_3) + E_SQ_4*(-9.0*Fcal_4 - 9.0*Fcal_SQ_4) + E_lnv_8*Fcal_0*(-15.0*logv - 1.5)) + E_6*(-8.0*E_6*Fcal_0
      + E_SO_3*(-10.0*Fcal_3 - 10.0*Fcal_SO_3) - 12.0*E_SQ_4*Fcal_2) - 15.0*E_8*E_SQ_4*Fcal_0 +
      E_SO_3*(E_SO_3*(-3.125*Fcal_6 - 3.125*Fcal_SO_6 - 3.125*Fcal_lnv_6*logv) + E_SO_5*(-8.75*Fcal_4 - 8.75*Fcal_SQ_4)
      - 11.25*E_SO_7*Fcal_2 + E_SQ_4*(-7.5*Fcal_5 - 7.5*Fcal_SO_5)) + E_SO_5*(-6.125*E_SO_5*Fcal_2 - 15.75*E_SO_7*Fcal_0
      + E_SQ_4*(-10.5*Fcal_3 - 10.5*Fcal_SO_3)) + E_SQ_4*(E_SQ_4*(-4.5*Fcal_4 - 4.5*Fcal_SQ_4) +
      E_lnv_8*Fcal_0*(-15.0*logv - 1.5)))/m + ((E_2*(E_2*(E_2*(4.0*Fcal_6 + 4.0*Fcal_SO_6 + 4.0*Fcal_lnv_6*logv) +
      E_4*(18.0*Fcal_4 + 18.0*Fcal_SQ_4) + 24.0*E_6*Fcal_2 + 30.0*E_8*Fcal_0 + E_SO_3*(15.0*Fcal_5 + 15.0*Fcal_SO_5) +
      E_SO_5*(21.0*Fcal_3 + 21.0*Fcal_SO_3) + E_SQ_4*(18.0*Fcal_4 + 18.0*Fcal_SQ_4) + E_lnv_8*Fcal_0*(30.0*logv + 3.0))
      + E_4*(27.0*E_4*Fcal_2 + 72.0*E_6*Fcal_0 + E_SO_3*(45.0*Fcal_3 + 45.0*Fcal_SO_3) + 54.0*E_SQ_4*Fcal_2) +
      72.0*E_6*E_SQ_4*Fcal_0 + E_SO_3*(E_SO_3*(18.75*Fcal_4 + 18.75*Fcal_SQ_4) + 52.5*E_SO_5*Fcal_2 + 67.5*E_SO_7*Fcal_0
      + E_SQ_4*(45.0*Fcal_3 + 45.0*Fcal_SO_3)) + 36.75*pow(E_SO_5, 2)*Fcal_0 + 27.0*pow(E_SQ_4, 2)*Fcal_2) +
      E_4*(E_4*(13.5*E_4*Fcal_0 + 40.5*E_SQ_4*Fcal_0) + E_SO_3*(28.125*E_SO_3*Fcal_2 + 78.75*E_SO_5*Fcal_0) +
      40.5*pow(E_SQ_4, 2)*Fcal_0) + 37.5*E_6*pow(E_SO_3, 2)*Fcal_0 + E_SO_3*(E_SO_3*(E_SO_3*(7.8125*Fcal_3 +
      7.8125*Fcal_SO_3) + 28.125*E_SQ_4*Fcal_2) + 78.75*E_SO_5*E_SQ_4*Fcal_0) + 13.5*pow(E_SQ_4, 3)*Fcal_0)/m +
      ((E_2*(E_2*(E_2*(E_2*(-8.0*Fcal_4 - 8.0*Fcal_SQ_4) - 48.0*E_4*Fcal_2 - 64.0*E_6*Fcal_0 + E_SO_3*(-40.0*Fcal_3 -
      40.0*Fcal_SO_3) - 48.0*E_SQ_4*Fcal_2) + E_4*(-108.0*E_4*Fcal_0 - 216.0*E_SQ_4*Fcal_0) +
      E_SO_3*(-75.0*E_SO_3*Fcal_2 - 210.0*E_SO_5*Fcal_0) - 108.0*pow(E_SQ_4, 2)*Fcal_0) - 225.0*E_4*pow(E_SO_3,
      2)*Fcal_0 - 225.0*pow(E_SO_3, 2)*E_SQ_4*Fcal_0) - 19.53125*pow(E_SO_3, 4)*Fcal_0)/m + (pow(E_2,
      3)*(E_2*(16.0*E_2*Fcal_2 + 120.0*E_4*Fcal_0 + 120.0*E_SQ_4*Fcal_0) + 250.0*pow(E_SO_3, 2)*Fcal_0)/m -
      32.0*pow(E_2, 6)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0)/E_0)/E_0 + ((-0.5*Fcal_11 - 0.5*Fcal_lnv_11*logv)/m +
      ((3.25*E_11*Fcal_0 + E_2*(1.0*Fcal_9 + 1.0*Fcal_lnv_9*logv) + E_4*(1.5*Fcal_7 + 1.5*Fcal_SO_7) + E_6*(2.0*Fcal_5 +
      2.0*Fcal_SO_5) + E_8*(2.5*Fcal_3 + 2.5*Fcal_SO_3) + E_SO_3*(1.25*Fcal_8 + 1.25*Fcal_SO_8 + 1.25*Fcal_lnv_8*logv) +
      E_SO_5*(1.75*Fcal_6 + 1.75*Fcal_SO_6 + 1.75*Fcal_lnv_6*logv) + E_SO_7*(2.25*Fcal_4 + 2.25*Fcal_SQ_4) +
      E_SQ_4*(1.5*Fcal_7 + 1.5*Fcal_SO_7) + E_lnv_8*(Fcal_3*(2.5*logv + 0.25) + Fcal_SO_3*(2.5*logv + 0.25)))/m +
      ((E_2*(E_2*(-2.0*Fcal_7 - 2.0*Fcal_SO_7) + E_4*(-6.0*Fcal_5 - 6.0*Fcal_SO_5) + E_6*(-8.0*Fcal_3 - 8.0*Fcal_SO_3) +
      E_SO_3*(-5.0*Fcal_6 - 5.0*Fcal_SO_6 - 5.0*Fcal_lnv_6*logv) + E_SO_5*(-7.0*Fcal_4 - 7.0*Fcal_SQ_4) -
      9.0*E_SO_7*Fcal_2 + E_SQ_4*(-6.0*Fcal_5 - 6.0*Fcal_SO_5)) + E_4*(E_4*(-4.5*Fcal_3 - 4.5*Fcal_SO_3) +
      E_SO_3*(-7.5*Fcal_4 - 7.5*Fcal_SQ_4) - 10.5*E_SO_5*Fcal_2 - 13.5*E_SO_7*Fcal_0 + E_SQ_4*(-9.0*Fcal_3 -
      9.0*Fcal_SO_3)) + E_6*(-10.0*E_SO_3*Fcal_2 - 14.0*E_SO_5*Fcal_0) - 12.5*E_8*E_SO_3*Fcal_0 +
      E_SO_3*(E_SO_3*(-3.125*Fcal_5 - 3.125*Fcal_SO_5) + E_SO_5*(-8.75*Fcal_3 - 8.75*Fcal_SO_3) + E_SQ_4*(-7.5*Fcal_4 -
      7.5*Fcal_SQ_4) + E_lnv_8*Fcal_0*(-12.5*logv - 1.25)) - 10.5*E_SO_5*E_SQ_4*Fcal_2 - 13.5*E_SO_7*E_SQ_4*Fcal_0 +
      pow(E_SQ_4, 2)*(-4.5*Fcal_3 - 4.5*Fcal_SO_3))/m + ((E_2*(E_2*(E_2*(4.0*Fcal_5 + 4.0*Fcal_SO_5) + E_4*(18.0*Fcal_3
      + 18.0*Fcal_SO_3) + E_SO_3*(15.0*Fcal_4 + 15.0*Fcal_SQ_4) + 21.0*E_SO_5*Fcal_2 + 27.0*E_SO_7*Fcal_0 +
      E_SQ_4*(18.0*Fcal_3 + 18.0*Fcal_SO_3)) + E_4*(45.0*E_SO_3*Fcal_2 + 63.0*E_SO_5*Fcal_0) + 60.0*E_6*E_SO_3*Fcal_0 +
      E_SO_3*(E_SO_3*(18.75*Fcal_3 + 18.75*Fcal_SO_3) + 45.0*E_SQ_4*Fcal_2) + 63.0*E_SO_5*E_SQ_4*Fcal_0) +
      E_4*(33.75*E_4*E_SO_3*Fcal_0 + 67.5*E_SO_3*E_SQ_4*Fcal_0) + E_SO_3*(E_SO_3*(7.8125*E_SO_3*Fcal_2 +
      32.8125*E_SO_5*Fcal_0) + 33.75*pow(E_SQ_4, 2)*Fcal_0))/m + (E_2*(E_2*(E_2*(E_2*(-8.0*Fcal_3 - 8.0*Fcal_SO_3) -
      40.0*E_SO_3*Fcal_2 - 56.0*E_SO_5*Fcal_0) - 180.0*E_4*E_SO_3*Fcal_0 - 180.0*E_SO_3*E_SQ_4*Fcal_0) -
      62.5*pow(E_SO_3, 3)*Fcal_0)/m + 100.0*pow(E_2, 4)*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0)/E_0) + ((-0.5*Fcal_10
      - 0.5*Fcal_lnv_10*logv)/m + ((3.0*E_10*Fcal_0 + E_2*(1.0*Fcal_8 + 1.0*Fcal_SO_8 + 1.0*Fcal_lnv_8*logv) +
      E_4*(1.5*Fcal_6 + 1.5*Fcal_SO_6 + 1.5*Fcal_lnv_6*logv) + E_6*(2.0*Fcal_4 + 2.0*Fcal_SQ_4) + 2.5*E_8*Fcal_2 +
      E_SO_3*(1.25*Fcal_7 + 1.25*Fcal_SO_7) + E_SO_5*(1.75*Fcal_5 + 1.75*Fcal_SO_5) + E_SO_7*(2.25*Fcal_3 +
      2.25*Fcal_SO_3) + E_SQ_4*(1.5*Fcal_6 + 1.5*Fcal_SO_6 + 1.5*Fcal_lnv_6*logv) + E_lnv_10*Fcal_0*(3.0*logv + 0.25) +
      E_lnv_8*Fcal_2*(2.5*logv + 0.25))/m + ((E_2*(E_2*(-2.0*Fcal_6 - 2.0*Fcal_SO_6 - 2.0*Fcal_lnv_6*logv) +
      E_4*(-6.0*Fcal_4 - 6.0*Fcal_SQ_4) - 8.0*E_6*Fcal_2 - 10.0*E_8*Fcal_0 + E_SO_3*(-5.0*Fcal_5 - 5.0*Fcal_SO_5) +
      E_SO_5*(-7.0*Fcal_3 - 7.0*Fcal_SO_3) + E_SQ_4*(-6.0*Fcal_4 - 6.0*Fcal_SQ_4) + E_lnv_8*Fcal_0*(-10.0*logv - 1.0)) +
      E_4*(-4.5*E_4*Fcal_2 - 12.0*E_6*Fcal_0 + E_SO_3*(-7.5*Fcal_3 - 7.5*Fcal_SO_3) - 9.0*E_SQ_4*Fcal_2) -
      12.0*E_6*E_SQ_4*Fcal_0 + E_SO_3*(E_SO_3*(-3.125*Fcal_4 - 3.125*Fcal_SQ_4) - 8.75*E_SO_5*Fcal_2 -
      11.25*E_SO_7*Fcal_0 + E_SQ_4*(-7.5*Fcal_3 - 7.5*Fcal_SO_3)) - 6.125*pow(E_SO_5, 2)*Fcal_0 - 4.5*pow(E_SQ_4,
      2)*Fcal_2)/m + ((E_2*(E_2*(E_2*(4.0*Fcal_4 + 4.0*Fcal_SQ_4) + 18.0*E_4*Fcal_2 + 24.0*E_6*Fcal_0 +
      E_SO_3*(15.0*Fcal_3 + 15.0*Fcal_SO_3) + 18.0*E_SQ_4*Fcal_2) + E_4*(27.0*E_4*Fcal_0 + 54.0*E_SQ_4*Fcal_0) +
      E_SO_3*(18.75*E_SO_3*Fcal_2 + 52.5*E_SO_5*Fcal_0) + 27.0*pow(E_SQ_4, 2)*Fcal_0) + 28.125*E_4*pow(E_SO_3, 2)*Fcal_0
      + 28.125*pow(E_SO_3, 2)*E_SQ_4*Fcal_0)/m + (pow(E_2, 2)*(E_2*(-8.0*E_2*Fcal_2 - 48.0*E_4*Fcal_0 -
      48.0*E_SQ_4*Fcal_0) - 75.0*pow(E_SO_3, 2)*Fcal_0)/m + 16.0*pow(E_2, 5)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0)/E_0) +
      ((-0.5*Fcal_9 - 0.5*Fcal_lnv_9*logv)/m + ((E_2*(1.0*Fcal_7 + 1.0*Fcal_SO_7) + E_4*(1.5*Fcal_5 + 1.5*Fcal_SO_5) +
      E_6*(2.0*Fcal_3 + 2.0*Fcal_SO_3) + E_SO_3*(1.25*Fcal_6 + 1.25*Fcal_SO_6 + 1.25*Fcal_lnv_6*logv) +
      E_SO_5*(1.75*Fcal_4 + 1.75*Fcal_SQ_4) + 2.25*E_SO_7*Fcal_2 + E_SQ_4*(1.5*Fcal_5 + 1.5*Fcal_SO_5))/m +
      ((E_2*(E_2*(-2.0*Fcal_5 - 2.0*Fcal_SO_5) + E_4*(-6.0*Fcal_3 - 6.0*Fcal_SO_3) + E_SO_3*(-5.0*Fcal_4 -
      5.0*Fcal_SQ_4) - 7.0*E_SO_5*Fcal_2 - 9.0*E_SO_7*Fcal_0 + E_SQ_4*(-6.0*Fcal_3 - 6.0*Fcal_SO_3)) +
      E_4*(-7.5*E_SO_3*Fcal_2 - 10.5*E_SO_5*Fcal_0) - 10.0*E_6*E_SO_3*Fcal_0 + E_SO_3*(E_SO_3*(-3.125*Fcal_3 -
      3.125*Fcal_SO_3) - 7.5*E_SQ_4*Fcal_2) - 10.5*E_SO_5*E_SQ_4*Fcal_0)/m + ((E_2*(E_2*(E_2*(4.0*Fcal_3 +
      4.0*Fcal_SO_3) + 15.0*E_SO_3*Fcal_2 + 21.0*E_SO_5*Fcal_0) + 45.0*E_4*E_SO_3*Fcal_0 + 45.0*E_SO_3*E_SQ_4*Fcal_0) +
      7.8125*pow(E_SO_3, 3)*Fcal_0)/m - 40.0*pow(E_2, 3)*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0) + ((-0.5*Fcal_8 -
      0.5*Fcal_SO_8 - 0.5*Fcal_lnv_8*logv)/m + ((E_2*(1.0*Fcal_6 + 1.0*Fcal_SO_6 + 1.0*Fcal_lnv_6*logv) +
      E_4*(1.5*Fcal_4 + 1.5*Fcal_SQ_4) + 2.0*E_6*Fcal_2 + 2.5*E_8*Fcal_0 + E_SO_3*(1.25*Fcal_5 + 1.25*Fcal_SO_5) +
      E_SO_5*(1.75*Fcal_3 + 1.75*Fcal_SO_3) + E_SQ_4*(1.5*Fcal_4 + 1.5*Fcal_SQ_4) + E_lnv_8*Fcal_0*(2.5*logv + 0.25))/m
      + ((E_2*(E_2*(-2.0*Fcal_4 - 2.0*Fcal_SQ_4) - 6.0*E_4*Fcal_2 - 8.0*E_6*Fcal_0 + E_SO_3*(-5.0*Fcal_3 -
      5.0*Fcal_SO_3) - 6.0*E_SQ_4*Fcal_2) + E_4*(-4.5*E_4*Fcal_0 - 9.0*E_SQ_4*Fcal_0) + E_SO_3*(-3.125*E_SO_3*Fcal_2 -
      8.75*E_SO_5*Fcal_0) - 4.5*pow(E_SQ_4, 2)*Fcal_0)/m + (E_2*(E_2*(4.0*E_2*Fcal_2 + 18.0*E_4*Fcal_0 +
      18.0*E_SQ_4*Fcal_0) + 18.75*pow(E_SO_3, 2)*Fcal_0)/m - 8.0*pow(E_2, 4)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0)/E_0) +
      ((-0.5*Fcal_7 - 0.5*Fcal_SO_7)/m + ((E_2*(1.0*Fcal_5 + 1.0*Fcal_SO_5) + E_4*(1.5*Fcal_3 + 1.5*Fcal_SO_3) +
      E_SO_3*(1.25*Fcal_4 + 1.25*Fcal_SQ_4) + 1.75*E_SO_5*Fcal_2 + 2.25*E_SO_7*Fcal_0 + E_SQ_4*(1.5*Fcal_3 +
      1.5*Fcal_SO_3))/m + ((E_2*(E_2*(-2.0*Fcal_3 - 2.0*Fcal_SO_3) - 5.0*E_SO_3*Fcal_2 - 7.0*E_SO_5*Fcal_0) -
      7.5*E_4*E_SO_3*Fcal_0 - 7.5*E_SO_3*E_SQ_4*Fcal_0)/m + 15.0*pow(E_2, 2)*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0)/E_0) +
      ((-0.5*Fcal_6 - 0.5*Fcal_SO_6 - 0.5*Fcal_lnv_6*logv)/m + ((E_2*(1.0*Fcal_4 + 1.0*Fcal_SQ_4) + 1.5*E_4*Fcal_2 +
      2.0*E_6*Fcal_0 + E_SO_3*(1.25*Fcal_3 + 1.25*Fcal_SO_3) + 1.5*E_SQ_4*Fcal_2)/m + ((E_2*(-2.0*E_2*Fcal_2 -
      6.0*E_4*Fcal_0 - 6.0*E_SQ_4*Fcal_0) - 3.125*pow(E_SO_3, 2)*Fcal_0)/m + 4.0*pow(E_2,
      3)*Fcal_0/(E_0*m))/E_0)/E_0)/E_0) + ((-0.5*Fcal_5 - 0.5*Fcal_SO_5)/m + ((E_2*(1.0*Fcal_3 + 1.0*Fcal_SO_3) +
      1.25*E_SO_3*Fcal_2 + 1.75*E_SO_5*Fcal_0)/m - 5.0*E_2*E_SO_3*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_4 -
      0.5*Fcal_SQ_4)/m + ((1.0*E_2*Fcal_2 + 1.5*E_4*Fcal_0 + 1.5*E_SQ_4*Fcal_0)/m - 2.0*pow(E_2,
      2)*Fcal_0/(E_0*m))/E_0)/E_0) + ((-0.5*Fcal_3 - 0.5*Fcal_SO_3)/m + 1.25*E_SO_3*Fcal_0/(E_0*m))/E_0) +
      (-0.5*Fcal_2/m + 1.0*E_2*Fcal_0/(E_0*m))/E_0) - 0.5*Fcal_0/(E_0*m))/(nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = -0.5*nu*v*(-2.0*E_0*m/Fcal_0 + pow(v, 2)*(v*(v*(v*(v*(v*(v*(v*(v*(v*(m*v*(-14.0*E_12 +
      E_lnv_12*(-14.0*logv - 1.0) + (E_0*(2.0*Fcal_12 + 2.0*Fcal_lnv2_12*pow(logv, 2) + 2.0*Fcal_lnv_12*logv) +
      12.0*E_10*Fcal_2 + E_2*(4.0*Fcal_10 + 4.0*Fcal_lnv_10*logv) + E_4*(6.0*Fcal_8 + 6.0*Fcal_SO_8 +
      6.0*Fcal_lnv_8*logv) + E_6*(8.0*Fcal_6 + 8.0*Fcal_SO_6 + 8.0*Fcal_lnv_6*logv) + E_8*(10.0*Fcal_4 + 10.0*Fcal_SQ_4)
      + E_SO_3*(5.0*Fcal_9 + 5.0*Fcal_lnv_9*logv) + E_SO_5*(7.0*Fcal_7 + 7.0*Fcal_SO_7) + E_SO_7*(9.0*Fcal_5 +
      9.0*Fcal_SO_5) + E_SQ_4*(6.0*Fcal_8 + 6.0*Fcal_SO_8 + 6.0*Fcal_lnv_8*logv) + E_lnv_10*Fcal_2*(12.0*logv + 1.0) +
      E_lnv_8*(Fcal_4*(10.0*logv + 1.0) + Fcal_SQ_4*(10.0*logv + 1.0)) + (E_0*(Fcal_2*(-4.0*Fcal_10 -
      4.0*Fcal_lnv_10*logv) + Fcal_3*(-4.0*Fcal_9 - 4.0*Fcal_lnv_9*logv) + Fcal_4*(-4.0*Fcal_8 - 4.0*Fcal_SO_8 -
      4.0*Fcal_lnv_8*logv) + Fcal_5*(-4.0*Fcal_7 - 4.0*Fcal_SO_7) + Fcal_6*(-2.0*Fcal_6 - 4.0*Fcal_SO_6 -
      4.0*Fcal_lnv_6*logv) - 4.0*Fcal_7*Fcal_SO_5 - 4.0*Fcal_8*Fcal_SQ_4 - 4.0*Fcal_9*Fcal_SO_3 -
      4.0*Fcal_SO_3*Fcal_lnv_9*logv - 4.0*Fcal_SO_5*Fcal_SO_7 + Fcal_SO_6*(-2.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) -
      4.0*Fcal_SO_8*Fcal_SQ_4 - 4.0*Fcal_SQ_4*Fcal_lnv_8*logv - 2.0*pow(Fcal_lnv_6, 2)*pow(logv, 2)) +
      E_2*(Fcal_2*(-8.0*Fcal_8 - 8.0*Fcal_SO_8 - 8.0*Fcal_lnv_8*logv) + Fcal_3*(-8.0*Fcal_7 - 8.0*Fcal_SO_7) +
      Fcal_4*(-8.0*Fcal_6 - 8.0*Fcal_SO_6 - 8.0*Fcal_lnv_6*logv) + Fcal_5*(-4.0*Fcal_5 - 8.0*Fcal_SO_5) -
      8.0*Fcal_6*Fcal_SQ_4 - 8.0*Fcal_7*Fcal_SO_3 - 8.0*Fcal_SO_3*Fcal_SO_7 - 4.0*pow(Fcal_SO_5, 2) -
      8.0*Fcal_SO_6*Fcal_SQ_4 - 8.0*Fcal_SQ_4*Fcal_lnv_6*logv) + E_4*(Fcal_2*(-12.0*Fcal_6 - 12.0*Fcal_SO_6 -
      12.0*Fcal_lnv_6*logv) + Fcal_3*(-12.0*Fcal_5 - 12.0*Fcal_SO_5) + Fcal_4*(-6.0*Fcal_4 - 12.0*Fcal_SQ_4) -
      12.0*Fcal_5*Fcal_SO_3 - 12.0*Fcal_SO_3*Fcal_SO_5 - 6.0*pow(Fcal_SQ_4, 2)) + E_6*(Fcal_2*(-16.0*Fcal_4 -
      16.0*Fcal_SQ_4) + Fcal_3*(-8.0*Fcal_3 - 16.0*Fcal_SO_3) - 8.0*pow(Fcal_SO_3, 2)) - 10.0*E_8*pow(Fcal_2, 2) +
      E_SO_3*(Fcal_2*(-10.0*Fcal_7 - 10.0*Fcal_SO_7) + Fcal_3*(-10.0*Fcal_6 - 10.0*Fcal_SO_6 - 10.0*Fcal_lnv_6*logv) +
      Fcal_4*(-10.0*Fcal_5 - 10.0*Fcal_SO_5) - 10.0*Fcal_5*Fcal_SQ_4 - 10.0*Fcal_6*Fcal_SO_3 +
      Fcal_SO_3*(-10.0*Fcal_SO_6 - 10.0*Fcal_lnv_6*logv) - 10.0*Fcal_SO_5*Fcal_SQ_4) + E_SO_5*(Fcal_2*(-14.0*Fcal_5 -
      14.0*Fcal_SO_5) + Fcal_3*(-14.0*Fcal_4 - 14.0*Fcal_SQ_4) - 14.0*Fcal_4*Fcal_SO_3 - 14.0*Fcal_SO_3*Fcal_SQ_4) +
      E_SO_7*Fcal_2*(-18.0*Fcal_3 - 18.0*Fcal_SO_3) + E_SQ_4*(Fcal_2*(-12.0*Fcal_6 - 12.0*Fcal_SO_6 -
      12.0*Fcal_lnv_6*logv) + Fcal_3*(-12.0*Fcal_5 - 12.0*Fcal_SO_5) + Fcal_4*(-6.0*Fcal_4 - 12.0*Fcal_SQ_4) -
      12.0*Fcal_5*Fcal_SO_3 - 12.0*Fcal_SO_3*Fcal_SO_5 - 6.0*pow(Fcal_SQ_4, 2)) + E_lnv_8*pow(Fcal_2, 2)*(-10.0*logv -
      1.0) + (E_0*(Fcal_2*(Fcal_2*(6.0*Fcal_8 + 6.0*Fcal_SO_8 + 6.0*Fcal_lnv_8*logv) + Fcal_3*(12.0*Fcal_7 +
      12.0*Fcal_SO_7) + Fcal_4*(12.0*Fcal_6 + 12.0*Fcal_SO_6 + 12.0*Fcal_lnv_6*logv) + Fcal_5*(6.0*Fcal_5 +
      12.0*Fcal_SO_5) + 12.0*Fcal_6*Fcal_SQ_4 + 12.0*Fcal_7*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SO_7 + 6.0*pow(Fcal_SO_5, 2)
      + 12.0*Fcal_SO_6*Fcal_SQ_4 + 12.0*Fcal_SQ_4*Fcal_lnv_6*logv) + Fcal_3*(Fcal_3*(6.0*Fcal_6 + 6.0*Fcal_SO_6 +
      6.0*Fcal_lnv_6*logv) + Fcal_4*(12.0*Fcal_5 + 12.0*Fcal_SO_5) + 12.0*Fcal_5*Fcal_SQ_4 + 12.0*Fcal_6*Fcal_SO_3 +
      Fcal_SO_3*(12.0*Fcal_SO_6 + 12.0*Fcal_lnv_6*logv) + 12.0*Fcal_SO_5*Fcal_SQ_4) + Fcal_4*(Fcal_4*(2.0*Fcal_4 +
      6.0*Fcal_SQ_4) + 12.0*Fcal_5*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SO_5 + 6.0*pow(Fcal_SQ_4, 2)) +
      12.0*Fcal_5*Fcal_SO_3*Fcal_SQ_4 + 6.0*Fcal_6*pow(Fcal_SO_3, 2) + Fcal_SO_3*(Fcal_SO_3*(6.0*Fcal_SO_6 +
      6.0*Fcal_lnv_6*logv) + 12.0*Fcal_SO_5*Fcal_SQ_4) + 2.0*pow(Fcal_SQ_4, 3)) + E_2*(Fcal_2*(Fcal_2*(12.0*Fcal_6 +
      12.0*Fcal_SO_6 + 12.0*Fcal_lnv_6*logv) + Fcal_3*(24.0*Fcal_5 + 24.0*Fcal_SO_5) + Fcal_4*(12.0*Fcal_4 +
      24.0*Fcal_SQ_4) + 24.0*Fcal_5*Fcal_SO_3 + 24.0*Fcal_SO_3*Fcal_SO_5 + 12.0*pow(Fcal_SQ_4, 2)) +
      Fcal_3*(Fcal_3*(12.0*Fcal_4 + 12.0*Fcal_SQ_4) + 24.0*Fcal_4*Fcal_SO_3 + 24.0*Fcal_SO_3*Fcal_SQ_4) +
      12.0*Fcal_4*pow(Fcal_SO_3, 2) + 12.0*pow(Fcal_SO_3, 2)*Fcal_SQ_4) + E_4*Fcal_2*(Fcal_2*(18.0*Fcal_4 +
      18.0*Fcal_SQ_4) + Fcal_3*(18.0*Fcal_3 + 36.0*Fcal_SO_3) + 18.0*pow(Fcal_SO_3, 2)) + 8.0*E_6*pow(Fcal_2, 3) +
      E_SO_3*(Fcal_2*(Fcal_2*(15.0*Fcal_5 + 15.0*Fcal_SO_5) + Fcal_3*(30.0*Fcal_4 + 30.0*Fcal_SQ_4) +
      30.0*Fcal_4*Fcal_SO_3 + 30.0*Fcal_SO_3*Fcal_SQ_4) + Fcal_3*(Fcal_3*(5.0*Fcal_3 + 15.0*Fcal_SO_3) +
      15.0*pow(Fcal_SO_3, 2)) + 5.0*pow(Fcal_SO_3, 3)) + E_SO_5*pow(Fcal_2, 2)*(21.0*Fcal_3 + 21.0*Fcal_SO_3) +
      E_SQ_4*Fcal_2*(Fcal_2*(18.0*Fcal_4 + 18.0*Fcal_SQ_4) + Fcal_3*(18.0*Fcal_3 + 36.0*Fcal_SO_3) + 18.0*pow(Fcal_SO_3,
      2)) + (E_0*(Fcal_2*(Fcal_2*(Fcal_2*(-8.0*Fcal_6 - 8.0*Fcal_SO_6 - 8.0*Fcal_lnv_6*logv) + Fcal_3*(-24.0*Fcal_5 -
      24.0*Fcal_SO_5) + Fcal_4*(-12.0*Fcal_4 - 24.0*Fcal_SQ_4) - 24.0*Fcal_5*Fcal_SO_3 - 24.0*Fcal_SO_3*Fcal_SO_5 -
      12.0*pow(Fcal_SQ_4, 2)) + Fcal_3*(Fcal_3*(-24.0*Fcal_4 - 24.0*Fcal_SQ_4) - 48.0*Fcal_4*Fcal_SO_3 -
      48.0*Fcal_SO_3*Fcal_SQ_4) - 24.0*Fcal_4*pow(Fcal_SO_3, 2) - 24.0*pow(Fcal_SO_3, 2)*Fcal_SQ_4) +
      Fcal_3*(Fcal_3*(Fcal_3*(-2.0*Fcal_3 - 8.0*Fcal_SO_3) - 12.0*pow(Fcal_SO_3, 2)) - 8.0*pow(Fcal_SO_3, 3)) -
      2.0*pow(Fcal_SO_3, 4)) + E_2*pow(Fcal_2, 2)*(Fcal_2*(-16.0*Fcal_4 - 16.0*Fcal_SQ_4) + Fcal_3*(-24.0*Fcal_3 -
      48.0*Fcal_SO_3) - 24.0*pow(Fcal_SO_3, 2)) - 6.0*E_4*pow(Fcal_2, 4) + E_SO_3*pow(Fcal_2, 3)*(-20.0*Fcal_3 -
      20.0*Fcal_SO_3) - 6.0*E_SQ_4*pow(Fcal_2, 4) + (E_0*pow(Fcal_2, 3)*(Fcal_2*(10.0*Fcal_4 + 10.0*Fcal_SQ_4) +
      Fcal_3*(20.0*Fcal_3 + 40.0*Fcal_SO_3) + 20.0*pow(Fcal_SO_3, 2)) - 2.0*E_0*pow(Fcal_2, 6)/Fcal_0 +
      4.0*E_2*pow(Fcal_2, 5))/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0 + m*(-13.0*E_11 + (E_0*(2.0*Fcal_11 +
      2.0*Fcal_lnv_11*logv) + E_2*(4.0*Fcal_9 + 4.0*Fcal_lnv_9*logv) + E_4*(6.0*Fcal_7 + 6.0*Fcal_SO_7) +
      E_6*(8.0*Fcal_5 + 8.0*Fcal_SO_5) + E_8*(10.0*Fcal_3 + 10.0*Fcal_SO_3) + E_SO_3*(5.0*Fcal_8 + 5.0*Fcal_SO_8 +
      5.0*Fcal_lnv_8*logv) + E_SO_5*(7.0*Fcal_6 + 7.0*Fcal_SO_6 + 7.0*Fcal_lnv_6*logv) + E_SO_7*(9.0*Fcal_4 +
      9.0*Fcal_SQ_4) + E_SQ_4*(6.0*Fcal_7 + 6.0*Fcal_SO_7) + E_lnv_8*(Fcal_3*(10.0*logv + 1.0) + Fcal_SO_3*(10.0*logv +
      1.0)) + (E_0*(Fcal_2*(-4.0*Fcal_9 - 4.0*Fcal_lnv_9*logv) + Fcal_3*(-4.0*Fcal_8 - 4.0*Fcal_SO_8 -
      4.0*Fcal_lnv_8*logv) + Fcal_4*(-4.0*Fcal_7 - 4.0*Fcal_SO_7) + Fcal_5*(-4.0*Fcal_6 - 4.0*Fcal_SO_6 -
      4.0*Fcal_lnv_6*logv) - 4.0*Fcal_6*Fcal_SO_5 - 4.0*Fcal_7*Fcal_SQ_4 - 4.0*Fcal_8*Fcal_SO_3 +
      Fcal_SO_3*(-4.0*Fcal_SO_8 - 4.0*Fcal_lnv_8*logv) + Fcal_SO_5*(-4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) -
      4.0*Fcal_SO_7*Fcal_SQ_4) + E_2*(Fcal_2*(-8.0*Fcal_7 - 8.0*Fcal_SO_7) + Fcal_3*(-8.0*Fcal_6 - 8.0*Fcal_SO_6 -
      8.0*Fcal_lnv_6*logv) + Fcal_4*(-8.0*Fcal_5 - 8.0*Fcal_SO_5) - 8.0*Fcal_5*Fcal_SQ_4 - 8.0*Fcal_6*Fcal_SO_3 +
      Fcal_SO_3*(-8.0*Fcal_SO_6 - 8.0*Fcal_lnv_6*logv) - 8.0*Fcal_SO_5*Fcal_SQ_4) + E_4*(Fcal_2*(-12.0*Fcal_5 -
      12.0*Fcal_SO_5) + Fcal_3*(-12.0*Fcal_4 - 12.0*Fcal_SQ_4) - 12.0*Fcal_4*Fcal_SO_3 - 12.0*Fcal_SO_3*Fcal_SQ_4) +
      E_6*Fcal_2*(-16.0*Fcal_3 - 16.0*Fcal_SO_3) + E_SO_3*(Fcal_2*(-10.0*Fcal_6 - 10.0*Fcal_SO_6 - 10.0*Fcal_lnv_6*logv)
      + Fcal_3*(-10.0*Fcal_5 - 10.0*Fcal_SO_5) + Fcal_4*(-5.0*Fcal_4 - 10.0*Fcal_SQ_4) - 10.0*Fcal_5*Fcal_SO_3 -
      10.0*Fcal_SO_3*Fcal_SO_5 - 5.0*pow(Fcal_SQ_4, 2)) + E_SO_5*(Fcal_2*(-14.0*Fcal_4 - 14.0*Fcal_SQ_4) +
      Fcal_3*(-7.0*Fcal_3 - 14.0*Fcal_SO_3) - 7.0*pow(Fcal_SO_3, 2)) - 9.0*E_SO_7*pow(Fcal_2, 2) +
      E_SQ_4*(Fcal_2*(-12.0*Fcal_5 - 12.0*Fcal_SO_5) + Fcal_3*(-12.0*Fcal_4 - 12.0*Fcal_SQ_4) - 12.0*Fcal_4*Fcal_SO_3 -
      12.0*Fcal_SO_3*Fcal_SQ_4) + (E_0*(Fcal_2*(Fcal_2*(6.0*Fcal_7 + 6.0*Fcal_SO_7) + Fcal_3*(12.0*Fcal_6 +
      12.0*Fcal_SO_6 + 12.0*Fcal_lnv_6*logv) + Fcal_4*(12.0*Fcal_5 + 12.0*Fcal_SO_5) + 12.0*Fcal_5*Fcal_SQ_4 +
      12.0*Fcal_6*Fcal_SO_3 + Fcal_SO_3*(12.0*Fcal_SO_6 + 12.0*Fcal_lnv_6*logv) + 12.0*Fcal_SO_5*Fcal_SQ_4) +
      Fcal_3*(Fcal_3*(6.0*Fcal_5 + 6.0*Fcal_SO_5) + Fcal_4*(6.0*Fcal_4 + 12.0*Fcal_SQ_4) + 12.0*Fcal_5*Fcal_SO_3 +
      12.0*Fcal_SO_3*Fcal_SO_5 + 6.0*pow(Fcal_SQ_4, 2)) + Fcal_4*(6.0*Fcal_4*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SQ_4) +
      6.0*Fcal_5*pow(Fcal_SO_3, 2) + Fcal_SO_3*(6.0*Fcal_SO_3*Fcal_SO_5 + 6.0*pow(Fcal_SQ_4, 2))) +
      E_2*(Fcal_2*(Fcal_2*(12.0*Fcal_5 + 12.0*Fcal_SO_5) + Fcal_3*(24.0*Fcal_4 + 24.0*Fcal_SQ_4) + 24.0*Fcal_4*Fcal_SO_3
      + 24.0*Fcal_SO_3*Fcal_SQ_4) + Fcal_3*(Fcal_3*(4.0*Fcal_3 + 12.0*Fcal_SO_3) + 12.0*pow(Fcal_SO_3, 2)) +
      4.0*pow(Fcal_SO_3, 3)) + E_4*pow(Fcal_2, 2)*(18.0*Fcal_3 + 18.0*Fcal_SO_3) + E_SO_3*Fcal_2*(Fcal_2*(15.0*Fcal_4 +
      15.0*Fcal_SQ_4) + Fcal_3*(15.0*Fcal_3 + 30.0*Fcal_SO_3) + 15.0*pow(Fcal_SO_3, 2)) + 7.0*E_SO_5*pow(Fcal_2, 3) +
      E_SQ_4*pow(Fcal_2, 2)*(18.0*Fcal_3 + 18.0*Fcal_SO_3) + (E_0*Fcal_2*(Fcal_2*(Fcal_2*(-8.0*Fcal_5 - 8.0*Fcal_SO_5) +
      Fcal_3*(-24.0*Fcal_4 - 24.0*Fcal_SQ_4) - 24.0*Fcal_4*Fcal_SO_3 - 24.0*Fcal_SO_3*Fcal_SQ_4) +
      Fcal_3*(Fcal_3*(-8.0*Fcal_3 - 24.0*Fcal_SO_3) - 24.0*pow(Fcal_SO_3, 2)) - 8.0*pow(Fcal_SO_3, 3)) + E_0*pow(Fcal_2,
      4)*(10.0*Fcal_3 + 10.0*Fcal_SO_3)/Fcal_0 + E_2*pow(Fcal_2, 3)*(-16.0*Fcal_3 - 16.0*Fcal_SO_3) -
      5.0*E_SO_3*pow(Fcal_2, 4))/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0) + m*(-12.0*E_10 + E_lnv_10*(-12.0*logv - 1.0) +
      (E_0*(2.0*Fcal_10 + 2.0*Fcal_lnv_10*logv) + E_2*(4.0*Fcal_8 + 4.0*Fcal_SO_8 + 4.0*Fcal_lnv_8*logv) +
      E_4*(6.0*Fcal_6 + 6.0*Fcal_SO_6 + 6.0*Fcal_lnv_6*logv) + E_6*(8.0*Fcal_4 + 8.0*Fcal_SQ_4) + 10.0*E_8*Fcal_2 +
      E_SO_3*(5.0*Fcal_7 + 5.0*Fcal_SO_7) + E_SO_5*(7.0*Fcal_5 + 7.0*Fcal_SO_5) + E_SO_7*(9.0*Fcal_3 + 9.0*Fcal_SO_3) +
      E_SQ_4*(6.0*Fcal_6 + 6.0*Fcal_SO_6 + 6.0*Fcal_lnv_6*logv) + E_lnv_8*Fcal_2*(10.0*logv + 1.0) +
      (E_0*(Fcal_2*(-4.0*Fcal_8 - 4.0*Fcal_SO_8 - 4.0*Fcal_lnv_8*logv) + Fcal_3*(-4.0*Fcal_7 - 4.0*Fcal_SO_7) +
      Fcal_4*(-4.0*Fcal_6 - 4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) + Fcal_5*(-2.0*Fcal_5 - 4.0*Fcal_SO_5) -
      4.0*Fcal_6*Fcal_SQ_4 - 4.0*Fcal_7*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SO_7 - 2.0*pow(Fcal_SO_5, 2) -
      4.0*Fcal_SO_6*Fcal_SQ_4 - 4.0*Fcal_SQ_4*Fcal_lnv_6*logv) + E_2*(Fcal_2*(-8.0*Fcal_6 - 8.0*Fcal_SO_6 -
      8.0*Fcal_lnv_6*logv) + Fcal_3*(-8.0*Fcal_5 - 8.0*Fcal_SO_5) + Fcal_4*(-4.0*Fcal_4 - 8.0*Fcal_SQ_4) -
      8.0*Fcal_5*Fcal_SO_3 - 8.0*Fcal_SO_3*Fcal_SO_5 - 4.0*pow(Fcal_SQ_4, 2)) + E_4*(Fcal_2*(-12.0*Fcal_4 -
      12.0*Fcal_SQ_4) + Fcal_3*(-6.0*Fcal_3 - 12.0*Fcal_SO_3) - 6.0*pow(Fcal_SO_3, 2)) - 8.0*E_6*pow(Fcal_2, 2) +
      E_SO_3*(Fcal_2*(-10.0*Fcal_5 - 10.0*Fcal_SO_5) + Fcal_3*(-10.0*Fcal_4 - 10.0*Fcal_SQ_4) - 10.0*Fcal_4*Fcal_SO_3 -
      10.0*Fcal_SO_3*Fcal_SQ_4) + E_SO_5*Fcal_2*(-14.0*Fcal_3 - 14.0*Fcal_SO_3) + E_SQ_4*(Fcal_2*(-12.0*Fcal_4 -
      12.0*Fcal_SQ_4) + Fcal_3*(-6.0*Fcal_3 - 12.0*Fcal_SO_3) - 6.0*pow(Fcal_SO_3, 2)) +
      (E_0*(Fcal_2*(Fcal_2*(6.0*Fcal_6 + 6.0*Fcal_SO_6 + 6.0*Fcal_lnv_6*logv) + Fcal_3*(12.0*Fcal_5 + 12.0*Fcal_SO_5) +
      Fcal_4*(6.0*Fcal_4 + 12.0*Fcal_SQ_4) + 12.0*Fcal_5*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SO_5 + 6.0*pow(Fcal_SQ_4, 2)) +
      Fcal_3*(Fcal_3*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + 12.0*Fcal_4*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SQ_4) +
      6.0*Fcal_4*pow(Fcal_SO_3, 2) + 6.0*pow(Fcal_SO_3, 2)*Fcal_SQ_4) + E_2*Fcal_2*(Fcal_2*(12.0*Fcal_4 +
      12.0*Fcal_SQ_4) + Fcal_3*(12.0*Fcal_3 + 24.0*Fcal_SO_3) + 12.0*pow(Fcal_SO_3, 2)) + 6.0*E_4*pow(Fcal_2, 3) +
      E_SO_3*pow(Fcal_2, 2)*(15.0*Fcal_3 + 15.0*Fcal_SO_3) + 6.0*E_SQ_4*pow(Fcal_2, 3) + (E_0*pow(Fcal_2,
      2)*(Fcal_2*(-8.0*Fcal_4 - 8.0*Fcal_SQ_4) + Fcal_3*(-12.0*Fcal_3 - 24.0*Fcal_SO_3) - 12.0*pow(Fcal_SO_3, 2)) +
      2.0*E_0*pow(Fcal_2, 5)/Fcal_0 - 4.0*E_2*pow(Fcal_2, 4))/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0) +
      m*(E_0*(2.0*Fcal_9 + 2.0*Fcal_lnv_9*logv) + E_2*(4.0*Fcal_7 + 4.0*Fcal_SO_7) + E_4*(6.0*Fcal_5 + 6.0*Fcal_SO_5) +
      E_6*(8.0*Fcal_3 + 8.0*Fcal_SO_3) + E_SO_3*(5.0*Fcal_6 + 5.0*Fcal_SO_6 + 5.0*Fcal_lnv_6*logv) + E_SO_5*(7.0*Fcal_4
      + 7.0*Fcal_SQ_4) + 9.0*E_SO_7*Fcal_2 + E_SQ_4*(6.0*Fcal_5 + 6.0*Fcal_SO_5) + (E_0*(Fcal_2*(-4.0*Fcal_7 -
      4.0*Fcal_SO_7) + Fcal_3*(-4.0*Fcal_6 - 4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) + Fcal_4*(-4.0*Fcal_5 - 4.0*Fcal_SO_5)
      - 4.0*Fcal_5*Fcal_SQ_4 - 4.0*Fcal_6*Fcal_SO_3 + Fcal_SO_3*(-4.0*Fcal_SO_6 - 4.0*Fcal_lnv_6*logv) -
      4.0*Fcal_SO_5*Fcal_SQ_4) + E_2*(Fcal_2*(-8.0*Fcal_5 - 8.0*Fcal_SO_5) + Fcal_3*(-8.0*Fcal_4 - 8.0*Fcal_SQ_4) -
      8.0*Fcal_4*Fcal_SO_3 - 8.0*Fcal_SO_3*Fcal_SQ_4) + E_4*Fcal_2*(-12.0*Fcal_3 - 12.0*Fcal_SO_3) +
      E_SO_3*(Fcal_2*(-10.0*Fcal_4 - 10.0*Fcal_SQ_4) + Fcal_3*(-5.0*Fcal_3 - 10.0*Fcal_SO_3) - 5.0*pow(Fcal_SO_3, 2)) -
      7.0*E_SO_5*pow(Fcal_2, 2) + E_SQ_4*Fcal_2*(-12.0*Fcal_3 - 12.0*Fcal_SO_3) + (E_0*(Fcal_2*(Fcal_2*(6.0*Fcal_5 +
      6.0*Fcal_SO_5) + Fcal_3*(12.0*Fcal_4 + 12.0*Fcal_SQ_4) + 12.0*Fcal_4*Fcal_SO_3 + 12.0*Fcal_SO_3*Fcal_SQ_4) +
      Fcal_3*(Fcal_3*(2.0*Fcal_3 + 6.0*Fcal_SO_3) + 6.0*pow(Fcal_SO_3, 2)) + 2.0*pow(Fcal_SO_3, 3)) + E_0*pow(Fcal_2,
      3)*(-8.0*Fcal_3 - 8.0*Fcal_SO_3)/Fcal_0 + E_2*pow(Fcal_2, 2)*(12.0*Fcal_3 + 12.0*Fcal_SO_3) +
      5.0*E_SO_3*pow(Fcal_2, 3))/Fcal_0)/Fcal_0)/pow(Fcal_0, 2)) + m*(-10.0*E_8 + E_lnv_8*(-10.0*logv - 1.0) +
      (E_0*(2.0*Fcal_8 + 2.0*Fcal_SO_8 + 2.0*Fcal_lnv_8*logv) + E_2*(4.0*Fcal_6 + 4.0*Fcal_SO_6 + 4.0*Fcal_lnv_6*logv) +
      E_4*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + 8.0*E_6*Fcal_2 + E_SO_3*(5.0*Fcal_5 + 5.0*Fcal_SO_5) + E_SO_5*(7.0*Fcal_3 +
      7.0*Fcal_SO_3) + E_SQ_4*(6.0*Fcal_4 + 6.0*Fcal_SQ_4) + (E_0*(Fcal_2*(-4.0*Fcal_6 - 4.0*Fcal_SO_6 -
      4.0*Fcal_lnv_6*logv) + Fcal_3*(-4.0*Fcal_5 - 4.0*Fcal_SO_5) + Fcal_4*(-2.0*Fcal_4 - 4.0*Fcal_SQ_4) -
      4.0*Fcal_5*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SO_5 - 2.0*pow(Fcal_SQ_4, 2)) + E_2*(Fcal_2*(-8.0*Fcal_4 -
      8.0*Fcal_SQ_4) + Fcal_3*(-4.0*Fcal_3 - 8.0*Fcal_SO_3) - 4.0*pow(Fcal_SO_3, 2)) - 6.0*E_4*pow(Fcal_2, 2) +
      E_SO_3*Fcal_2*(-10.0*Fcal_3 - 10.0*Fcal_SO_3) - 6.0*E_SQ_4*pow(Fcal_2, 2) + (E_0*Fcal_2*(Fcal_2*(6.0*Fcal_4 +
      6.0*Fcal_SQ_4) + Fcal_3*(6.0*Fcal_3 + 12.0*Fcal_SO_3) + 6.0*pow(Fcal_SO_3, 2)) - 2.0*E_0*pow(Fcal_2, 4)/Fcal_0 +
      4.0*E_2*pow(Fcal_2, 3))/Fcal_0)/Fcal_0)/Fcal_0)/Fcal_0) + m*(-9.0*E_SO_7 + (E_0*(2.0*Fcal_7 + 2.0*Fcal_SO_7) +
      E_2*(4.0*Fcal_5 + 4.0*Fcal_SO_5) + E_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) + E_SO_3*(5.0*Fcal_4 + 5.0*Fcal_SQ_4) +
      7.0*E_SO_5*Fcal_2 + E_SQ_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) + (E_0*(Fcal_2*(-4.0*Fcal_5 - 4.0*Fcal_SO_5) +
      Fcal_3*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) - 4.0*Fcal_4*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SQ_4) + E_0*pow(Fcal_2,
      2)*(6.0*Fcal_3 + 6.0*Fcal_SO_3)/Fcal_0 + E_2*Fcal_2*(-8.0*Fcal_3 - 8.0*Fcal_SO_3) - 5.0*E_SO_3*pow(Fcal_2,
      2))/Fcal_0)/Fcal_0)/Fcal_0) + m*(-8.0*E_6 + (E_0*(2.0*Fcal_6 + 2.0*Fcal_SO_6 + 2.0*Fcal_lnv_6*logv) +
      E_2*(4.0*Fcal_4 + 4.0*Fcal_SQ_4) + 6.0*E_4*Fcal_2 + E_SO_3*(5.0*Fcal_3 + 5.0*Fcal_SO_3) + 6.0*E_SQ_4*Fcal_2 +
      (E_0*(Fcal_2*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) + Fcal_3*(-2.0*Fcal_3 - 4.0*Fcal_SO_3) - 2.0*pow(Fcal_SO_3, 2)) +
      2.0*E_0*pow(Fcal_2, 3)/Fcal_0 - 4.0*E_2*pow(Fcal_2, 2))/Fcal_0)/Fcal_0)/Fcal_0) + m*(-7.0*E_SO_5 +
      (E_0*(2.0*Fcal_5 + 2.0*Fcal_SO_5) + E_0*Fcal_2*(-4.0*Fcal_3 - 4.0*Fcal_SO_3)/Fcal_0 + E_2*(4.0*Fcal_3 +
      4.0*Fcal_SO_3) + 5.0*E_SO_3*Fcal_2)/Fcal_0)/Fcal_0) + m*(-6.0*E_4 - 6.0*E_SQ_4 + (E_0*(2.0*Fcal_4 + 2.0*Fcal_SQ_4)
      - 2.0*E_0*pow(Fcal_2, 2)/Fcal_0 + 4.0*E_2*Fcal_2)/Fcal_0)/Fcal_0) + m*(E_0*(2.0*Fcal_3 + 2.0*Fcal_SO_3)/Fcal_0 -
      5.0*E_SO_3)/Fcal_0) + m*(2.0*E_0*Fcal_2/Fcal_0 - 4.0*E_2)/Fcal_0))/Fcal_coeff;
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ellHat = FrameFromAngularVelocity_Integrand(rfrak_ellHat, OmegaVec().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], OmegaVec_chiVec_1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], OmegaVec_chiVec_2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ellHat[0];
    dydt[6] = rfrakdot_ellHat[1];
    dydt[7] = rfrakdot_ellHat[2];
    dydt[8] = v*v*v/m;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
}; // class TaylorTn_6p0PN_Q : public TaylorTn_Q
