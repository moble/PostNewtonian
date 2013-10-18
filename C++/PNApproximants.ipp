// File produced automatically by OrbitalEvolution.ipynb

void cross(const double a[3], const double b[3], double c[3]) {
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = -a[0]*b[2]+a[2]*b[0];
  c[2] = a[0]*b[1]-a[1]*b[0];
  return;
}

class TaylorTn {
private:
  const double m1;
  double v, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, ellHat_x, ellHat_y, ellHat_z, nHat_x, nHat_y, nHat_z;
  const double m2, delta, nu;
  std::vector<double> ellHat, nHat;
  const double chi1chi1, chi2chi2;
  double chi1_l, chi1_n, chi1_la, chi2_l, chi2_n, chi2_la, S_l, S_n, Sigma_l, Sigma_n, logv, Fcal_coeff;
  const double Fcal_0, Fcal_2, Fcal_3, Fcal_4, Fcal_5, Fcal_6, Fcal_lnv_6, Fcal_7;
  double Fcal_SQ_4, Fcal_SO_3, Fcal_SO_5, Fcal_SO_6, Fcal_SO_7;
  const double E_0, E_2, E_4, E_6;
  double E_SQ_4, E_SO_3, E_SO_5, E_SO_7, MDot_Alvi_5;
  double Phi, gamma;

public:
  TaylorTn(const double m1_i, const double v_i, const double chi1_x_i, const double chi1_y_i, const double chi1_z_i,
           const double chi2_x_i, const double chi2_y_i, const double chi2_z_i, const double ellHat_x_i, const double
           ellHat_y_i, const double ellHat_z_i, const double nHat_x_i, const double nHat_y_i, const double nHat_z_i) :
    m1(m1_i), v(v_i), chi1_x(chi1_x_i), chi1_y(chi1_y_i), chi1_z(chi1_z_i), chi2_x(chi2_x_i), chi2_y(chi2_y_i),
    chi2_z(chi2_z_i), ellHat_x(ellHat_x_i), ellHat_y(ellHat_y_i), ellHat_z(ellHat_z_i), nHat_x(nHat_x_i),
    nHat_y(nHat_y_i), nHat_z(nHat_z_i), m2(-m1 + 1.0), delta(m1 - m2), nu(m1*m2), ellHat(3), nHat(3),
    chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z, 2)), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z,
    2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x + chi1_y*nHat_y +
    chi1_z*nHat_z), chi1_la(chi1_x*(ellHat_y*nHat_z - ellHat_z*nHat_y) + chi1_y*(-ellHat_x*nHat_z + ellHat_z*nHat_x) +
    chi1_z*(ellHat_x*nHat_y - ellHat_y*nHat_x)), chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z),
    chi2_n(chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z), chi2_la(chi2_x*(ellHat_y*nHat_z - ellHat_z*nHat_y) +
    chi2_y*(-ellHat_x*nHat_z + ellHat_z*nHat_x) + chi2_z*(ellHat_x*nHat_y - ellHat_y*nHat_x)), S_l(chi1_l*pow(m1, 2) +
    chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2)), Sigma_l(-chi1_l*m1 + chi2_l*m2), Sigma_n(-chi1_n*m1
    + chi2_n*m2), logv(log(v)), Fcal_coeff(6.4*pow(nu, 2)*pow(v, 10)), Fcal_0(1.00000000000000),
    Fcal_2(-2.91666666666667*nu - 3.71130952380952), Fcal_3(12.5663706143592), Fcal_4(3.61111111111111*pow(nu, 2) +
    18.3948412698413*nu - 4.92846119929453), Fcal_5(-76.3145215434521*nu - 38.2928354546934),
    Fcal_6(-2.39197530864198*pow(nu, 3) - 31.2179232804233*pow(nu, 2) - 8.87205344238227*nu + 115.731716675611),
    Fcal_lnv_6(-16.3047619047619), Fcal_7(200.905057974359*pow(nu, 2) + 390.417427312002*nu - 101.509595959742),
    Fcal_SQ_4(chi1_l*(chi1_l*(1.03125*delta - 2.0625*nu + 1.03125) + 3.875*chi2_l*nu) +
    chi1_la*(chi1_la*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
    2.14583333333333*chi2_la*nu) + chi1_n*(chi1_n*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667)
    - 2.14583333333333*chi2_n*nu) + pow(chi2_l, 2)*(-1.03125*delta - 2.0625*nu + 1.03125) + pow(chi2_la,
    2)*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) + pow(chi2_n, 2)*(0.463541666666667*delta +
    0.927083333333333*nu - 0.463541666666667)), Fcal_SO_3(-4.0*S_l - 1.25*Sigma_l*delta),
    Fcal_SO_5(S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125)), Fcal_SO_6(-50.2654824574367*S_l -
    16.2315620435473*Sigma_l*delta), Fcal_SO_7(S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu +
    70.053644914756) + Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762)),
    E_0(1.00000000000000), E_2(-0.0833333333333333*nu - 0.75), E_4(-0.0416666666666667*pow(nu, 2) + 2.375*nu - 3.375),
    E_6(-0.00675154320987654*pow(nu, 3) - 1.61458333333333*pow(nu, 2) + 38.7246294907293*nu - 10.546875),
    E_SQ_4(-0.5*pow(chi1_l, 2)*delta + pow(chi1_l, 2)*nu - 0.5*pow(chi1_l, 2) - 2.0*chi1_l*chi2_l*nu + 0.25*pow(chi1_la,
    2)*delta - 0.5*pow(chi1_la, 2)*nu + 0.25*pow(chi1_la, 2) + chi1_la*chi2_la*nu + 0.25*pow(chi1_n, 2)*delta -
    0.5*pow(chi1_n, 2)*nu + 0.25*pow(chi1_n, 2) + chi1_n*chi2_n*nu + 0.5*pow(chi2_l, 2)*delta + pow(chi2_l, 2)*nu -
    0.5*pow(chi2_l, 2) - 0.25*pow(chi2_la, 2)*delta - 0.5*pow(chi2_la, 2)*nu + 0.25*pow(chi2_la, 2) - 0.25*pow(chi2_n,
    2)*delta - 0.5*pow(chi2_n, 2)*nu + 0.25*pow(chi2_n, 2)), E_SO_3(4.66666666666667*S_l + 2.0*Sigma_l*delta),
    E_SO_5(S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0)),
    E_SO_7(S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75) + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu + 6.75)),
    MDot_Alvi_5(-0.25*chi1_l*pow(m1, 3)*(3.0*chi1chi1 + 1.0) - 0.25*chi2_l*pow(m2, 3)*(3.0*chi2chi2 + 1.0)), Phi(0.0)
  {
    ellHat[0] = ellHat_x;
    ellHat[1] = ellHat_y;
    ellHat[2] = ellHat_z;
    nHat[0] = nHat_x;
    nHat[1] = nHat_y;
    nHat[2] = nHat_z;
  }

  void Recalculate(double t, const double* y) {
    v = y[0];
    chi1_x = y[1];
    chi1_y = y[2];
    chi1_z = y[3];
    chi2_x = y[4];
    chi2_y = y[5];
    chi2_z = y[6];
    ellHat_x = y[7];
    ellHat_y = y[8];
    ellHat_z = y[9];
    Phi = y[10];
    gamma = y[11];

    const Quaternion Rax =
        sqrtOfRotor(-normalized(Quaternion(0., ellHat_x, ellHat_y, ellHat_z))*zHat);
    const Quaternion R = Rax * exp(((gamma+Phi)/2.)*zHat);
    const Quaternion nHatQ = R*xHat*R.conjugate();
    nHat_x = nHatQ[1];
    nHat_y = nHatQ[2];
    nHat_z = nHatQ[3];

    ellHat[0] = ellHat_x;
    ellHat[1] = ellHat_y;
    ellHat[2] = ellHat_z;
    nHat[0] = nHat_x;
    nHat[1] = nHat_y;
    nHat[2] = nHat_z;
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi1_la = chi1_x*(ellHat_y*nHat_z - ellHat_z*nHat_y) + chi1_y*(-ellHat_x*nHat_z + ellHat_z*nHat_x) +
      chi1_z*(ellHat_x*nHat_y - ellHat_y*nHat_x);
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    chi2_la = chi2_x*(ellHat_y*nHat_z - ellHat_z*nHat_y) + chi2_y*(-ellHat_x*nHat_z + ellHat_z*nHat_x) +
      chi2_z*(ellHat_x*nHat_y - ellHat_y*nHat_x);
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_n*pow(m2, 2);
    Sigma_l = -chi1_l*m1 + chi2_l*m2;
    Sigma_n = -chi1_n*m1 + chi2_n*m2;
    logv = log(v);
    Fcal_coeff = 6.4*pow(nu, 2)*pow(v, 10);
    Fcal_SQ_4 = chi1_l*(chi1_l*(1.03125*delta - 2.0625*nu + 1.03125) + 3.875*chi2_l*nu) +
      chi1_la*(chi1_la*(-0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) -
      2.14583333333333*chi2_la*nu) + chi1_n*(chi1_n*(-0.463541666666667*delta + 0.927083333333333*nu -
      0.463541666666667) - 2.14583333333333*chi2_n*nu) + pow(chi2_l, 2)*(-1.03125*delta - 2.0625*nu + 1.03125) +
      pow(chi2_la, 2)*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667) + pow(chi2_n,
      2)*(0.463541666666667*delta + 0.927083333333333*nu - 0.463541666666667);
    Fcal_SO_3 = -4.0*S_l - 1.25*Sigma_l*delta;
    Fcal_SO_5 = S_l*(30.2222222222222*nu - 4.5) + Sigma_l*delta*(10.75*nu - 0.8125);
    Fcal_SO_6 = -50.2654824574367*S_l - 16.2315620435473*Sigma_l*delta;
    Fcal_SO_7 = S_l*(-104.074074074074*pow(nu, 2) + 32.6560846560847*nu + 70.053644914756) +
      Sigma_l*delta*(-41.6944444444444*pow(nu, 2) + 14.6746031746032*nu + 28.3779761904762);
    E_SQ_4 = -0.5*pow(chi1_l, 2)*delta + pow(chi1_l, 2)*nu - 0.5*pow(chi1_l, 2) - 2.0*chi1_l*chi2_l*nu +
      0.25*pow(chi1_la, 2)*delta - 0.5*pow(chi1_la, 2)*nu + 0.25*pow(chi1_la, 2) + chi1_la*chi2_la*nu + 0.25*pow(chi1_n,
      2)*delta - 0.5*pow(chi1_n, 2)*nu + 0.25*pow(chi1_n, 2) + chi1_n*chi2_n*nu + 0.5*pow(chi2_l, 2)*delta + pow(chi2_l,
      2)*nu - 0.5*pow(chi2_l, 2) - 0.25*pow(chi2_la, 2)*delta - 0.5*pow(chi2_la, 2)*nu + 0.25*pow(chi2_la, 2) -
      0.25*pow(chi2_n, 2)*delta - 0.5*pow(chi2_n, 2)*nu + 0.25*pow(chi2_n, 2);
    E_SO_3 = 4.66666666666667*S_l + 2.0*Sigma_l*delta;
    E_SO_5 = S_l*(-6.77777777777778*nu + 11.0) + Sigma_l*delta*(-3.33333333333333*nu + 3.0);
    E_SO_7 = S_l*(2.41666666666667*pow(nu, 2) - 91.75*nu + 33.75) + Sigma_l*delta*(1.25*pow(nu, 2) - 39.0*nu + 6.75);
    MDot_Alvi_5 = -0.25*chi1_l*pow(m1, 3)*(3.0*chi1chi1 + 1.0) - 0.25*chi2_l*pow(m2, 3)*(3.0*chi2chi2 + 1.0);
  }

  std::vector<double> OmegaVec_chiVec_1() {
    double OmegaVec1_2 = delta*(0.625*nu - 0.5625) - 0.0416666666666667*pow(nu, 2) + 1.25*nu + 0.5625;
    double OmegaVec1_0 = -0.75*delta + 0.5*nu + 0.75;
    double OmegaVec1_4 = delta*(-0.15625*pow(nu, 2) + 4.875*nu - 0.84375) - 0.0208333333333333*pow(nu, 3) -
      3.28125*pow(nu, 2) + 0.1875*nu + 0.84375;
    double OmegaVec1_coeff = pow(v, 5);
    return OmegaVec1_coeff*ellHat*(OmegaVec1_0 + pow(v, 2)*(OmegaVec1_2 + OmegaVec1_4*pow(v, 2)));
  }
  std::vector<double> OmegaVec_chiVec_2() {
    double OmegaVec2_0 = 0.75*delta + 0.5*nu + 0.75;
    double OmegaVec2_coeff = pow(v, 5);
    double OmegaVec2_2 = -delta*(0.625*nu - 0.5625) - 0.0416666666666667*pow(nu, 2) + 1.25*nu + 0.5625;
    double OmegaVec2_4 = -delta*(-0.15625*pow(nu, 2) + 4.875*nu - 0.84375) - 0.0208333333333333*pow(nu, 3) -
      3.28125*pow(nu, 2) + 0.1875*nu + 0.84375;
    return OmegaVec2_coeff*ellHat*(OmegaVec2_0 + pow(v, 2)*(OmegaVec2_2 + OmegaVec2_4*pow(v, 2)));
  }
  std::vector<double> OmegaVec_ellHat() {
    double OmegaVec_ellHat_4 = S_n*(9.0*pow(nu, 2) - 29.5*nu - 1.5) + Sigma_n*delta*(4.33333333333333*pow(nu, 2) -
      9.625*nu - 1.5);
    double OmegaVec_ellHat_coeff = pow(v, 6);
    double OmegaVec_ellHat_2 = S_n*(-12.0*nu - 3.0) + Sigma_n*delta*(-5.5*nu - 3.0);
    double OmegaVec_ellHat_0 = 7.0*S_n + 3.0*Sigma_n*delta;
    return OmegaVec_ellHat_coeff*nHat*(OmegaVec_ellHat_0 + pow(v, 2)*(OmegaVec_ellHat_2 + OmegaVec_ellHat_4*pow(v, 2)));
  }

  int TaylorT1_3p5PN(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = Fcal_coeff*(Fcal_0 + pow(v, 2)*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 +
      v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7)))))));
    const double dEdv = -0.5*nu*v*(2.0*E_0 + pow(v, 2)*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 +
      v*(7.0*E_SO_5 + v*(8.0*E_6 + 9.0*E_SO_7*v))))));
    const double Absorption = Fcal_coeff*MDot_Alvi_5*pow(v, 5);
    const double dvdt_T1 = (-Absorption - Flux)/dEdv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, y, dydt);
  }

  int TaylorT4_3p5PN(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = -2.0*Fcal_coeff*(pow(v, 2)*(v*(v*(v*(v*(v*(-0.5*Fcal_7 - 0.5*Fcal_SO_7 + (E_2*(1.0*Fcal_5 +
      1.0*Fcal_SO_5 + 1.0*MDot_Alvi_5) + E_4*(1.5*Fcal_3 + 1.5*Fcal_SO_3) + E_SO_3*(1.25*Fcal_4 + 1.25*Fcal_SQ_4) +
      1.75*E_SO_5*Fcal_2 + 2.25*E_SO_7*Fcal_0 + E_SQ_4*(1.5*Fcal_3 + 1.5*Fcal_SO_3) + (E_2*(E_2*(-2.0*Fcal_3 -
      2.0*Fcal_SO_3) - 5.0*E_SO_3*Fcal_2 - 7.0*E_SO_5*Fcal_0) - 7.5*E_4*E_SO_3*Fcal_0 - 7.5*E_SO_3*E_SQ_4*Fcal_0 +
      15.0*pow(E_2, 2)*E_SO_3*Fcal_0/E_0)/E_0)/E_0)/E_0 + (-0.5*Fcal_6 - 0.5*Fcal_SO_6 - 0.5*Fcal_lnv_6*logv +
      (E_2*(1.0*Fcal_4 + 1.0*Fcal_SQ_4) + 1.5*E_4*Fcal_2 + 2.0*E_6*Fcal_0 + E_SO_3*(1.25*Fcal_3 + 1.25*Fcal_SO_3) +
      1.5*E_SQ_4*Fcal_2 + (E_2*(-2.0*E_2*Fcal_2 - 6.0*E_4*Fcal_0 - 6.0*E_SQ_4*Fcal_0) - 3.125*pow(E_SO_3, 2)*Fcal_0 +
      4.0*pow(E_2, 3)*Fcal_0/E_0)/E_0)/E_0)/E_0) + (-0.5*Fcal_5 - 0.5*Fcal_SO_5 - 0.5*MDot_Alvi_5 + (E_2*(1.0*Fcal_3 +
      1.0*Fcal_SO_3) + 1.25*E_SO_3*Fcal_2 + 1.75*E_SO_5*Fcal_0 - 5.0*E_2*E_SO_3*Fcal_0/E_0)/E_0)/E_0) + (-0.5*Fcal_4 -
      0.5*Fcal_SQ_4 + (1.0*E_2*Fcal_2 + 1.5*E_4*Fcal_0 + 1.5*E_SQ_4*Fcal_0 - 2.0*pow(E_2, 2)*Fcal_0/E_0)/E_0)/E_0) +
      (-0.5*Fcal_3 - 0.5*Fcal_SO_3 + 1.25*E_SO_3*Fcal_0/E_0)/E_0) + (-0.5*Fcal_2 + 1.0*E_2*Fcal_0/E_0)/E_0) -
      0.5*Fcal_0/E_0)/(nu*v);
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, y, dydt);
  }

  int TaylorT5_3p5PN(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = -0.5*nu*v*(-2.0*E_0/Fcal_0 + pow(v, 2)*(v*(v*(v*(v*(v*(-9.0*E_SO_7 + (E_0*(2.0*Fcal_7 +
      2.0*Fcal_SO_7) + E_2*(4.0*Fcal_5 + 4.0*Fcal_SO_5 + 4.0*MDot_Alvi_5) + E_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) +
      E_SO_3*(5.0*Fcal_4 + 5.0*Fcal_SQ_4) + 7.0*E_SO_5*Fcal_2 + E_SQ_4*(6.0*Fcal_3 + 6.0*Fcal_SO_3) +
      (E_0*(Fcal_2*(-4.0*Fcal_5 - 4.0*Fcal_SO_5 - 4.0*MDot_Alvi_5) + Fcal_3*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) -
      4.0*Fcal_4*Fcal_SO_3 - 4.0*Fcal_SO_3*Fcal_SQ_4) + E_0*pow(Fcal_2, 2)*(6.0*Fcal_3 + 6.0*Fcal_SO_3)/Fcal_0 +
      E_2*Fcal_2*(-8.0*Fcal_3 - 8.0*Fcal_SO_3) - 5.0*E_SO_3*pow(Fcal_2, 2))/Fcal_0)/Fcal_0)/Fcal_0 + (-8.0*E_6 +
      (E_0*(2.0*Fcal_6 + 2.0*Fcal_SO_6 + 2.0*Fcal_lnv_6*logv) + E_2*(4.0*Fcal_4 + 4.0*Fcal_SQ_4) + 6.0*E_4*Fcal_2 +
      E_SO_3*(5.0*Fcal_3 + 5.0*Fcal_SO_3) + 6.0*E_SQ_4*Fcal_2 + (E_0*(Fcal_2*(-4.0*Fcal_4 - 4.0*Fcal_SQ_4) +
      Fcal_3*(-2.0*Fcal_3 - 4.0*Fcal_SO_3) - 2.0*pow(Fcal_SO_3, 2)) + 2.0*E_0*pow(Fcal_2, 3)/Fcal_0 -
      4.0*E_2*pow(Fcal_2, 2))/Fcal_0)/Fcal_0)/Fcal_0) + (-7.0*E_SO_5 + (E_0*(2.0*Fcal_5 + 2.0*Fcal_SO_5 +
      2.0*MDot_Alvi_5) + E_0*Fcal_2*(-4.0*Fcal_3 - 4.0*Fcal_SO_3)/Fcal_0 + E_2*(4.0*Fcal_3 + 4.0*Fcal_SO_3) +
      5.0*E_SO_3*Fcal_2)/Fcal_0)/Fcal_0) + (-6.0*E_4 - 6.0*E_SQ_4 + (E_0*(2.0*Fcal_4 + 2.0*Fcal_SQ_4) -
      2.0*E_0*pow(Fcal_2, 2)/Fcal_0 + 4.0*E_2*Fcal_2)/Fcal_0)/Fcal_0) + (E_0*(2.0*Fcal_3 + 2.0*Fcal_SO_3)/Fcal_0 -
      5.0*E_SO_3)/Fcal_0) + (2.0*E_0*Fcal_2/Fcal_0 - 4.0*E_2)/Fcal_0))/Fcal_coeff;
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, y, dydt);
  }

  int CommonRHS(const double dvdt, const double* y, double* dydt) {
    dydt[0] = dvdt;
    cross(&OmegaVec_chiVec_1()[0], &y[1], &dydt[1]);
    cross(&OmegaVec_chiVec_2()[0], &y[4], &dydt[4]);
    cross(&OmegaVec_ellHat()[0], &y[7], &dydt[7]);
    dydt[10] = v*v*v;
    const Quaternion adot(0., dydt[7], dydt[8], dydt[9]);
    const Quaternion Rax =
        sqrtOfRotor(-Quaternion(0., ellHat_x, ellHat_y, ellHat_z).normalized()*zHat);
    const Quaternion Raxdot = ( (-1.0/std::sqrt(2+2*y[9]))*adot*zHat - (dydt[9]/(2+2*y[9]))*Rax );
    dydt[11] = -2*(Rax.conjugate() * Raxdot)[3];

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
};
