// File produced automatically by OrbitalEvolution.ipynb

class TaylorTn {
private:
  const Quaternion xHat, yHat, zHat;
  const double m1;
  double v;
  const double chi1Mag, chi2Mag;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z;
  const double m2, delta, nu, nu__2, nu__3;
  Quaternion R, nHat;
  double nHat_x, nHat_y, nHat_z;
  Quaternion lambdaHat;
  double lambdaHat_x, lambdaHat_y, lambdaHat_z;
  Quaternion ellHat;
  double ellHat_x, ellHat_y, ellHat_z;
  Quaternion R_S1, R_S2, chi1;
  double chi1_x, chi1_y, chi1_z;
  Quaternion chi2;
  double chi2_x, chi2_y, chi2_z;
  const double chi1chi1, chi2chi2;
  double chi1_l, chi1_n, chi1_la, chi2_l, chi2_n, chi2_la;
  const double sqrt1Mchi1chi1, sqrt1Mchi2chi2;
  double S_l, S_n, Sigma_l, Sigma_n, logv, Omega_chi1_ellHat, Omega_chi2_ellHat, Omega_ellHat_nHat;
  double Phi;

public:
  TaylorTn(const Quaternion xHat_i, const Quaternion yHat_i, const Quaternion zHat_i, const double m1_i, const double
           v_i, const double chi1Mag_i, const double chi2Mag_i, const double rfrak_chi1_x_i, const double
           rfrak_chi1_y_i, const double rfrak_chi2_x_i, const double rfrak_chi2_y_i, const double rfrak_ell_x_i, const
           double rfrak_ell_y_i, const double rfrak_ell_z_i) :
    xHat(xHat_i), yHat(yHat_i), zHat(zHat_i), m1(m1_i), v(v_i), chi1Mag(chi1Mag_i), chi2Mag(chi2Mag_i),
    rfrak_chi1_x(rfrak_chi1_x_i), rfrak_chi1_y(rfrak_chi1_y_i), rfrak_chi2_x(rfrak_chi2_x_i),
    rfrak_chi2_y(rfrak_chi2_y_i), rfrak_ell_x(rfrak_ell_x_i), rfrak_ell_y(rfrak_ell_y_i), rfrak_ell_z(rfrak_ell_z_i),
    m2(-1.0*m1 + 1.0), delta(m1 - m2), nu(m1*m2), nu__2(pow(m1, 2)*pow(m2, 2)), nu__3(pow(m1, 3)*pow(m2, 3)),
    R(exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat)), nHat(R*xHat*conjugate(R)), nHat_x(nHat[1]),
    nHat_y(nHat[2]), nHat_z(nHat[3]), lambdaHat(R*yHat*conjugate(R)), lambdaHat_x(lambdaHat[1]),
    lambdaHat_y(lambdaHat[2]), lambdaHat_z(lambdaHat[3]), ellHat(R*zHat*conjugate(R)), ellHat_x(ellHat[1]),
    ellHat_y(ellHat[2]), ellHat_z(ellHat[3]), R_S1(exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)),
    R_S2(exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)), chi1(chi1Mag*R_S1*zHat*conjugate(R_S1)), chi1_x(chi1[1]),
    chi1_y(chi1[2]), chi1_z(chi1[3]), chi2(chi2Mag*R_S2*zHat*conjugate(R_S2)), chi2_x(chi2[1]), chi2_y(chi2[2]),
    chi2_z(chi2[3]), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z, 2)), chi2chi2(pow(chi2_x, 2) + pow(chi2_y,
    2) + pow(chi2_z, 2)), chi1_l(chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z), chi1_n(chi1_x*nHat_x +
    chi1_y*nHat_y + chi1_z*nHat_z), chi1_la(chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z),
    chi2_l(chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z), chi2_n(chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z),
    chi2_la(chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z), sqrt1Mchi1chi1(sqrt(-chi1chi1 + 1.0)),
    sqrt1Mchi2chi2(sqrt(-chi2chi2 + 1.0)), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) +
    chi2_l*pow(m2, 2)), Sigma_l(-chi1_l*m1 + chi2_l*m2), Sigma_n(-chi1_n*m1 + chi2_n*m2), logv(log(v)),
    Omega_chi1_ellHat(pow(v, 5)*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu - 0.5625) +
    nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) +
    nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75)), Omega_chi2_ellHat(pow(v,
    5)*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu + 0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v,
    2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) +
    0.5625) + 0.75)), Omega_ellHat_nHat(pow(v, 6)*(7.0*S_n + 3.0*Sigma_n*delta + pow(v, 2)*(S_n*(-12.0*nu - 3.0) +
    Sigma_n*delta*(-5.5*nu - 3.0) + v*(S_l*(11.6666666666667*S_n + 5.0*Sigma_n*delta) + 7.0*S_n*Sigma_l*delta +
    3.0*Sigma_l*Sigma_n*pow(delta, 2) + v*(S_n*(nu*(9.0*nu - 29.5) - 1.5) + Sigma_n*delta*(nu*(4.33333333333333*nu -
    9.625) - 1.5) + v*(S_l*(S_n*(-9.88888888888889*nu + 6.66666666666667) - 4.83333333333333*Sigma_n*delta*nu) +
    S_n*Sigma_l*delta*(-9.66666666666667*nu + 4.0) - 4.5*Sigma_l*Sigma_n*pow(delta, 2)*nu +
    v*(S_n*(nu*(nu*(-1.83950617283951*nu + 97.75) + 37.8775531435105) - 1.5) +
    Sigma_n*delta*(nu*(nu*(-0.907407407407407*nu + 43.25) + 27.6796656329331) - 1.5) +
    v*(S_l*(S_n*(nu*(-40.962962962963*nu - 90.6111111111111) + 4.16666666666667) +
    Sigma_n*delta*(nu*(-17.2777777777778*nu - 36.875) - 2.5)) + S_n*Sigma_l*delta*(nu*(-12.8888888888889*nu - 75.75) +
    2.5) + Sigma_l*Sigma_n*pow(delta, 2)*(nu*(-5.16666666666667*nu - 30.375) - 1.5) +
    v*(S_n*(nu*(nu*(nu*(-0.119341563786008*nu - 92.9104938271605) - 108.905311483896) + 26.7761145568897) - 8.5) +
    Sigma_n*delta*(nu*(nu*(nu*(-0.0555555555555556*nu - 44.0462962962963) - 71.3424151160663) + 14.3906687341338) - 4.5)
    + v*(S_l*(S_n*(nu*(nu*(63.1358024691358*nu + 194.675925925926) + 108.0) - 45.0) +
    Sigma_n*delta*(nu*(nu*(29.5185185185185*nu + 101.180555555556) + 72.75) - 25.0)) +
    S_n*Sigma_l*delta*(nu*(nu*(25.7777777777778*nu + 136.5) + 102.166666666667) - 27.0) + Sigma_l*Sigma_n*pow(delta,
    2)*(nu*(nu*(12.0*nu + 67.4166666666667) + 65.75) - 15.0) + v*(S_n*(nu*(nu*(nu*(nu*(0.0713305898491084*nu +
    36.9351851851852) + 76.6168560337971) - 28.6617828603012) + 10.2773328164665) + 1.5) +
    Sigma_n*delta*(nu*(nu*(nu*(nu*(0.0349794238683128*nu + 18.1358024691358) + 49.6152860607331) - 14.8337253664952) +
    4.65233281646654) + 1.5) + v*(S_l*(S_n*(nu*(nu*(nu*(-34.6666666666667*nu - 149.648148148148) - 136.215277777778) +
    57.875) + 7.5) + Sigma_n*delta*(nu*(nu*(nu*(-17.0*nu - 84.7361111111111) - 91.40625) + 29.75) + 7.5)) +
    S_n*Sigma_l*delta*(nu*(nu*(nu*(-15.4074074074074*nu - 98.0740740740741) - 136.625) + 29.0) + 4.5) +
    Sigma_l*Sigma_n*pow(delta, 2)*(nu*(nu*(nu*(-7.55555555555556*nu - 53.1388888888889) - 88.2708333333333) + 12.125) +
    4.5))))))))))))), Phi(0.0)
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    rfrak_chi1_x = y[1];
    rfrak_chi1_y = y[2];
    rfrak_chi2_x = y[3];
    rfrak_chi2_y = y[4];
    rfrak_ell_x = y[5];
    rfrak_ell_y = y[6];
    rfrak_ell_x = y[7];
    Phi = y[8];

    R = exp(rfrak_ell_x*xHat + rfrak_ell_y*yHat + rfrak_ell_z*zHat);
    nHat = R*xHat*conjugate(R);
    nHat_x = nHat[1];
    nHat_y = nHat[2];
    nHat_z = nHat[3];
    lambdaHat = R*yHat*conjugate(R);
    lambdaHat_x = lambdaHat[1];
    lambdaHat_y = lambdaHat[2];
    lambdaHat_z = lambdaHat[3];
    ellHat = R*zHat*conjugate(R);
    ellHat_x = ellHat[1];
    ellHat_y = ellHat[2];
    ellHat_z = ellHat[3];
    R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat);
    R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat);
    chi1 = chi1Mag*R_S1*zHat*conjugate(R_S1);
    chi1_x = chi1[1];
    chi1_y = chi1[2];
    chi1_z = chi1[3];
    chi2 = chi2Mag*R_S2*zHat*conjugate(R_S2);
    chi2_x = chi2[1];
    chi2_y = chi2[2];
    chi2_z = chi2[3];
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi1_la = chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    chi2_la = chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z;
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_l*pow(m2, 2);
    Sigma_l = -chi1_l*m1 + chi2_l*m2;
    Sigma_n = -chi1_n*m1 + chi2_n*m2;
    logv = log(v);
    Omega_chi1_ellHat = pow(v, 5)*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu - 0.5625) +
      nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75);
    Omega_chi2_ellHat = pow(v, 5)*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu + 0.5625) +
      nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) +
      nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75);
    Omega_ellHat_nHat = pow(v, 6)*(7.0*S_n + 3.0*Sigma_n*delta + pow(v, 2)*(S_n*(-12.0*nu - 3.0) +
      Sigma_n*delta*(-5.5*nu - 3.0) + v*(S_l*(11.6666666666667*S_n + 5.0*Sigma_n*delta) + 7.0*S_n*Sigma_l*delta +
      3.0*Sigma_l*Sigma_n*pow(delta, 2) + v*(S_n*(nu*(9.0*nu - 29.5) - 1.5) + Sigma_n*delta*(nu*(4.33333333333333*nu -
      9.625) - 1.5) + v*(S_l*(S_n*(-9.88888888888889*nu + 6.66666666666667) - 4.83333333333333*Sigma_n*delta*nu) +
      S_n*Sigma_l*delta*(-9.66666666666667*nu + 4.0) - 4.5*Sigma_l*Sigma_n*pow(delta, 2)*nu +
      v*(S_n*(nu*(nu*(-1.83950617283951*nu + 97.75) + 37.8775531435105) - 1.5) +
      Sigma_n*delta*(nu*(nu*(-0.907407407407407*nu + 43.25) + 27.6796656329331) - 1.5) +
      v*(S_l*(S_n*(nu*(-40.962962962963*nu - 90.6111111111111) + 4.16666666666667) +
      Sigma_n*delta*(nu*(-17.2777777777778*nu - 36.875) - 2.5)) + S_n*Sigma_l*delta*(nu*(-12.8888888888889*nu - 75.75) +
      2.5) + Sigma_l*Sigma_n*pow(delta, 2)*(nu*(-5.16666666666667*nu - 30.375) - 1.5) +
      v*(S_n*(nu*(nu*(nu*(-0.119341563786008*nu - 92.9104938271605) - 108.905311483896) + 26.7761145568897) - 8.5) +
      Sigma_n*delta*(nu*(nu*(nu*(-0.0555555555555556*nu - 44.0462962962963) - 71.3424151160663) + 14.3906687341338) -
      4.5) + v*(S_l*(S_n*(nu*(nu*(63.1358024691358*nu + 194.675925925926) + 108.0) - 45.0) +
      Sigma_n*delta*(nu*(nu*(29.5185185185185*nu + 101.180555555556) + 72.75) - 25.0)) +
      S_n*Sigma_l*delta*(nu*(nu*(25.7777777777778*nu + 136.5) + 102.166666666667) - 27.0) + Sigma_l*Sigma_n*pow(delta,
      2)*(nu*(nu*(12.0*nu + 67.4166666666667) + 65.75) - 15.0) + v*(S_n*(nu*(nu*(nu*(nu*(0.0713305898491084*nu +
      36.9351851851852) + 76.6168560337971) - 28.6617828603012) + 10.2773328164665) + 1.5) +
      Sigma_n*delta*(nu*(nu*(nu*(nu*(0.0349794238683128*nu + 18.1358024691358) + 49.6152860607331) - 14.8337253664952) +
      4.65233281646654) + 1.5) + v*(S_l*(S_n*(nu*(nu*(nu*(-34.6666666666667*nu - 149.648148148148) - 136.215277777778) +
      57.875) + 7.5) + Sigma_n*delta*(nu*(nu*(nu*(-17.0*nu - 84.7361111111111) - 91.40625) + 29.75) + 7.5)) +
      S_n*Sigma_l*delta*(nu*(nu*(nu*(-15.4074074074074*nu - 98.0740740740741) - 136.625) + 29.0) + 4.5) +
      Sigma_l*Sigma_n*pow(delta, 2)*(nu*(nu*(nu*(-7.55555555555556*nu - 53.1388888888889) - 88.2708333333333) + 12.125)
      + 4.5))))))))))));
  }

  Quaternion Omega_chi1() {
    return Omega_chi1_ellHat*ellHat;
  }
  Quaternion Omega_chi2() {
    return Omega_chi2_ellHat*ellHat;
  }
  Quaternion Omega_ellHat() {
    return Omega_ellHat_nHat*nHat;
  }

  int TaylorT1_3p5PN(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = pow(v, 10)*(6.4*nu__2 + pow(v, 2)*(-18.6666666666667*nu*nu__2 - 23.752380952381*nu__2 +
      v*(-25.6*S_l*nu__2 - 8.0*Sigma_l*delta*nu__2 + 80.4247719318987*nu__2 + v*(chi1_l*(chi1_l*(6.6*delta*nu__2 -
      13.2*nu*nu__2 + 6.6*nu__2) + 24.8*chi2_l*nu*nu__2) + chi1_la*(chi1_la*(-2.96666666666667*delta*nu__2 +
      5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) - 13.7333333333333*chi2_la*nu*nu__2) +
      chi1_n*(chi1_n*(-2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) -
      13.7333333333333*chi2_n*nu*nu__2) + pow(chi2_l, 2)*(-6.6*delta*nu__2 - 13.2*nu*nu__2 + 6.6*nu__2) + pow(chi2_la,
      2)*(2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) + pow(chi2_n,
      2)*(2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) + 117.726984126984*nu*nu__2
      + nu__2*(23.1111111111111*nu__2 - 31.542151675485) + v*(S_l*(193.422222222222*nu*nu__2 - 28.8*nu__2) +
      Sigma_l*delta*(68.8*nu*nu__2 - 5.2*nu__2) - 488.412937878093*nu*nu__2 - 245.074146910038*nu__2 +
      v*(-321.699087727595*S_l*nu__2 - 103.881997078702*Sigma_l*delta*nu__2 - 104.350476190476*logv*nu__2 -
      56.7811420312465*nu*nu__2 + nu__2*(-199.794708994709*nu__2 - 15.3086419753086*nu__3 + 740.682986723913) +
      v*(S_l*(208.998941798942*nu*nu__2 + nu__2*(-666.074074074074*nu__2 + 448.343327454439)) +
      Sigma_l*delta*(93.9174603174603*nu*nu__2 + nu__2*(-266.844444444444*nu__2 + 181.619047619048)) +
      2498.67153479682*nu*nu__2 + nu__2*(1285.7923710359*nu__2 - 649.661414142347) + v*(S_l*(3875.74795014869*nu*nu__2 -
      729.896693184029*nu__2) + Sigma_l*delta*(1302.34474121815*nu*nu__2 - 214.316458834892*nu__2) +
      337.555736961451*logv*nu__2 - 752.028100625135*nu__2 + v*(-1311.30675759439*logv*nu__2 + 4602.42139029395*nu__2 +
      v*(746.4952102023*logv*nu__2 - 7788.20474442907*nu__2 + v*(3031.19666031508*logv*nu__2 + 6137.18380876523*nu__2 +
      v*(logv*(850.70483446712*logv*nu__2 - 12164.6676812038*nu__2) + 13022.6558856344*nu__2))))))))))));
    const double dEnergydv = v*(-1.0*nu + pow(v, 2)*(nu*(0.166666666666667*nu + 1.5) + v*(-11.6666666666667*S_l*nu -
      5.0*Sigma_l*delta*nu + v*(chi1_l*(chi1_l*(1.5*delta*nu + nu*(-3.0*nu + 1.5)) + 6.0*chi2_l*nu__2) +
      chi1_la*(chi1_la*(-0.75*delta*nu + nu*(1.5*nu - 0.75)) - 3.0*chi2_la*nu__2) + chi1_n*(chi1_n*(-0.75*delta*nu +
      nu*(1.5*nu - 0.75)) - 3.0*chi2_n*nu__2) + pow(chi2_l, 2)*(-1.5*delta*nu + nu*(-3.0*nu + 1.5)) + pow(chi2_la,
      2)*(0.75*delta*nu + nu*(1.5*nu - 0.75)) + pow(chi2_n, 2)*(0.75*delta*nu + nu*(1.5*nu - 0.75)) + nu*(nu*(0.125*nu -
      7.125) + 10.125) + v*(S_l*nu*(23.7222222222222*nu - 38.5) + Sigma_l*delta*nu*(11.6666666666667*nu - 10.5) +
      v*(nu*(nu*(nu*(0.0270061728395062*nu + 6.45833333333333) - 154.898517962917) + 42.1875) +
      v*(S_l*nu*(nu*(-10.875*nu + 412.875) - 151.875) + Sigma_l*delta*nu*(nu*(-5.625*nu + 175.5) - 30.375) +
      v*(-298.666666666667*logv*nu__2 + nu*(nu*(nu*(nu*(-0.012377829218107*nu - 0.870949074074074) + 450.663995131026) -
      769.4015) + 155.0390625) + pow(v, 2)*(logv*(1523.2*nu*nu__2 + 1710.17142857143*nu__2) + nu*(330.78*nu +
      538.20703125) + pow(v, 2)*(16016.0*logv*nu__2 + nu*(-4116.0*nu + 1808.9736328125)))))))))));
    const double Absorption = pow(v, 15)*(chi1_l*pow(m1, 3)*(-4.8*chi1chi1*nu__2 - 1.6*nu__2) + chi2_l*pow(m2,
      3)*(-4.8*chi2chi2*nu__2 - 1.6*nu__2) + pow(v, 3)*(pow(m1, 4)*(chi1chi1*nu__2*(9.6*sqrt1Mchi1chi1 + 9.6) +
      nu__2*(3.2*sqrt1Mchi1chi1 + 3.2)) + pow(m2, 4)*(chi2chi2*nu__2*(9.6*sqrt1Mchi2chi2 + 9.6) +
      nu__2*(3.2*sqrt1Mchi2chi2 + 3.2))));
    const double dvdt_T1 = (-Absorption - Flux)/dEnergydv;
    if(dvdt_T1<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T1, t, y, dydt);
  }

  int TaylorT4_3p5PN(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt_T4 = pow(v, 9)*(6.4*nu + pow(v, 2)*(nu*(-17.6*nu - 14.152380952381) + v*(-100.266666666667*S_l*nu
      - 40.0*Sigma_l*delta*nu + 80.4247719318987*nu + v*(chi1_l*(chi1_l*(16.2*delta*nu + nu*(-32.4*nu + 16.2)) +
      63.2*chi2_l*pow(nu, 2)) + chi1_la*(chi1_la*(-7.76666666666667*delta*nu + nu*(15.5333333333333*nu -
      7.76666666666667)) - 32.9333333333333*chi2_la*pow(nu, 2)) + chi1_n*(chi1_n*(-7.76666666666667*delta*nu +
      nu*(15.5333333333333*nu - 7.76666666666667)) - 32.9333333333333*chi2_n*pow(nu, 2)) + pow(chi2_l,
      2)*(-16.2*delta*nu + nu*(-32.4*nu + 16.2)) + pow(chi2_la, 2)*(7.76666666666667*delta*nu + nu*(15.5333333333333*nu
      - 7.76666666666667)) + pow(chi2_n, 2)*(7.76666666666667*delta*nu + nu*(15.5333333333333*nu - 7.76666666666667)) +
      nu*(nu*(20.9777777777778*nu + 43.368253968254) + 12.0292768959436) + v*(S_l*nu*(533.866666666667*nu -
      260.488888888889) + Sigma_l*delta*nu*(224.8*nu - 61.6380952380952) + chi1_l*pow(m1, 3)*(-4.8*chi1chi1*nu - 1.6*nu)
      + chi2_l*pow(m2, 3)*(-4.8*chi2chi2*nu - 1.6*nu) + nu*(-475.008809222777*nu - 124.43698901219) +
      v*(S_l*(1169.77777777778*S_l*nu + 968.0*Sigma_l*delta*nu - 1259.98809359975*nu) +
      Sigma_l*(200.0*Sigma_l*pow(delta, 2)*nu - 506.005856738196*delta*nu) + chi1_l*(chi1_l*(delta*nu*(-23.7*nu +
      3.07142857142857) + nu*(nu*(47.4*nu - 29.8428571428571) + 3.07142857142857)) + chi2_l*pow(nu,
      2)*(-95.0666666666667*nu + 9.88571428571429)) + chi1_la*(chi1_la*(delta*nu*(11.9055555555556*nu -
      1.03571428571429) + nu*(nu*(-23.8111111111111*nu + 13.9769841269841) - 1.03571428571429)) + chi2_la*pow(nu,
      2)*(47.3111111111111*nu - 6.94285714285714)) + chi1_n*(chi1_n*(delta*nu*(11.9055555555556*nu - 1.03571428571429) +
      nu*(nu*(-23.8111111111111*nu + 13.9769841269841) - 1.03571428571429)) + chi2_n*pow(nu, 2)*(47.3111111111111*nu -
      6.94285714285714)) + pow(chi2_l, 2)*(delta*nu*(23.7*nu - 3.07142857142857) + nu*(nu*(47.4*nu - 29.8428571428571) +
      3.07142857142857)) + pow(chi2_la, 2)*(delta*nu*(-11.9055555555556*nu + 1.03571428571429) +
      nu*(nu*(-23.8111111111111*nu + 13.9769841269841) - 1.03571428571429)) + pow(chi2_n,
      2)*(delta*nu*(-11.9055555555556*nu + 1.03571428571429) + nu*(nu*(-23.8111111111111*nu + 13.9769841269841) -
      1.03571428571429)) - 104.350476190476*logv*nu + nu*(nu*(nu*(-13.8395061728395*nu + 3.8642857142857) -
      1058.43868227316) + 885.434044924971) + v*(S_l*(chi1_l*(chi1_l*(-339.4*delta*nu + nu*(678.8*nu - 339.4)) -
      1338.93333333333*chi2_l*pow(nu, 2)) + chi1_la*(chi1_la*(165.811111111111*delta*nu + nu*(-331.622222222222*nu +
      165.811111111111)) + 685.022222222222*chi2_la*pow(nu, 2)) + chi1_n*(chi1_n*(165.811111111111*delta*nu +
      nu*(-331.622222222222*nu + 165.811111111111)) + 685.022222222222*chi2_n*pow(nu, 2)) + pow(chi2_l,
      2)*(339.4*delta*nu + nu*(678.8*nu - 339.4)) + pow(chi2_la, 2)*(-165.811111111111*delta*nu +
      nu*(-331.622222222222*nu + 165.811111111111)) + pow(chi2_n, 2)*(-165.811111111111*delta*nu +
      nu*(-331.622222222222*nu + 165.811111111111)) + nu*(nu*(-1321.48148148148*nu + 4159.09523809524) -
      1525.06490299824)) + Sigma_l*(chi1_l*(chi1_l*delta*(-141.0*delta*nu + nu*(282.0*nu - 141.0)) -
      556.0*chi2_l*delta*pow(nu, 2)) + chi1_la*(chi1_la*delta*(68.8333333333333*delta*nu + nu*(-137.666666666667*nu +
      68.8333333333333)) + 284.666666666667*chi2_la*delta*pow(nu, 2)) + chi1_n*(chi1_n*delta*(68.8333333333333*delta*nu
      + nu*(-137.666666666667*nu + 68.8333333333333)) + 284.666666666667*chi2_n*delta*pow(nu, 2)) + pow(chi2_l,
      2)*delta*(141.0*delta*nu + nu*(282.0*nu - 141.0)) + pow(chi2_la, 2)*delta*(-68.8333333333333*delta*nu +
      nu*(-137.666666666667*nu + 68.8333333333333)) + pow(chi2_n, 2)*delta*(-68.8333333333333*delta*nu +
      nu*(-137.666666666667*nu + 68.8333333333333)) + delta*nu*(nu*(-580.6*nu + 1631.89206349206) - 421.784479717813)) +
      chi1_l*pow(m1, 3)*(chi1chi1*nu*(-0.8*nu - 7.2) + nu*(-0.266666666666667*nu - 2.4)) +
      chi1_l*(chi1_l*(120.637157897848*delta*nu + nu*(-241.274315795696*nu + 120.637157897848)) +
      482.548631591392*chi2_l*pow(nu, 2)) + chi1_la*(chi1_la*(-60.318578948924*delta*nu + nu*(120.637157897848*nu -
      60.318578948924)) - 241.274315795696*chi2_la*pow(nu, 2)) + chi1_n*(chi1_n*(-60.318578948924*delta*nu +
      nu*(120.637157897848*nu - 60.318578948924)) - 241.274315795696*chi2_n*pow(nu, 2)) + pow(chi2_l,
      2)*(-120.637157897848*delta*nu + nu*(-241.274315795696*nu + 120.637157897848)) + chi2_l*pow(m2,
      3)*(chi2chi2*nu*(-0.8*nu - 7.2) + nu*(-0.266666666666667*nu - 2.4)) + pow(chi2_la, 2)*(60.318578948924*delta*nu +
      nu*(120.637157897848*nu - 60.318578948924)) + pow(chi2_n, 2)*(60.318578948924*delta*nu + nu*(120.637157897848*nu -
      60.318578948924)) + nu*(nu*(1216.67733265692*nu + 1192.39232277917) - 22.0160818501571))))))));
    if(dvdt_T4<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T4, t, y, dydt);
  }

  int TaylorT5_3p5PN(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = pow(m1, 3)*(chi1_l*(0.6640625*chi1chi1 + 0.221354166666667)/pow(v, 2) +
      (chi1_l*(0.694056919642857*chi1chi1 + 0.231352306547619) + chi1_l*(0.1171875*chi1chi1 + 0.0390625)/pow(v,
      2))/(nu*pow(v, 2))) + pow(m2, 3)*(chi2_l*(0.6640625*chi2chi2 + 0.221354166666667)/pow(v, 2) +
      (chi2_l*(0.694056919642857*chi2chi2 + 0.231352306547619) + chi2_l*(0.1171875*chi2chi2 + 0.0390625)/pow(v,
      2))/(nu*pow(v, 2))) + (S_l*(chi1_l*(8.212890625*chi1_l - 15.6575520833333*chi2_l) +
      chi1_la*(-3.78634982638889*chi1_la + 8.46896701388889*chi2_la) + chi1_n*(-3.78634982638889*chi1_n +
      8.46896701388889*chi2_n) + 8.212890625*pow(chi2_l, 2) - 3.78634982638889*pow(chi2_la, 2) -
      3.78634982638889*pow(chi2_n, 2) + 0.0661892361111116*nu - 68.065837880291) +
      Sigma_l*(chi1_l*(3.0029296875*chi1_l*delta - 5.712890625*chi2_l*delta) + chi1_la*(-1.37939453125*chi1_la*delta +
      3.1005859375*chi2_la*delta) + chi1_n*(-1.37939453125*chi1_n*delta + 3.1005859375*chi2_n*delta) +
      3.0029296875*pow(chi2_l, 2)*delta - 1.37939453125*pow(chi2_la, 2)*delta - 1.37939453125*pow(chi2_n, 2)*delta +
      delta*(-0.256890190972224*nu - 33.4405808221726)) + chi1_l*(-13.989904785517*chi1_l + 26.9980618667873*chi2_l) +
      chi1_la*(6.58589084932235*chi1_la - 14.317154020266*chi2_la) + chi1_n*(6.58589084932235*chi1_n -
      14.317154020266*chi2_n) - 13.989904785517*pow(chi2_l, 2) + 6.58589084932235*pow(chi2_la, 2) +
      6.58589084932235*pow(chi2_n, 2) + 2.40388587172471*nu - 6.14428294102153 + (chi1_l*(chi1_l*(-1.5966796875*delta +
      3.193359375*nu + 2.05165318080357) + chi2_l*(-6.16536458333333*nu - 7.06531343005952)) +
      chi1_la*(chi1_la*(0.752224392361111*delta - 1.50444878472222*nu - 0.975545247395833) +
      chi2_la*(3.26714409722222*nu + 3.72545030381944)) + chi1_n*(chi1_n*(0.752224392361111*delta - 1.50444878472222*nu
      - 0.975545247395833) + chi2_n*(3.26714409722222*nu + 3.72545030381944)) + pow(chi2_l, 2)*(1.5966796875*delta +
      3.193359375*nu + 2.05165318080357) + pow(chi2_la, 2)*(-0.752224392361111*delta - 1.50444878472222*nu -
      0.975545247395833) + pow(chi2_n, 2)*(-0.752224392361111*delta - 1.50444878472222*nu - 0.975545247395833) +
      nu*(0.770550009645062*nu - 0.343853985821765) + 25.8462542787179 + (0.429687499999999*S_l -
      0.1171875*Sigma_l*delta + 0.797670009700533 + (chi1_l*(0.791015625*chi1_l - 1.54296875*chi2_l) +
      chi1_la*(-0.379231770833333*chi1_la + 0.804036458333333*chi2_la) + chi1_n*(-0.379231770833333*chi1_n +
      0.804036458333333*chi2_n) + 0.791015625*pow(chi2_l, 2) - 0.379231770833333*pow(chi2_la, 2) -
      0.379231770833333*pow(chi2_n, 2) + 0.669487847222222*nu + 0.841548859126984 + 0.4296875/pow(v, 2))/v)/v)/v)/pow(v,
      2) + (S_l*(pow(chi1_l, 2)*(-4.1064453125*delta - 4.1064453125) + pow(chi1_la, 2)*(1.89317491319444*delta +
      1.89317491319444) + pow(chi1_n, 2)*(1.89317491319444*delta + 1.89317491319444) + pow(chi2_l,
      2)*(4.1064453125*delta - 4.1064453125) + pow(chi2_la, 2)*(-1.89317491319444*delta + 1.89317491319444) +
      pow(chi2_n, 2)*(-1.89317491319444*delta + 1.89317491319444) + 92.0671273639423) + Sigma_l*(pow(chi1_l,
      2)*delta*(-1.50146484375*delta - 1.50146484375) + pow(chi1_la, 2)*delta*(0.689697265625*delta + 0.689697265625) +
      pow(chi1_n, 2)*delta*(0.689697265625*delta + 0.689697265625) + pow(chi2_l, 2)*delta*(1.50146484375*delta -
      1.50146484375) + pow(chi2_la, 2)*delta*(-0.689697265625*delta + 0.689697265625) + pow(chi2_n,
      2)*delta*(-0.689697265625*delta + 0.689697265625) + 27.6075954538174*delta) + pow(chi1_l,
      2)*(6.99495239275852*delta + 6.99495239275852) + pow(chi1_la, 2)*(-3.29294542466118*delta - 3.29294542466118) +
      pow(chi1_n, 2)*(-3.29294542466118*delta - 3.29294542466118) + pow(chi2_l, 2)*(-6.99495239275852*delta +
      6.99495239275852) + pow(chi2_la, 2)*(3.29294542466118*delta - 3.29294542466118) + pow(chi2_n,
      2)*(3.29294542466118*delta - 3.29294542466118) - 7.44928308515137 + (S_l*(9.79166666666667*S_l +
      6.96614583333333*Sigma_l*delta - 30.7614280664001) + Sigma_l*(1.220703125*Sigma_l*pow(delta, 2) -
      12.1900339943979*delta) + pow(chi1_l, 2)*(-1.82416643415179*delta - 1.82416643415179) + pow(chi1_la,
      2)*(0.863884819878472*delta + 0.863884819878472) + pow(chi1_n, 2)*(0.863884819878472*delta + 0.863884819878472) +
      pow(chi2_l, 2)*(1.82416643415179*delta - 1.82416643415179) + pow(chi2_la, 2)*(-0.863884819878472*delta +
      0.863884819878472) + pow(chi2_n, 2)*(-0.863884819878472*delta + 0.863884819878472) + 2.54761904761905*logv +
      3.44765918029945 + (17.185794890873*S_l + 5.82380022321429*Sigma_l*delta - 5.64577976646101 + (pow(chi1_l,
      2)*(-0.3955078125*delta - 0.3955078125) + pow(chi1_la, 2)*(0.189615885416667*delta + 0.189615885416667) +
      pow(chi1_n, 2)*(0.189615885416667*delta + 0.189615885416667) + pow(chi2_l, 2)*(0.3955078125*delta - 0.3955078125)
      + pow(chi2_la, 2)*(-0.189615885416667*delta + 0.189615885416667) + pow(chi2_n, 2)*(-0.189615885416667*delta +
      0.189615885416667) + 0.4703617648593 + (2.44791666666667*S_l + 0.9765625*Sigma_l*delta - 1.96349540849362 +
      (0.345517113095238 + 0.15625/pow(v, 2))/v)/v)/v)/v)/v)/(nu*pow(v, 2));
    const double dvdt_T5 = 1.0/dtdv;
    if(dvdt_T5<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt_T5, t, y, dydt);
  }

  int CommonRHS(const double dvdt, double t, const double* y, double* dydt) {
    std::vector<double> rfrak_ellHat(3);
    rfrak_ellHat[0] = y[5];
    rfrak_ellHat[1] = y[6];
    rfrak_ellHat[2] = y[7];
    const std::vector<double> rfrakdot_ell = FrameFromAngularVelocity_Integrand(rfrak_ellHat, Omega_ellHat().vec());
    dydt[0] = dvdt;
    FrameFromAngularVelocity_2D_Integrand(y[1], y[2], Omega_chi1().vec(), dydt[1], dydt[2]);
    FrameFromAngularVelocity_2D_Integrand(y[3], y[4], Omega_chi2().vec(), dydt[3], dydt[4]);
    dydt[5] = rfrakdot_ell[0];
    dydt[6] = rfrakdot_ell[1];
    dydt[7] = rfrakdot_ell[2];
    dydt[8] = v*v*v;

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
};
