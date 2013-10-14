// File produced automatically by OrbitalEvolution.ipynb

void AngularVelocityIntegrand(const double r0, const double r1, const double r2,
                              const double omega0, const double omega1, const double omega2,
                              double& rdot0, double& rdot1, double& rdot2)
{
  const double rmag = sqrt(r0*r0+r1*r1+r2*r2);

  // If rmag is basically zero, return an answer assuming exactly zero
  if(rmag<Quaternion_Epsilon) {
    rdot0 = omega0/2.0;
    rdot1 = omega1/2.0;
    rdot2 = omega2/2.0;
    return;
  }

  // Otherwise, actually do the calculation
  const double dotTerm = (r0*omega0+r1*omega1+r2*omega2)/(rmag*rmag);
  const double cotTerm = rmag/(2*tan(rmag));
  rdot0 = (omega0 - r0*dotTerm)*cotTerm + r0*dotTerm + 0.5*(omega1*r2 - omega2*r1);
  rdot1 = (omega1 - r1*dotTerm)*cotTerm + r1*dotTerm + 0.5*(omega2*r0 - omega0*r2);
  rdot2 = (omega2 - r2*dotTerm)*cotTerm + r2*dotTerm + 0.5*(omega0*r1 - omega1*r0);
  return;
}

void AngularVelocityIntegrand(const double r0, const double r1,
                              const double omega0, const double omega1, double omega2,
                              double& rdot0, double& rdot1)
{
  const double rmag = sqrt(r0*r0+r1*r1);

  // If rmag is basically zero, return an answer assuming exactly zero
  if(rmag<Quaternion_Epsilon) {
    rdot0 = omega0/2.0;
    rdot1 = omega1/2.0;
    return;
  }

  // Otherwise, actually do the calculation
  const double dotTerm = (r0*omega0+r1*omega1)/(rmag*rmag);
  const double cotTerm = rmag/(2*tan(rmag));
  rdot0 = (omega0 - r0*dotTerm)*cotTerm + r0*dotTerm - 0.5*omega2*r1;
  rdot1 = (omega1 - r1*dotTerm)*cotTerm + r1*dotTerm + 0.5*omega2*r0;
  return;
}

const Quaternions::Quaternion xHat(0,1.0,0,0);
const Quaternions::Quaternion yHat(0,0,1.0,0);
const Quaternions::Quaternion zHat(0,0,0,1.0);

class TaylorTn {
private:
  const double m1, chi1Mag, chi2Mag, m2, delta, nu, nu__2, nu__3, chi1chi1, chi2chi2, sqrt1Mchi1chi1, sqrt1Mchi2chi2;
  const Quaternion xHat;
  const Quaternion yHat;
  const Quaternion zHat;
  double v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_ell_x, rfrak_ell_y, rfrak_ell_z, nHat_x,
         nHat_y, nHat_z, lambdaHat_x, lambdaHat_y, lambdaHat_z, ellHat_x, ellHat_y, ellHat_z, chi1_x, chi1_y, chi1_z,
         chi2_x, chi2_y, chi2_z, chi1_l, chi1_n, chi1_la, chi2_l, chi2_n, chi2_la, S_l, Sigma_l, logv;
  Quaternion R;
  Quaternion nHat;
  Quaternion lambdaHat;
  Quaternion ellHat;
  Quaternion R_S1;
  Quaternion R_S2;
  Quaternion chi1;
  Quaternion chi2;
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
    sqrt1Mchi2chi2(sqrt(-chi2chi2 + 1.0)), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), Sigma_l(-chi1_l*m1 + chi2_l*m2),
    logv(log(v)), Phi(0.0)
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

    m2 = -1.0*m1 + 1.0;
    delta = m1 - m2;
    nu = m1*m2;
    nu__2 = pow(m1, 2)*pow(m2, 2);
    nu__3 = pow(m1, 3)*pow(m2, 3);
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
    chi1chi1 = pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z, 2);
    chi2chi2 = pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z, 2);
    chi1_l = chi1_x*ellHat_x + chi1_y*ellHat_y + chi1_z*ellHat_z;
    chi1_n = chi1_x*nHat_x + chi1_y*nHat_y + chi1_z*nHat_z;
    chi1_la = chi1_x*lambdaHat_x + chi1_y*lambdaHat_y + chi1_z*lambdaHat_z;
    chi2_l = chi2_x*ellHat_x + chi2_y*ellHat_y + chi2_z*ellHat_z;
    chi2_n = chi2_x*nHat_x + chi2_y*nHat_y + chi2_z*nHat_z;
    chi2_la = chi2_x*lambdaHat_x + chi2_y*lambdaHat_y + chi2_z*lambdaHat_z;
    sqrt1Mchi1chi1 = sqrt(-chi1chi1 + 1.0);
    sqrt1Mchi2chi2 = sqrt(-chi2chi2 + 1.0);
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    Sigma_l = -chi1_l*m1 + chi2_l*m2;
    logv = log(v);
  }

  int TaylorT1RHS(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double Flux = pow(v, 10)*(6.4*nu__2 + pow(v, 2)*(-18.6666666666667*nu*nu__2 - 23.752380952381*nu__2 +
    v*(-25.6*S_l*nu__2 - 8.0*Sigma_l*delta*nu__2 + 80.4247719318987*nu__2 + v*(chi1_l*(chi1_l*(6.6*delta*nu__2 -
    13.2*nu*nu__2 + 6.6*nu__2) + 24.8*chi2_l*nu*nu__2) + chi1_la*(chi1_la*(-2.96666666666667*delta*nu__2 +
    5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) - 13.7333333333333*chi2_la*nu*nu__2) +
    chi1_n*(chi1_n*(-2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) -
    13.7333333333333*chi2_n*nu*nu__2) + pow(chi2_l, 2)*(-6.6*delta*nu__2 - 13.2*nu*nu__2 + 6.6*nu__2) + pow(chi2_la,
    2)*(2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) + pow(chi2_n,
    2)*(2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) + 117.726984126984*nu*nu__2 +
    nu__2*(23.1111111111111*nu__2 - 31.542151675485) + v*(S_l*(193.422222222222*nu*nu__2 - 28.8*nu__2) +
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
    v*(nu*(nu*(nu*(0.0270061728395062*nu + 6.45833333333333) - 154.898517962917) + 42.1875) + v*(S_l*nu*(nu*(-10.875*nu
    + 412.875) - 151.875) + Sigma_l*delta*nu*(nu*(-5.625*nu + 175.5) - 30.375) + v*(-298.666666666667*logv*nu__2 +
    nu*(nu*(nu*(nu*(-0.012377829218107*nu - 0.870949074074074) + 450.663995131026) - 769.4015) + 155.0390625) + pow(v,
    2)*(logv*(1523.2*nu*nu__2 + 1710.17142857143*nu__2) + nu*(330.78*nu + 538.20703125) + pow(v, 2)*(16016.0*logv*nu__2
    + nu*(-4116.0*nu + 1808.9736328125)))))))))));
    const double Absorption = pow(v, 15)*(chi1_l*pow(m1, 3)*(-4.8*chi1chi1*nu__2 - 1.6*nu__2) + chi2_l*pow(m2,
    3)*(-4.8*chi2chi2*nu__2 - 1.6*nu__2) + pow(v, 3)*(pow(m1, 4)*(chi1chi1*nu__2*(9.6*sqrt1Mchi1chi1 + 9.6) +
    nu__2*(3.2*sqrt1Mchi1chi1 + 3.2)) + pow(m2, 4)*(chi2chi2*nu__2*(9.6*sqrt1Mchi2chi2 + 9.6) +
    nu__2*(3.2*sqrt1Mchi2chi2 + 3.2))));
    const double dvdt = (-Absorption - Flux)/dEnergydv;
    if(dvdt<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt, t, y, dydt);
  }

  int TaylorT4RHS(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt = (-Absorption - Flux)/dEnergydv;
    if(dvdt<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt, t, y, dydt);
  }

  int TaylorT5RHS(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = chi1_l*pow(m1, 3)*(0.1171875*chi1chi1 + 0.0390625)/(nu*pow(v, 4)) + chi2_l*pow(m2,
    3)*(0.1171875*chi2chi2 + 0.0390625)/(nu*pow(v, 4)) + (chi1_l*(chi1_l*(-1.5966796875*delta + 3.193359375*nu +
    2.05165318080357) + chi2_l*(-6.16536458333333*nu - 7.06531343005952)) + chi1_la*(chi1_la*(0.752224392361111*delta -
    1.50444878472222*nu - 0.975545247395833) + chi2_la*(3.26714409722222*nu + 3.72545030381944)) +
    chi1_n*(chi1_n*(0.752224392361111*delta - 1.50444878472222*nu - 0.975545247395833) + chi2_n*(3.26714409722222*nu +
    3.72545030381944)) + pow(chi2_l, 2)*(1.5966796875*delta + 3.193359375*nu + 2.05165318080357) + pow(chi2_la,
    2)*(-0.752224392361111*delta - 1.50444878472222*nu - 0.975545247395833) + pow(chi2_n, 2)*(-0.752224392361111*delta -
    1.50444878472222*nu - 0.975545247395833) + nu*(0.770550009645062*nu - 0.343853985821765) + 25.8462542787179 +
    (0.429687499999999*S_l - 0.1171875*Sigma_l*delta + 0.797670009700533 + (chi1_l*(0.791015625*chi1_l -
    1.54296875*chi2_l) + chi1_la*(-0.379231770833333*chi1_la + 0.804036458333333*chi2_la) +
    chi1_n*(-0.379231770833333*chi1_n + 0.804036458333333*chi2_n) + 0.791015625*pow(chi2_l, 2) -
    0.379231770833333*pow(chi2_la, 2) - 0.379231770833333*pow(chi2_n, 2) + 0.669487847222222*nu + 0.841548859126984 +
    0.4296875/pow(v, 2))/v)/v)/pow(v, 3) + (S_l*(9.79166666666667*S_l + 6.96614583333333*Sigma_l*delta -
    30.7614280664001) + Sigma_l*(1.220703125*Sigma_l*pow(delta, 2) - 12.1900339943979*delta) + pow(chi1_l,
    2)*(-1.82416643415179*delta - 1.82416643415179) + pow(chi1_la, 2)*(0.863884819878472*delta + 0.863884819878472) +
    pow(chi1_n, 2)*(0.863884819878472*delta + 0.863884819878472) + pow(chi2_l, 2)*(1.82416643415179*delta -
    1.82416643415179) + pow(chi2_la, 2)*(-0.863884819878472*delta + 0.863884819878472) + pow(chi2_n,
    2)*(-0.863884819878472*delta + 0.863884819878472) + 2.54761904761905*logv + 3.44765918029945 + (17.185794890873*S_l
    + 5.82380022321429*Sigma_l*delta - 5.64577976646101 + (pow(chi1_l, 2)*(-0.3955078125*delta - 0.3955078125) +
    pow(chi1_la, 2)*(0.189615885416667*delta + 0.189615885416667) + pow(chi1_n, 2)*(0.189615885416667*delta +
    0.189615885416667) + pow(chi2_l, 2)*(0.3955078125*delta - 0.3955078125) + pow(chi2_la, 2)*(-0.189615885416667*delta
    + 0.189615885416667) + pow(chi2_n, 2)*(-0.189615885416667*delta + 0.189615885416667) + 0.4703617648593 +
    (2.44791666666667*S_l + 0.9765625*Sigma_l*delta - 1.96349540849362 + (0.345517113095238 + 0.15625/pow(v,
    2))/v)/v)/v)/v)/(nu*pow(v, 3));
    const double dvdt = (-Absorption - Flux)/dEnergydv;
    if(dvdt<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt, t, y, dydt);
  }

  int CommonRHS(const double dvdt, double t, const double* y, double* dydt) {
    dydt[0] = dvdt;
    AngularVelocityIntegrand(rfrak_chi1_x, rfrak_chi1_y,
                             Omega1Mag*ellHat_x, Omega1Mag*ellHat_y, Omega1Mag*ellHat_z,
                             dydt[1], dydt[2]);
    AngularVelocityIntegrand(rfrak_chi2_x, rfrak_chi2_y,
                             Omega2Mag*ellHat_x, Omega2Mag*ellHat_y, Omega2Mag*ellHat_z,
                             dydt[3], dydt[4]);
    AngularVelocityIntegrand(rfrak_ell_x, rfrak_ell_y, rfrak_ell_z,
                             OmegaLMag*nHat_x,  OmegaLMag*nHat_y,  OmegaLMag*nHat_z,
                             dydt[5], dydt[6], dydt[7]);
    dydt[8] = v*v*v; // Phi

    return GSL_SUCCESS; // GSL expects this if everything went well
  }
};
