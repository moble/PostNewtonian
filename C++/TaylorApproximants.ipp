
// File produced automatically by OrbitalEvolution.ipynb

inline void cross(const double* a, const double* b, double* c) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
  return;
}

const Quaternions::Quaternion xHat(0,1.0,0,0);
const Quaternions::Quaternion yHat(0,0,1.0,0);
const Quaternions::Quaternion zHat(0,0,0,1.0);

class TaylorTn {
private:
  const double m1;
  double v, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, Lhat_Nx, Lhat_Ny, Lhat_Nz, nhat_x, nhat_y, nhat_z;
  const double m2, delta, nu, nu__2, nu__3, chi1chi1, chi2chi2, sqrt1Mchi1chi1, sqrt1Mchi2chi2;
  double chi1_l, chi1_n, chi1_la, chi2_l, chi2_n, chi2_la, S_l, S_n, Sigma_l, Sigma_n, logv, gamma_L, a_L_overvcubed,
         OmegaLMag, Omega1Mag, Omega2Mag;
  double Phi, gamma;

public:
  TaylorTn(const double m1_in, double v_0, double chi1_x_0, double chi1_y_0, double chi1_z_0, double chi2_x_0, double chi2_y_0,
            double chi2_z_0, double Lhat_Nx_0, double Lhat_Ny_0, double Lhat_Nz_0, double nhat_x_0, double nhat_y_0,
            double nhat_z_0) :
    m1(m1_in), v(v_0), chi1_x(chi1_x_0), chi1_y(chi1_y_0), chi1_z(chi1_z_0), chi2_x(chi2_x_0), chi2_y(chi2_y_0),
    chi2_z(chi2_z_0), Lhat_Nx(Lhat_Nx_0), Lhat_Ny(Lhat_Ny_0), Lhat_Nz(Lhat_Nz_0), nhat_x(nhat_x_0), nhat_y(nhat_y_0),
    nhat_z(nhat_z_0), m2(-m1 + 1.0), delta(m1 - m2), nu(m1*m2), nu__2(pow(m1, 2)*pow(m2, 2)), nu__3(pow(m1, 3)*pow(m2,
    3)), chi1chi1(pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z, 2)), chi2chi2(pow(chi2_x, 2) + pow(chi2_y, 2) +
    pow(chi2_z, 2)), sqrt1Mchi1chi1(sqrt(-chi1chi1 + 1.0)), sqrt1Mchi2chi2(sqrt(-chi2chi2 + 1.0)), chi1_l(Lhat_Nx*chi1_x
    + Lhat_Ny*chi1_y + Lhat_Nz*chi1_z), chi1_n(chi1_x*nhat_x + chi1_y*nhat_y + chi1_z*nhat_z),
    chi1_la(Lhat_Nx*(-chi1_y*nhat_z + chi1_z*nhat_y) + Lhat_Ny*(chi1_x*nhat_z - chi1_z*nhat_x) + Lhat_Nz*(-chi1_x*nhat_y
    + chi1_y*nhat_x)), chi2_l(Lhat_Nx*chi2_x + Lhat_Ny*chi2_y + Lhat_Nz*chi2_z), chi2_n(chi2_x*nhat_x + chi2_y*nhat_y +
    chi2_z*nhat_z), chi2_la(Lhat_Nx*(-chi2_y*nhat_z + chi2_z*nhat_y) + Lhat_Ny*(chi2_x*nhat_z - chi2_z*nhat_x) +
    Lhat_Nz*(-chi2_x*nhat_y + chi2_y*nhat_x)), S_l(chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2)), S_n(chi1_n*pow(m1, 2) +
    chi2_l*pow(m2, 2)), Sigma_l(-chi1_l*m1 + chi2_l*m2), Sigma_n(-chi1_n*m1 + chi2_n*m2), logv(log(v)), gamma_L(pow(v,
    2)*(pow(v, 2)*(-0.333333333333333*nu + v*(1.66666666666667*S_l + Sigma_l*delta + v*(-5.41666666666667*nu +
    v*(S_l*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_l*delta + v*(nu*(nu*(0.0123456790123457*nu +
    6.36111111111111) - 2.98177812235564) + v*(S_l*(nu*(-6.0*nu - 10.5833333333333) + 5.0) +
    Sigma_l*delta*(nu*(-2.66666666666667*nu - 10.1666666666667) + 3.0)) + 1.0)) + 1.0)) + 1.0) + 1.0)),
    a_L_overvcubed(pow(v, 4)*(7.0*S_n + 3.0*Sigma_n*delta + pow(v, 2)*(S_n*(-9.66666666666667*nu - 10.0) +
    Sigma_n*delta*(-4.5*nu - 6.0) + pow(v, 2)*(S_n*(nu*(5.77777777777778*nu + 14.75) + 1.5) +
    Sigma_n*delta*(nu*(2.83333333333333*nu + 9.125) + 1.5))))), OmegaLMag(a_L_overvcubed*gamma_L), Omega1Mag(pow(v,
    5)*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu - 0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v,
    2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) +
    0.5625) + 0.75)), Omega2Mag(pow(v, 5)*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu + 0.5625) +
    nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) +
    nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75))
  { }

  void Recalculate(double t, const double* y) {
    v = y[0];
    Phi = y[1];
    chi1_x = y[2];
    chi1_y = y[3];
    chi1_z = y[4];
    chi2_x = y[5];
    chi2_y = y[6];
    chi2_z = y[7];
    Lhat_Nx = y[8];
    Lhat_Ny = y[9];
    Lhat_Nz = y[10];
    gamma = y[11];
    const Quaternions::Quaternion Rax =
        Quaternions::sqrtOfRotor(-Quaternions::normalized(Quaternions::Quaternion(0., Lhat_Nx, Lhat_Ny, Lhat_Nz))*zHat);
    const Quaternions::Quaternion R = Rax * Quaternions::exp(((gamma+Phi)/2.)*zHat);
    const Quaternions::Quaternion nhatQ = R*xHat*R.conjugate();
    nhat_x = nhatQ[1];
    nhat_y = nhatQ[2];
    nhat_z = nhatQ[3];
    chi1_l = Lhat_Nx*chi1_x + Lhat_Ny*chi1_y + Lhat_Nz*chi1_z;
    chi1_n = chi1_x*nhat_x + chi1_y*nhat_y + chi1_z*nhat_z;
    chi1_la = Lhat_Nx*(-chi1_y*nhat_z + chi1_z*nhat_y) + Lhat_Ny*(chi1_x*nhat_z - chi1_z*nhat_x) +
        Lhat_Nz*(-chi1_x*nhat_y + chi1_y*nhat_x);
    chi2_l = Lhat_Nx*chi2_x + Lhat_Ny*chi2_y + Lhat_Nz*chi2_z;
    chi2_n = chi2_x*nhat_x + chi2_y*nhat_y + chi2_z*nhat_z;
    chi2_la = Lhat_Nx*(-chi2_y*nhat_z + chi2_z*nhat_y) + Lhat_Ny*(chi2_x*nhat_z - chi2_z*nhat_x) +
        Lhat_Nz*(-chi2_x*nhat_y + chi2_y*nhat_x);
    S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    S_n = chi1_n*pow(m1, 2) + chi2_l*pow(m2, 2);
    Sigma_l = -chi1_l*m1 + chi2_l*m2;
    Sigma_n = -chi1_n*m1 + chi2_n*m2;
    logv = log(v);
    gamma_L = pow(v, 2)*(pow(v, 2)*(-0.333333333333333*nu + v*(1.66666666666667*S_l + Sigma_l*delta +
        v*(-5.41666666666667*nu + v*(S_l*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_l*delta +
        v*(nu*(nu*(0.0123456790123457*nu + 6.36111111111111) - 2.98177812235564) + v*(S_l*(nu*(-6.0*nu -
        10.5833333333333) + 5.0) + Sigma_l*delta*(nu*(-2.66666666666667*nu - 10.1666666666667) + 3.0)) + 1.0)) + 1.0)) +
        1.0) + 1.0);
    a_L_overvcubed = pow(v, 4)*(7.0*S_n + 3.0*Sigma_n*delta + pow(v, 2)*(S_n*(-9.66666666666667*nu - 10.0) +
        Sigma_n*delta*(-4.5*nu - 6.0) + pow(v, 2)*(S_n*(nu*(5.77777777777778*nu + 14.75) + 1.5) +
        Sigma_n*delta*(nu*(2.83333333333333*nu + 9.125) + 1.5))));
    OmegaLMag = a_L_overvcubed*gamma_L;
    Omega1Mag = pow(v, 5)*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu - 0.5625) + nu*(-0.0416666666666667*nu +
        1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) +
        0.1875) + 0.84375) + 0.5625) + 0.75);
    Omega2Mag = pow(v, 5)*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu + 0.5625) + nu*(-0.0416666666666667*nu +
        1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) +
        0.1875) + 0.84375) + 0.5625) + 0.75);
  }

  int TaylorT1RHS(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dEnergydv = v*(-1.0*nu + pow(v, 2)*(nu*(0.166666666666667*nu + 1.5) + v*(-11.6666666666667*S_l*nu -
        5.0*Sigma_l*delta*nu + v*(chi1_l*(chi1_l*(1.5*delta*nu + nu*(-3.0*nu + 1.5)) + 6.0*chi2_l*pow(nu, 2)) +
        chi1_la*(chi1_la*(-0.75*delta*nu + nu*(1.5*nu - 0.75)) - 3.0*chi2_la*pow(nu, 2)) +
        chi1_n*(chi1_n*(-0.75*delta*nu + nu*(1.5*nu - 0.75)) - 3.0*chi2_n*pow(nu, 2)) + pow(chi2_l, 2)*(-1.5*delta*nu +
        nu*(-3.0*nu + 1.5)) + pow(chi2_la, 2)*(0.75*delta*nu + nu*(1.5*nu - 0.75)) + pow(chi2_n, 2)*(0.75*delta*nu +
        nu*(1.5*nu - 0.75)) + nu*(-7.125*nu + 0.125*nu__2 + 10.125) + v*(S_l*nu*(23.7222222222222*nu - 38.5) +
        Sigma_l*delta*nu*(11.6666666666667*nu - 10.5) + v*(nu*(-154.898517962917*nu + 6.45833333333333*nu__2 +
        0.0270061728395062*nu__3 + 42.1875) + v*(S_l*nu*(412.875*nu - 10.875*nu__2 - 151.875) +
        Sigma_l*delta*nu*(175.5*nu - 5.625*nu__2 - 30.375) + v*(-298.666666666667*logv*pow(nu, 2) + nu*(-769.4015*nu +
        nu__2*(-0.012377829218107*nu__2 + 450.663995131026) - 0.870949074074074*nu__3 + 155.0390625) + pow(v,
        2)*(logv*pow(nu, 2)*(1523.2*nu + 1710.17142857143) + nu*(330.78*nu + 538.20703125) + pow(v,
        2)*(16016.0*logv*pow(nu, 2) + nu*(-4116.0*nu + 1808.9736328125)))))))))));
    const double Flux = pow(v, 10)*(6.4*nu__2 + pow(v, 2)*(-18.6666666666667*nu*nu__2 - 23.752380952381*nu__2 +
        v*(25.6*S_l*nu__2 - 8.0*Sigma_l*delta*nu__2 + 80.4247719318987*nu__2 + v*(chi1_l*(chi1_l*(6.6*delta*nu__2 -
        13.2*nu*nu__2 + 6.6*nu__2) + 24.8*chi2_l*nu*nu__2) + chi1_la*(chi1_la*(-2.96666666666667*delta*nu__2 +
        5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) - 13.7333333333333*chi2_la*nu*nu__2) +
        chi1_n*(chi1_n*(-2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) -
        13.7333333333333*chi2_n*nu*nu__2) + pow(chi2_l, 2)*(-6.6*delta*nu__2 - 13.2*nu*nu__2 + 6.6*nu__2) + pow(chi2_la,
        2)*(2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) + pow(chi2_n,
        2)*(2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) +
        117.726984126984*nu*nu__2 + nu__2*(23.1111111111111*nu__2 - 31.542151675485) + v*(S_l*(193.422222222222*nu*nu__2
        - 28.8*nu__2) + Sigma_l*delta*(68.8*nu*nu__2 - 5.2*nu__2) - 488.412937878093*nu*nu__2 - 245.074146910038*nu__2 +
        v*(-321.699087727595*S_l*nu__2 - 103.881997078702*Sigma_l*delta*nu__2 - 104.350476190476*logv*nu__2 -
        56.7811420312465*nu*nu__2 + nu__2*(-199.794708994709*nu__2 - 15.3086419753086*nu__3 + 945.576192944022) -
        204.89320622011*nu__2 + v*(S_l*(208.998941798942*nu*nu__2 + nu__2*(-666.074074074074*nu__2 + 448.343327454439))
        + Sigma_l*delta*(93.9174603174603*nu*nu__2 + nu__2*(-266.844444444444*nu__2 + 181.619047619048)) +
        2498.67153479682*nu*nu__2 + nu__2*(1285.7923710359*nu__2 - 649.661414142346) + v*(S_l*(3875.74795014869*nu*nu__2
        - 729.896693184029*nu__2) + Sigma_l*delta*(1302.34474121815*nu*nu__2 - 214.316458834892*nu__2) +
        337.555736961451*logv*nu__2 - 752.028100625135*nu__2 + v*(-1311.30675759439*logv*nu__2 + 4602.42139029395*nu__2
        + v*(746.4952102023*logv*nu__2 - 7788.20474442907*nu__2 + v*(3031.19666031508*logv*nu__2 +
        6137.18380876523*nu__2 + v*(3340.73494332217*logv*nu__2 + logv*(850.70483446712*logv*nu__2 -
        15505.402624526*nu__2) + 13022.6558856344*nu__2))))))))))));
    const double Absorption = pow(v, 15)*(chi1_l*pow(m1, 3)*(-4.8*chi1chi1*nu__2 - 1.6*nu__2) + chi2_l*pow(m2,
        3)*(-4.8*chi2chi2*nu__2 - 1.6*nu__2) + pow(v, 3)*(pow(m1, 4)*(chi1chi1*nu__2*(9.6*sqrt1Mchi1chi1 + 9.6) +
        nu__2*(3.2*sqrt1Mchi1chi1 + 3.2)) + pow(m2, 4)*(chi2chi2*nu__2*(9.6*sqrt1Mchi2chi2 + 9.6) +
        nu__2*(3.2*sqrt1Mchi2chi2 + 3.2))));
    const double dvdt = - (Flux + Absorption) / dEnergydv;
    if(dvdt<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt, t, y, dydt);
  }

  int TaylorT4RHS(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dvdt = pow(v, 9)*(pow(v, 2)*(-17.6*nu__2 + v*(v*(chi1_l*(-32.4*chi1_l*nu__2 + 63.2*chi2_l*nu__2) +
        chi1_la*(15.5333333333333*chi1_la*nu__2 - 32.9333333333333*chi2_la*nu__2) +
        chi1_n*(15.5333333333333*chi1_n*nu__2 - 32.9333333333333*chi2_n*nu__2) - 32.4*pow(chi2_l, 2)*nu__2 +
        15.5333333333333*pow(chi2_la, 2)*nu__2 + 15.5333333333333*pow(chi2_n, 2)*nu__2 - 2.93333333333333*nu*nu__2 +
        43.368253968254*nu__2 + v*(542.4*S_l*nu__2 + 224.8*Sigma_l*delta*nu__2 + chi1_l*pow(m1, 3)*(-4.8*chi1chi1*nu__2
        - 1.6*nu__2)/nu + chi2_l*pow(m2, 3)*(-4.8*chi2chi2*nu__2 - 1.6*nu__2)/nu - 475.008809222777*nu__2 +
        v*(chi1_l*(chi1_l*(-23.7*delta*nu__2 + 47.4*nu*nu__2 - 29.8428571428571*nu__2) +
        chi2_l*(-95.0666666666667*nu*nu__2 + 9.88571428571429*nu__2)) + chi1_la*(chi1_la*(11.9055555555556*delta*nu__2 -
        23.8111111111111*nu*nu__2 + 13.9769841269841*nu__2) + chi2_la*(47.3111111111111*nu*nu__2 -
        6.94285714285714*nu__2)) + chi1_n*(chi1_n*(11.9055555555556*delta*nu__2 - 23.8111111111111*nu*nu__2 +
        13.9769841269841*nu__2) + chi2_n*(47.3111111111111*nu*nu__2 - 6.94285714285714*nu__2)) + pow(chi2_l,
        2)*(23.7*delta*nu__2 + 47.4*nu*nu__2 - 29.8428571428571*nu__2) + pow(chi2_la, 2)*(-11.9055555555556*delta*nu__2
        - 23.8111111111111*nu*nu__2 + 13.9769841269841*nu__2) + pow(chi2_n, 2)*(-11.9055555555556*delta*nu__2 -
        23.8111111111111*nu*nu__2 + 13.9769841269841*nu__2) + nu*(-0.488888888888889*nu*nu__2 + 128.228042328042*nu__2)
        + nu__2*(1.78518518518519*nu__2 - 1058.43868227317) + (S_l*(572.444444444444*S_l*nu__2 +
        712.0*Sigma_l*delta*nu__2 - 1259.98809359975*nu__2) + Sigma_l*(200.0*Sigma_l*pow(delta, 2)*nu__2 -
        506.005856738196*delta*nu__2) + pow(chi1_l, 2)*(3.07142857142857*delta*nu__2 + 3.07142857142857*nu__2) +
        pow(chi1_la, 2)*(-1.03571428571429*delta*nu__2 - 1.03571428571429*nu__2) + pow(chi1_n,
        2)*(-1.03571428571429*delta*nu__2 - 1.03571428571429*nu__2) + pow(chi2_l, 2)*(-3.07142857142857*delta*nu__2 +
        3.07142857142857*nu__2) + pow(chi2_la, 2)*(1.03571428571429*delta*nu__2 - 1.03571428571429*nu__2) + pow(chi2_n,
        2)*(1.03571428571429*delta*nu__2 - 1.03571428571429*nu__2) - 104.350476190476*logv*nu__2 +
        nu__2*(-124.363756613757*nu__2 - 15.1358024691358*nu__3 + 1090.32725114508) - 204.89320622011*nu__2)/nu) +
        (-183.688888888889*S_l*nu__2 - 61.6380952380952*Sigma_l*delta*nu__2 - 124.43698901219*nu__2)/nu) + (pow(chi1_l,
        2)*(16.2*delta*nu__2 + 16.2*nu__2) + pow(chi1_la, 2)*(-7.76666666666667*delta*nu__2 - 7.76666666666667*nu__2) +
        pow(chi1_n, 2)*(-7.76666666666667*delta*nu__2 - 7.76666666666667*nu__2) + pow(chi2_l, 2)*(-16.2*delta*nu__2 +
        16.2*nu__2) + pow(chi2_la, 2)*(7.76666666666667*delta*nu__2 - 7.76666666666667*nu__2) + pow(chi2_n,
        2)*(7.76666666666667*delta*nu__2 - 7.76666666666667*nu__2) + nu__2*(23.9111111111111*nu__2 +
        12.0292768959436))/nu) + (-49.0666666666667*S_l*nu__2 - 40.0*Sigma_l*delta*nu__2 + 80.4247719318987*nu__2)/nu) -
        14.152380952381*nu__2/nu) + 6.4*nu__2/nu);
    if(dvdt<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt, t, y, dydt);
  }

  int TaylorT5RHS(double t, const double* y, double* dydt) {
    Recalculate(t, y);
    if(v>=1.0) { return GSL_EDOM; } // Beyond domain of PN validity
    const double dtdv = chi1_l*pow(m1, 3)*(0.1171875*chi1chi1*nu + 0.0390625*nu)/(nu__2*pow(v, 4)) + chi2_l*pow(m2,
        3)*(0.1171875*chi2chi2*nu + 0.0390625*nu)/(nu__2*pow(v, 4)) + (nu*(-3.2543041087963*nu + 0.454443876074735) -
        0.583767361111111*nu/pow(v, 2))/pow(v, 3) + (S_l*(-4.79166666666667*S_l*nu - 2.40885416666667*Sigma_l*delta*nu +
        0.654498469497875*nu) + Sigma_l*(1.220703125*Sigma_l*pow(delta, 2)*nu - 12.1900339943979*delta*nu) +
        chi1_l*(chi1_l*(delta*nu*(-1.5966796875*nu - 1.82416643415179) + nu*(nu*(3.193359375*nu + 2.05165318080357) -
        1.82416643415179)) + chi2_l*pow(nu, 2)*(-6.16536458333333*nu - 7.06531343005952)) +
        chi1_la*(chi1_la*(delta*nu*(0.752224392361111*nu + 0.863884819878472) + nu*(nu*(-1.50444878472222*nu -
        0.975545247395833) + 0.863884819878472)) + chi2_la*pow(nu, 2)*(3.26714409722222*nu + 3.72545030381944)) +
        chi1_n*(chi1_n*(delta*nu*(0.752224392361111*nu + 0.863884819878472) + nu*(nu*(-1.50444878472222*nu -
        0.975545247395833) + 0.863884819878472)) + chi2_n*pow(nu, 2)*(3.26714409722222*nu + 3.72545030381944)) +
        pow(chi2_l, 2)*(delta*nu*(1.5966796875*nu + 1.82416643415179) + nu*(nu*(3.193359375*nu + 2.05165318080357) -
        1.82416643415179)) + pow(chi2_la, 2)*(delta*nu*(-0.752224392361111*nu - 0.863884819878472) +
        nu*(nu*(-1.50444878472222*nu - 0.975545247395833) + 0.863884819878472)) + pow(chi2_n,
        2)*(delta*nu*(-0.752224392361111*nu - 0.863884819878472) + nu*(nu*(-1.50444878472222*nu - 0.975545247395833) +
        0.863884819878472)) + 2.54761904761905*logv*nu + nu*(nu*(nu*(3.65532769097222*nu - 0.798297861896495) +
        25.8462542787179) + 0.369526427469136*nu__3 - 1.55461636218371) + 5.00227554248315*nu +
        (S_l*nu*(-6.65364583333333*nu + 9.78252108134921) + Sigma_l*delta*nu*(-0.1171875*nu + 5.82380022321429) +
        nu*(0.797670009700533*nu - 5.64577976646101) + (chi1_l*(chi1_l*(-0.3955078125*delta*nu + nu*(0.791015625*nu -
        0.3955078125)) - 1.54296875*chi2_l*pow(nu, 2)) + chi1_la*(chi1_la*(0.189615885416667*delta*nu +
        nu*(-0.379231770833333*nu + 0.189615885416667)) + 0.804036458333333*chi2_la*pow(nu, 2)) +
        chi1_n*(chi1_n*(0.189615885416667*delta*nu + nu*(-0.379231770833333*nu + 0.189615885416667)) +
        0.804036458333333*chi2_n*pow(nu, 2)) + pow(chi2_l, 2)*(0.3955078125*delta*nu + nu*(0.791015625*nu -
        0.3955078125)) + pow(chi2_la, 2)*(-0.189615885416667*delta*nu + nu*(-0.379231770833333*nu + 0.189615885416667))
        + pow(chi2_n, 2)*(-0.189615885416667*delta*nu + nu*(-0.379231770833333*nu + 0.189615885416667)) +
        nu*(nu*(1.25325520833333*nu + 0.841548859126984) + 0.4703617648593) + (1.19791666666667*S_l*nu +
        0.9765625*Sigma_l*delta*nu - 1.96349540849362*nu + (nu*(0.4296875*nu + 0.345517113095238) + 0.15625*nu/pow(v,
        2))/v)/v)/v)/v)/(nu__2*pow(v, 3));
    const double dvdt = 1.0/dtdv;
    if(dvdt<0.0) { return GSL_EDIVERGE; } // v is decreasing
    return CommonRHS(dvdt, t, y, dydt);
  }

  int CommonRHS(const double dvdt, double t, const double* y, double* dydt) {
    dydt[0] = dvdt;
    dydt[1] = v*v*v;
    const double Omega1[3] = {Omega1Mag*Lhat_Nx, Omega1Mag*Lhat_Ny, Omega1Mag*Lhat_Nz};
    const double Omega2[3] = {Omega2Mag*Lhat_Nx, Omega2Mag*Lhat_Ny, Omega2Mag*Lhat_Nz};
    const double OmegaL[3]  = {OmegaLMag*nhat_x,  OmegaLMag*nhat_y,  OmegaLMag*nhat_z};
    cross(&Omega1[0], &y[2], &dydt[2]);
    cross(&Omega2[0], &y[5], &dydt[5]);
    cross(&OmegaL[0], &y[8], &dydt[8]);
    const Quaternions::Quaternion adot(0., dydt[8], dydt[9], dydt[10]);
    const Quaternions::Quaternion Rax =
      Quaternions::sqrtOfRotor(-Quaternions::normalized(Quaternions::Quaternion(0., OmegaL[0], OmegaL[1], OmegaL[2]))*zHat);
    const Quaternions::Quaternion Raxdot = ( (-1.0/std::sqrt(2+2*y[10]))*adot*zHat - (dydt[10]/(2+2*y[10]))*Rax );
    const Quaternions::Quaternion dgammadt = 2*(Rax.conjugate() * Raxdot * zHat);
    dydt[11] = dgammadt[0];
    return GSL_SUCCESS; // GSL expects this if everything went well
  }
};
