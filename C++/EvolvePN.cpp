#include <vector>
#include <iostream>
#include <iomanip>
#include <string>

#include "PNEvolution.hpp"


int main() {
  const std::string Approximant = "TaylorT4";
  const double PNOrder = 3.5;
  const double v0 = 0.3073533404944769;
  const double v_i = v0;
  const double m1 = 6./7.;
  std::vector<double> chi1_i(3);
  std::vector<double> chi2_i(3);
  const Quaternions::Quaternion R_frame_i(1,0,0,0);
  std::vector<double> t;
  std::vector<double> v;
  std::vector<std::vector<double> > chi1;
  std::vector<std::vector<double> > chi2;
  std::vector<Quaternions::Quaternion> R_frame;
  std::vector<double> Phi;

  chi1_i[0] = 7.393069129077587398e-01;
  chi1_i[1] = 1.901752366053321708e-01;
  chi1_i[2] = -4.953378721708571741e-01;
  chi2_i[0] = -1.896513685161916041e-01;
  chi2_i[1] = 5.104595303989780536e-02;
  chi2_i[2] = -2.268760863004861961e-01;

  PostNewtonian::EvolvePN(Approximant, PNOrder,
			  v0, v_i,
			  m1,
			  chi1_i, chi2_i,
			  R_frame_i,
			  t, v,
			  chi1, chi2,
			  R_frame,
			  Phi
			  );

  std::cerr << "NTimes = " << t.size()
	    << "\nPrinting results to std::cout... "
	    << std::flush;

  std::cout << std::setprecision(12)
	    << "## const std::string Approximant = \"" << Approximant << "\";\n"
	    << "## const double PNOrder = " << PNOrder << ";\n"
	    << "## const double v0 = " << v0 << ";\n"
	    << "## const double v_i = " << v_i << ";\n"
	    << "## const double m1 = " << m1 << ";\n"
	    << "## std::vector<double> chi1_i(3);\n"
	    << "## std::vector<double> chi2_i(3);\n"
	    << "## const Quaternions::Quaternion R_frame_i("
	    << R_frame_i[0] << ", " << R_frame_i[1] << ", " << R_frame_i[2] << ", " << R_frame_i[3] << ");\n"
	    << "## chi1_i[0] = " << chi1_i[0] << ";\n"
	    << "## chi1_i[1] = " << chi1_i[1] << ";\n"
	    << "## chi1_i[2] = " << chi1_i[2] << ";\n"
	    << "## chi2_i[0] = " << chi2_i[0] << ";\n"
	    << "## chi2_i[1] = " << chi2_i[1] << ";\n"
	    << "## chi2_i[2] = " << chi2_i[2] << ";\n"
	    << "# [1] = t/m\n"
	    << "# [2] = v\n"
	    << "# [3] = chi1_x\n"
	    << "# [4] = chi1_y\n"
	    << "# [5] = chi1_z\n"
	    << "# [6] = chi2_x\n"
	    << "# [7] = chi2_y\n"
	    << "# [8] = chi2_z\n"
	    << "# [9] = R_frame[0]\n"
	    << "# [10] = R_frame[1]\n"
	    << "# [11] = R_frame[2]\n"
	    << "# [12] = R_frame[3]\n"
	    << "# [13] = Phi\n";
  const unsigned int tsize = t.size();
  for(unsigned int i=0; i<tsize; ++i) {
    std::cout << t[i] << " " << v[i]
	      << " " << chi1[i][0] << " " << chi1[i][1] << " " << chi1[i][2]
	      << " " << chi2[i][0] << " " << chi2[i][1] << " " << chi2[i][2]
	      << " " << R_frame[i][0] << " " << R_frame[i][1] << " " << R_frame[i][2] << " " << R_frame[i][3]
	      << " " << Phi[i] << std::endl;
  }

  std::cerr << "succeeded." << std::endl;

  return 0;
}
