#include <vector>
#include <iostream>
#include <iomanip>
#include <string>

#include "PNEvolution.hpp"


int main() {
  const std::string Approximant = "TaylorT1";
  const double PNOrder = 3.5;
  const double v0 = 0.05;
  const double v_i = v0;
  const double m1 = 0.5;
  std::vector<double> chi1_i(3);
  std::vector<double> chi2_i(3);
  const Quaternions::Quaternion R_frame_i(1,0,0,0);
  std::vector<double> t;
  std::vector<double> v;
  std::vector<std::vector<double> > chi1;
  std::vector<std::vector<double> > chi2;
  std::vector<Quaternions::Quaternion> R_frame;
  std::vector<double> Phi;

  chi1_i[0] = 0.0;
  chi1_i[1] = 0.0;
  chi1_i[2] = 0.7;
  chi2_i[0] = 0.0;
  chi2_i[1] = 0.0;
  chi2_i[2] = 0.7;

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
	    << "\nPrinting results to std::cout...\n"
	    << std::endl;

  std::cout << "# [1] = t/m\n"
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
	    << "# [13] = Phi\n"
	    << std::setprecision(12);
  const unsigned int tsize = t.size();
  for(unsigned int i=0; i<tsize; ++i) {
    std::cout << t[i] << " " << v[i]
	      << " " << chi1[i][0] << " " << chi1[i][1] << " " << chi1[i][2]
	      << " " << chi2[i][0] << " " << chi2[i][1] << " " << chi2[i][2]
	      << " " << R_frame[i][0] << " " << R_frame[i][1] << " " << R_frame[i][2] << " " << R_frame[i][3]
	      << " " << Phi[i] << std::endl;
  }


  return 0;
}
