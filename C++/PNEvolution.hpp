#ifndef PNEVOLUTION_HPP
#define PNEVOLUTION_HPP

#include <vector>
#include "Quaternions.hpp"

namespace PostNewtonian {

  /// These functions evolve the PN orbital system, returning (by
  /// reference) the orbital speed, spin vectors, and frame
  /// information, as well as the completely redundant but
  /// occassionally useful phase information `Phi`.  Additional
  /// functions are given to calculate the basis vectors of the
  /// system's frame given the quaternion information on the frame.

  void EvolvePN(const std::string& Approximant,
		const double v_i, const double m1, const double m2,
		const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
		std::vector<double>& t, std::vector<double>& v,
		std::vector<std::vector<double> >& chi1, std::vector<std::vector<double> >& chi2,
		std::vector<Quaternions::Quaternion>& R_frame,
		std::vector<double>& Phi, std::vector<std::vector<double> >& L
		);

  void EvolvePN(const std::string& Approximant, const double PNOrder,
		const double v0, const double v_i,
		const double m1, const double m2,
		const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
		const Quaternions::Quaternion& R_frame_i,
		std::vector<double>& t, std::vector<double>& v,
		std::vector<std::vector<double> >& chi1, std::vector<std::vector<double> >& chi2,
		std::vector<Quaternions::Quaternion>& R_frame,
		std::vector<double>& Phi, std::vector<std::vector<double> >& L,
		const bool ForwardInTime=true
		);

  void EvolvePN_Q(const std::string& Approximant,
		  const double v_i, const double m1, const double m2,
		  const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
		  std::vector<double>& t, std::vector<double>& v,
		  std::vector<std::vector<double> >& chi1, std::vector<std::vector<double> >& chi2,
		  std::vector<Quaternions::Quaternion>& R_frame,
		  std::vector<double>& Phi, std::vector<std::vector<double> >& L
		  );

  void EvolvePN_Q(const std::string& Approximant, const double PNOrder,
		  const double v0, const double v_i,
		  const double m1, const double m2,
		  const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
		  const Quaternions::Quaternion& R_frame_i,
		  std::vector<double>& t, std::vector<double>& v,
		  std::vector<std::vector<double> >& chi1, std::vector<std::vector<double> >& chi2,
		  std::vector<Quaternions::Quaternion>& R_frame,
		  std::vector<double>& Phi, std::vector<std::vector<double> >& L,
		  const bool ForwardInTime=true
		  );

  std::vector<std::vector<double> > ellHat(const std::vector<Quaternions::Quaternion>& R);
  std::vector<std::vector<double> > nHat(const std::vector<Quaternions::Quaternion>& R);
  std::vector<std::vector<double> > lambdaHat(const std::vector<Quaternions::Quaternion>& R);

};



#endif // PNEVOLUTION_HPP
