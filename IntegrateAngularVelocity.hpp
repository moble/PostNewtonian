// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#ifndef INTEGRATEANGULARVELOCITY_HPP
#define INTEGRATEANGULARVELOCITY_HPP

#include <vector>
#include "Quaternions.hpp"

#define Quaternion_Epsilon 1.0e-14

namespace Quaternions {

  std::vector<Quaternion> FrameFromAngularVelocity(const std::vector<Quaternion>& Omega, const std::vector<double>& T);

} // namespace Quaternions

#endif // INTEGRATEANGULARVELOCITY_HPP
