// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_EIGEN_SOLVER_HPP
#define PALACE_DRIVERS_EIGEN_SOLVER_HPP

#include <complex>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "utils/configfile.hpp"

namespace palace
{

class ErrorIndicator;
class LumpedPortOperator;
class Mesh;
class PostOperator;

//
// Driver class for eigenmode simulations.
//
class EigenSolver : public BaseSolver
{
private:
  void Postprocess(const PostOperator &post_op, const LumpedPortOperator &lumped_port_op,
                   int i, std::complex<double> omega, double error_bkwd, double error_abs,
                   int num_conv, double E_elec, double E_mag,
                   const ErrorIndicator *indicator) const;

  void PostprocessEigen(int i, std::complex<double> omega, double error_bkwd,
                        double error_abs, int num_conv) const;

  void PostprocessPorts(const PostOperator &post_op,
                        const LumpedPortOperator &lumped_port_op, int i) const;

  void PostprocessEPR(const PostOperator &post_op, const LumpedPortOperator &lumped_port_op,
                      int i, std::complex<double> omega, double E_m) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

  // Store Josephson configuration
  JosephsonConfig josephson_config;
  // Store previous field values for convergence check
  mutable std::vector<std::complex<double>> previous_field_values;

public:
  using BaseSolver::BaseSolver;

  // Add setter for Josephson configuration
  void SetJosephsonElements(const std::vector<mfem::Vector>& locations, 
                          double conv_threshold = 1e-6);
  // Add field convergence check
  bool CheckFieldConvergence() const;
  // Add field value getter
  std::vector<std::complex<double>> GetFieldValuesAtJosephsonElements() const;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_EIGEN_SOLVER_HPP
