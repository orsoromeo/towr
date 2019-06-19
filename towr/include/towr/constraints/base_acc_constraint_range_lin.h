/******************************************************************************
Copyright (c) 2018, Alexander W. Winkler. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#ifndef TOWR_CONSTRAINTS_BASE_ACC_CONSTRAINT_RANGE_LIN_H_
#define TOWR_CONSTRAINTS_BASE_ACC_CONSTRAINT_RANGE_LIN_H_

#include <ifopt/problem.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/spline.h>
#include <towr/variables/nodes_variables_phase_based.h>
#include <towr/terrain/height_map.h>
#include <towr/models/dynamic_model.h>
#include <ifopt/constraint_set.h>
#include <ifopt/composite.h>
#include "time_discretization_constraint.h"
#include "geometry.h"
#include "derivative.h"
#include <towr/constraints/total_duration_constraint.h>
#include <iostream>
namespace towr {

/**
 * @brief Keeps the base acc in a specified range.
 *
 * 
 *
 * @ingroup Constraints
 */
class BaseAccConstraintRangeLin : public TimeDiscretizationConstraint {
public:
  /**
   * @brief Links the base variables and sets hardcoded bounds on the state.
   * @param T  The total time of the optimization horizon.
   * @param dt The discretization interval of the constraints.
   * @param spline_holder  Holds pointers to the base variables.
   */
  BaseAccConstraintRangeLin (const DynamicModel::Ptr &model,
                             double T,
                             double dt,
                             const NodeSpline::Ptr& spline,
                             std::string name,
                             const HeightMap::Ptr &terrain,
                             const SplineHolder& spline_holder,
                             int numberoflegs,
                             const std::vector<std::vector<double>> &phase_durations,
                             Geometry  geom,
                             Derivative  der);
  virtual ~BaseAccConstraintRangeLin () = default;

  void UpdateConstraintAtInstance (double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance (double t, int k, VecBound&) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string, Jacobian&) const override;
  std::vector<NodeSpline::Ptr> ee_motion_in_touch_;
  using Jac = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  Eigen::VectorXd ComputeWrench (State com, double t) const;
  VectorXd FillConstraint (State com, double t) const;
  NodeSpline::Jacobian FillJacobianLinWrenchWrtEENodes(double t, int ee) const;
  NodeSpline::Jacobian FillJacobianAngWrenchWrtEENodes(double t, int ee) const;
    NodeSpline::Jacobian FillJacobianEdgesWrtLambda (double t, int k) const;
private:
  std::vector<std::vector<double>> phase_durations_;
  int numberofleg_;
  std::vector<VecTimes> duration_;
  ifopt::Composite::Ptr variables_;
  int NumberNodes_;
  mutable DynamicModel::Ptr model_;
  double m_;
  double g_;
  Eigen::SparseMatrix<double, Eigen::RowMajor> I_b;
  NodeSpline::Ptr base_linear_;
  EulerConverter base_angular_;
  NodeSpline::Ptr lambda_;
  NodeSpline::Ptr spline_;        ///< a spline comprised of polynomials
  VecBound node_bounds_;
  std::string node_variables_id_;
  std::vector<NodeSpline::Ptr> ee_motion_;
  NodesVariablesPhaseBased::Ptr ee_force_;
  double mu_;
  Eigen::MatrixXd  LinearEdges_;
  Eigen::MatrixXd  AngularEdges_;
  Eigen::Vector3d base_;
  Geometry geom_;
  Derivative der_;
  int GetRow (int node, int dim) const;
  int GetNumberOfFeetInTouch (double t);
  NodeSpline::Jacobian FillJacobian(NodeSpline::Ptr spline_, double t) const;
  Eigen::Vector3d ComputeLinearWrench (State com) const;
  Eigen::Vector3d ComputeAngularWrench (Eigen::VectorXd acc, Jac I_w, State com) const;

  NodeSpline::Jacobian FillJacobianLinWrenchWrtLin(double t) const;
  NodeSpline::Jacobian FillJacobianAngWrenchWrtLin(double t, int k ) const;
  BaseAccConstraintRangeLin::Jac DerivativeOfrxma(double t) const;
  NodeSpline::Jacobian FillJacobianAngWrenchWrtAng(double t) const;
  bool IsInTouch(double t, int ee) const;

};

} /* namespace towr */

#endif /* TOWR_CONSTRAINTS_BASE_MOTION_CONSTRAINT_H_ */

